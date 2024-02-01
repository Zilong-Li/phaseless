#include "fastphase.hpp"

#include "common.hpp"
#include "io.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>
#include <stdexcept>

using namespace std;

void FastPhaseK2::initRecombination(const Int2D & pos, std::string rfile, int B, double Ne)
{
    nGen = 4 * Ne / C;
    int nchunks = pos.size();
    pos_chunk.resize(nchunks + 1);
    int i{0}, ss{0};
    dist.reserve(M);
    for(i = 0; i < nchunks; i++)
    {
        pos_chunk[i] = ss;
        auto tmp = calc_position_distance(pos[i]);
        dist.insert(dist.end(), tmp.begin(), tmp.end());
        R.middleCols(ss, pos[i].size()) = calc_transRate_diploid(tmp, nGen);
        ss += pos[i].size();
    }
    pos_chunk[nchunks] = ss; // add sentinel
    if(!rfile.empty()) load_csv(R, rfile, true);
    er = R.row(0).sqrt();
    protect_er(er);
    R = er2R(er);
}

void FastPhaseK2::setFlags(double tol_p, double tol_f, double tol_q, bool debug_, bool nQ, bool nP, bool nF, bool nR)
{
    alleleEmitThreshold = tol_p;
    clusterFreqThreshold = tol_f;
    admixtureThreshold = tol_q;
    debug = debug_;
    NQ = nQ;
    NP = nP;
    NF = nF;
    NR = nR;
}

void FastPhaseK2::refillHaps()
{
    int s{0};
    for(int c = 0; c < C; c++)
    {
        // bin hapsum per 100 snps ?
        for(int m = 0; m < M; m++)
        {
            if(HapSum(c, m) >= minHapfreq) continue;
            MyArr1D h = HapSum.col(m);
            h(c) = 0; // do not re-sample current
            h /= h.sum();
            MyFloat1D p(h.data(), h.data() + h.size());
            std::discrete_distribution<int> distribution{p.begin(), p.end()};
            int choice = distribution(rng);
            assert(choice != c);
            h.maxCoeff(&choice); // if no binning, this may be better
            F(m, c) = F(m, choice);
            s++;
        }
    }
    cao.warn("refill ", 100 * s / (C * M), "% infrequently used haps");
}

void FastPhaseK2::initIteration()
{
    // initial temp variables
    Ezj.setZero(C, M); // reset post(Z,j)
    Ezg1.setZero(C, M); // reset pos(Z,g)
    Ezg2.setZero(C, M); // reset pos(Z,g)
    HapSum.setZero(C, M); // reset post(Z,j)
}

void FastPhaseK2::updateIteration()
{
    // update R
    if(!NR) er = 1.0 - Ezj.colwise().sum() / N;
    // update F
    if(!NP) F = (Ezg2 / (Ezg1 + Ezg2)).transpose();
    // update PI
    if(!NF)
    {
        PI = Ezj;
        PI.rowwise() /= PI.colwise().sum();
    }
    protectPars();
}

void FastPhaseK2::protectPars()
{
    // protect F
    if(!NP)
    {
        if(F.isNaN().any())
        {
            cao.warn("NaN in F in FastPhaseK2 model. will fill it with AF");
            if(AF.size() == 0) cao.error("AF is not assigned!\n");
            for(int i = 0; i < M; i++) F.row(i) = F.row(i).isNaN().select(AF(i), F.row(i));
        }
        // map F to domain but no normalization
        F = (F < alleleEmitThreshold).select(alleleEmitThreshold, F); // lower bound
        F = (F > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, F); // upper bound
    }
    // protect R
    if(!NR)
    {
        for(int i = 0; i < er.size(); i++)
        {
            const double miner = std::exp(-nGen * maxRate * dist[i] / 100 / 1e6);
            const double maxer = std::exp(-nGen * minRate * dist[i] / 100 / 1e6);
            er(i) = er(i) < miner ? miner : er(i);
            er(i) = er(i) > maxer ? maxer : er(i);
        }
        protect_er(er);
        R = er2R(er);
    }
    // protect PI
    if(!NF)
    {
        if(PI.isNaN().any()) cao.warn("NaN in PI. reset cluster frequency to ", clusterFreqThreshold);
        PI = (PI < clusterFreqThreshold).select(clusterFreqThreshold, PI);
        PI = (PI > 1 - clusterFreqThreshold).select(1 - clusterFreqThreshold, PI);
        // re-normalize F per site. hope should work okay. otherwise do the complicated.
        PI.rowwise() /= PI.colwise().sum();
    }
    // norm HapSum
    HapSum.rowwise() /= HapSum.colwise().sum();
}

/*
** @param GL        genotype likelihood of all individuals in snp major form
** @param ind       current individual i
** @param finalIter boolean, call genotype of not
** @return individual total likelihood
*/
double FastPhaseK2::hmmIterWithJumps(const MyFloat1D & GL, const int ic, const int ind, bool finalIter)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * S * 3, S, 3);
    MyArr2D emit = get_emission_by_gl(gli, F.middleRows(pos_chunk[ic], S)).transpose(); // CC x S
    const auto [alpha, beta, cs] =
        forward_backwards_diploid(emit, R.middleCols(pos_chunk[ic], S), PI.middleCols(pos_chunk[ic], S));
    if(!((1 - ((alpha * beta).colwise().sum())).abs() < 1e-9).all())
        cao.error((alpha * beta).colwise().sum(), "\ngamma sum is not 1.0!\n");
    // now get posterios
    MyArr2D ind_post_zg1(C, S), ind_post_zg2(C, S), ind_post_zj(C, S), gammaC(C, S);
    MyArr1D gamma_div_emit(CC), beta_mult_emit(CC);
    MyArr1D alphatmp(C);
    int z1, m, s;
    for(s = 0; s < S; s++)
    {
        m = s + pos_chunk[ic];
        gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // C2
        gammaC.col(s) = (alpha.col(s) * beta.col(s)).reshaped(C, C).colwise().sum();
        for(z1 = 0; z1 < C; z1++)
        {
            ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(m, z1))
                                   * (gli(s, 0) * (1 - F.row(m)) + gli(s, 1) * F.row(m)).transpose())
                                      .sum();
            ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(m, z1))
                                   * (gli(s, 1) * (1 - F.row(m)) + gli(s, 2) * F.row(m)).transpose())
                                      .sum();
            if(s > 0) alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum() * R(1, m);
            if(s == 0) ind_post_zj(z1, s) = (alpha.col(0) * beta.col(0)).segment(z1 * C, C).sum();
            if(finalIter) callGenoLoopC(z1, m, ind, gli.row(s), F.row(m), gamma_div_emit);
        }
        if(s == 0) continue;
        alphatmp += PI.col(m) * R(2, m) * 1.0; // inner alpha.col(s-1).sum == 1
        beta_mult_emit = emit.col(s) * beta.col(s); // C2
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1, s) = cs(s) * (PI(z1, m) * alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
    }
    { // sum over all samples for updates
        std::scoped_lock<std::mutex> lock(mutex_it);
        Ezj.middleCols(pos_chunk[ic], S) += ind_post_zj;
        Ezg1.middleCols(pos_chunk[ic], S) += ind_post_zg1;
        Ezg2.middleCols(pos_chunk[ic], S) += ind_post_zg2;
        HapSum.middleCols(pos_chunk[ic], S) += gammaC;
    }

    return (1 / cs).log().sum();
}

void FastPhaseK2::callGenoLoopC(int z1,
                                int m,
                                int ind,
                                const MyArr1D & gli,
                                const MyArr1D & p,
                                const MyArr1D & gamma_div_emit)
{
    MyArr1D tmp_zg(4);
    for(int z2 = 0; z2 < C; z2++)
    {
        int z12 = z1 * C + z2;
        tmp_zg(0) = gli(0) * (1 - p(z1)) * (1 - p(z2));
        tmp_zg(1) = gli(1) * (1 - p(z1)) * p(z2);
        tmp_zg(2) = gli(1) * p(z1) * (1 - p(z2));
        tmp_zg(3) = gli(2) * p(z1) * p(z2);
        GP(3 * m + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
        GP(3 * m + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
        GP(3 * m + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
    }
}

double FastPhaseK2::runAllChunks(const MyFloat2D & GL, const int ind, bool finalIter)
{
    if(pos_chunk.size() == 0) cao.error("please run initRecombination first");
    double loglike{0};
    for(size_t ic = 0; ic < GL.size(); ic++)
    {
        loglike += hmmIterWithJumps(GL[ic], ic, ind, finalIter);
    }
    return loglike;
}

int run_impute_main(Options & opts)
{
    cao.cao.open(opts.out + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running fastphase");
    int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < allthreads ? opts.nthreads : allthreads;
    cao.print(tim.date(), allthreads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool pool(opts.nthreads);

    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    init_bigass(genome, opts);
    Eigen::IOFormat fmt(6, Eigen::DontAlignCols, " ", "\n");
    vector<future<double>> res;
    FastPhaseK2 faith(genome->nsamples, genome->nsnps, opts.C, opts.seed);
    faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.nP, opts.nF, opts.nR);
    faith.initRecombination(genome->pos, opts.in_rfile);
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; SIG_COND && it <= opts.nimpute; it++)
    {
        tim.clock();
        if(it > 4 && it < opts.nimpute / 2 && it % 4 == 1) faith.refillHaps();
        faith.initIteration();
        for(int i = 0; i < faith.N; i++)
            res.emplace_back(
                pool.enqueue(&FastPhaseK2::runAllChunks, &faith, std::ref(genome->gls), i, it == opts.nimpute));
        loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        faith.updateIteration();
        diff = it ? loglike - prevlike : NAN;
        prevlike = loglike;
        cao.print(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff =", diff, ", time",
                  tim.reltime(), " sec");
    }
    // reuse Ezj for AE
    if(opts.eHap)
    {
        cao.print("use hapsum");
        faith.Ezj.setZero(faith.CC, faith.HapSum.cols());
        for(int m = 0; m < faith.HapSum.cols(); m++)
            faith.Ezj.col(m) =
                (faith.HapSum.col(m).matrix() * faith.HapSum.col(m).transpose().matrix()).reshaped().array();
    }
    else
    {
        faith.Ezj = get_cluster_frequency(faith.R, faith.PI);
    }
    auto bw = make_bcfwriter(opts.out + ".vcf.gz", genome->chrs, genome->sampleids);
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        const int S = faith.pos_chunk[ic + 1] - faith.pos_chunk[ic];
        MyArr2D out = faith.Ezj.middleCols(faith.pos_chunk[ic], S);
        genome->AE.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.R.middleCols(faith.pos_chunk[ic], S);
        genome->R.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.PI.middleCols(faith.pos_chunk[ic], S);
        genome->PI.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.F.middleRows(faith.pos_chunk[ic], S);
        genome->F.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.GP.middleRows(faith.pos_chunk[ic], S * 3);
        write_bigass_to_bcf(bw, out.data(), genome->chrs[ic], genome->pos[ic]);
    }
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::ofstream ofs(opts.out + ".pars.bin", std::ios::out | std::ios::binary);
    auto bytes_written = alpaca::serialize<OPTIONS, BigAss>(*genome, ofs);
    ofs.close();
    assert(std::filesystem::file_size(opts.out + ".pars.bin") == bytes_written);
    cao.done(tim.date(), "imputation done and outputting.", bytes_written, " bytes written to file");
    std::ofstream orecomb(opts.out + ".recomb");
    orecomb << faith.R.transpose().format(fmt) << "\n";
    std::ofstream opi(opts.out + ".pi");
    opi << faith.PI.transpose().format(fmt) << "\n";
    std::ofstream ohap(opts.out + ".hapfreq");
    ohap << faith.HapSum.transpose().format(fmt) << "\n";
    std::ofstream oae(opts.out + ".ae");
    oae << faith.Ezj.transpose().format(fmt) << "\n";
    return 0;
}
