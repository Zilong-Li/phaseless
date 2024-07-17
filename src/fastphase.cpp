#include "fastphase.hpp"

#include "common.hpp"
#include "io.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>
#include <stdexcept>

using namespace std;

void FastPhaseK2::initRecombination(const Int2D & pos, std::string rfile, int B_, double Ne)
{
    nGen = 4 * Ne / C;
    B = B_;
    if(B == 2) cao.error("-B can not be 2");
    int nchunks = pos.size();
    pos_chunk.resize(nchunks + 1);
    grid_chunk.resize(nchunks + 1);
    int i{0}, ss{0}, sg{0};
    // how many grids in total
    for(i = 0, G = 0; i < nchunks; i++) G += B > 1 ? (pos[i].size() + B - 1) / B : pos[i].size();
    R = MyArr2D(3, G);
    PI = MyArr2D::Ones(C, G);
    PI.rowwise() /= PI.colwise().sum(); // normalize it per site
    dist.reserve(G);
    Int1D tmpdist;
    for(i = 0, ss = 0, sg = 0; i < nchunks; i++)
    {
        pos_chunk[i] = ss;
        grid_chunk[i] = sg;
        collapse.segment(ss, pos[i].size()) = find_grid_to_collapse(pos[i], B);
        tmpdist = calc_grid_distance(pos[i], collapse.segment(ss, pos[i].size()));
        dist.insert(dist.end(), tmpdist.begin(), tmpdist.end());
        R.middleCols(sg, tmpdist.size()) = calc_transRate_diploid(tmpdist, nGen);
        ss += pos[i].size();
        sg += tmpdist.size();
    }
    pos_chunk[nchunks] = ss; // add sentinel
    grid_chunk[nchunks] = sg; // add sentinel
    if(!rfile.empty()) load_csv(R, rfile, true);
    er = R.row(0).sqrt();
    protect_er(er);
    R = er2R(er);
}

void FastPhaseK2::setFlags(double tol_p,
                           double tol_f,
                           double tol_q,
                           bool debug_,
                           bool nQ,
                           bool nP,
                           bool nF,
                           bool nR)
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

void FastPhaseK2::refillHaps(int strategy)
{
    int s{0}, c{0}, ic{0}, g{0}, i{0}, sg{0};
    int nchunks = pos_chunk.size() - 1;
    // bin hapsum per 100 snps ?
    for(c = 0; c < C; c++)
    {
        for(ic = 0, sg = 0; ic < nchunks; ic++)
        {
            const int S = pos_chunk[ic + 1] - pos_chunk[ic];
            const auto se = find_grid_start_end(collapse.segment(pos_chunk[ic], S));
            for(g = 0; g < (int)se.size(); g++, sg++)
            {
                if(HapSum(c, sg) >= minHapfreq) continue;
                MyArr1D h = HapSum.col(sg);
                h(c) = 0; // do not re-sample current
                h /= h.sum();
                MyFloat1D p(h.data(), h.data() + h.size());
                std::discrete_distribution<int> distribution{p.begin(), p.end()};
                int choice = distribution(rng);
                assert(choice != c);
                // now go through all sites in the grid
                for(i = se[g][0]; i <= se[g][1]; i++)
                {
                    if(strategy == 1)
                    {
                        P(i + pos_chunk[ic], c) = alleleEmitThreshold;
                    }
                    else if(strategy == 2)
                    {
                        h.maxCoeff(&choice); // if no binning, this may be better
                        P(i + pos_chunk[ic], c) = P(i + pos_chunk[ic], choice);
                    }
                    else
                    {
                        P(i + pos_chunk[ic], c) = P(i + pos_chunk[ic], choice);
                    }
                    s++;
                }
            }
        }
    }
    cao.warn("refill ", 100 * s / (C * M), "% infrequently used haps");
}

void FastPhaseK2::initIteration()
{
    // initial temp variables
    Ezg1.setZero(C, M); // reset pos(Z,g)
    Ezg2.setZero(C, M); // reset pos(Z,g)
    Ezj.setZero(C, G); // reset post(Z,j)
    HapSum.setZero(C, G); // reset post(Z,j)
}

void FastPhaseK2::updateIteration()
{
    // update R
    if(!NR) er = 1.0 - Ezj.colwise().sum() / N;
    // update F
    if(!NP) P = (Ezg2 / (Ezg1 + Ezg2)).transpose();
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
        if(P.isNaN().any())
        {
            cao.warn("NaN in F in FastPhaseK2 model. will fill it with AF");
            if(AF.size() == 0) cao.error("AF is not assigned!\n");
            for(int i = 0; i < M; i++) P.row(i) = P.row(i).isNaN().select(AF(i), P.row(i));
        }
        // map F to domain but no normalization
        P = (P < alleleEmitThreshold).select(alleleEmitThreshold, P); // lower bound
        P = (P > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, P); // upper bound
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
    const int nGrids = grid_chunk[ic + 1] - grid_chunk[ic];
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * S * 3, S, 3);
    int start = pos_chunk[ic], nsize = S;
    if(nGrids != S)
    {
        start = grid_chunk[ic];
        nsize = nGrids;
    }
    MyArr2D emit_grid =
        get_emission_by_grid(gli, P.middleRows(pos_chunk[ic], S), collapse.segment(pos_chunk[ic], S));
    MyArr2D emit = get_emission_by_gl(gli, P.middleRows(pos_chunk[ic], S));
    const auto [alpha, beta, cs] =
        forward_backwards_diploid(emit_grid, R.middleCols(start, nsize), PI.middleCols(start, nsize));
    if(!((1 - ((alpha * beta).colwise().sum())).abs() < 1e-9).all())
        cao.error((alpha * beta).colwise().sum(), "\ngamma sum is not 1.0!\n");
    // now get posterios
    MyArr2D ind_post_zg1(C, S), ind_post_zg2(C, S), ind_post_zj(C, nGrids), gammaC(C, nGrids);
    MyArr1D gamma_div_emit(CC), beta_mult_emit(CC), igamma(CC);
    MyArr1D alphatmp(C);
    int z1, m, s, g{0}, gg{0};
    const auto se = find_grid_start_end(collapse.segment(pos_chunk[ic], S));
    for(g = 0; g < nGrids; g++)
    {
        gg = g + grid_chunk[ic];
        gammaC.col(g) = (alpha.col(g) * beta.col(g)).reshaped(C, C).colwise().sum();
        igamma = alpha.col(g) * beta.col(g);
        if(B == 1) gamma_div_emit = igamma / emit_grid.col(g);
        for(z1 = 0; z1 < C; z1++)
        {
            for(s = se[g][0]; s <= se[g][1]; s++)
            {
                m = s + pos_chunk[ic];
                // if(B > 1) gamma_div_emit = igamma / get_emission_by_site(gli.row(s), P.row(m));
                if(B > 1) gamma_div_emit = igamma / emit.col(s);
                ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - P(m, z1))
                                       * (gli(s, 0) * (1 - P.row(m)) + gli(s, 1) * P.row(m)).transpose())
                                          .sum();
                ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (P(m, z1))
                                       * (gli(s, 1) * (1 - P.row(m)) + gli(s, 2) * P.row(m)).transpose())
                                          .sum();
                if(finalIter) callGenoLoopC(z1, m, ind, gli.row(s), P.row(m), gamma_div_emit);
            }
            if(g == 0) ind_post_zj(z1, g) = (alpha.col(g) * beta.col(g)).segment(z1 * C, C).sum();
            if(g > 0) alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), g - 1).sum() * R(1, gg);
        }
        if(g == 0) continue;
        alphatmp += PI.col(gg) * R(2, gg) * 1.0; // inner alpha.col(s-1).sum == 1
        beta_mult_emit = emit_grid.col(g) * beta.col(g); // C2
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1, g) =
                cs(g) * (PI(z1, gg) * alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
    }
    { // sum over all samples for updates
        std::scoped_lock<std::mutex> lock(mutex_it);
        Ezg1.middleCols(pos_chunk[ic], S) += ind_post_zg1;
        Ezg2.middleCols(pos_chunk[ic], S) += ind_post_zg2;
        Ezj.middleCols(grid_chunk[ic], nGrids) += ind_post_zj;
        HapSum.middleCols(grid_chunk[ic], nGrids) += gammaC;
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
    cao.print(tim.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
              ", M =", genome->nsnps, ", nchunks =", genome->nchunks, ", B =", opts.gridsize,
              ", G =", genome->G, ", seed =", opts.seed);
    // estimate the RAM
    // input data : N x M x 3 x float
    // output data : N x M x 3 x float
    // model :C x M x 4 +  (CxCxMx3 + C x M x 4) x threads
    double ram = (double)(genome->nsamples * genome->nsnps * 6
                          + genome->C * genome->chunksize * 4 * (opts.nthreads + 1)
                          + genome->C * genome->C * genome->chunksize * 3 * opts.nthreads)
                 * sizeof(MyFloat) * 8 / (1024 * 1024 * 1024);
    cao.print(tim.date(), "roughly estimated RAM usage would be ", ram, " Gbs");
    vector<future<double>> res;
    FastPhaseK2 faith(genome->nsamples, genome->nsnps, opts.C, opts.seed);
    faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.nP, opts.nF, opts.nR);
    faith.initRecombination(genome->pos, opts.in_rfile, opts.gridsize);
    genome->G = faith.G;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; SIG_COND && it <= opts.nimpute; it++)
    {
        tim.clock();
        if(opts.refillHaps && it > 4 && it < opts.nimpute / 2 && it % 4 == 1)
            faith.refillHaps(opts.refillHaps);
        faith.initIteration();
        for(int i = 0; i < faith.N; i++)
            res.emplace_back(pool.enqueue(&FastPhaseK2::runAllChunks, &faith, std::ref(genome->gls), i,
                                          it == opts.nimpute));
        loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        faith.updateIteration();
        diff = it ? loglike - prevlike : NAN;
        prevlike = loglike;
        cao.print(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff =", diff,
                  ", time", tim.reltime(), " sec");
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
    genome->collapse = Char1D(faith.collapse.data(), faith.collapse.data() + faith.collapse.size());
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        const int S = faith.pos_chunk[ic + 1] - faith.pos_chunk[ic];
        const int G = faith.grid_chunk[ic + 1] - faith.grid_chunk[ic];
        MyArr2D out = faith.Ezj.middleCols(faith.grid_chunk[ic], G);
        genome->AE.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.R.middleCols(faith.grid_chunk[ic], G);
        genome->R.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.PI.middleCols(faith.grid_chunk[ic], G);
        genome->PI.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.P.middleRows(faith.pos_chunk[ic], S);
        genome->P.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
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
    orecomb << faith.R.transpose().format(fmt6) << "\n";
    std::ofstream opi(opts.out + ".pi");
    opi << faith.PI.transpose().format(fmt6) << "\n";
    std::ofstream ohap(opts.out + ".hapfreq");
    ohap << faith.HapSum.transpose().format(fmt6) << "\n";
    std::ofstream oae(opts.out + ".ae");
    oae << faith.Ezj.transpose().format(fmt6) << "\n";
    return 0;
}
