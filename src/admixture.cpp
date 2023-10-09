#include "admixture.hpp"

#include "common.hpp"
#include "io.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>
#include <fstream>

using namespace std;

// the complexity of this version should be O(CC + 2CK)
double Admixture::runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D kapa, Ekg;
    MyArr2D alpha, beta, ae;
    MyArr1D iQ = MyArr1D::Zero(K);
    MyArr1D Hz(C);
    double norm = 0, llike = 0, tmp = 0;
    int c1, k1, s, c2, c12;
    for(int ic = 0, m = 0; ic < genome->nchunks; ic++)
    {
        const int nsnps = genome->pos[ic].size();
        const int nGrids = genome->B > 1 ? (nsnps + genome->B - 1) / genome->B : nsnps;
        alpha.setZero(C * C, nGrids);
        beta.setZero(C * C, nGrids);
        get_cluster_likelihood(ind, nsnps, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                               genome->F[ic]); // return gamma
        ae.setZero(C * C, nGrids);
        get_cluster_pairs_probabity(ae, genome->R[ic], genome->PI[ic]);
        kapa.setZero(C * K, nGrids); // C x K x M layout
        Ekg.setZero(K, nGrids);
        for(s = 0; s < nGrids; s++, m++)
        {
            for(c1 = 0; c1 < C; c1++) Hz(c1) = (Q.col(ind) * F(Eigen::seqN(c1, K, C), m)).sum();
            for(norm = 0, c1 = 0; c1 < C; c1++)
            {
                for(tmp = 0, c2 = 0; c2 < C; c2++)
                {
                    c12 = c1 * C + c2;
                    double xz = alpha(c12, s) * beta(c12, s) / ae(c12, s);
                    double zy = Hz(c1) * Hz(c2);
                    tmp += xz * zy;
                }
                norm += tmp;
                kapa(Eigen::seqN(c1, K, C), s) = (Q.col(ind) * F(Eigen::seqN(c1, K, C), m)) * tmp / Hz(c1);
            }
            llike += log(norm);
            kapa.col(s) /= kapa.col(s).sum();
            for(k1 = 0; k1 < K; k1++) Ekg(k1, s) = 2 * kapa.middleRows(k1 * C, C).col(s).sum();
        }
        iQ += Ekg.rowwise().sum();
        { // for update F
            std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(m - nGrids, nGrids) += 2 * kapa;
            NormF.middleCols(m - nGrids, nGrids) += Ekg;
        }
    }
    // update Q, iQ.sum() should be 2M
    if(!nonewQ) Q.col(ind) = iQ / (2 * M);

    return llike;
}

// the complexity of this version is O(CCKK)
double Admixture::runNativeWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D w((C * C + C) / 2, K * K);
    MyArr2D Ekg, iEkc, alpha, beta, ae;
    double norm = 0, llike = 0;
    int c1, c2, c12, cc;
    int k1, k2, k12, s;
    MyArr1D iQ = MyArr1D::Zero(K);
    for(int ic = 0, m = 0; ic < genome->nchunks; ic++)
    {
        const int nsnps = genome->pos[ic].size();
        const int nGrids = genome->B > 1 ? (nsnps + genome->B - 1) / genome->B : nsnps;
        alpha.setZero(C * C, nGrids);
        beta.setZero(C * C, nGrids);
        get_cluster_likelihood(ind, nsnps, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic], genome->F[ic]);
        ae.setZero(C * C, nGrids);
        get_cluster_pairs_probabity(ae, genome->R[ic], genome->PI[ic]);
        iEkc.setZero(C * K, nGrids);
        Ekg.setZero(K, nGrids);
        for(s = 0; s < nGrids; s++, m++)
        {
            for(norm = 0, cc = 0, c1 = 0; c1 < C; c1++)
            {
                for(c2 = c1; c2 < C; c2++)
                {
                    c12 = c1 * C + c2;
                    double xz = alpha(c12, s) * beta(c12, s) / ae(c12, s);
                    for(k1 = 0; k1 < K; k1++)
                    {
                        for(k2 = 0; k2 < K; k2++)
                        {
                            k12 = k1 * K + k2;
                            w(cc, k12) = xz * F(k1 * C + c1, m) * Q(k1, ind) * F(k2 * C + c2, m) * Q(k2, ind);
                            if(c1 != c2) w(cc, k12) *= 2;
                            norm += w(cc, k12);
                        }
                    }
                    ++cc;
                }
            }
            llike += log(norm);
            w /= norm;
            for(cc = 0, c1 = 0; c1 < C; c1++)
            {
                for(c2 = c1; c2 < C; c2++)
                {
                    for(k1 = 0; k1 < K; k1++)
                    {
                        for(k2 = 0; k2 < K; k2++)
                        {
                            k12 = k1 * K + k2;
                            Ekg(k1, s) += w(cc, k12);
                            Ekg(k2, s) += w(cc, k12);
                            iEkc(k1 * C + c1, s) += w(cc, k12);
                            iEkc(k2 * C + c2, s) += w(cc, k12);
                        }
                    }
                    ++cc;
                }
            }
        }
        iQ += Ekg.rowwise().sum();
        {
            std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(m - nGrids, nGrids) += iEkc;
            NormF.middleCols(m - nGrids, nGrids) += Ekg;
        }
    }
    // update Q, iQ.sum() should be 2M
    if(!nonewQ) Q.col(ind) = iQ / (2 * M);

    return llike;
}

void Admixture::initIteration()
{
    Ekc.setZero(C * K, M);
    NormF.setZero(K, M);
}

void Admixture::updateIteration()
{
    for(int k = 0; k < K; k++) F.middleRows(k * C, C) = Ekc.middleRows(k * C, C).rowwise() / NormF.row(k);
    protectPars();
}

void Admixture::protectPars()
{
    if(Q.isNaN().any()) cao.error("NaN in Q\n");
    Q = (Q < admixtureThreshold).select(admixtureThreshold, Q); // lower bound
    Q = (Q > 1 - admixtureThreshold).select(1 - admixtureThreshold, Q); // upper bound
    Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual

    if(F.isNaN().any()) cao.error("NaN in F\n");
    F = (F < clusterFreqThreshold).select(clusterFreqThreshold, F); // lower bound
    F = (F > 1 - clusterFreqThreshold).select(1 - clusterFreqThreshold, F); // upper bound
    normalizeF();
}

void Admixture::normalizeF()
{
    for(int k = 0; k < K; k++) // normalize F per snp per k
        F.middleRows(k * C, C).rowwise() /= F.middleRows(k * C, C).colwise().sum();
}

void Admixture::setStartPoint(std::string qfile)
{
    F = RandomUniform<MyArr2D, std::default_random_engine>(C * K, M, rng, clusterFreqThreshold,
                                                           1 - clusterFreqThreshold);
    normalizeF();
    Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, admixtureThreshold, 1 - admixtureThreshold);
    Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
    if(!qfile.empty()) load_csv(qfile, Q);
}

void Admixture::writeQ(std::string out)
{
    std::ofstream ofs(out);
    if(!ofs) cao.error(out, strerror(errno));
    Q = (Q * 1e6).round() / 1e6;
    ofs << std::fixed << Q.transpose() << "\n";
    ofs.close();
}

void Admixture::setFlags(bool debug_, bool nonewQ_)
{
    debug = debug_;
    nonewQ = nonewQ_;
}

int run_admix_main(Options & opts)
{
    cao.cao.open(opts.out.string() + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running admixture");
    int avail_threads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < avail_threads ? opts.nthreads : avail_threads;
    cao.print(tim.date(), avail_threads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool poolit(opts.nthreads);
    // Deserialize from file
    auto filesize = std::filesystem::file_size(opts.in_bin);
    std::error_code ec;
    std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    cao.done(tim.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
    cao.print(tim.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks, ", B =", genome->B, ", G =", genome->G);
    assert(opts.K < genome->C);

    cao.warn(tim.date(), "-> running admixture with seed =", opts.seed);
    Admixture admixer(genome->nsamples, genome->G, genome->C, opts.K, opts.seed);
    admixer.setFlags(opts.debug, opts.nonewQ);
    admixer.setStartPoint(opts.in_qfile);
    vector<future<double>> llike;
    if(!opts.noaccel)
    {
        MyArr2D F0, Q0, F1, Q1;
        const int istep{4};
        double alpha{std::numeric_limits<double>::lowest()}, qdiff, ldiff, stepMax{4}, alphaMax{1280};
        double prevlike{std::numeric_limits<double>::lowest()};
        for(int it = 0; it < opts.nadmix / 3; it++)
        {
            // first accel iteration
            admixer.initIteration();
            F0 = admixer.F;
            Q0 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            for(auto && ll : llike) ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            // second accel iteration
            F1 = admixer.F;
            Q1 = admixer.Q;
            qdiff = (Q1 - Q0).square().sum();
            tim.clock();
            admixer.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            double loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            ldiff = loglike - prevlike;
            prevlike = loglike;
            cao.print(tim.date(), "SqS3 iteration", it * 3 + 1, ", diff(Q) =", std::scientific, qdiff,
                      ", alpha=", alpha, ", likelihoods =", std::fixed, loglike, ", diff(likelihoods)=", ldiff,
                      ", elapsed", tim.reltime(), " sec");
            if(qdiff < opts.qtol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff(Q) =", std::scientific, qdiff, " <", opts.qtol);
                break;
            }
            if(ldiff < opts.ltol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff(loglikelihood) =", std::scientific, ldiff, " <",
                          opts.ltol);
                break;
            }
            // accel iteration with steplen
            alpha = ((F1 - F0).square().sum() + (Q1 - Q0).square().sum())
                    / ((admixer.F - 2 * F1 + F0).square().sum() + (admixer.Q - 2 * Q1 + Q0).square().sum());
            // alpha = ((Q1 - Q0).square().sum()) / ((admixer.Q - 2 * Q1 + Q0).square().sum());
            alpha = max(1.0, sqrt(alpha));
            if(alpha >= stepMax)
            {
                alpha = min(stepMax, alphaMax);
                stepMax = min(stepMax * istep, alphaMax);
            }
            admixer.F = F0 + 2 * alpha * (F1 - F0) + alpha * alpha * (admixer.F - 2 * F1 + F0);
            admixer.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (admixer.Q - 2 * Q1 + Q0);
            admixer.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            for(auto && ll : llike) ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
        }
    }
    else
    {
        MyArr2D Q0;
        double qdiff, ldiff, loglike;
        double prevlike{std::numeric_limits<double>::lowest()};
        for(int it = 0; it < opts.nadmix; it++)
        {
            tim.clock();
            admixer.initIteration();
            Q0 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));

            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            ldiff = loglike - prevlike;
            prevlike = loglike;
            qdiff = (admixer.Q - Q0).square().sum();
            cao.print(tim.date(), "normal iteration", it, ", diff(Q) =", std::scientific, qdiff,
                      ", likelihoods =", std::fixed, loglike, ", diff(likelihoods)=", ldiff, ", elapsed", tim.reltime(),
                      " sec");
            if(qdiff < opts.qtol) break;
            if(ldiff < opts.ltol) break;
        }
    }
    cao.done(tim.date(), "admixture done and outputting");
    admixer.writeQ(opts.out.string() + ".Q");
    cao.done(tim.date(), "-> good job. have a nice day, bye!");

    return 0;
}
