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
    MyArr1D iQ = MyArr1D::Zero(K);
    MyArr1D Hz(C);
    double norm = 0, llike = 0, tmp = 0;
    int c1, k1, s, c2, c12, ic, g, ss, ng;
    for(ic = 0, g = 0, ss = 0, ng = 0; ic < genome->nchunks; ic++)
    {
        const int S = genome->pos[ic].size();
        const int nGrids = grids[ic];
        Eigen::Map<const MyArr2D> gli(genome->gls[ic].data() + ind * S * 3, S, 3);
        Eigen::Map<const MyArr2D> P(genome->P[ic].data(), S, C);
        Eigen::Map<const MyArr2D> PI(genome->PI[ic].data(), C, nGrids);
        Eigen::Map<const MyArr2D> R(genome->R[ic].data(), 3, nGrids);
        Eigen::Map<const MyArr2D> AE(genome->AE[ic].data(), C * C, nGrids);
        const auto cl = get_cluster_likelihoods(gli, P, R, PI, AE, collapse.segment(ss, S));
        ss += S;
        kapa.setZero(C * K, nGrids); // C x K x M layout
        Ekg.setZero(K, nGrids);
        for(s = 0; s < nGrids; s++, g++)
        {
            if(magicTol > 0 && pi.col(g).maxCoeff() < magicTol)
            {
                kapa.col(s).fill(1.0);
            }
            else
            {
                ng++;
                for(c1 = 0; c1 < C; c1++) Hz(c1) = (Q.col(ind) * F(Eigen::seqN(c1, K, C), g)).sum();
                for(norm = 0, c1 = 0; c1 < C; c1++)
                {
                    for(tmp = 0, c2 = 0; c2 < C; c2++)
                    {
                        c12 = c1 * C + c2;
                        double xz = cl(c12, s);
                        double zy = Hz(c1) * Hz(c2);
                        tmp += xz * zy;
                    }
                    norm += tmp;
                    kapa(Eigen::seqN(c1, K, C), s) = (Q.col(ind) * F(Eigen::seqN(c1, K, C), g)) * tmp / Hz(c1);
                }
                llike += log(norm);
            }
            kapa.col(s) /= kapa.col(s).sum();
            for(k1 = 0; k1 < K; k1++) Ekg(k1, s) = 2 * kapa.middleRows(k1 * C, C).col(s).sum();
        }
        iQ += Ekg.rowwise().sum();
        { // for update F
            std::scoped_lock<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(g - nGrids, nGrids) += 2 * kapa;
            NormF.middleCols(g - nGrids, nGrids) += Ekg;
        }
    }
    // update Q, iQ.sum() should be 2M
    if(!nonewQ) Q.col(ind) = iQ / (2 * ng);

    return llike;
}

// the complexity of this version is O(CCKK)
double Admixture::runNativeWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D w((C * C + C) / 2, K * K);
    MyArr2D Ekg, iEkc;
    double norm = 0, llike = 0;
    int c1, c2, c12, cc;
    int k1, k2, k12, s;
    MyArr1D iQ = MyArr1D::Zero(K);
    for(int ic = 0, g = 0, ss = 0; ic < genome->nchunks; ic++)
    {
        const int S = genome->pos[ic].size();
        const int nGrids = grids[ic];
        Eigen::Map<const MyArr2D> gli(genome->gls[ic].data() + ind * S * 3, S, 3);
        Eigen::Map<const MyArr2D> P(genome->P[ic].data(), S, C);
        Eigen::Map<const MyArr2D> PI(genome->PI[ic].data(), C, nGrids);
        Eigen::Map<const MyArr2D> R(genome->R[ic].data(), 3, nGrids);
        Eigen::Map<const MyArr2D> AE(genome->AE[ic].data(), C * C, nGrids);
        const auto cl = get_cluster_likelihoods(gli, P, R, PI, AE, collapse.segment(ss, S));
        ss += S;
        iEkc.setZero(C * K, nGrids);
        Ekg.setZero(K, nGrids);
        for(s = 0; s < nGrids; s++, g++)
        {
            for(norm = 0, cc = 0, c1 = 0; c1 < C; c1++)
            {
                for(c2 = c1; c2 < C; c2++)
                {
                    c12 = c1 * C + c2;
                    double xz = cl(c12, s);
                    if(magicTol > 0 && pi.col(g).maxCoeff() < magicTol) xz = 0.0;
                    for(k1 = 0; k1 < K; k1++)
                    {
                        for(k2 = 0; k2 < K; k2++)
                        {
                            k12 = k1 * K + k2;
                            w(cc, k12) = xz * F(k1 * C + c1, g) * Q(k1, ind) * F(k2 * C + c2, g) * Q(k2, ind);
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
            std::scoped_lock<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(g - nGrids, nGrids) += iEkc;
            NormF.middleCols(g - nGrids, nGrids) += Ekg;
        }
    }
    // update Q, iQ.sum() should be 2M
    if(!nonewQ) Q.col(ind) = iQ / (2 * G);

    return llike;
}

void Admixture::initIteration()
{
    Ekc.setZero(C * K, G);
    NormF.setZero(K, G);
}

void Admixture::updateIteration()
{
    for(int k = 0; k < K; k++) F.middleRows(k * C, C) = Ekc.middleRows(k * C, C).rowwise() / NormF.row(k);
    protectPars();
}

void Admixture::protectPars()
{
    if(!nonewQ)
    {
        if(Q.isNaN().any()) cao.error("NaN in Q\n", Q);
        Q = (Q < admixtureThreshold).select(admixtureThreshold, Q); // lower bound
        Q = (Q > 1 - admixtureThreshold).select(1 - admixtureThreshold, Q); // upper bound
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
    }

    if(F.isNaN().any()) cao.error("NaN in F\n");
    if(cF)
    {
        constrainF();
    }
    else
    {
        F = (F < clusterFreqThreshold).select(clusterFreqThreshold, F); // lower bound
        F = (F > 1 - clusterFreqThreshold).select(1 - clusterFreqThreshold, F); // upper bound
        for(int k = 0; k < K; k++) F.middleRows(k * C, C).rowwise() /= F.middleRows(k * C, C).colwise().sum();
    }
}

// need to fix the cluster ordering
void Admixture::constrainF()
{
    for(int k = 0; k < K; k++)
    {
        for(int c = 0; c < C; c++)
            for(int g = 0; g < G; g++)
                if(F(k * C + c, g) < P(c, g)) F(k * C + c, g) = P(c, g);
        F.middleRows(k * C, C).rowwise() /= F.middleRows(k * C, C).colwise().sum();
    }
}

void Admixture::setStartPoint(const std::unique_ptr<BigAss> & genome, std::string qfile, std::string pifile)
{
    if(!pifile.empty())
    {
        pi.setZero(C, G);
        load_csv(pi, pifile);
    }
    P = MyArr2D(C, G);
    collapse = Bool1D::Constant(genome->nsnps, false);
    int ic{0}, sg{0}, ss{0};
    for(auto c : genome->collapse) collapse(ic++) = (c == 1);
    grids.resize(genome->nchunks);
    for(ic = 0, sg = 0, ss = 0; ic < genome->nchunks; ic++)
    {
        const int S = genome->pos[ic].size();
        const auto se = find_grid_start_end(collapse.segment(ss, S));
        const int iG = se.size();
        grids[ic] = iG;
        Eigen::Map<const MyArr2D> AE(genome->AE[ic].data(), C * C, iG);
        for(int g = 0; g < iG; g++) P.col(sg + g) = AE.col(g).reshaped(C, C).colwise().sum();
        sg += iG;
        ss += S;
    }
    assert(sg == G);
    if(!qfile.empty()) load_csv(Q, qfile);
}

void Admixture::setFlags(double cftol, double Ftol, double Qtol, bool debug_, bool nonewQ_, bool cF_)
{
    magicTol = cftol;
    clusterFreqThreshold = Ftol;
    admixtureThreshold = Qtol;
    debug = debug_;
    nonewQ = nonewQ_;
    cF = cF_;
}

int run_admix_main(Options & opts)
{
    cao.cao.open(opts.out + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running admixture");
    int avail_threads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < avail_threads ? opts.nthreads : avail_threads;
    cao.print(tim.date(), avail_threads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool poolit(opts.nthreads);
    // Deserialize from file
    auto filesize = std::filesystem::file_size(std::filesystem::path(opts.in_bin));
    std::error_code ec;
    std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    cao.done(tim.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
    cao.print(tim.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks, ", B =", genome->B, ", G =", genome->G);

    Admixture admixer(genome->nsamples, genome->G, genome->C, opts.K, opts.seed);
    cao.warn(tim.date(), "-> running admixture with seed =", opts.seed);
    admixer.setFlags(opts.tol_pi, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.cF);
    admixer.setStartPoint(genome, opts.in_qfile, opts.pi_file);
    vector<future<double>> llike;
    if(!opts.noaccel)
    {
        MyArr2D F0, Q0, F1, Q1, F2, Q2, Ft, Qt;
        const double istep{4};
        double alpha{std::numeric_limits<double>::lowest()}, qdiff, ldiff, stepMax{4}, alphaMax{1280};
        double prevlike{std::numeric_limits<double>::lowest()}, logcheck{0}, loglike{0};
        for(int it = 0; SIG_COND && (it < opts.nadmix / 3); it++)
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
            tim.clock();
            admixer.initIteration();
            F1 = admixer.F;
            Q1 = admixer.Q;
            qdiff = (Q1 - Q0).square().sum();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            ldiff = it ? loglike - prevlike : NAN;
            prevlike = loglike;
            cao.print(tim.date(), "SqS3 iteration", it * 3 + 1, ", diff(Q) =", std::scientific, qdiff,
                      ", alpha=", alpha, ", likelihoods =", std::fixed, loglike, ", diff(likelihoods)=", ldiff,
                      ", elapsed", tim.reltime(), " sec");
            if(ldiff < opts.ltol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff(loglikelihood) =", std::scientific, ldiff, " <",
                          opts.ltol);
                break;
            }
            if(!opts.force)
            {
                // save for later comparison
                Ft = admixer.F;
                Qt = admixer.Q;
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
            // third accel iter
            // update Q and F using the second em iter
            admixer.F = F0 + 2 * alpha * (F1 - F0) + alpha * alpha * (admixer.F - 2 * F1 + F0);
            admixer.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (admixer.Q - 2 * Q1 + Q0);
            admixer.protectPars();
            admixer.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            if(!opts.force)
            { // save current pars
                F2 = admixer.F;
                Q2 = admixer.Q;
                // check if normal third iter is better
                admixer.Q = Qt;
                admixer.F = Ft;
                admixer.initIteration();
                for(int i = 0; i < genome->nsamples; i++)
                    llike.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
                logcheck = 0;
                for(auto && ll : llike) logcheck += ll.get();
                llike.clear(); // clear future and renew
                admixer.updateIteration();
                if(logcheck - loglike > 0.1)
                {
                    stepMax = istep;
                    cao.warn(tim.date(),
                             "reset stepMax to 4, normal EM yields better likelihoods than the accelerated EM.",
                             logcheck, " -", loglike, " > 0.1");
                }
                else
                {
                    admixer.Q = Q2;
                    admixer.F = F2;
                }
            }
        }
    }
    else
    {
        MyArr2D Q0;
        double qdiff, ldiff, loglike;
        double prevlike{std::numeric_limits<double>::lowest()};
        for(int it = 0; SIG_COND && it < opts.nadmix; it++)
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
            ldiff = it ? loglike - prevlike : NAN;
            prevlike = loglike;
            qdiff = (admixer.Q - Q0).square().sum();
            cao.print(tim.date(), "normal iteration", it, ", diff(Q) =", std::scientific, qdiff,
                      ", likelihoods =", std::fixed, loglike, ", diff(likelihoods)=", ldiff, ", elapsed", tim.reltime(),
                      " sec");
            if(ldiff < opts.ltol) break;
        }
    }
    cao.done(tim.date(), "admixture done and outputting");
    std::ofstream oq(opts.out + ".Q");
    oq << std::fixed << admixer.Q.transpose().format(fmt10) << "\n";
    if(opts.oF)
    {
        std::ofstream of(opts.out + ".F");
        of << admixer.F.transpose().format(fmt6) << "\n";
    }
    cao.done(tim.date(), "-> good job. have a nice day, bye!");

    return 0;
}
