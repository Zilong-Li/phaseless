/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/phaseless.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "phaseless.hpp"

#include "common.hpp"
#include "io.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>

using namespace std;

void Phaseless::initRecombination(const Int1D & pos, std::string rfile, double Ne, int B)
{
    nGen = 4 * Ne / C;
    dist = calc_position_distance(pos);
    if(rfile.empty())
        R = calc_transRate_diploid(dist, nGen);
    else
        load_csv(rfile, R);
    er = R.row(0).sqrt();
}

void Phaseless::initRecombination(const Int2D & pos, std::string rfile, double Ne, int B)
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
        dist.insert(dist.end(),tmp.begin(), tmp.end());
        R.middleCols(ss, pos[i].size()) = calc_transRate_diploid(dist, nGen);
        ss += pos[i].size();
    }
    pos_chunk[nchunks] = ss; // add sentinel
    if(!rfile.empty()) load_csv(rfile, R);
    er = R.row(0).sqrt();
}

void Phaseless::setFlags(double tol_p, double tol_f, double tol_q, bool debug_, bool nQ, bool nP, bool nF, bool nR)
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

void Phaseless::setStartPoint(std::string qfile, std::string pfile)
{
    if(!qfile.empty()) load_csv(qfile, Q);
    if(!pfile.empty()) load_csv2(pfile, P);
}

void Phaseless::setStartPoint(const std::unique_ptr<Pars> & par)
{
    er = Eigen::Map<MyArr1D>(par->er.data(), M);
    R = er2R(er);
    Q = Eigen::Map<MyArr2D>(par->Q.data(), K, N);
    P = Eigen::Map<MyArr2D>(par->P.data(), M, C);
    for(int i = 0; i < K; i++) F[i] = Eigen::Map<MyArr2D>(par->F[i].data(), C, M);
}

void Phaseless::protectPars()
{
    // if we accelerate pars, protect them!
    if(!NQ)
    { // protect Q
        if(debug && Q.isNaN().any()) cao.warn("NaN in Q in Phaseless model. reset it to the threshold");
        Q = (Q < admixtureThreshold).select(admixtureThreshold, Q); // lower bound
        Q = (Q > 1 - admixtureThreshold).select(1 - admixtureThreshold, Q); // upper bound
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
        if(debug && !(1.0 - Q.colwise().sum() == 0).any()) cao.warn("Q colsum is not 1.0");
    }
    if(!NP)
    { // protect P
        if(P.isNaN().any()) cao.warn("NaN in P in Phaseless model. will fill it with AF");
        P = (P < alleleEmitThreshold).select(alleleEmitThreshold, P); // lower bound
        P = (P > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, P); // upper bound
    }
    if(!NF)
    { // protect F
        for(int k = 0; k < K; k++)
        {
            // could cluster jump be zero?
            if(F[k].isNaN().any()) cao.warn("NaN in F in Phaseless model. reset it to the threshold. k =", k);
            F[k] = (F[k] < clusterFreqThreshold).select(clusterFreqThreshold, F[k]);
            F[k] = (F[k] > 1 - clusterFreqThreshold).select(1 - clusterFreqThreshold, F[k]);
            // re-normalize F per site. hope should work well. otherwise do the complicated.
            F[k].rowwise() /= F[k].colwise().sum();
        }
    }
    if(!NR)
    {
        for(int i = 0; i < er.size(); i++)
        {
            double miner = std::exp(-nGen * maxRate * dist[i] / 100 / 1e6);
            double maxer = std::exp(-nGen * minRate * dist[i] / 100 / 1e6);
            er(i) = er(i) < miner ? miner : er(i);
            er(i) = er(i) > maxer ? maxer : er(i);
        }
        R = er2R(er);
    }
}

void Phaseless::initIteration()
{
    EclusterK.setZero(C * K, M);
    EclusterA1.setZero(C, M);
    EclusterA2.setZero(C, M);
    Eancestry.setZero(K, N);
}

void Phaseless::updateIteration()
{
    // update P
    if(!NP) P = (EclusterA2 / (EclusterA1 + EclusterA2)).transpose();
    if(!NQ)
    { // update Q
        Q = Eancestry;
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
    }
    if(!NF)
    { // update F
        for(int k = 0; k < K; k++)
        {
            F[k] = EclusterK.middleRows(k * C, C); // C x M
            F[k].rowwise() /= F[k].colwise().sum(); // normalize F per site per K
        }
    }
    if(!NR) er = 1.0 - EclusterK.colwise().sum() / N;
    protectPars();
}

void Phaseless::callGenoLoopC(int ind, int s, int z1, const MyArr2D & gli, const MyArr1D & gamma_div_emit)
{
    MyArr1D tmp_zg(4);
    for(int z2 = 0; z2 < C; z2++)
    {
        int z12 = z1 * C + z2;
        tmp_zg(0) = gli(s, 0) * (1 - P(s, z1)) * (1 - P(s, z2));
        tmp_zg(1) = gli(s, 1) * (1 - P(s, z1)) * P(s, z2);
        tmp_zg(2) = gli(s, 1) * P(s, z1) * (1 - P(s, z2));
        tmp_zg(3) = gli(s, 2) * P(s, z1) * P(s, z2);
        GP(3 * s + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
        GP(3 * s + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
        GP(3 * s + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
    }
}

void Phaseless::getPosterios(const int ind,
                             const int ic,
                             const MyArr2D & gli,
                             const MyArr2D & emit,
                             const MyArr2D & H,
                             const MyArr1D & cs,
                             const MyArr2D & alpha,
                             const MyArr2D & beta,
                             bool finalIter)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    int m{0}, s{0}, z1{0}, z2{0}, y1{0}, zz{0};
    MyArr2D ind_post_zg1(C, S), ind_post_zg2(C, S);
    MyArr2D ind_post_zy(C * K, S);
    MyArr1D gamma_div_emit(CC);
    MyArr2D ind_post_y = MyArr2D::Zero(K, S);
    for(s = 0; s < S; s++)
    {
        m = s + pos_chunk[ic];
        gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // what if emit is 0
        for(z1 = 0; z1 < C; z1++)
        {
            if(finalIter) callGenoLoopC(ind, m, z1, gli, gamma_div_emit);
            ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - P(m, z1))
                                   * (gli(s, 0) * (1 - P.row(m)) + gli(s, 1) * P.row(m)).transpose())
                                      .sum();
            ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (P(m, z1))
                                   * (gli(s, 1) * (1 - P.row(m)) + gli(s, 2) * P.row(m)).transpose())
                                      .sum();
            if(s == 0)
            {
                auto tmp = (alpha.col(0) * beta.col(0)).segment(z1 * C, C).sum();
                for(y1 = 0; y1 < K; y1++)
                {
                    ind_post_zy(y1 * C + z1, 0) = tmp * Q(y1, ind);
                    ind_post_y(y1, 0) += ind_post_zy(y1 * C + z1, 0);
                }
            }
        }
        if(s == 0) continue;
        MyArr1D alphaprev(C); // previous alpha colsums
        for(z1 = 0; z1 < C; z1++) alphaprev(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum();
        for(z1 = 0; z1 < C; z1++)
        {
            double tmp{0};
            for(z2 = 0; z2 < C; z2++)
            {
                zz = z1 * C + z2;
                double eb = emit(zz, s) * beta(zz, s);
                tmp += eb * (R(1, m) * alphaprev(z2) + R(2, m) * H(z2, s));
                for(y1 = 0; y1 < K; y1++)
                {
                    ind_post_y(y1, s) += eb * cs(s) * Q(y1, ind)
                                         * (R(0, m) * alpha(zz, s - 1)
                                            + R(1, m) * (alphaprev(z2) * F[y1](z1, m) + alphaprev(z1) * H(z2, s))
                                            + R(2, m) * F[y1](z1, m) * H(z2, s));
                }
            }
            for(y1 = 0; y1 < K; y1++) ind_post_zy(y1 * C + z1, s) = tmp * Q(y1, ind) * F[y1](z1, m) * cs(s);
        }
    }
    Eancestry.col(ind) += ind_post_y.rowwise().sum(); // need to take care
    { // sum over all samples for updates
        std::scoped_lock<std::mutex> lock(mutex_it);
        EclusterA1.middleCols(pos_chunk[ic], S) += ind_post_zg1;
        EclusterA2.middleCols(pos_chunk[ic], S) += ind_post_zg2;
        EclusterK.middleCols(pos_chunk[ic], S) += ind_post_zy;
    }
}

double Phaseless::runForwardBackwards(const int ind, const int ic, const MyFloat1D & GL, bool finalIter)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * S * 3, S, 3);
    MyArr2D emit = get_emission_by_gl(gli, P.middleRows(pos_chunk[ic], S)).transpose(); // CC x S
    MyArr2D alpha(CC, S), beta(CC, S);
    // first get H ie old PI in fastphase
    MyArr2D H = MyArr2D::Zero(C, S);
    int z1, y1, s;
    for(s = 0; s < S; s++)
        for(z1 = 0; z1 < C; z1++)
            for(y1 = 0; y1 < K; y1++) H(z1, s) += Q(y1, ind) * F[y1](z1, s + pos_chunk[ic]);
    // cs is 1 / colsum(alpha)
    auto cs = forward_backwards_diploid(alpha, beta, emit, R.middleCols(pos_chunk[ic], S), H);
    // get posterios
    getPosterios(ind, ic, gli, emit, H, cs, alpha, beta, finalIter);
    return (1 / cs).log().sum();
}

double Phaseless::runBigass(int ind, const MyFloat2D & GL, bool finalIter)
{
    if(pos_chunk.size() == 0) cao.error("please run initRecombination first");
    int nchunks = GL.size();
    double loglike{0};
    for(int ic = 0; ic < nchunks; ic++)
    {
        loglike += runForwardBackwards(ind, ic, GL[ic], finalIter);
    }
    return loglike;
}

int run_phaseless_main(Options & opts)
{
    cao.cao.open(opts.out + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running phaseless");
    int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < allthreads ? opts.nthreads : allthreads;
    cao.print(tim.date(), allthreads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool pool(opts.nthreads);

    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    init_bigass(genome, opts);
    Eigen::IOFormat fmt(6, Eigen::DontAlignCols, " ", "\n");
    vector<future<double>> res;
    std::ofstream oanc(opts.out + ".Q");
    std::ofstream op(opts.out + ".P");
    Phaseless faith(opts.K, opts.C, genome->nsamples, genome->nsnps, opts.seed);
    faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.nP, opts.nF, opts.nR);
    faith.setStartPoint(opts.in_qfile, opts.in_pfile);
    faith.initRecombination(genome->pos, opts.in_rfile);
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    if(opts.noaccel)
    {
        for(int it = 0; SIG_COND && it <= opts.nimpute; it++)
        {
            tim.clock();
            faith.initIteration();
            for(int i = 0; i < faith.N; i++)
                res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            diff = it ? loglike - prevlike : NAN;
            prevlike = loglike;
            cao.print(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff =", diff,
                      ", time", tim.reltime(), " sec");
            if(diff < opts.ltol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff =", std::scientific, diff, " <", opts.ltol);
                break;
            }
        }
    }
    else
    {
        MyArr2D Q0, Q1, Q2, Qt;
        MyArr2D F0, F1, F2, Ft;
        const int istep{4};
        double alpha{0}, stepMax{4}, alphaMax{1280}, logcheck{0};
        for(int it = 0; SIG_COND && (it < opts.nimpute / 4); it++)
        {
            // first normal iter
            faith.initIteration();
            Q0 = faith.Q;
            F0 = cat_stdvec_of_eigen(faith.F);
            for(int i = 0; i < faith.N; i++)
                res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            // second normal iter
            tim.clock();
            faith.initIteration();
            Q1 = faith.Q;
            F1 = cat_stdvec_of_eigen(faith.F);
            for(int i = 0; i < faith.N; i++)
                res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            diff = it ? loglike - prevlike : NAN;
            prevlike = loglike;
            cao.print(tim.date(), "SqS3 iteration", it * 4 + 1, ", alpha=", alpha, ", likelihoods =", std::fixed,
                      loglike, ", diff =", diff, ", time", tim.reltime(), " sec");
            if(diff < opts.ltol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff =", std::scientific, diff, " <", opts.ltol);
                break;
            }
            // save for later comparison
            Qt = faith.Q;
            Ft = cat_stdvec_of_eigen(faith.F);
            // calculate alpha based on first two pars
            // alpha = ((Q1 - Q0).square().sum()) / ((faith.Q - 2 * Q1 + Q0).square().sum());
            alpha = ((F1 - F0).square().sum() + (Q1 - Q0).square().sum())
                    / ((cat_stdvec_of_eigen(faith.F) - 2 * F1 + F0).square().sum()
                       + (faith.Q - 2 * Q1 + Q0).square().sum());
            alpha = max(1.0, sqrt(alpha));
            if(alpha >= stepMax)
            {
                alpha = min(stepMax, alphaMax);
                stepMax = min(stepMax * istep, alphaMax);
            }
            // third accel iter
            // update Q and F using the second em iter
            faith.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (faith.Q - 2 * Q1 + Q0);
            for(int k = 0; k < faith.K; k++)
                faith.F[k] = F0.middleRows(k * faith.C, faith.C)
                             + 2 * alpha * (F1.middleRows(k * faith.C, faith.C) - F0.middleRows(k * faith.C, faith.C))
                             + alpha * alpha
                                   * (faith.F[k] - 2 * F1.middleRows(k * faith.C, faith.C)
                                      + F0.middleRows(k * faith.C, faith.C));
            faith.protectPars();
            faith.initIteration();
            for(int i = 0; i < faith.N; i++)
                res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            // save current pars
            Q2 = faith.Q;
            F2 = cat_stdvec_of_eigen(faith.F);
            // check if normal third iter is better
            faith.Q = Qt;
            for(int k = 0; k < faith.K; k++) faith.F[k] = Ft.middleRows(k * faith.C, faith.C);
            faith.initIteration();
            for(int i = 0; i < faith.N; i++)
                res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            logcheck = 0;
            for(auto && ll : res) logcheck += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            if(logcheck - loglike > 0.1)
            {
                stepMax = istep;
                cao.warn(tim.date(), "reset stepMax to 4, normal EM yields better likelihoods than the accelerated EM.",
                         logcheck, " -", loglike, "> 0.1");
            }
            else
            {
                faith.Q = Q2;
                for(int k = 0; k < faith.K; k++) faith.F[k] = F2.middleRows(k * faith.C, faith.C);
            }
        }
    }
    cao.done(tim.date(), "joint model done. run one more iteration to output vcf.");
    faith.initIteration();
    faith.GP.setZero(faith.M * 3, faith.N);
    for(int i = 0; i < faith.N; i++)
        res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), true));
    loglike = 0;
    for(auto && ll : res) loglike += ll.get();
    res.clear(); // clear future and renew
    faith.updateIteration();
    auto bw = make_bcfwriter(opts.out + ".vcf.gz", genome->chrs, genome->sampleids);
    for(int ic = 0; ic < genome->nchunks; ic++)
        write_bigass_to_bcf(bw, faith.GP.data(), genome->chrs[ic], genome->pos[ic]);
    faith.Q = (faith.Q * 1e6).round() / 1e6;
    oanc << std::fixed << faith.Q.transpose().format(fmt) << "\n";
    op << std::fixed << faith.P.format(fmt) << "\n";
    std::unique_ptr<Pars> par = std::make_unique<Pars>();
    par->init(faith.K, faith.C, faith.M, faith.N, faith.er, faith.P, faith.Q, faith.F);
    par->pos = genome->pos;
    par->gls = genome->gls;
    std::ofstream opar(opts.out + ".pars.bin", std::ios::out | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    auto bytes_written = alpaca::serialize<OPTIONS, Pars>(*par, opar);
    opar.close();
    assert(std::filesystem::file_size(opts.out + ".pars.bin") == bytes_written);
    cao.done(tim.date(), "joint model done and outputting.", bytes_written, " bytes written to file");

    return 0;
}
