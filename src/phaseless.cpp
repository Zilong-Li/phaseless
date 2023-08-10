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

void Phaseless::initRecombination(const Int1D & pos, double Ne = 20000, int B = 1)
{
    dist = calc_position_distance(pos);
    nGen = 4 * Ne / C;
    er = calc_distRate(pos, C, 1.0);
    et = calc_distRate(pos, C, 0.05);
}

void Phaseless::initRecombination(const Int2D & pos, double Ne = 20000, int B = 1)
{
    int nchunks = pos.size();
    pos_chunk.resize(nchunks + 1);
    er = MyArr1D::Ones(M);
    et = MyArr1D::Ones(M);
    int i{0}, ss{0};
    for(i = 0; i < nchunks; i++)
    {
        pos_chunk[i] = ss;
        dist = calc_position_distance(pos[i]);
        er.segment(ss, pos[i].size()) = calc_distRate(pos[i], C, 1.0);
        et.segment(ss, pos[i].size()) = calc_distRate(pos[i], C, 0.05);
        ss += pos[i].size();
    }
    pos_chunk[nchunks] = ss; // add sentinel
}

void Phaseless::protectPars()
{
    // if we accelerate pars, protect them!
    // protect Q
    Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
    if(debug && P.isNaN().any()) cao.warn("NaN in P in Phaseless model. will fill it with AF");
    // protect P
    P = (P < alleleEmitThreshold).select(alleleEmitThreshold, P); // lower bound
    P = (P > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, P); // upper bound
    // protect F
    for(int k = 0; k < K; k++)
    {
        F[k].rowwise() /= F[k].colwise().sum(); // normalize F per site
        // could cluster jump be zero?
        if(F[k].isNaN().any() || (F[k] < clusterFreqThreshold).any())
        {
            if(debug && F[k].isNaN().any()) cao.warn("NaN in F in Phaseless model. reset it to the threshold. k =", k);
            F[k] = (F[k] < clusterFreqThreshold).select(0, F[k]); // reset to 0 first
            for(auto m = 0; m < F[k].cols(); m++)
            {
                if(auto c = (F[k].col(m) == 0).count() > 0)
                {
                    // for columns with an entry below 0
                    // each 0 entry becomes threshold
                    // then rest re-scaled so whole thing has sum 1
                    double xsum = 1 - c * clusterFreqThreshold;
                    double csum = F[k].col(m).sum();
                    F[k].col(m) = (F[k].col(m) > 0).select(F[k].col(m) * xsum / csum, clusterFreqThreshold);
                }
            }
        }
    }
}

void Phaseless::initIteration()
{
    PclusterK.setZero(C * K, M);
    EclusterA1.setZero(C, M);
    EclusterA2.setZero(C, M);
    Eancestry.setZero(K, N);
}

void Phaseless::updateIteration()
{
    // update Q
    Q = Eancestry;
    // update P
    P = (EclusterA2 / (EclusterA1 + EclusterA2)).transpose();
    // update F
    for(int k = 0; k < K; k++) F[k] = PclusterK.middleRows(k * C, C);
    protectPars();
}

void Phaseless::getForwardPrevSums(const MyArr2D & alpha,
                                   MyArr2D & prevsum_z,
                                   MyArr2D & prevsum_zz,
                                   MyArr2D & prevsum_zy,
                                   MyArr1D & prevsum_zzy)
{
    int z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    prevsum_z.setZero();
    prevsum_zz.setZero();
    prevsum_zy.setZero();
    prevsum_zzy.setZero();
    for(z1 = 0; z1 < C; z1++)
        for(z2 = 0; z2 < C; z2++)
            for(y1 = 0; y1 < K; y1++)
                for(y2 = 0; y2 < K; y2++)
                {
                    zz = z1 * C + z2;
                    yy = y1 * K + y2;
                    auto aa = alpha(zz, yy);
                    prevsum_z(z2, yy) += aa;
                    prevsum_zz(y1, y2) += aa;
                    prevsum_zy(z1, y1) += aa;
                    prevsum_zzy(y1) += aa;
                }
}

void Phaseless::getBackwardPrevSums(const MyArr2D & beta_mult_emit,
                                    MyArr2D & prevsum_z,
                                    MyArr2D & prevsum_zz,
                                    MyArr2D & prevsum_zy,
                                    MyArr1D & prevsum_zzy,
                                    double & prevsum_zzyy,
                                    int ind,
                                    int s)
{
    int z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    prevsum_z.setZero();
    prevsum_zz.setZero();
    prevsum_zy.setZero();
    prevsum_zzy.setZero();
    prevsum_zzyy = 0.0;
    for(z1 = 0; z1 < C; z1++)
        for(z2 = 0; z2 < C; z2++)
            for(y1 = 0; y1 < K; y1++)
                for(y2 = 0; y2 < K; y2++)
                {
                    zz = z1 * C + z2;
                    yy = y1 * K + y2;
                    prevsum_z(z1, yy) += beta_mult_emit(zz, yy) * F[y2](z2, s);
                    prevsum_zz(y1, y2) += beta_mult_emit(zz, yy) * F[y1](z1, s) * F[y2](z2, s);
                    prevsum_zy(z1, y1) += beta_mult_emit(zz, yy) * F[y2](z2, s) * Q(y2, ind);
                    prevsum_zzy(y1) += beta_mult_emit(zz, yy) * F[y1](z1, s) * F[y2](z2, s) * Q(y2, ind);
                    prevsum_zzyy += beta_mult_emit(zz, yy) * F[y1](z1, s) * F[y2](z2, s) * Q(y1, ind) * Q(y2, ind);
                }
}

void Phaseless::moveForward(int ind,
                            int s,
                            MyArr2D & curr,
                            const MyArr1D & emit,
                            const MyArr2D & prev,
                            const MyArr2D & prevsum_z,
                            const MyArr2D & prevsum_zz,
                            const MyArr2D & prevsum_zy,
                            const MyArr1D & prevsum_zzy)
{
    int z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    for(z1 = 0; z1 < C; z1++)
        for(z2 = 0; z2 < C; z2++)
            for(y1 = 0; y1 < K; y1++)
                for(y2 = 0; y2 < K; y2++)
                {
                    zz = z1 * C + z2;
                    yy = y1 * K + y2;
                    int yy2 = y2 * K + y1;
                    double ss = et(s) * et(s)
                                * (er(s) * er(s) * prev(zz, yy)
                                   + er(s) * (1 - er(s))
                                         * (F[y2](z2, s) * prevsum_z(z1, yy2) + F[y1](z1, s) * prevsum_z(z2, yy))
                                   + (1 - er(s)) * (1 - er(s)) * (F[y1](z1, s) * F[y2](z2, s) * prevsum_zz(y1, y2)));
                    double sj = et(s) * (1 - et(s))
                                * (Q(y2, ind) * F[y2](z2, s)
                                       * (er(s) * prevsum_zy(z1, y1) + (1 - er(s)) * prevsum_zzy(y1) * F[y1](z1, s))
                                   + Q(y1, ind) * F[y1](z1, s)
                                         * (er(s) * prevsum_zy(z2, y2) + (1 - er(s)) * prevsum_zzy(y2) * F[y2](z2, s)));
                    double jj = (1 - et(s)) * (1 - et(s)) * Q(y1, ind) * F[y1](z1, s) * Q(y2, ind) * F[y2](z2, s);
                    curr(zz, yy) = emit(zz) * (ss + sj + jj);
                }
}

void Phaseless::moveBackward(int ind,
                             int s,
                             MyFloat cs,
                             MyArr2D & curr,
                             const MyArr2D & beta_mult_emit,
                             const MyArr2D & prevsum_z,
                             const MyArr2D & prevsum_zz,
                             const MyArr2D & prevsum_zy,
                             const MyArr1D & prevsum_zzy,
                             const double prevsum_zzyy)
{
    int z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    for(z1 = 0; z1 < C; z1++)
        for(z2 = 0; z2 < C; z2++)
            for(y1 = 0; y1 < K; y1++)
                for(y2 = 0; y2 < K; y2++)
                {
                    zz = z1 * C + z2;
                    yy = y1 * K + y2;
                    int yy2 = y2 * K + y1;
                    double ss = et(s) * et(s)
                                * (er(s) * er(s) * beta_mult_emit(zz, yy)
                                   + er(s) * (1 - er(s)) * (prevsum_z(z1, yy) + prevsum_z(z2, yy2))
                                   + (1 - er(s)) * (1 - er(s)) * prevsum_zz(y1, y2));
                    double sj = (1 - et(s)) * et(s)
                                * (er(s) * (prevsum_zy(z1, y1) + prevsum_zy(z2, y2))
                                   + (1 - er(s)) * (prevsum_zzy(y1) + prevsum_zzy(y2)));
                    double jj = (1 - et(s)) * (1 - et(s)) * prevsum_zzyy;
                    curr(zz, yy) = (ss + sj + jj) / cs;
                }
}

void Phaseless::getPosterios(const int ind,
                             const int ic,
                             const MyArr2D & gli,
                             const MyArr2D & emit,
                             const MyArr1D & cs,
                             const std::vector<MyArr2D> & alpha,
                             const std::vector<MyArr2D> & beta)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    assert(S == alpha.size());
    int m{0}, s{0}, z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    MyArr2D ind_post_zg1(C, S), ind_post_zg2(C, S), ind_post_zz(CC, S);
    MyArr2D ind_post_y = MyArr2D::Zero(K, S);
    MyArr2D ind_post_zy = MyArr2D::Zero(C * K, S);
    MyArr2D prevsum_z(C, KK), prevsum_zz(K, K), prevsum_zy(C, K);
    MyArr1D prevsum_zzy(K);
    for(s = 0; s < S; s++)
    {
        m = s + pos_chunk[ic];
        MyArr2D gamma = alpha[s] * beta[s]; // post(z, y)
        if(debug && !(std::abs(1.0 - gamma.sum()) < 1e-6))
            cao.error(s, ", gamma sum is not approx 1.0 (tol = 1e-6)\n", gamma.sum());
        ind_post_zz.col(s) = gamma.rowwise().sum();
        if(debug && !(std::abs(1.0 - ind_post_zz.col(s).sum()) < 1e-6))
            cao.error(s, ", zpost is not approx 1.0 (tol = 1e-6)\n", ind_post_zz.col(s).sum());

        auto gamma_div_emit = ind_post_zz.col(s) / emit.col(s);
        for(z1 = 0; z1 < C; z1++)
        {
            ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - P(m, z1))
                                   * (gli(s, 0) * (1 - P.row(m)) + gli(s, 1) * P.row(m)).transpose())
                                      .sum();
            ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (P(m, z1))
                                   * (gli(s, 1) * (1 - P.row(m)) + gli(s, 2) * P.row(m)).transpose())
                                      .sum();
            if(s == 0)
            {
                for(y1 = 0; y1 < K; y1++)
                    ind_post_zy(y1 * C + z1, s) = gamma(Eigen::seqN(z1, C, C), Eigen::seqN(y1, K, K)).sum();
            }
        }
        if(debug && !(std::abs(1.0 - (ind_post_zg1.col(s) + ind_post_zg2.col(s)).sum()) < 1e-6))
            cao.error(s, ", clusterallele is not approx 1.0 (tol = 1e-6)\n", gamma.sum());
        if(s == 0) ind_post_y.col(s) = gamma.colwise().sum().reshaped(K, K).colwise().sum();
        if(s == 0) continue;
        getForwardPrevSums(alpha[s - 1], prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy);
        for(z1 = 0; z1 < C; z1++)
            for(z2 = 0; z2 < C; z2++)
                for(y1 = 0; y1 < K; y1++)
                    for(y2 = 0; y2 < K; y2++)
                    {
                        zz = z1 * C + z2;
                        yy = y1 * K + y2;
                        double ebc = emit(zz, s) * beta[s](zz, yy) / cs(s);
                        ind_post_y(y1, s) +=
                            ebc
                            * ((1 - et(m)) * (1 - et(m)) * F[y1](z1, m) * Q(y1, ind) * F[y2](z2, m) * Q(y2, ind)
                               + et(m) * (1 - et(m)) * Q(y1, ind) * F[y1](z1, m)
                                     * (er(m) * prevsum_zy(z2, y2) + (1 - er(m)) * F[y2](z2, m) * prevsum_zzy(y2)));
                        ind_post_zy(y1 * C + z1, s) +=
                            ebc
                            * ((1 - et(m)) * Q(y1, ind) * F[y1](z1, m)
                                   * ((1 - et(m)) * F[y2](z2, m) * Q(y2, ind)
                                      + et(m) * (1 - er(m)) * F[y2](z2, m) * prevsum_zzy(y2)
                                      + et(m) * er(m) * prevsum_zy(z2, y2))
                               + et(m) * (1 - er(m)) * F[y1](z1, m)
                                     * ((1 - et(m)) * F[y2](z2, m) * Q(y2, ind) * prevsum_zzy(y1)
                                        + et(m) * (1 - er(m)) * F[y2](z2, m) * prevsum_zz(y1, y2)
                                        + et(m) * er(m) * prevsum_z(z2, yy)));
                    }
    }
    { // sum over all samples for updates
        std::lock_guard<std::mutex> lock(mutex_it);
        EclusterA1.middleCols(pos_chunk[ic], S) += ind_post_zg1;
        EclusterA2.middleCols(pos_chunk[ic], S) += ind_post_zg2;
        PclusterK.middleCols(pos_chunk[ic], S) += ind_post_zy; // cluster ancestry allele jump
        Eancestry.col(ind) += ind_post_y.rowwise().sum(); // need to take care
    }
}

double Phaseless::runForwardBackwards(const int ind, const int ic, const MyFloat1D & GL, bool finalIter)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * S * 3, S, 3);
    MyArr2D emit = get_emission_by_gl(gli, P.middleRows(pos_chunk[ic], S)).transpose(); // CC x S
    std::vector<MyArr2D> alpha(S, MyArr2D(CC, KK));
    std::vector<MyArr2D> beta(S, MyArr2D(CC, KK));
    MyArr2D prevsum_z(C, KK), prevsum_zz(K, K), prevsum_zy(C, K);
    MyArr1D cs(S), prevsum_zzy(K);
    // init
    int m{0}, s{0}, z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    m = s + pos_chunk[ic];
    for(z1 = 0; z1 < C; z1++)
        for(z2 = 0; z2 < C; z2++)
            for(y1 = 0; y1 < K; y1++)
                for(y2 = 0; y2 < K; y2++)
                {
                    zz = z1 * C + z2;
                    yy = y1 * K + y2;
                    alpha[s](zz, yy) = emit(zz, s) * Q(y1, ind) * F[y1](z1, m) * Q(y2, ind) * F[y2](z2, m);
                }
    cs(s) = alpha[s].sum();
    alpha[s] /= cs(s); // norm alpha
    // forward recursion
    for(s = 1; s < S; s++)
    {
        getForwardPrevSums(alpha[s - 1], prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy);
        m = s + pos_chunk[ic];
        moveForward(ind, m, alpha[s], emit.col(s), alpha[s - 1], prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy);
        cs(s) = alpha[s].sum();
        alpha[s] /= cs(s); // norm alpha
    }
    // backwards recursion
    s = S - 1;
    beta[s].setConstant(1.0);
    MyArr2D beta_mult_emit(CC, KK);
    double prevsum_zzyy{0};
    for(s = S - 2; s >= 0; s--)
    {
        for(int yy = 0; yy < KK; yy++) beta_mult_emit.col(yy) = beta[s + 1].col(yy) * emit.col(s + 1);
        m = s + 1 + pos_chunk[ic];
        getBackwardPrevSums(beta_mult_emit, prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy, prevsum_zzyy, ind, m);
        moveBackward(ind, m, cs(s + 1), beta[s], beta_mult_emit, prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy,
                     prevsum_zzyy);
    }
    // get posterios
    getPosterios(ind, ic, gli, emit, cs, alpha, beta);
    return cs.log().sum();
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
    cao.cao.open(opts.out.string() + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running fastphase");
    int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < allthreads ? opts.nthreads : allthreads;
    cao.print(tim.date(), allthreads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool pool(opts.nthreads);

    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    init_bigass(genome, opts);
    Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
    vector<future<double>> res;
    std::ofstream oanc(opts.out.string() + ".Q");
    Phaseless faith(opts.K, opts.C, genome->nsamples, genome->nsnps, opts.seed);
    faith.initRecombination(genome->pos);
    faith.debug = opts.debug;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    if(opts.noaccel)
    {
        for(int it = 0; it <= opts.nimpute; it++)
        {
            tim.clock();
            faith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == opts.nimpute)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), true));
                else
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            }
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
        MyArr2D Q0, Q1;
        MyArr2D F0, F1;
        const int istep{4};
        double alpha{0}, stepMax{4}, alphaMax{1280};
        for(int it = 0; it < opts.nimpute / 3; it++)
        {
            // first normal iter
            faith.initIteration();
            Q0 = faith.Q;
            F0 = cat_stdvec_of_eigen(faith.F);
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == opts.nimpute)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), true));
                else
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            }
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            // second normal iter
            tim.clock();
            faith.initIteration();
            Q1 = faith.Q;
            F1 = cat_stdvec_of_eigen(faith.F);
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == opts.nimpute)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), true));
                else
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            }
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
            diff = it ? loglike - prevlike : NAN;
            prevlike = loglike;
            cao.print(tim.date(), "SqS3 iteration", it * 3 + 1, ", alpha=", alpha, ", likelihoods =", std::fixed,
                      loglike, ", diff =", diff, ", time", tim.reltime(), " sec");
            if(diff < opts.ltol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff =", std::scientific, diff, " <", opts.ltol);
                break;
            }
            // third accel iter
            faith.initIteration();
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
            faith.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (faith.Q - 2 * Q1 + Q0);
            for(int k = 0; k < opts.K; k++)
                faith.F[k] =
                    F0.middleRows(k * opts.C, opts.C)
                    + 2 * alpha * (F1.middleRows(k * opts.C, opts.C) - F0.middleRows(k * opts.C, opts.C))
                    + alpha * alpha
                          * (faith.F[k] - 2 * F1.middleRows(k * opts.C, opts.C) + F0.middleRows(k * opts.C, opts.C));
            faith.protectPars();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == opts.nimpute)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), true));
                else
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
            }
            loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            faith.updateIteration();
        }
    }
    faith.Q = (faith.Q * 1e6).round() / 1e6;
    oanc << std::fixed << faith.Q.transpose().format(fmt) << "\n";

    return 0;
}
