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

void Phaseless::initIteration()
{
    PclusterK.setZero(C * K, M);
    EclusterA1.setZero(C, M);
    EclusterA2.setZero(C, M);
}

void Phaseless::updateIteration()
{
    P = (EclusterA2 / (EclusterA1 + EclusterA2)).transpose();
    if(P.isNaN().any())
    {
        if(debug) cao.warn("NaN in P in Phaseless model. will fill it with AF");
    }
    // map F to domain but no normalization
    P = (P < alleleEmitThreshold).select(alleleEmitThreshold, P); // lower bound
    P = (P > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, P); // upper bound

    for(int i = 0; i < K; i++)
    {
        F[i] = PclusterK.middleRows(i * C, C);
        F[i].rowwise() /= F[i].colwise().sum(); // normalize F per site
    }
    Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
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
                            const MyArr2D & prev,
                            const MyArr2D & emit,
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
                    curr(zz, yy) = emit(zz, s) * (ss + sj + jj);
                }
}

void Phaseless::moveBackward(int ind,
                             int s,
                             MyArr2D & curr,
                             const MyArr1D & cs,
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
                    curr(zz, yy) = (ss + sj + jj) / cs(s);
                }
}

double Phaseless::runForwardBackwards(int ind, const MyFloat1D & GL, bool finalIter)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D emit = get_emission_by_gl(gli, P).transpose(); // CC x M
    std::vector<MyArr2D> alpha(M, MyArr2D(CC, KK));
    std::vector<MyArr2D> beta(M, MyArr2D(CC, KK));
    MyArr2D prevsum_z(C, KK), prevsum_zz(K, K), prevsum_zy(C, K);
    MyArr1D cs(M), prevsum_zzy(K);
    // init
    int s{0}, z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    for(z1 = 0; z1 < C; z1++)
        for(z2 = 0; z2 < C; z2++)
            for(y1 = 0; y1 < K; y1++)
                for(y2 = 0; y2 < K; y2++)
                {
                    zz = z1 * C + z2;
                    yy = y1 * K + y2;
                    alpha[s](zz, yy) = emit(zz, s) * Q(y1, ind) * F[y1](z1, s) * Q(y2, ind) * F[y2](z2, s);
                }
    cs(s) = alpha[s].sum();
    alpha[s] /= cs(s); // norm alpha
    // forward recursion
    for(s = 1; s < M; s++)
    {
        getForwardPrevSums(alpha[s - 1], prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy);
        moveForward(ind, s, alpha[s], alpha[s - 1], emit, prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy);
        cs(s) = alpha[s].sum();
        alpha[s] /= cs(s); // norm alpha
    }
    // backwards recursion
    s = M - 1;
    beta[s].setConstant(1.0);
    MyArr2D beta_mult_emit(CC, KK);
    double prevsum_zzyy{0};
    for(s = M - 2; s >= 0; s--)
    {
        for(int yy = 0; yy < KK; yy++) beta_mult_emit.col(yy) = beta[s + 1].col(yy) * emit.col(s + 1);
        getBackwardPrevSums(beta_mult_emit, prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy, prevsum_zzyy, ind, s + 1);
        moveBackward(ind, s + 1, beta[s], cs, beta_mult_emit, prevsum_z, prevsum_zz, prevsum_zy, prevsum_zzy,
                     prevsum_zzyy);
    }
    // get posterios
    getPosterios(ind, gli, emit, cs, alpha, beta);
    return (1 / cs).log().sum();
}

void Phaseless::getPosterios(int ind,
                             const MyArr2D & gli,
                             const MyArr2D & emit,
                             const MyArr1D & cs,
                             const std::vector<MyArr2D> & alpha,
                             const std::vector<MyArr2D> & beta)
{
    int s{0}, z1{0}, z2{0}, y1{0}, y2{0}, zz{0}, yy{0};
    MyArr2D ind_post_zg1(C, M), ind_post_zg2(C, M), ind_post_zz(CC, M);
    MyArr2D ind_post_y = MyArr2D::Zero(K, M);
    MyArr2D ind_post_zy = MyArr2D::Zero(C * K, M);
    MyArr2D prevsum_z(C, KK), prevsum_zz(K, K), prevsum_zy(C, K);
    MyArr1D prevsum_zzy(K);
    for(s = 0; s < M; s++)
    {
        MyArr2D gamma = alpha[s] * beta[s]; // post(z, y)
        if(debug && !(std::abs(1.0 - gamma.sum()) < 1e-6))
            cao.error(s, ", gamma sum is not approx 1.0 (tol = 1e-6)\n", gamma.sum());
        ind_post_zz.col(s) = gamma.rowwise().sum();
        auto ind_post_yy = gamma.colwise().sum();
        ind_post_y.col(s) = ind_post_yy.reshaped(K, K).colwise().sum();
        auto gamma_div_emit = ind_post_zz.col(s) / emit.col(s);
        for(z1 = 0; z1 < C; z1++)
        {
            ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - P(s, z1))
                                   * (gli(s, 0) * (1 - P.row(s)) + gli(s, 1) * P.row(s)).transpose())
                                      .sum();
            ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (P(s, z1))
                                   * (gli(s, 1) * (1 - P.row(s)) + gli(s, 2) * P.row(s)).transpose())
                                      .sum();
            for(y1 = 0; y1 < K; y1++)
                ind_post_zy(y1 * C + z1, s) = gamma(Eigen::seqN(z1, C, C), Eigen::seqN(y1, K, K)).sum();
        }
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
                            * ((1 - et(s)) * (1 - et(s)) * F[y1](z1, s) * Q(y1, ind) * F[y2](z2, s) * Q(y2, ind)
                               + et(s) * (1 - et(s)) * Q(y1, ind) * F[y1](z1, s)
                                     * (er(s) * prevsum_zy(z2, y2) + (1 - er(s)) * F[y2](z2, s) * prevsum_zzy(y2)));
                        ind_post_zy(y1 * C + z1, s) +=
                            ebc
                            * ((1 - et(s)) * Q(y1, ind) * F[y1](z1, s)
                                   * ((1 - et(s)) * F[y2](z2, s) * Q(y2, ind)
                                      + et(s) * (1 - er(s)) * F[y2](z2, s) * prevsum_zzy(y2)
                                      + et(s) * er(s) * prevsum_zy(z2, y2))
                               + et(s) * (1 - er(s)) * F[y1](z1, s)
                                     * ((1 - et(s)) * F[y2](z2, s) * Q(y2, ind) * prevsum_zzy(y1)
                                        + et(s) * (1 - er(s)) * F[y2](z2, s) * prevsum_zz(y1, y2)
                                        + et(s) * er(s) * prevsum_z(z2, yy)));
                    }
    }
    { // sum over all samples for updates
        std::lock_guard<std::mutex> lock(mutex_it);
        EclusterA1 += ind_post_zg1;
        EclusterA2 += ind_post_zg2;
        PclusterK += ind_post_zy; // cluster ancestry allele jump
        Q.col(ind) = ind_post_y.rowwise().sum();
    }
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
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        Phaseless faith(opts.K, opts.C, genome->nsamples, genome->pos[ic].size(), opts.seed);
        faith.initRecombination(genome->pos[ic]);
        faith.debug = opts.debug;
        for(int it = 0; it <= opts.nimpute; it++)
        {
            tim.clock();
            faith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == opts.nimpute)
                    res.emplace_back(
                        pool.enqueue(&Phaseless::runForwardBackwards, &faith, i, std::ref(genome->gls[ic]), true));
                else
                    res.emplace_back(
                        pool.enqueue(&Phaseless::runForwardBackwards, &faith, i, std::ref(genome->gls[ic]), false));
            }
            double loglike = 0;
            for(auto && ll : res) loglike += ll.get();
            res.clear(); // clear future and renew
            cao.print(tim.date(), "run single chunk", ic, ", iteration", it, ", likelihoods =", loglike, ",",
                      tim.reltime(), " sec");
            if(it != opts.nimpute) faith.updateIteration();
        }
        oanc << faith.Q.transpose().format(fmt) << "\n";
    }

    return 0;
}
