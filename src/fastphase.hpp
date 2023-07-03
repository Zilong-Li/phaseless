/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/fastphase.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef FASTPHASE_H_
#define FASTPHASE_H_

#include "common.hpp"
#include <cmath>
#include <mutex>

class FastPhaseK2
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK2(const Int1D & pos, int n, int c, int seed);
    ~FastPhaseK2();

    // BOUNDING
    const double alphaMatThreshold = 1e-4; // Ezj
    const double emissionThreshold = 1e-4;

    // FLAGS
    bool debug = true;

    // SHARED VARIBALES
    const int M, N, C, C2; // C2 = C x C
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr1D pi; // C, PI in first SNP
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr2D Ezj; // C x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D Ezg1, Ezg2; // C x M

    void initIteration();
    void updateIteration();
    double runWithOneThread(int, const MyFloat1D &);
    double forwardAndBackwardsLowRam(int, const MyFloat1D &, bool);
    auto forwardAndBackwardsHighRam(int, const MyFloat1D &, bool);
};

inline FastPhaseK2::FastPhaseK2(const Int1D & pos, int n, int c, int seed)
: M(pos.size()), N(n), C(c), C2(c * c)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, emissionThreshold,
                                                           1 - emissionThreshold);
    PI = MyArr2D::Ones(C, M);
    PI.rowwise() /= PI.colwise().sum(); // normalize it per site
    GP.setZero(M * 3, N);
    R = calc_transRate(pos, C);
}

inline FastPhaseK2::~FastPhaseK2() {}

inline void FastPhaseK2::initIteration()
{
    std::cerr << R << std::endl;
    // initial temp variables
    pi.setZero(C); // reset pi at first SNP
    Ezj.setZero(C, M); // reset post(Z,j)
    Ezg1.setZero(C, M); // reset pos(Z,g)
    Ezg2.setZero(C, M); // reset pos(Z,g)
}

inline void FastPhaseK2::updateIteration()
{
    // morgans per SNP assuming T=100, 0.5 cM/Mb
    // x1 <- exp(-nGen * minRate * dl/100/1000000) # lower
    // x2 <- exp(-nGen * maxRate * dl/100/1000000) # upper
    // update R, e^-r = 1 - Ezj / N
    MyArr1D er = 1.0 - Ezj.colwise().sum() / N;
    const double miner = 0.9;
    const double maxer = std::exp(-1.0e-9);
    er = (er < miner).select(miner, er);
    er = (er > maxer).select(maxer, er);
    R.row(0) = er.square();
    R.row(1) = (1 - er) * er;
    R.row(2) = (1 - er).square();

    // update F
    F = (Ezg2 / (Ezg1 + Ezg2)).transpose();
    // map F to domain but no normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in F\n");
    // F = (F < 1e-6).select(1e-6, F); // lower bound
    // F = (F > 1 - 1e-6).select(1 - 1e-6, F); // upper bound

    // update PI(C, M) except the first snp
    // first we normalize Ezj so that each col sum to 1
    Ezj.col(0) = pi / pi.sum(); // now update the first SNP
    Ezj.rowwise() /= Ezj.colwise().sum();
    // what if Ezj.colwise().sum() = 0 ? of course, the first col is nan
    // Ezj.col(0) = pi / pi.sum(); // now update the first SNP
    // Ezj = (Ezj < alphaMatThreshold).select(alphaMatThreshold, Ezj); // reset to 0 first
    // Ezj.rowwise() /= Ezj.colwise().sum();

    if(Ezj.isNaN().any() || (Ezj < alphaMatThreshold).any())
    {
        // std::cerr << "reset values below threshold\n";
        Ezj = (Ezj < alphaMatThreshold).select(0, Ezj); // reset to 0 first
        for(int i = 1; i < M; i++)
        {
            // for columns with an entry below 0
            // each 0 entry becomes threshold
            // then rest re-scaled so whole thing has sum 1
            if(auto c = (Ezj.col(i) == 0).count() > 0)
            {
                double xsum = 1 - c * alphaMatThreshold;
                double csum = Ezj.col(i).sum();
                Ezj.col(i) = (Ezj.col(i) > 0).select(Ezj.col(i) * xsum / csum, alphaMatThreshold);
            }
        }
    }
    // Ezj.col(0) = pi / pi.sum(); // now update the first SNP

    if(Ezj.isNaN().any())
    {
        std::cerr << Ezj << "\n";
        throw std::runtime_error("NaN in PI\n");
    }
    PI = Ezj;
    if(debug && !((1 - PI.colwise().sum()).abs() < 1e-3).all())
    {
        std::cerr << PI.colwise().sum() << "\n";
        throw std::runtime_error("colsum of PI is not 1.0!\n");
    }
}

/*
** @param niters    number of iterations
** @param GL        genotype likelihood of all individuals in snp major form
** @param pos       SNP position
** @return likelihood difference between last two iters
*/
inline double FastPhaseK2::runWithOneThread(int niters, const MyFloat1D & GL)
{
    double loglike, diff, prevlike;
    for(int it = 0; it <= niters; it++)
    {
        initIteration();
        loglike = 0;
        for(int i = 0; i < N; i++)
        {
            if(it == niters)
                loglike += forwardAndBackwardsLowRam(i, GL, true);
            else
                loglike += forwardAndBackwardsLowRam(i, GL, false);
        }
        diff = it ? loglike - prevlike : 0;
        prevlike = loglike;
        if(it != niters) updateIteration();
    }
    return diff;
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param call_geno boolean, call genotype of not
** @return individual total likelihood
*/
inline double FastPhaseK2::forwardAndBackwardsLowRam(int ind, const MyFloat1D & GL, bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D alpha(C2, M), beta(C2, M);
    const MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
    const MyArr1D cs = forward_backwards_diploid(alpha, beta, emit, R, F, PI);
    if(debug && !((1 - ((alpha * beta).colwise().sum())).abs() < 1e-4).all())
    {
        std::cerr << (alpha * beta).colwise().sum() / cs.transpose() << "\n";
        throw std::runtime_error("gamma sum is not 1.0!\n");
    }
    const MyArr1D ind_gamma_one = alpha.col(0) * beta.col(0); //  gamma in first snp
    const MyArr1D gamma1_sum = ind_gamma_one.reshaped(C, C).colwise().sum();
    MyArr1D ind_post_zj(C);
    MyArr1D ind_post_zg1(C);
    MyArr1D ind_post_zg2(C);
    MyArr1D tmp_zg(4);
    int z1, z2, s = 0;
    // now get expectation of post(Z,J)
    const MyArr1D gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // C2
    for(z1 = 0; z1 < C; z1++)
    {

        ind_post_zg1(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                            * (gli(s, 0) * (1 - F.row(s)) + gli(s, 1) * F.row(s)).transpose())
                               .sum();
        ind_post_zg2(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                            * (gli(s, 1) * (1 - F.row(s)) + gli(s, 2) * F.row(s)).transpose())
                               .sum();
        if(call_geno)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                int z12 = z1 * C + z2;
                tmp_zg(0) = gli(s, 0) * (1 - F(s, z1)) * (1 - F(s, z2));
                tmp_zg(1) = gli(s, 1) * (1 - F(s, z1)) * F(s, z2);
                tmp_zg(2) = gli(s, 1) * F(s, z1) * (1 - F(s, z2));
                tmp_zg(3) = gli(s, 2) * F(s, z1) * F(s, z2);
                GP(3 * s + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                GP(3 * s + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                GP(3 * s + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
            }
        }
    }
    {
        std::lock_guard<std::mutex> lock(mutex_it);
        Ezg1.col(s) += ind_post_zg1;
        Ezg2.col(s) += ind_post_zg2;
        pi += gamma1_sum;
    }
    for(s = 1; s < M; s++)
    {
        const MyArr1D beta_mult_emit = emit.col(s) * beta.col(s); // C2
        const MyArr1D gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // C2
        MyArr1D alphatmp(C);
        for(z1 = 0; z1 < C; z1++)
        {
            ind_post_zg1(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                                * (gli(s, 0) * (1 - F.row(s)) + gli(s, 1) * F.row(s)).transpose())
                                   .sum();
            ind_post_zg2(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                                * (gli(s, 1) * (1 - F.row(s)) + gli(s, 2) * F.row(s)).transpose())
                                   .sum();
            alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum() * R(1, s);
            if(call_geno)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    int z12 = z1 * C + z2;
                    tmp_zg(0) = gli(s, 0) * (1 - F(s, z1)) * (1 - F(s, z2));
                    tmp_zg(1) = gli(s, 1) * (1 - F(s, z1)) * F(s, z2);
                    tmp_zg(2) = gli(s, 1) * F(s, z1) * (1 - F(s, z2));
                    tmp_zg(3) = gli(s, 2) * F(s, z1) * F(s, z2);
                    GP(3 * s + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                    GP(3 * s + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                    GP(3 * s + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
                }
            }
        }
        alphatmp += PI.col(s) * R(2, s) * 1.0; // inner alpha.col(s-1).sum == 1
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1) = cs(s) * PI(z1, s) * (alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
        { // sum over all samples for updates
            std::lock_guard<std::mutex> lock(mutex_it);
            Ezj.col(s) += ind_post_zj;
            Ezg1.col(s) += ind_post_zg1;
            Ezg2.col(s) += ind_post_zg2;
        }
    }

    // return (1 / cs).log().sum();
    return log((1 / cs).sum());
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param call_geno boolean, call genotype of not
** @return individual total likelihood
*/
inline auto FastPhaseK2::forwardAndBackwardsHighRam(int ind, const MyFloat1D & GL, bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D alpha(C2, M), beta(C2, M);
    MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
    auto cs = forward_backwards_diploid(alpha, beta, emit, R, F, PI);
    if(debug && !((1 - ((alpha * beta).colwise().sum())).abs() < 1e-4).all())
    {
        std::cerr << (alpha * beta).colwise().sum() << "\n";
        throw std::runtime_error("gamma sum is not 1.0!\n");
    }
    const MyArr1D ind_gamma_one = alpha.col(0) * beta.col(0); //  gamma in first snp
    MyArr1D gamma1 = ind_gamma_one.reshaped(C, C).colwise().sum();
    MyArr2D ind_post_zj = MyArr2D::Zero(C, M);
    MyArr2D ind_post_zg1 = MyArr2D::Zero(C, M);
    MyArr2D ind_post_zg2 = MyArr2D::Zero(C, M);
    MyArr1D tmp_zg(4);

    int z1, z2, s = 0;
    // now get expectation of post(Z,J)
    const MyArr1D gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // C2
    for(z1 = 0; z1 < C; z1++)
    {

        ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                               * (gli(s, 0) * (1 - F.row(s)) + gli(s, 1) * F.row(s)).transpose())
                                  .sum();
        ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                               * (gli(s, 1) * (1 - F.row(s)) + gli(s, 2) * F.row(s)).transpose())
                                  .sum();
        if(call_geno)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                int z12 = z1 * C + z2;
                tmp_zg(0) = gli(s, 0) * (1 - F(s, z1)) * (1 - F(s, z2));
                tmp_zg(1) = gli(s, 1) * (1 - F(s, z1)) * F(s, z2);
                tmp_zg(2) = gli(s, 1) * F(s, z1) * (1 - F(s, z2));
                tmp_zg(3) = gli(s, 2) * F(s, z1) * F(s, z2);
                GP(3 * s + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                GP(3 * s + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                GP(3 * s + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
            }
        }
    }
    for(s = 1; s < M; s++)
    {
        const MyArr1D beta_mult_emit = emit.col(s) * beta.col(s); // C2
        const MyArr1D gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // C2
        MyArr1D alphatmp(C);
        for(z1 = 0; z1 < C; z1++)
        {
            alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum() * R(1, s);
            ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                                   * (gli(s, 0) * (1 - F.row(s)) + gli(s, 1) * F.row(s)).transpose())
                                      .sum();
            ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                                   * (gli(s, 1) * (1 - F.row(s)) + gli(s, 2) * F.row(s)).transpose())
                                      .sum();
            if(call_geno)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    int z12 = z1 * C + z2;
                    tmp_zg(0) = gli(s, 0) * (1 - F(s, z1)) * (1 - F(s, z2));
                    tmp_zg(1) = gli(s, 1) * (1 - F(s, z1)) * F(s, z2);
                    tmp_zg(2) = gli(s, 1) * F(s, z1) * (1 - F(s, z2));
                    tmp_zg(3) = gli(s, 2) * F(s, z1) * F(s, z2);
                    GP(3 * s + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                    GP(3 * s + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                    GP(3 * s + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
                }
            }
        }
        alphatmp += PI.col(s) * R(2, s) * 1.0;
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1, s) = cs(s) * (PI(z1, s) * alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
    }

    // if(debug && call_geno) std::cerr << GP.col(ind).sum() << "\n";

    double indLogLike = (1 / cs).log().sum();
    return std::tuple(indLogLike, ind_post_zj, ind_post_zg1, ind_post_zg2, gamma1);
}

#endif // FASTPHASE_H_
