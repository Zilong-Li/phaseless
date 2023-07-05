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
    FastPhaseK2(int m, int n, int c, int seed);
    ~FastPhaseK2();

    // BOUNDING
    double minRate{0.1}, maxRate{100}; // threshold for R
    double clusterFreqThreshold{1e-4}; // threshold for PI
    double alleleEmitThreshold{1e-4}; // threshold for F

    // FLAGS
    bool debug{0};

    // SHARED VARIBALES
    const int M, N, C, C2; // C2 = C x C
    int nGrids, B;
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr1D pi; // C, PI in first SNP
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr2D Ezj; // C x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D Ezg1, Ezg2; // C x M
    Int1D dist; // physical position distance between two markers
    MyArr1D AF;
    double nGen;

    void initRecombination(const Int1D & pos, int B, double Ne);
    void initIteration();
    void updateIteration();
    double runWithOneThread(int, const MyFloat1D &);
    double forwardAndBackwardsLowRam(int, const MyFloat1D &, bool);
    auto forwardAndBackwardsHighRam(int, const MyFloat1D &, bool);
};

inline FastPhaseK2::FastPhaseK2(int m, int n, int c, int seed) : M(m), N(n), C(c), C2(c * c)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, alleleEmitThreshold,
                                                           1 - alleleEmitThreshold);
    GP.setZero(M * 3, N);
}

inline FastPhaseK2::~FastPhaseK2() {}

inline void FastPhaseK2::initRecombination(const Int1D & pos, int B_ = 1, double Ne = 20000)
{
    B = B_;
    nGrids = B > 1 ? (pos.size() + B - 1) / B : M;
    PI = MyArr2D::Ones(C, nGrids);
    PI.rowwise() /= PI.colwise().sum(); // normalize it per site
    if(B > 1)
    {
        auto [grid_pos, central] = turn_pos_into_grid(pos, B);
        dist = calc_grid_distance(grid_pos, central);
    }
    else
    {
        dist = calc_position_distance(pos);
    }
    nGen = 4 * Ne / C;
    R = calc_transRate_diploid(dist, nGen);
}

inline void FastPhaseK2::initIteration()
{
    if(debug) std::cerr << R << std::endl;
    // initial temp variables
    pi.setZero(C); // reset pi at first SNP
    Ezj.setZero(C, nGrids); // reset post(Z,j)
    Ezg1.setZero(C, M); // reset pos(Z,g)
    Ezg2.setZero(C, M); // reset pos(Z,g)
}

inline void FastPhaseK2::updateIteration()
{
    MyArr1D er = 1.0 - Ezj.colwise().sum() / N;
    // morgans per SNP assuming T=100, 0.5 cM/Mb
    // x1 <- exp(-nGen * minRate * dl/100/1000000) # lower
    // x2 <- exp(-nGen * maxRate * dl/100/1000000) # upper
    for(int i = 0; i < er.size(); i++)
    {
        double miner = std::exp(-nGen * maxRate * dist[i] / 1e8);
        double maxer = std::exp(-nGen * minRate * dist[i] / 1e8);
        er(i) = er(i) < miner ? miner : er(i);
        er(i) = er(i) > maxer ? maxer : er(i);
    }
    R.row(0) = er.square();
    R.row(1) = (1 - er) * er;
    R.row(2) = (1 - er).square();

    // update F
    F = (Ezg2 / (Ezg1 + Ezg2)).transpose();
    if(F.isNaN().any())
    {
        if(debug) cao.warn("NaN in F in FastPhaseK2 model. will fill it with AF");
        if(AF.size() == 0) cao.error("AF is not assigned!\n");
        for(int i = 0; i < M; i++) F.row(i) = F.row(i).isNaN().select(AF(i), F.row(i));
    }
    // map F to domain but no normalization
    F = (F < alleleEmitThreshold).select(alleleEmitThreshold, F); // lower bound
    F = (F > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, F); // upper bound

    // update PI(C, M) except the first snp
    // first we normalize Ezj so that each col sum to 1
    Ezj.col(0) = pi / pi.sum(); // now update the first SNP
    Ezj.rowwise() /= Ezj.colwise().sum();
    if(Ezj.isNaN().any() || (Ezj < clusterFreqThreshold).any())
    {
        // std::cerr << "reset values below threshold\n";
        Ezj = (Ezj < clusterFreqThreshold).select(0, Ezj); // reset to 0 first
        for(int i = 0; i < nGrids; i++)
        {
            // for columns with an entry below 0
            // each 0 entry becomes threshold
            // then rest re-scaled so whole thing has sum 1
            if(auto c = (Ezj.col(i) == 0).count() > 0)
            {
                double xsum = 1 - c * clusterFreqThreshold;
                double csum = Ezj.col(i).sum();
                Ezj.col(i) = (Ezj.col(i) > 0).select(Ezj.col(i) * xsum / csum, clusterFreqThreshold);
            }
        }
    }
    PI = Ezj;

    if(Ezj.isNaN().any()) cao.error(Ezj, "NaN in PI from FastPhaseK2\n");
    if(debug && !((1 - PI.colwise().sum()).abs() < 1e-3).all())
        cao.error(PI.colwise().sum(), "\ncolsum of PI is not 1.0!\n");
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
    const MyArr1D cs = forward_backwards_diploid(alpha, beta, emit, R, PI);
    if(debug && !((1 - ((alpha * beta).colwise().sum())).abs() < 1e-4).all())
        cao.error((alpha * beta).colwise().sum() / cs.transpose(), "\ngamma sum is not 1.0!\n");
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

    return (1 / cs).log().sum();
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
    if(nGrids < M)
    {
        MyArr2D alpha(C2, nGrids), beta(C2, nGrids);
        MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
        MyArr2D emitGrids = collapse_emission_by_grid(emit, B, nGrids);
        auto cs = forward_backwards_diploid(alpha, beta, emitGrids, R, PI);
        if(debug && !((1 - ((alpha * beta).colwise().sum())).abs() < 1e-4).all())
            cao.cerr((alpha * beta).colwise().sum() / cs.transpose(), "\ngamma sum is not 1.0!\n");

        MyArr1D gamma = alpha.col(0) * beta.col(0); //  gamma in first snp
        MyArr1D gamma1 = gamma.reshaped(C, C).colwise().sum();
        MyArr2D ind_post_zj = MyArr2D::Zero(C, nGrids);
        MyArr2D ind_post_zg1 = MyArr2D::Zero(C, M);
        MyArr2D ind_post_zg2 = MyArr2D::Zero(C, M);
        MyArr1D tmp_zg(4);

        int z1, z2, s, e, i, g = 0;
        // now get expectation of post(Z,J)
        s = g * B;
        e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
        for(z1 = 0; z1 < C; z1++)
        {
            for(i = s; i <= e; i++)
            {
                auto gamma_div_emit = gamma / emit.col(i);
                ind_post_zg1(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(i, z1))
                                       * (gli(i, 0) * (1 - F.row(i)) + gli(i, 1) * F.row(i)).transpose())
                                          .sum();
                ind_post_zg2(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(i, z1))
                                       * (gli(i, 1) * (1 - F.row(i)) + gli(i, 2) * F.row(i)).transpose())
                                          .sum();
                if(call_geno)
                {
                    for(z2 = 0; z2 < C; z2++)
                    {
                        int z12 = z1 * C + z2;
                        tmp_zg(0) = gli(i, 0) * (1 - F(i, z1)) * (1 - F(i, z2));
                        tmp_zg(1) = gli(i, 1) * (1 - F(i, z1)) * F(i, z2);
                        tmp_zg(2) = gli(i, 1) * F(i, z1) * (1 - F(i, z2));
                        tmp_zg(3) = gli(i, 2) * F(i, z1) * F(i, z2);
                        GP(3 * i + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                        GP(3 * i + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                        GP(3 * i + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
                    }
                }
            }
        }
        for(g = 1; g < nGrids; g++)
        {
            const MyArr1D beta_mult_emit = emitGrids.col(g) * beta.col(g); // C2
            gamma = (alpha.col(g) * beta.col(g)); // C2
            MyArr1D alphatmp(C);
            s = g * B;
            e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
            for(z1 = 0; z1 < C; z1++)
            {
                alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), g - 1).sum() * R(1, g);
                for(i = s; i <= e; i++)
                {
                    auto gamma_div_emit = gamma / emit.col(i);
                    ind_post_zg1(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(i, z1))
                                           * (gli(i, 0) * (1 - F.row(i)) + gli(i, 1) * F.row(i)).transpose())
                                              .sum();
                    ind_post_zg2(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(i, z1))
                                           * (gli(i, 1) * (1 - F.row(i)) + gli(i, 2) * F.row(i)).transpose())
                                              .sum();
                    if(call_geno)
                    {
                        for(z2 = 0; z2 < C; z2++)
                        {
                            int z12 = z1 * C + z2;
                            tmp_zg(0) = gli(i, 0) * (1 - F(i, z1)) * (1 - F(i, z2));
                            tmp_zg(1) = gli(i, 1) * (1 - F(i, z1)) * F(i, z2);
                            tmp_zg(2) = gli(i, 1) * F(i, z1) * (1 - F(i, z2));
                            tmp_zg(3) = gli(i, 2) * F(i, z1) * F(i, z2);
                            GP(3 * i + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                            GP(3 * i + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                            GP(3 * i + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
                        }
                    }
                }
            }
            alphatmp += PI.col(g) * R(2, g) * 1.0;
            for(z1 = 0; z1 < C; z1++)
                ind_post_zj(z1, g) =
                    cs(g) * (PI(z1, g) * alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
        }
        double indLogLike = (1 / cs).log().sum();
        return std::tuple(indLogLike, ind_post_zj, ind_post_zg1, ind_post_zg2, gamma1);
    }
    else
    {

        MyArr2D alpha(C2, M), beta(C2, M);
        MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
        auto cs = forward_backwards_diploid(alpha, beta, emit, R, PI);
        if(debug && !((1 - ((alpha * beta).colwise().sum())).abs() < 1e-4).all())
            cao.error((alpha * beta).colwise().sum() / cs.transpose(), "\ngamma sum is not 1.0!\n");
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
                ind_post_zj(z1, s) =
                    cs(s) * (PI(z1, s) * alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
        }
        double indLogLike = (1 / cs).log().sum();
        return std::tuple(indLogLike, ind_post_zj, ind_post_zg1, ind_post_zg2, gamma1);
    }
}

#endif // FASTPHASE_H_
