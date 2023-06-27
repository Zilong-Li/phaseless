/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/fastphase.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef FASTPHASE_H_
#define FASTPHASE_H_

#include "common.hpp"
#include <mutex>

class FastPhaseK2
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK2(const Int1D & pos, int n, int c, int seed);
    ~FastPhaseK2();

    // SHARED VARIBALES
    const int M, N, C, C2; // C2 = C x C
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr1D pi; // C, PI in first SNP
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr2D Ezj; // C x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D Ezg1, Ezg2; // C x M
    MyArr2D GZP1, GZP2; // M x C

    void initIteration(double tol = 1e-6);
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
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.001, 0.999);
    PI = MyArr2D::Ones(C, M);
    PI.rowwise() /= PI.colwise().sum(); // normalize it per site
    GP.setZero(M * 3, N);
    GZP1.setZero(M, C);
    GZP2.setZero(M, C);
    R = calc_transRate(pos, C);
}

inline FastPhaseK2::~FastPhaseK2() {}

inline void FastPhaseK2::initIteration(double tol)
{
    // map PI to domain with normalization
    if(PI.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    PI = (PI < tol).select(tol, PI); // lower bound
    PI = (PI > 1 - tol).select(1 - tol, PI); // upper bound
    PI.rowwise() /= PI.colwise().sum(); // normalize it per site
    // map F to domain but no normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in F\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
    // initial temp variables
    pi.setZero(C); // reset pi at first SNP
    Ezj.setZero(C, M); // reset post(Z,j)
    Ezg1.setZero(C, M); // reset pos(Z,g)
    Ezg2.setZero(C, M); // reset pos(Z,g)
}

inline void FastPhaseK2::updateIteration()
{
    // update R = e^-r = 1 - Ezj / N
    auto er = Ezj.colwise().sum() / N;
    R.row(0) = (1 - er).square();
    R.row(1) = (1 - er) * er;
    R.row(2) = er.square();
    // update F
    F = (Ezg2 / (Ezg1 + Ezg2)).transpose();
    // update PI(C, M) except the first snp
    PI = Ezj.colwise() / Ezj.rowwise().sum(); // make sure the first col of Ezj is 0s
    // updat PI in first SNP,
    PI.col(0) = pi / pi.sum();
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
        updateIteration();
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
    MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
    double indLogLike = forward_backwards_diploid(alpha, beta, emit, R, F, PI);
    MyArr1D ind_gamma_one = alpha.col(0) * beta.col(0); //  gamma in first snp
    MyArr1D gamma1 = ind_gamma_one.reshaped(C, C).colwise().sum();
    MyArr1D ind_post_zj(C);
    MyArr1D ind_post_zg1(C);
    MyArr1D ind_post_zg2(C);
    int z1, s = 0;
    // now get expectation of post(Z,J)
    auto gamma_div_emit = ind_gamma_one / emit.col(s); // C2
    for(z1 = 0; z1 < C; z1++)
    {

        ind_post_zg1(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                            * (gli(s, 0) * (1 - F(s, z1)) + gli(s, 1) * F(s, z1)))
                               .sum();
        ind_post_zg2(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                            * (gli(s, 1) * (1 - F(s, z1)) + gli(s, 2) * F(s, z1)))
                               .sum();
    }
    {
        std::lock_guard<std::mutex> lock(mutex_it);
        Ezg1.col(s) += ind_post_zg1;
        Ezg2.col(s) += ind_post_zg2;
        pi += gamma1;
    }
    for(s = 1; s < M; s++)
    {
        MyArr1D beta_mult_emit = emit.col(s) * beta.col(s); // C2
        auto gamma_div_emit = alpha.col(s) * beta.col(s) / emit.col(s); // C2
        MyArr1D alphatmp(C);
        for(z1 = 0; z1 < C; z1++)
        {
            ind_post_zg1(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                                * (gli(s, 0) * (1 - F(s, z1)) + gli(s, 1) * F(s, z1)))
                                   .sum();
            ind_post_zg2(z1) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                                * (gli(s, 1) * (1 - F(s, z1)) + gli(s, 2) * F(s, z1)))
                                   .sum();
            alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum() * R(1, s);
        }
        alphatmp += PI.col(s) * R(2, s);
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1) = PI(z1, s) * (alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
        { // sum over all samples, for update PI and F
            std::lock_guard<std::mutex> lock(mutex_it);
            Ezj.col(s) += ind_post_zj;
            Ezg1.col(s) += ind_post_zg1;
            Ezg2.col(s) += ind_post_zg2;
        }
    }

    return indLogLike;
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
    MyArr2D alpha(C2, M); // likelihood of forward recursion for ind i, not log
    MyArr2D beta(C2, M); // likelihood of backward recursion for ind i, not log
    MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
    double indLogLike = forward_backwards_diploid(alpha, beta, emit, R, F, PI);
    MyArr1D ind_gamma_one = alpha.col(0) * beta.col(0); //  gamma in first snp
    MyArr1D gamma1 = ind_gamma_one.reshaped(C, C).colwise().sum();
    MyArr2D ind_post_zj(C, M);
    MyArr2D ind_post_zg1(C, M);
    MyArr2D ind_post_zg2(C, M);
    int z1, s = 0;
    // now get expectation of post(Z,J)
    auto gamma_div_emit = ind_gamma_one / emit.col(s); // C2
    for(z1 = 0; z1 < C; z1++)
    {

        ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                               * (gli(s, 0) * (1 - F(s, z1)) + gli(s, 1) * F(s, z1)))
                                  .sum();
        ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                               * (gli(s, 1) * (1 - F(s, z1)) + gli(s, 2) * F(s, z1)))
                                  .sum();
    }
    for(s = 1; s < M; s++)
    {
        MyArr1D beta_mult_emit = emit.col(s) * beta.col(s); // C2
        auto gamma_div_emit = alpha.col(s) * beta.col(s) / emit.col(s); // C2
        MyArr1D alphatmp(C);
        for(z1 = 0; z1 < C; z1++)
        {
            ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(s, z1))
                                   * (gli(s, 0) * (1 - F(s, z1)) + gli(s, 1) * F(s, z1)))
                                      .sum();
            ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(s, z1))
                                   * (gli(s, 1) * (1 - F(s, z1)) + gli(s, 2) * F(s, z1)))
                                      .sum();
            alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum() * R(1, s);
        }
        alphatmp += PI.col(s) * R(2, s);
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1, s) = PI(z1, s) * (alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
    }

    return std::tuple(indLogLike, ind_post_zj, ind_post_zg1, ind_post_zg2, gamma1);
}

class FastPhaseK4
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK4(int n, int m, int c, int seed);
    ~FastPhaseK4();

    // SHARED VARIBALES
    std::ofstream ofs;
    const int N, M, C, C2; // C2 = C x C
    MyArr2D GP; // genotype probabilies for all individuals, N x (M x 3)
    MyArr2D PI; // nsnps x C
    MyArr2D F; // nsnps x C
    MyArr2D transHap; // nsnps x (C x C)
    MyArr2D transDip; // nsnps x (C x C x C x C)

    void updateClusterFreqPI(const MyArr2D &, double);
    void updateAlleleFreqWithinCluster(const MyArr2D &, double);
    void transitionCurIter(const MyArr1D &);
    auto callGenotypeInd(const MyArr2D &);
    void openClusterFile(std::string out);
    double forwardAndBackwards(int, const MyFloat1D &, MyArr2D &, MyArr2D &, bool);
};

inline FastPhaseK4::FastPhaseK4(int n, int m, int c, int seed)
: N(n), M(m), C(c), C2(c * c), GP(M * 3, N), transHap(M, C2), transDip(M, C2 * C2)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = PI.colwise() / PI.rowwise().sum(); // normalize it
}

inline FastPhaseK4::~FastPhaseK4()
{
    if(ofs.is_open()) ofs.close();
}

inline void FastPhaseK4::openClusterFile(std::string out)
{
    ofs.open(out, std::ios::binary);
    if(!ofs) throw std::runtime_error(out + ": " + strerror(errno));
    ofs.write((char *)&N, 4);
    ofs.write((char *)&M, 4);
    ofs.write((char *)&C, 4);
}

inline auto FastPhaseK4::callGenotypeInd(const MyArr2D & indPostProbsZandG)
{
    MyArr1D geno = MyArr1D::Zero(M * 3);
    int z1, z2, z12, g1, g2, g12, g3;
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    g3 = g1 + g2;
                    geno(Eigen::seqN(g3, M, 3)) += indPostProbsZandG.col(z12 * 4 + g12);
                }
            }
        }
    }
    return geno;
}

/*
** @param ind current individual i
** @param GL  genotype likelihood of all individuals in snp major form
** @return individual total likelihood
*/
inline double FastPhaseK4::forwardAndBackwards(int ind,
                                               const MyFloat1D & GL,
                                               MyArr2D & postProbsZ,
                                               MyArr2D & postProbsZandG,
                                               bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D emitDip(M, C2); // emission probabilies, nsnps x (C x C)
    int k1, k2, g1, g2, g12;
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.col(g1 + g2)
                                                * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                                * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                }
            }
        }
    }
    emitDip = emitDip.log();
    MyArr2D logLikeForwardInd(C2, M); // log likelihood of forward recursion for ind i
    MyArr2D logLikeBackwardInd(C2, M); // log likelihood of backward recursion for ind i

    // ======== forward recursion ===========
    int z1, z2, z12, z12s, z12e;
    double protect_me = 0; // for underflow protection
    // first site
    int s{0};
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            logLikeForwardInd(z12, s) = emitDip(s, z12) + PI(s, z1) * PI(s, z2);
        }
    }

    // do forward recursion
    for(s = 1; s < M; s++)
    {
        // MyArr1D likeForwardTmp;
        auto likeForwardTmp = Eigen::exp(logLikeForwardInd.col(s - 1) - protect_me);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                z12s = (z1 * C + z2) * C2; // fist index
                z12e = z12s + C2 - 1; // last index
                logLikeForwardInd(z12, s) =
                    protect_me + emitDip(s, z12)
                    + log((likeForwardTmp * transDip(s, Eigen::seq(z12s, z12e)).transpose()).sum());
            }
        }
        // protect_me = logLikeForwardInd.row(s).maxCoeff();
        protect_me = logLikeForwardInd(C2 - 1, s);
    }
    // total likelhoods of the individual
    double indLogLikeForwardAll = protect_me + log((logLikeForwardInd.col(M - 1) - protect_me).exp().sum());

    // ======== backward recursion ===========
    // set last site
    logLikeBackwardInd.col(M - 1).setZero();
    // site M-2 to 0
    protect_me = 0;
    for(s = M - 2; s >= 0; s--)
    {
        auto likeBackwardTmp =
            Eigen::exp(emitDip.row(s + 1).transpose() + logLikeBackwardInd.col(s + 1) - protect_me);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                logLikeBackwardInd(z12, s) =
                    protect_me
                    + log((likeBackwardTmp * transDip(s + 1, Eigen::seqN(z12, C2, C2)).transpose()).sum());
            }
        }
        // protect_me = logLikeBackwardInd.row(s).maxCoeff();
        protect_me = logLikeBackwardInd(C2 - 1, s);
    }

    std::lock_guard<std::mutex> lock(mutex_it); // lock here if RAM cost really matters

    // ======== post decoding get p(Z|X, G) ===========
    // MyArr2D indPostProbsZ; // post probabilities for ind i, M x (C x C)
    MyArr2D indPostProbsZ = (logLikeBackwardInd + logLikeForwardInd - indLogLikeForwardAll).exp().transpose();
    // ======== post decoding get p(Z,G|X, theta) ===========
    // MyArr2D indPostProbsZandG, M x (C x C x 2 x 2)
    MyArr2D indPostProbsZandG(M, C2 * 4);
    MyArr1D tmpSum(M);
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            tmpSum.setZero();
            z12 = z1 * C + z2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    indPostProbsZandG.col(z12 * 4 + g12) = gli.col(g1 + g2)
                                                           * (g1 * F.col(z1) + (1 - g1) * (1 - F.col(z1)))
                                                           * (g2 * F.col(z2) + (1 - g2) * (1 - F.col(z2)));
                    tmpSum += indPostProbsZandG.col(z12 * 4 + g12);
                }
            }
            indPostProbsZandG.middleCols(z12 * 4, 4).colwise() *= indPostProbsZ.col(z12);
            indPostProbsZandG.middleCols(z12 * 4, 4).colwise() /= tmpSum;
        }
    }

    if(call_geno)
    {
        // std::lock_guard<std::mutex> lock(mutex_it);
        GP.col(ind) = callGenotypeInd(indPostProbsZandG);
        // output likelihood of each cluster
        MyArr2D likeCluster = (logLikeBackwardInd + logLikeForwardInd).exp().transpose();
        ofs.write((char *)likeCluster.data(), M * C2 * 8);
    }
    postProbsZ += indPostProbsZ;
    postProbsZandG += indPostProbsZandG;
    // std::cout << std::this_thread::get_id() << ": " << ind << '\n';

    return indLogLikeForwardAll;
}

/*
** @param distRate distance or recombination rate between two markers
*/
inline void FastPhaseK4::transitionCurIter(const MyArr1D & distRate)
{
    int to1, to2, from1, from2;
    int k4d, z2d1, z2d2;
    for(from1 = 0; from1 < C; from1++)
    {
        for(to1 = 0; to1 < C; to1++)
        {
            z2d1 = from1 * C + to1;
            if(from1 == to1)
                transHap.col(z2d1) = distRate + (1 - distRate) * PI.col(to1);
            else
                transHap.col(z2d1) = (1 - distRate) * PI.col(to1);
        }
    }

    for(from1 = 0; from1 < C; from1++)
    {
        for(from2 = 0; from2 < C; from2++)
        {
            for(to1 = 0; to1 < C; to1++)
            {
                for(to2 = 0; to2 < C; to2++)
                {
                    k4d = to2 + (to1 + (from2 + from1 * C) * C) * C;
                    z2d1 = from1 * C + to1;
                    z2d2 = from2 * C + to2;
                    transDip.col(k4d) = transHap.col(z2d1) * transHap.col(z2d2);
                }
            }
        }
    }
}

inline void FastPhaseK4::updateClusterFreqPI(const MyArr2D & postProbsZ, double tol)
{
    int z1, z2, z12;
    PI.setZero();
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            PI.col(z1) += postProbsZ.col(z12);
            PI.col(z2) += postProbsZ.col(z12);
        }
    }
    PI /= 2 * N;
    // map to domain
    if(PI.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    PI = (PI < tol).select(tol, PI); // lower bound
    PI = (PI > 1 - tol).select(1 - tol, PI); // upper bound
    // normalize it now
    PI = PI.colwise() / PI.rowwise().sum();
}

inline void FastPhaseK4::updateAlleleFreqWithinCluster(const MyArr2D & postProbsZandG, double tol)
{
    int z1, z2, z12, g1, g2, g12;
    MyArr2D Ekg = MyArr2D::Zero(M, C * 2);
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    Ekg.col(z1 * 2 + g1) += postProbsZandG.col(z12 * 4 + g12);
                    Ekg.col(z2 * 2 + g2) += postProbsZandG.col(z12 * 4 + g12);
                }
            }
        }
    }
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2))
        / (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
    // map to domain but normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
}

#endif // FASTPHASE_H_
