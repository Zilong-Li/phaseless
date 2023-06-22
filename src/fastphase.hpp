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
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr2D J; // 3 x M, jumping / recombination rate
    MyArr2D Ek, Ekg; // M x C, M x C x 2
    MyArr2D GZP1, GZP2; // M x C

    void initIteration(double tol = 1e-6);
    void updateIteration();
    auto forwardAndBackwardsHighRam(int, const MyFloat1D &, bool);
    double forwardAndBackwardsLowRam(int, const MyFloat1D &, bool);
    double runWithOneThread(int, const MyFloat1D &);
};

inline FastPhaseK2::FastPhaseK2(const Int1D & pos, int n, int c, int seed) : M(pos.size()), N(n), C(c), C2(c * c)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.001, 0.999);
    PI = MyArr2D::Ones(C, M);
    PI.rowwise() /= PI.colwise().sum(); // normalize it per site
    GP.setZero(M * 3, N);
    GZP1.setZero(M, C);
    GZP2.setZero(M, C);
    J = calc_transRate(pos, C);
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
    Ek.setZero(M, C);
    Ekg.setZero(M, C * 2);
}

inline void FastPhaseK2::updateIteration()
{
    PI = Ek.transpose() / (2 * N);
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2))
        / (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
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
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @return individual total likelihood
*/
inline double FastPhaseK2::forwardAndBackwardsLowRam(int ind, const MyFloat1D & GL, bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D LikeForwardInd(C2, M); // likelihood of forward recursion for ind i, not log
    MyArr2D LikeBackwardInd(C2, M); // likelihood of backward recursion for ind i, not log
    double indLogLikeForwardAll = getClusterLikelihoods(LikeForwardInd, LikeBackwardInd, gli, J, PI, F, true);
    MyArr1D ind_post_z_col(M); // col of indPostProbsZ
    MyArr2D ind_post_z_g(M, 4); // cols of indPostProbsZandG
    LikeForwardInd *= LikeBackwardInd;
    LikeForwardInd.rowwise() /= LikeForwardInd.colwise().sum(); // normlize it so that colwise.sum()==1
    int g1, g2, g3, g12, z1, z2, z12;
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            ind_post_z_col = LikeForwardInd.row(z12).transpose();
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    ind_post_z_g.col(g12) = gli.col(g1 + g2) * (g1 * F.col(z1) + (1 - g1) * (1 - F.col(z1)))
                                            * (g2 * F.col(z2) + (1 - g2) * (1 - F.col(z2)));
                }
            }
            ind_post_z_g.colwise() *= ind_post_z_col / ind_post_z_g.rowwise().sum();
            if(call_geno)
            {
                for(g1 = 0; g1 < 2; g1++)
                {
                    for(g2 = 0; g2 < 2; g2++)
                    {
                        g12 = g1 * 2 + g2;
                        g3 = g1 + g2;
                        GP(Eigen::seqN(g3, M, 3), ind) += ind_post_z_g.col(g12);
                    }
                }
            }
            { // sum over all samples, for update PI and F
                std::lock_guard<std::mutex> lock(mutex_it);
                Ek.col(z1) += ind_post_z_col;
                Ek.col(z2) += ind_post_z_col;
                for(g1 = 0; g1 < 2; g1++)
                {
                    for(g2 = 0; g2 < 2; g2++)
                    {
                        g12 = g1 * 2 + g2;
                        Ekg.col(z1 * 2 + g1) += ind_post_z_g.col(g12);
                        Ekg.col(z2 * 2 + g2) += ind_post_z_g.col(g12);
                    }
                }
                if(call_geno)
                {
                    GZP1.col(z1) += ind_post_z_g.col(1) + ind_post_z_g.col(2);
                    GZP1.col(z2) += ind_post_z_g.col(1) + ind_post_z_g.col(2);
                    GZP2.col(z1) += ind_post_z_g.col(3);
                    GZP2.col(z2) += ind_post_z_g.col(3);
                }
            }
        }
    }

    return indLogLikeForwardAll;
}
/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @return individual total likelihood
*/
inline auto FastPhaseK2::forwardAndBackwardsHighRam(int ind, const MyFloat1D & GL, bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D LikeForwardInd(C2, M); // likelihood of forward recursion for ind i, not log
    MyArr2D LikeBackwardInd(C2, M); // likelihood of backward recursion for ind i, not log
    double indLogLikeForwardAll = getClusterLikelihoods(LikeForwardInd, LikeBackwardInd, gli, J, PI, F, true);
    MyArr1D ind_post_z_col(M); // col of indPostProbsZ
    MyArr2D ind_post_z_g(M, 4); // cols of indPostProbsZandG
    LikeForwardInd *= LikeBackwardInd;
    LikeForwardInd.rowwise() /= LikeForwardInd.colwise().sum(); // normlize it so that colwise.sum()==1
    MyArr2D iEk = MyArr2D::Zero(M, C), iEkg = MyArr2D::Zero(M, C * 2); // MxC, MxCx2
    int g1, g2, g3, g12, z1, z2, z12;
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            ind_post_z_col = LikeForwardInd.row(z12).transpose();
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    ind_post_z_g.col(g12) = gli.col(g1 + g2) * (g1 * F.col(z1) + (1 - g1) * (1 - F.col(z1)))
                                            * (g2 * F.col(z2) + (1 - g2) * (1 - F.col(z2)));
                }
            }
            ind_post_z_g.colwise() *= ind_post_z_col / ind_post_z_g.rowwise().sum();
            if(call_geno)
            {
                for(g1 = 0; g1 < 2; g1++)
                {
                    for(g2 = 0; g2 < 2; g2++)
                    {
                        g12 = g1 * 2 + g2;
                        g3 = g1 + g2;
                        GP(Eigen::seqN(g3, M, 3), ind) += ind_post_z_g.col(g12);
                    }
                }
            }
            { // for update PI and F
                iEk.col(z1) += ind_post_z_col;
                iEk.col(z2) += ind_post_z_col;
                for(g1 = 0; g1 < 2; g1++)
                {
                    for(g2 = 0; g2 < 2; g2++)
                    {
                        g12 = g1 * 2 + g2;
                        iEkg.col(z1 * 2 + g1) += ind_post_z_g.col(g12);
                        iEkg.col(z2 * 2 + g2) += ind_post_z_g.col(g12);
                    }
                }
            }
        }
    }

    return std::tuple(indLogLikeForwardAll, iEk, iEkg);
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
    auto emitDip = emissionCurIterInd(gli, F, true);
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
        auto likeBackwardTmp = Eigen::exp(emitDip.row(s + 1).transpose() + logLikeBackwardInd.col(s + 1) - protect_me);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                logLikeBackwardInd(z12, s) =
                    protect_me + log((likeBackwardTmp * transDip(s + 1, Eigen::seqN(z12, C2, C2)).transpose()).sum());
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
    int g1, g2, g12;
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
