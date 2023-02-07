#ifndef FASTPHASE_H_
#define FASTPHASE_H_

#include "common.hpp"
#include <fstream>
#include <iostream>
#include <mutex>

class FastPhaseK2
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK2(int n, int m, int c, int seed, const std::string & out);
    ~FastPhaseK2();

    // SHARED VARIBALES
    std::ofstream ofs;
    const int N, M, C, C2; // C2 = C x C
    ArrDouble2D GP; // genotype probabilies for all individuals, N x (M x 3)
    ArrDouble2D PI, Ek; // nsnps x C
    ArrDouble2D F; // M x C
    ArrDouble2D Ekg; // M x C x 2

    void initIteration();
    ArrDouble2D getClusterLikelihoods(int ind, const DoubleVec1D & GL, const ArrDouble2D & transRate);
    void updateClusterFreqPI(double tol);
    void updateAlleleFreqWithinCluster(double tol);
    double forwardAndBackwards(int ind, const DoubleVec1D & GL, const ArrDouble2D & transRate, bool call_geno);
};

inline FastPhaseK2::FastPhaseK2(int n, int m, int c, int seed, const std::string & out)
: N(n), M(m), C(c), C2(c * c), GP(M * 3, N), Ek(M, C), Ekg(M, C * 2)
{
    ofs.open(out, std::ios::binary);
    if(!ofs) throw std::runtime_error(out + ": " + strerror(errno));
    ofs.write((char *)&N, 4);
    ofs.write((char *)&M, 4);
    ofs.write((char *)&C, 4);
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<ArrDouble2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = RandomUniform<ArrDouble2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = PI.colwise() / PI.rowwise().sum(); // normalize it
}

inline FastPhaseK2::~FastPhaseK2()
{
    ofs.close();
}

inline void FastPhaseK2::initIteration()
{
    Ek.setZero();
    Ekg.setZero();
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @return cluster likelihoods
*/
inline ArrDouble2D FastPhaseK2::getClusterLikelihoods(int ind, const DoubleVec1D & GL, const ArrDouble2D & transRate)
{
    Eigen::Map<const ArrDouble2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, F, false);
    ArrDouble2D LikeForwardInd(C2, M); // likelihood of forward recursion for ind i, not log
    ArrDouble2D LikeBackwardInd(C2, M); // likelihood of backward recursion for ind i, not log
    ArrDouble1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    ArrDouble1D cs(M);
    double constTmp;

    // ======== forward recursion ===========
    int k1, k2, k12;
    // first site
    int s{0};
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            LikeForwardInd(k12, s) = emitDip(s, k12) * PI(s, k1) * PI(s, k2);
        }
    }
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it

    // do forward recursion
    for(s = 1; s < M; s++)
    {
        // precompute outside of internal double loop
        sumTmp1.setZero();
        sumTmp2.setZero();
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += LikeForwardInd(k12, s - 1);
                sumTmp2(k2) += LikeForwardInd(k12, s - 1);
            }
        }
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        sumTmp1 *= transRate(1, s);
        sumTmp2 *= transRate(1, s);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeForwardInd(k12, s) = emitDip(s, k12)
                                         * (LikeForwardInd(k12, s - 1) * transRate(0, s) + PI(s, k1) * sumTmp1(k2)
                                            + PI(s, k2) * sumTmp2(k1) + PI(s, k1) * PI(s, k2) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }

    // ======== backward recursion ===========
    // set last site
    s = M - 1;
    LikeBackwardInd.col(s).setOnes(); // not log scale
    LikeBackwardInd.col(s) *= cs(s);
    // site M-2 to 0
    for(s = M - 2; s >= 0; s--)
    {
        // precompute outside of internal double loop
        sumTmp1.setZero();
        sumTmp2.setZero();
        constTmp = 0;
        auto beta_mult_emit = emitDip.row(s + 1).transpose() * LikeBackwardInd.col(s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI(s + 1, k2);
                sumTmp2(k2) += beta_mult_emit(k12) * PI(s + 1, k1);
                constTmp += beta_mult_emit(k12) * PI(s + 1, k1) * PI(s + 1, k2);
            }
        }
        // add transRate to constants
        sumTmp1 *= transRate(1, s + 1);
        sumTmp2 *= transRate(1, s + 1);
        constTmp *= transRate(2, s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeBackwardInd(k12, s) =
                    beta_mult_emit(k12) * transRate(0, s + 1) + sumTmp1(k1) + sumTmp2(k2) + constTmp;
            }
        }
        // apply scaling
        LikeBackwardInd.col(s) *= cs(s);
    }

    auto clusterLike = LikeForwardInd * LikeBackwardInd;
    return clusterLike;
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @return individual total likelihood
*/
inline double FastPhaseK2::forwardAndBackwards(int ind,
                                               const DoubleVec1D & GL,
                                               const ArrDouble2D & transRate,
                                               bool call_geno)
{
    Eigen::Map<const ArrDouble2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, F, false);
    ArrDouble2D LikeForwardInd(C2, M); // likelihood of forward recursion for ind i, not log
    ArrDouble2D LikeBackwardInd(C2, M); // likelihood of backward recursion for ind i, not log
    ArrDouble1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    ArrDouble1D cs(M);
    double constTmp;

    // ======== forward recursion ===========
    int k1, k2, k12;
    // first site
    int s{0};
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            LikeForwardInd(k12, s) = emitDip(s, k12) * PI(s, k1) * PI(s, k2);
        }
    }
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it

    // do forward recursion
    for(s = 1; s < M; s++)
    {
        // precompute outside of internal double loop
        sumTmp1.setZero();
        sumTmp2.setZero();
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += LikeForwardInd(k12, s - 1);
                sumTmp2(k2) += LikeForwardInd(k12, s - 1);
            }
        }
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        sumTmp1 *= transRate(1, s);
        sumTmp2 *= transRate(1, s);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeForwardInd(k12, s) = emitDip(s, k12)
                                         * (LikeForwardInd(k12, s - 1) * transRate(0, s) + PI(s, k1) * sumTmp1(k2)
                                            + PI(s, k2) * sumTmp2(k1) + PI(s, k1) * PI(s, k2) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    // total likelhoods of the individual
    double indLogLikeForwardAll = log((LikeForwardInd.col(M - 1) / cs(M - 1)).sum());

    // ======== backward recursion ===========
    // set last site
    s = M - 1;
    LikeBackwardInd.col(s).setOnes(); // not log scale
    LikeBackwardInd.col(s) *= cs(s);
    // site M-2 to 0
    for(s = M - 2; s >= 0; s--)
    {
        // precompute outside of internal double loop
        sumTmp1.setZero();
        sumTmp2.setZero();
        constTmp = 0;
        auto beta_mult_emit = emitDip.row(s + 1).transpose() * LikeBackwardInd.col(s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI(s + 1, k2);
                sumTmp2(k2) += beta_mult_emit(k12) * PI(s + 1, k1);
                constTmp += beta_mult_emit(k12) * PI(s + 1, k1) * PI(s + 1, k2);
            }
        }
        // add transRate to constants
        sumTmp1 *= transRate(1, s + 1);
        sumTmp2 *= transRate(1, s + 1);
        constTmp *= transRate(2, s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeBackwardInd(k12, s) =
                    beta_mult_emit(k12) * transRate(0, s + 1) + sumTmp1(k1) + sumTmp2(k2) + constTmp;
            }
        }
        // apply scaling
        LikeBackwardInd.col(s) *= cs(s);
    }

    cs *= LikeForwardInd.col(M - 1).sum(); // get back last forward likelihood

    // std::lock_guard<std::mutex> lock(mutex_it); // lock here if RAM cost really matters
    // ======== post decoding get p(Z|X, G),  M x (C x C) ===========
    // ======== post decoding get p(Z,G|X, theta), M x (C x C x 2 x 2) =======
    ArrDouble1D ind_post_z_col(M); // col of indPostProbsZ
    ArrDouble2D ind_post_z_g(M, 4); // cols of indPostProbsZandG
    ArrDouble1D tmpSum(M);
    ArrDouble1D geno;
    if(call_geno) geno.setZero(M * 3);
    int g1, g2, g12, g3;
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            tmpSum.setZero();
            k12 = k1 * C + k2;
            ind_post_z_col = (LikeForwardInd.row(k12) * LikeBackwardInd.row(k12)).transpose() / cs;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    ind_post_z_g.col(g12) = gli.row(g1 + g2).transpose() * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                            * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                    tmpSum += ind_post_z_g.col(g12);
                }
            }
            ind_post_z_g.colwise() *= ind_post_z_col;
            ind_post_z_g.colwise() /= tmpSum;
            if(call_geno)
            {
                for(g1 = 0; g1 < 2; g1++)
                {
                    for(g2 = 0; g2 < 2; g2++)
                    {
                        g12 = g1 * 2 + g2;
                        g3 = g1 + g2;
                        geno(Eigen::seqN(g3, M, 3)) += ind_post_z_g.col(g12);
                    }
                }
            }
            // for update PI and F
            {
                std::lock_guard<std::mutex> lock(mutex_it);
                Ek.col(k1) += ind_post_z_col;
                Ek.col(k2) += ind_post_z_col;
                for(g1 = 0; g1 < 2; g1++)
                {
                    for(g2 = 0; g2 < 2; g2++)
                    {
                        g12 = g1 * 2 + g2;
                        Ekg.col(k1 * 2 + g1) += ind_post_z_g.col(g12);
                        Ekg.col(k2 * 2 + g2) += ind_post_z_g.col(g12);
                    }
                }
            }
        }
    }

    if(call_geno)
    {
        std::lock_guard<std::mutex> lock(mutex_it);
        GP.col(ind) = geno;
        // output likelihood of each cluster
        ArrDouble2D likeCluster = (LikeBackwardInd * LikeForwardInd).transpose();
        ofs.write((char *)likeCluster.data(), M * C2 * 8);
    }

    return indLogLikeForwardAll;
}

inline void FastPhaseK2::updateClusterFreqPI(double tol)
{
    PI = Ek / (2 * N);
    // map to domain
    if(PI.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    PI = (PI < tol).select(tol, PI); // lower bound
    PI = (PI > 1 - tol).select(1 - tol, PI); // upper bound
    // normalize it now
    PI = PI.colwise() / PI.rowwise().sum();
}

inline void FastPhaseK2::updateAlleleFreqWithinCluster(double tol)
{
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2))
        / (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
    // map to domain but no normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
}

class FastPhaseK4
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK4(int n, int m, int c, int seed, const std::string & out);
    ~FastPhaseK4();

    // SHARED VARIBALES
    std::ofstream ofs;
    const int N, M, C, C2; // C2 = C x C
    ArrDouble2D GP; // genotype probabilies for all individuals, N x (M x 3)
    ArrDouble2D PI; // nsnps x C
    ArrDouble2D F; // nsnps x C
    ArrDouble2D transHap; // nsnps x (C x C)
    ArrDouble2D transDip; // nsnps x (C x C x C x C)

    void updateClusterFreqPI(const ArrDouble2D & postProbsZ, double tol);
    void updateAlleleFreqWithinCluster(const ArrDouble2D & postProbsZandG, double tol);
    void transitionCurIter(const ArrDouble1D & distRate);
    ArrDouble2D callGenotypeInd(const ArrDouble2D & indPostProbsZandG);
    double forwardAndBackwards(int ind,
                               const DoubleVec1D & GL,
                               ArrDouble2D & postProbsZ,
                               ArrDouble2D & postProbsZandG,
                               bool call_geno = false);
};

inline FastPhaseK4::FastPhaseK4(int n, int m, int c, int seed, const std::string & out)
: N(n), M(m), C(c), C2(c * c), GP(M * 3, N), transHap(M, C2), transDip(M, C2 * C2)
{
    ofs.open(out, std::ios::binary);
    if(!ofs) throw std::runtime_error(out + ": " + strerror(errno));
    ofs.write((char *)&N, 4);
    ofs.write((char *)&M, 4);
    ofs.write((char *)&C, 4);
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<ArrDouble2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = RandomUniform<ArrDouble2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = PI.colwise() / PI.rowwise().sum(); // normalize it
}

inline FastPhaseK4::~FastPhaseK4()
{
    ofs.close();
}

/*
** @param ind current individual i
** @param GL  genotype likelihood of all individuals in snp major form
** @return individual total likelihood
*/
inline double FastPhaseK4::forwardAndBackwards(int ind,
                                               const DoubleVec1D & GL,
                                               ArrDouble2D & postProbsZ,
                                               ArrDouble2D & postProbsZandG,
                                               bool call_geno)
{
    Eigen::Map<const ArrDouble2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, F, true);
    ArrDouble2D logLikeForwardInd(C2, M); // log likelihood of forward recursion for ind i
    ArrDouble2D logLikeBackwardInd(C2, M); // log likelihood of backward recursion for ind i

    // ======== forward recursion ===========
    int k1, k2, k12, k12s, k12e;
    double protect_me = 0; // for underflow protection
    // first site
    int s{0};
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            logLikeForwardInd(k12, s) = emitDip(s, k12) + PI(s, k1) * PI(s, k2);
        }
    }

    // do forward recursion
    for(s = 1; s < M; s++)
    {
        // ArrDouble1D likeForwardTmp;
        auto likeForwardTmp = Eigen::exp(logLikeForwardInd.col(s - 1) - protect_me);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                k12s = (k1 * C + k2) * C2; // fist index
                k12e = k12s + C2 - 1; // last index
                logLikeForwardInd(k12, s) =
                    protect_me + emitDip(s, k12)
                    + log((likeForwardTmp * transDip(s, Eigen::seq(k12s, k12e)).transpose()).sum());
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
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                logLikeBackwardInd(k12, s) =
                    protect_me + log((likeBackwardTmp * transDip(s + 1, Eigen::seqN(k12, C2, C2)).transpose()).sum());
            }
        }
        // protect_me = logLikeBackwardInd.row(s).maxCoeff();
        protect_me = logLikeBackwardInd(C2 - 1, s);
    }

    std::lock_guard<std::mutex> lock(mutex_it); // lock here if RAM cost really matters

    // ======== post decoding get p(Z|X, G) ===========
    // ArrDouble2D indPostProbsZ; // post probabilities for ind i, M x (C x C)
    ArrDouble2D indPostProbsZ = (logLikeBackwardInd + logLikeForwardInd - indLogLikeForwardAll).exp().transpose();
    // ======== post decoding get p(Z,G|X, theta) ===========
    // ArrDouble2D indPostProbsZandG, M x (C x C x 2 x 2)
    ArrDouble2D indPostProbsZandG(M, C2 * 4);
    ArrDouble1D tmpSum(M);
    int g1, g2, g12;
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            tmpSum.setZero();
            k12 = k1 * C + k2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    indPostProbsZandG.col(k12 * 4 + g12) = gli.row(g1 + g2).transpose()
                                                           * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                                           * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                    tmpSum += indPostProbsZandG.col(k12 * 4 + g12);
                }
            }
            indPostProbsZandG.middleCols(k12 * 4, 4).colwise() *= indPostProbsZ.col(k12);
            indPostProbsZandG.middleCols(k12 * 4, 4).colwise() /= tmpSum;
        }
    }

    if(call_geno)
    {
        // std::lock_guard<std::mutex> lock(mutex_it);
        GP.col(ind) = callGenotypeInd(indPostProbsZandG);
        // output likelihood of each cluster
        ArrDouble2D likeCluster = (logLikeBackwardInd + logLikeForwardInd).exp().transpose();
        ofs.write((char *)likeCluster.data(), M * C2 * 8);
    }
    postProbsZ += indPostProbsZ;
    postProbsZandG += indPostProbsZandG;
    // std::cout << std::this_thread::get_id() << ": " << ind << '\n';

    return indLogLikeForwardAll;
}

inline ArrDouble2D FastPhaseK4::callGenotypeInd(const ArrDouble2D & indPostProbsZandG)
{
    ArrDouble1D geno = ArrDouble1D::Zero(M * 3);
    int k1, k2, k12, g1, g2, g12, g3;
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    g3 = g1 + g2;
                    geno(Eigen::seqN(g3, M, 3)) += indPostProbsZandG.col(k12 * 4 + g12);
                }
            }
        }
    }
    return geno;
}

/*
** @param distRate distance or recombination rate between two markers
*/
inline void FastPhaseK4::transitionCurIter(const ArrDouble1D & distRate)
{
    int to1, to2, from1, from2;
    int k4d, k2d1, k2d2;
    for(from1 = 0; from1 < C; from1++)
    {
        for(to1 = 0; to1 < C; to1++)
        {
            k2d1 = from1 * C + to1;
            if(from1 == to1)
                transHap.col(k2d1) = Eigen::exp(-distRate) + (1 - Eigen::exp(-distRate)) * PI.col(to1);
            else
                transHap.col(k2d1) = (1 - Eigen::exp(-distRate)) * PI.col(to1);
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
                    k2d1 = from1 * C + to1;
                    k2d2 = from2 * C + to2;
                    transDip.col(k4d) = transHap.col(k2d1) * transHap.col(k2d2);
                }
            }
        }
    }
}

inline void FastPhaseK4::updateClusterFreqPI(const ArrDouble2D & postProbsZ, double tol)
{
    int k1, k2, k12;
    PI.setZero();
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            PI.col(k1) += postProbsZ.col(k12);
            PI.col(k2) += postProbsZ.col(k12);
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

inline void FastPhaseK4::updateAlleleFreqWithinCluster(const ArrDouble2D & postProbsZandG, double tol)
{
    int k1, k2, k12, g1, g2, g12;
    ArrDouble2D Ekg = ArrDouble2D::Zero(M, C * 2);
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    Ekg.col(k1 * 2 + g1) += postProbsZandG.col(k12 * 4 + g12);
                    Ekg.col(k2 * 2 + g2) += postProbsZandG.col(k12 * 4 + g12);
                }
            }
        }
    }
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2))
        / (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
    // std::cout << F.rows() << "," << F.cols() << std::endl;
    // std::cout << F.topRows(3) << std::endl;
    // map to domain but normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
    // std::cout << "updateAlleleFreqWithinCluster" << std::endl;
}

#endif // FASTPHASE_H_
