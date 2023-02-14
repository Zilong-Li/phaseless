/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/fastphase.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef FASTPHASE_H_
#define FASTPHASE_H_

#include "common.hpp"
#include <fstream>
#include <mutex>

class FastPhaseK2
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK2(int n, int m, int c, int seed);
    ~FastPhaseK2();

    // SHARED VARIBALES
    std::ofstream ofs;
    const int N, M, C, C2; // C2 = C x C
    MyArr2D GP; // genotype probabilies for all individuals, N x (M x 3)
    MyArr2D PI, Ek; // nsnps x C
    MyArr2D F; // M x C
    MyArr2D Ekg; // M x C x 2

    void openClusterFile(std::string out);
    void initIteration(double tol = 1e-6);
    void updateIteration();
    double forwardAndBackwards(int ind, const MyFloat1D & GL, const MyArr2D & transRate, bool call_geno);
};

inline FastPhaseK2::FastPhaseK2(int n, int m, int c, int seed)
: N(n), M(m), C(c), C2(c * c), GP(M * 3, N), Ek(M, C), Ekg(M, C * 2)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
}

inline FastPhaseK2::~FastPhaseK2()
{
    if(ofs.is_open()) ofs.close();
}

inline void FastPhaseK2::openClusterFile(std::string out)
{
    ofs.open(out, std::ios::binary);
    if(!ofs) throw std::runtime_error(out + ": " + strerror(errno));
    ofs.write((char *)&N, 4);
    ofs.write((char *)&M, 4);
    ofs.write((char *)&C, 4);
}

inline void FastPhaseK2::initIteration(double tol)
{
    // map PI to domain with normalization
    if(PI.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    PI = (PI < tol).select(tol, PI); // lower bound
    PI = (PI > 1 - tol).select(1 - tol, PI); // upper bound
    // normalize it now
    PI = PI.colwise() / PI.rowwise().sum();
    // map F to domain but no normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in F\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
    Ek.setZero();
    Ekg.setZero();
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @return individual total likelihood
*/
inline double FastPhaseK2::forwardAndBackwards(int ind,
                                               const MyFloat1D & GL,
                                               const MyArr2D & transRate,
                                               bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, 3, M);
    const double maxEmission = 1e-10;
    MyArr2D emitDip(C2, M);
    MyArr2D LikeForwardInd(C2, M); // likelihood of forward recursion for ind i, not log
    MyArr2D LikeBackwardInd(C2, M); // likelihood of backward recursion for ind i, not log
    MyArr1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    MyArr1D cs(M);
    double constTmp;

    // ======== forward recursion ===========
    int g1, g2, g3, g12, k1, k2, k12;
    int s{0};
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            emitDip(k12, s) = 0;
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip(k12, s) += gli(g1 + g1, s) * (g1 * F(s, k1) + (1 - g1) * (1 - F(s, k1)))
                                       * (g2 * F(s, k2) + (1 - g2) * (1 - F(s, k2)));
                }
            }
            if(emitDip(k12, s) < maxEmission) emitDip(k12, s) = maxEmission;
            LikeForwardInd(k12, s) = emitDip(k12, s) * PI(s, k1) * PI(s, k2);
        }
    }
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1.setZero();
        sumTmp2.setZero();
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += LikeForwardInd(k12, s - 1) * transRate(1, s);
                sumTmp2(k2) += LikeForwardInd(k12, s - 1) * transRate(1, s);
            }
        }
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                emitDip(k12, s) = 0;
                for(g1 = 0; g1 <= 1; g1++)
                {
                    for(g2 = 0; g2 <= 1; g2++)
                    {
                        emitDip(k12, s) += gli(g1 + g1, s) * (g1 * F(s, k1) + (1 - g1) * (1 - F(s, k1)))
                                           * (g2 * F(s, k2) + (1 - g2) * (1 - F(s, k2)));
                    }
                }
                if(emitDip(k12, s) < maxEmission) emitDip(k12, s) = maxEmission;
                LikeForwardInd(k12, s) =
                    emitDip(k12, s)
                    * (LikeForwardInd(k12, s - 1) * transRate(0, s) + PI(s, k1) * sumTmp1(k2)
                       + PI(s, k2) * sumTmp2(k1) + PI(s, k1) * PI(s, k2) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    // total likelhoods of the individual
    double indLogLikeForwardAll = log((LikeForwardInd.col(M - 1).sum() / cs).sum());

    // ======== backward recursion ===========
    s = M - 1; // set last site
    LikeBackwardInd.col(s).setConstant(cs(s)); // not log scale
    for(s = M - 2; s >= 0; s--)
    {
        sumTmp1.setZero();
        sumTmp2.setZero();
        constTmp = 0;
        auto beta_mult_emit = emitDip.col(s + 1) * LikeBackwardInd.col(s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI(s + 1, k2) * transRate(1, s + 1);
                sumTmp2(k2) += beta_mult_emit(k12) * PI(s + 1, k1) * transRate(1, s + 1);
                constTmp += beta_mult_emit(k12) * PI(s + 1, k1) * PI(s + 1, k2) * transRate(2, s + 1);
            }
        }
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                // apply scaling
                LikeBackwardInd(k12, s) =
                    (beta_mult_emit(k12) * transRate(0, s + 1) + sumTmp1(k1) + sumTmp2(k2) + constTmp)
                    * cs(s);
            }
        }
    }
    cs *= LikeForwardInd.col(M - 1).sum(); // get last forward likelihood back
    MyArr1D ind_post_z_col(M); // col of indPostProbsZ
    MyArr2D ind_post_z_g(M, 4); // cols of indPostProbsZandG
    MyArr1D tmpSum(M);
    MyArr1D geno;
    if(call_geno) geno.setZero(M * 3);
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
                    ind_post_z_g.col(g12) = gli.row(g1 + g2).transpose()
                                            * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
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
        if(ofs.is_open()) // output likelihood of each cluster
        {
            MyArr2D likeCluster = (LikeBackwardInd * LikeForwardInd).transpose();
            ofs.write((char *)likeCluster.data(), M * C2 * 8);
        }
    }

    return indLogLikeForwardAll;
}

inline void FastPhaseK2::updateIteration()
{
    PI = Ek / (2 * N);
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2))
        / (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
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

    void updateClusterFreqPI(const MyArr2D & postProbsZ, double tol);
    void updateAlleleFreqWithinCluster(const MyArr2D & postProbsZandG, double tol);
    void transitionCurIter(const MyArr1D & distRate);
    auto callGenotypeInd(const MyArr2D & indPostProbsZandG);
    void openClusterFile(std::string out);
    double forwardAndBackwards(int ind,
                               const MyFloat1D & GL,
                               MyArr2D & postProbsZ,
                               MyArr2D & postProbsZandG,
                               bool call_geno = false);
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
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, F, true);
    MyArr2D logLikeForwardInd(C2, M); // log likelihood of forward recursion for ind i
    MyArr2D logLikeBackwardInd(C2, M); // log likelihood of backward recursion for ind i

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
        // MyArr1D likeForwardTmp;
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
        auto likeBackwardTmp =
            Eigen::exp(emitDip.row(s + 1).transpose() + logLikeBackwardInd.col(s + 1) - protect_me);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                logLikeBackwardInd(k12, s) =
                    protect_me
                    + log((likeBackwardTmp * transDip(s + 1, Eigen::seqN(k12, C2, C2)).transpose()).sum());
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

inline void FastPhaseK4::updateClusterFreqPI(const MyArr2D & postProbsZ, double tol)
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

inline void FastPhaseK4::updateAlleleFreqWithinCluster(const MyArr2D & postProbsZandG, double tol)
{
    int k1, k2, k12, g1, g2, g12;
    MyArr2D Ekg = MyArr2D::Zero(M, C * 2);
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
