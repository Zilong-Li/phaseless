#ifndef FASTPHASE_H_
#define FASTPHASE_H_

#include <Eigen/Dense>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>


using MyMat2D = Eigen::MatrixXf;     // use MatrixXf if no accuracy drop
using MyArr2D = Eigen::ArrayXXf;     // use ArrayXf if no accuracy drop
using MatDouble2D = Eigen::MatrixXd; // use matrix for linear algebra operations
using ArrDouble2D = Eigen::ArrayXXd; // use array for element-wise operations
using ArrDouble1D = Eigen::ArrayXd;
using ArrFloat2D = Eigen::ArrayXXf;
using ArrFloat1D = Eigen::ArrayXf;
using DoubleVec1D = std::vector<double>;
using FloatVec1D = std::vector<float>;

template <typename MatrixType, typename RandomEngineType>
inline MatrixType RandomUniform(const Eigen::Index numRows, const Eigen::Index numCols,
                                RandomEngineType& engine)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{
        0.0, 1.0}; // or using 0.05, 0.95
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
};


class fastPhaseK4
{
private:
    std::mutex mutex_it; // in case of race condition

public:
    fastPhaseK4(int n, int m, int c, int seed);
    ~fastPhaseK4();

    // SHARED VARIBALES
    const int N, M, C;
    ArrDouble2D GP;       // genotype probabilies for all individuals, N x (M x 3)
    ArrDouble2D PI;       // nsnps x C
    ArrDouble2D F;        // nsnps x C
    ArrDouble2D transHap; // nsnps x (C x C)
    ArrDouble2D transDip; // nsnps x (C x C x C x C)

    ArrDouble2D emissionCurIterInd(const ArrDouble2D& gli, bool use_log = true);
    void transitionCurIter(const ArrDouble1D& distRate);
    void updateClusterFreqPI(const ArrDouble2D& postProbsZ, double tol);
    void updateAlleleFreqWithinCluster(const ArrDouble2D& postProbsZandG, double tol);
    double forwardAndBackwards(int ind, const DoubleVec1D& GL, ArrDouble2D& postProbsZ,
                               ArrDouble2D& postProbsZandG, bool call_geno = false);
    ArrDouble2D callGenotypeInd(const ArrDouble2D& indPostProbsZandG);
};

inline fastPhaseK4::fastPhaseK4(int n, int m, int c, int seed = 1)
    : N(n), M(m), C(c), GP(M * 3, N), transHap(M, C * C), transDip(M, C * C * C * C)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<ArrDouble2D, std::default_random_engine>(M, C, rng);
    PI = RandomUniform<ArrDouble2D, std::default_random_engine>(M, C, rng);
    PI = PI.colwise() / PI.rowwise().sum(); // normalize it
}

inline fastPhaseK4::~fastPhaseK4()
{
}


inline void fastPhaseK4::updateClusterFreqPI(const ArrDouble2D& postProbsZ, double tol)
{
    int k1, k2, k12;
    PI.setZero();
    for (k1 = 0; k1 < C; k1++)
    {
        for (k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            PI.col(k1) += postProbsZ.col(k12);
            PI.col(k2) += postProbsZ.col(k12);
        }
    }
    PI /= 2 * N;
    // map to domain
    if (PI.isNaN().any())
        throw std::runtime_error("NaN in PI\n");
    PI = (PI < tol).select(tol, PI);         // lower bound
    PI = (PI > 1 - tol).select(1 - tol, PI); // upper bound
    // normalize it now
    PI = PI.colwise() / PI.rowwise().sum();
}

inline void fastPhaseK4::updateAlleleFreqWithinCluster(const ArrDouble2D& postProbsZandG, double tol)
{
    int k1, k2, k12, g1, g2, g12;
    ArrDouble2D Ekg = ArrDouble2D::Zero(M, C * 2);
    for (k1 = 0; k1 < C; k1++)
    {
        for (k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            for (g1 = 0; g1 < 2; g1++)
            {
                for (g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    Ekg.col(k1 * 2 + g1) += postProbsZandG.col(k12 * 4 + g12);
                    Ekg.col(k2 * 2 + g2) += postProbsZandG.col(k12 * 4 + g12);
                }
            }
        }
    }
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) /
        (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
    // std::cout << F.rows() << "," << F.cols() << std::endl;
    // std::cout << F.topRows(3) << std::endl;
    // map to domain but normalization
    if (F.isNaN().any())
        throw std::runtime_error("NaN in PI\n");
    F = (F < tol).select(tol, F);         // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
    // std::cout << "updateAlleleFreqWithinCluster" << std::endl;
}


/*
** @param ind current individual i
** @param GL  genotype likelihood of all individuals in snp major form
** @return individual total likelihood
*/
inline double fastPhaseK4::forwardAndBackwards(int ind, const DoubleVec1D& GL, ArrDouble2D& postProbsZ,
                                               ArrDouble2D& postProbsZandG, bool call_geno)
{
    Eigen::Map<const ArrDouble2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli);
    ArrDouble2D logLikeForwardInd(M, C * C);  // log likelihood of forward recursion for ind i
    ArrDouble2D logLikeBackwardInd(M, C * C); // log likelihood of backward recursion for ind i

    // ======== forward recursion ===========
    int k1, k2, k12, k12s, k12e;
    double protect_me = 0; // for underflow protection
    // first site
    int s{0};
    for (k1 = 0; k1 < C; k1++)
    {
        for (k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            logLikeForwardInd(s, k12) = emitDip(s, k12) + PI(s, k1) * PI(s, k2);
        }
    }

    // do forward recursion
    for (s = 1; s < M; s++)
    {
        // ArrDouble1D likeForwardTmp;
        auto likeForwardTmp = Eigen::exp(logLikeForwardInd.row(s - 1) - protect_me);
        for (k1 = 0; k1 < C; k1++)
        {
            for (k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                k12s = (k1 * C + k2) * C * C; // fist index
                k12e = k12s + C * C - 1;      // last index
                logLikeForwardInd(s, k12) = protect_me + emitDip(s, k12) +
                                            log((likeForwardTmp * transDip(s, Eigen::seq(k12s, k12e))).sum());
            }
        }
        // protect_me = logLikeForwardInd.row(s).maxCoeff();
        protect_me = logLikeForwardInd(s, C * C - 1);
    }
    // total likelhoods of the individual
    double indLikeForwardAll = protect_me + log((logLikeForwardInd.row(M - 1) - protect_me).exp().sum());

    // ======== backward recursion ===========
    // set last site
    logLikeBackwardInd.row(M - 1).setZero();
    // site M-2 to 0
    protect_me = 0;
    for (s = M - 2; s >= 0; s--)
    {
        auto likeBackwardTmp = Eigen::exp(emitDip.row(s + 1) + logLikeBackwardInd.row(s + 1) - protect_me);
        for (k1 = 0; k1 < C; k1++)
        {
            for (k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                logLikeBackwardInd(s, k12) =
                    protect_me +
                    log((likeBackwardTmp * transDip(s + 1, Eigen::seqN(k12, C * C, C * C))).sum());
            }
        }
        // protect_me = logLikeBackwardInd.row(s).maxCoeff();
        protect_me = logLikeBackwardInd(s, C * C - 1);
    }

    std::lock_guard<std::mutex> lock(mutex_it); // lock here if RAM cost really matters

    // ======== post decoding get p(Z|X, G) ===========
    // ArrDouble2D indPostProbsZ; // post probabilities for ind i, M x (C x C)
    ArrDouble2D indPostProbsZ = (logLikeBackwardInd + logLikeForwardInd - indLikeForwardAll).exp();

    // ======== post decoding get p(Z,G|X, theta) ===========
    // ArrDouble2D indPostProbsZandG, M x (C x C x 2 x 2)
    ArrDouble2D indPostProbsZandG(M, C * C * 4);
    ArrDouble1D tmpSum(M);
    int g1, g2, g12;
    for (k1 = 0; k1 < C; k1++)
    {
        for (k2 = 0; k2 < C; k2++)
        {
            tmpSum.setZero();
            k12 = k1 * C + k2;
            for (g1 = 0; g1 < 2; g1++)
            {
                for (g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    indPostProbsZandG.col(k12 * 4 + g12) = gli.row(g1 + g2).transpose() *
                                                           (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1))) *
                                                           (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                    tmpSum += indPostProbsZandG.col(k12 * 4 + g12);
                }
            }
            indPostProbsZandG.middleCols(k12 * 4, 4).colwise() *= indPostProbsZ.col(k12);
            indPostProbsZandG.middleCols(k12 * 4, 4).colwise() /= tmpSum;
        }
    }

    if (call_geno)
    {
        // std::lock_guard<std::mutex> lock(mutex_it);
        GP.col(ind) = callGenotypeInd(indPostProbsZandG);
    }
    postProbsZ += indPostProbsZ;
    postProbsZandG += indPostProbsZandG;
    // std::cout << std::this_thread::get_id() << ": " << ind << '\n';

    return indLikeForwardAll;
}

inline ArrDouble2D fastPhaseK4::callGenotypeInd(const ArrDouble2D& indPostProbsZandG)
{
    ArrDouble1D geno = ArrDouble1D::Zero(M * 3);
    int k1, k2, k12, g1, g2, g12, g3;
    for (k1 = 0; k1 < C; k1++)
    {
        for (k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            for (g1 = 0; g1 < 2; g1++)
            {
                for (g2 = 0; g2 < 2; g2++)
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
inline void fastPhaseK4::transitionCurIter(const ArrDouble1D& distRate)
{
    int to1, to2, from1, from2;
    int k4d, k2d1, k2d2;
    for (from1 = 0; from1 < C; from1++)
    {
        for (to1 = 0; to1 < C; to1++)
        {
            k2d1 = from1 * C + to1;
            if (from1 == to1)
                transHap.col(k2d1) = Eigen::exp(-distRate) + (1 - Eigen::exp(-distRate)) * PI.col(to1);
            else
                transHap.col(k2d1) = (1 - Eigen::exp(-distRate)) * PI.col(to1);
        }
    }

    for (from1 = 0; from1 < C; from1++)
    {
        for (from2 = 0; from2 < C; from2++)
        {
            for (to1 = 0; to1 < C; to1++)
            {
                for (to2 = 0; to2 < C; to2++)
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

/*
** @param gli  genotype likelihoods of current individual i,3 x nsnps
*/
inline ArrDouble2D fastPhaseK4::emissionCurIterInd(const ArrDouble2D& gli, bool use_log)
{
    int k1, k2, g1, g2;
    ArrDouble2D emitDip(M, C * C); // emission probabilies, nsnps x (C x C)
    for (k1 = 0; k1 < C; k1++)
    {
        for (k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for (g1 = 0; g1 <= 1; g1++)
            {
                for (g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.row(g1 + g2).transpose() *
                                                (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1))) *
                                                (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                }
            }
        }
    }
    if (use_log)
        return emitDip.log();
    else
        return emitDip;
}

#endif // FASTPHASE_H_
