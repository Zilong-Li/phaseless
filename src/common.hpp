#ifndef COMMON_H_
#define COMMON_H_

#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <vector>

// using MyMat2D = Eigen::MatrixXd; // use MatrixXf if no accuracy drop
// using MyMat1D = Eigen::VectorXd; // use MatrixXf if no accuracy drop
// using MyArr2D = Eigen::ArrayXXd; // use ArrayXf if no accuracy drop
// using MyArr1D = Eigen::ArrayXd; // use ArrayXf if no accuracy drop
using MyMat2D = Eigen::MatrixXf;
using MyMat1D = Eigen::VectorXf;
using MyArr2D = Eigen::ArrayXXf;
using MyArr1D = Eigen::ArrayXf;
using MatDouble2D = Eigen::MatrixXd; // use matrix for linear algebra operations
using ArrDouble2D = Eigen::ArrayXXd; // use array for element-wise operations
using ArrDouble1D = Eigen::ArrayXd;
using ArrFloat2D = Eigen::ArrayXXf;
using ArrFloat1D = Eigen::ArrayXf;
using DoubleVec1D = std::vector<double>;
using FloatVec1D = std::vector<float>;
using IntVec1D = std::vector<int>;

template<typename MatrixType, typename RandomEngineType>
inline MatrixType RandomUniform(const Eigen::Index numRows,
                                const Eigen::Index numCols,
                                RandomEngineType & engine,
                                double a,
                                double b)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{
        a, b}; // or using 0.05, 0.95
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
};

// check initialize_sigmaCurrent_m in STITCH
inline auto calc_distRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    ArrDouble1D distRate(markers.size());
    distRate(0) = 1e20;
    // int nGen = 4 * Ne / C;
    for(size_t i = 1; i < markers.size(); i++) distRate(i) = (markers[i] - markers[i - 1]) / 1e6;
    // distRate(i) = (markers[i] - markers[i - 1]) * nGen * expRate / 1e8;
    return distRate;
}

inline auto calc_transRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    ArrDouble1D distRate(markers.size());
    distRate(0) = exp(-1e20);
    // int nGen = 4 * Ne / C;
    for(size_t i = 1; i < markers.size(); i++) distRate(i) = exp(-(markers[i] - markers[i - 1]) / 1e6);
    ArrDouble2D transRate(3, markers.size());
    transRate.row(0) = distRate.square();
    transRate.row(1) = distRate * (1 - distRate);
    transRate.row(2) = (1 - distRate).square();
    return transRate;
}

/*
** @param gli  genotype likelihoods of current individual i,3 x nsnps
*/
inline auto emissionCurIterInd(const ArrDouble2D & gli, const ArrDouble2D & F, bool use_log)
{
    int k1, k2, g1, g2;
    const int M = F.rows();
    const int C = F.cols();
    ArrDouble2D emitDip(M, C * C); // emission probabilies, nsnps x (C x C)
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.row(g1 + g2).transpose()
                                                * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                                * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                }
            }
        }
    }
    if(use_log)
    {
        emitDip = emitDip.log();
    }
    else
    {
        // be careful with underflow
        const double maxEmissionMatrixDifference = 1e-10;
        auto x = emitDip.rowwise().maxCoeff();
        emitDip = emitDip.colwise() / x;
        emitDip = (emitDip < maxEmissionMatrixDifference).select(maxEmissionMatrixDifference, emitDip);
    }
    return emitDip;
}

inline auto getClusterLikelihoods(int ind,
                                  const DoubleVec1D & GL,
                                  const ArrDouble2D & transRate,
                                  const ArrDouble2D & PI,
                                  const ArrDouble2D & F)
{
    int k1, k2, k12;
    const int M = F.rows();
    const int C = F.cols();
    // ======== forward and backward recursion ===========
    Eigen::Map<const ArrDouble2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, F, false);
    ArrDouble2D LikeForwardInd(C * C, M); // likelihood of forward recursion for ind i, not log
    ArrDouble2D LikeBackwardInd(C * C, M); // likelihood of backward recursion for ind i, not log
    ArrDouble1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    ArrDouble1D cs(M);
    double constTmp;

    // ======== forward recursion ===========
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
    for(s = 1; s < M; s++)
    {
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
                LikeForwardInd(k12, s) =
                    emitDip(s, k12)
                    * (LikeForwardInd(k12, s - 1) * transRate(0, s) + PI(s, k1) * sumTmp1(k2)
                       + PI(s, k2) * sumTmp2(k1) + PI(s, k1) * PI(s, k2) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    // ======== backward recursion ===========
    s = M - 1;
    LikeBackwardInd.col(s).setOnes(); // not log scale
    LikeBackwardInd.col(s) *= cs(s);
    for(s = M - 2; s >= 0; s--)
    {
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
    ArrDouble2D icluster = LikeForwardInd * LikeBackwardInd; // C x C x M
    return icluster;
}

#endif // COMMON_H_
