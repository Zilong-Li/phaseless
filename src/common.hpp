/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/common.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#ifndef COMMON_H_
#define COMMON_H_

#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <vector>

using MatDouble2D = Eigen::MatrixXd; // use matrix for linear algebra operations
using ArrDouble2D = Eigen::ArrayXXd; // use array for element-wise operations
using ArrDouble1D = Eigen::ArrayXd;
using ArrFloat2D = Eigen::ArrayXXf;
using ArrFloat1D = Eigen::ArrayXf;
using IntVec1D = std::vector<int>;
using IntVec2D = std::vector<IntVec1D>;
using FloatVec1D = std::vector<float>;
using FloatVec2D = std::vector<FloatVec1D>;
using DoubleVec1D = std::vector<double>;
using DoubleVec2D = std::vector<DoubleVec1D>;
using StringVec1D = std::vector<std::string>;
using MapStringInt1D = std::map<std::string, IntVec1D>;

// using MyMat2D = Eigen::MatrixXd; // use MatrixXf if no accuracy drop
// using MyMat1D = Eigen::VectorXd; // use MatrixXf if no accuracy drop
// using MyArr2D = Eigen::ArrayXXd; // use ArrayXf if no accuracy drop
// using MyArr1D = Eigen::ArrayXd; // use ArrayXf if no accuracy drop
using MyMat2D = Eigen::MatrixXf;
using MyMat1D = Eigen::VectorXf;
using MyArr2D = Eigen::ArrayXXf;
using MyArr1D = Eigen::ArrayXf;
using MyFloat1D = FloatVec1D;
using MyFloat2D = FloatVec2D;

template<typename MatrixType, typename RandomEngineType>
inline MatrixType RandomUniform(const Eigen::Index numRows,
                                const Eigen::Index numCols,
                                RandomEngineType & engine,
                                typename MatrixType::Scalar a,
                                typename MatrixType::Scalar b)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{
        a, b}; // or using 0.05, 0.95
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
};

// all the genome info I need from fastphase
struct BigAss
{
    IntVec2D pos; // store position of markers of each chunk
    MyFloat2D gls; // store gl of each chunk
    StringVec1D sampleids, chrs;
    int chunksize, nsamples, nsnps, nchunks;
    MyFloat2D PI, F, transRate;
    int C; // fastphase pars
};

// check initialize_sigmaCurrent_m in STITCH
inline auto calc_distRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    MyArr1D distRate(markers.size());
    distRate(0) = 1e20;
    // int nGen = 4 * Ne / C;
    for(size_t i = 1; i < markers.size(); i++) distRate(i) = (markers[i] - markers[i - 1]) / 1e6;
    // distRate(i) = (markers[i] - markers[i - 1]) * nGen * expRate / 1e8;
    return distRate;
}

inline auto calc_transRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    MyArr1D distRate(markers.size());
    distRate(0) = exp(-1e20);
    // int nGen = 4 * Ne / C;
    for(size_t i = 1; i < markers.size(); i++) distRate(i) = exp(-(markers[i] - markers[i - 1]) / 1e6);
    MyArr2D transRate(3, markers.size());
    transRate.row(0) = distRate.square();
    transRate.row(1) = distRate * (1 - distRate);
    transRate.row(2) = (1 - distRate).square();
    return transRate;
}

/*
** @param gli  genotype likelihoods of current individual i,3 x nsnps
*/
inline auto emissionCurIterInd(const MyArr2D & gli, const MyArr2D & F, bool use_log)
{
    int k1, k2, g1, g2;
    const int M = F.rows();
    const int C = F.cols();
    MyArr2D emitDip(M, C * C); // emission probabilies, nsnps x (C x C)
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
        emitDip = emitDip.colwise() / emitDip.rowwise().maxCoeff(); // normalize it
        const double maxEmissionMatrixDifference = 1e-10;
        emitDip = (emitDip < maxEmissionMatrixDifference).select(maxEmissionMatrixDifference, emitDip);
    }
    return emitDip;
}

inline auto getClusterLikelihoods(int ind,
                                  const MyFloat1D & GL,
                                  const MyArr2D & transRate,
                                  const MyArr2D & PI,
                                  const MyArr2D & F)
{
    int k1, k2, k12;
    const int M = F.rows();
    const int C = F.cols();
    // ======== forward and backward recursion ===========
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, F, false);
    MyArr2D LikeForwardInd(C * C, M); // likelihood of forward recursion for ind i, not log
    MyArr2D LikeBackwardInd(C * C, M); // likelihood of backward recursion for ind i, not log
    MyArr1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    MyArr1D cs(M);
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
    MyArr2D icluster = LikeForwardInd * LikeBackwardInd; // C x C x M
    return icluster;
}

inline auto getClusterLikelihoods(int ind,
                                  int M,
                                  int C,
                                  const MyFloat1D & GL,
                                  const MyFloat1D & transRate_,
                                  const MyFloat1D & PI_,
                                  const MyFloat1D & F_)
{
    int g1, g2, k1, k2, k12;
    int C2 = C * C;
    int igs = ind * M * 3;
    const double maxEmission = 1e-10;
    // ======== forward and backward recursion ===========
    MyArr2D emitDip(C2, M);
    MyArr2D LikeForwardInd(C2, M); // likelihood of forward recursion for ind i, not log
    MyArr2D LikeBackwardInd(C2, M); // likelihood of backward recursion for ind i, not log
    MyArr1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    MyArr1D cs(M);
    double constTmp;

    // ======== forward recursion ===========
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
                    emitDip(k12, s) += GL[igs + s * 3 + g1 + g2]
                                       * (g1 * F_[k1 * M + s] + (1 - g1) * (1 - F_[k1 * M + s]))
                                       * (g2 * F_[k2 * M + s] + (1 - g2) * (1 - F_[k2 * M + s]));
                }
            }
            if(emitDip(k12, s) < maxEmission) emitDip(k12, s) = maxEmission;
            LikeForwardInd(k12, s) = emitDip(k12, s) * PI_[k1 * M + s] * PI_[k2 * M + s];
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
                sumTmp1(k1) += LikeForwardInd(k12, s - 1) * transRate_[s * 3 + 1];
                sumTmp2(k2) += LikeForwardInd(k12, s - 1) * transRate_[s * 3 + 1];
            }
        }
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate_[s * 3 + 2];
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
                        emitDip(k12, s) += GL[igs + s * 3 + g1 + g2]
                                           * (g1 * F_[k1 * M + s] + (1 - g1) * (1 - F_[k1 * M + s]))
                                           * (g2 * F_[k2 * M + s] + (1 - g2) * (1 - F_[k2 * M + s]));
                    }
                }
                if(emitDip(k12, s) < maxEmission) emitDip(k12, s) = maxEmission;
                LikeForwardInd(k12, s) =
                    emitDip(k12, s)
                    * (LikeForwardInd(k12, s - 1) * transRate_[s * 3 + 0] + PI_[k1 * M + s] * sumTmp1(k2)
                       + PI_[k2 * M + s] * sumTmp2(k1) + PI_[k1 * M + s] * PI_[k2 * M + s] * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    // ======== backward recursion ===========
    s = M - 1;
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
                sumTmp1(k1) += beta_mult_emit(k12) * PI_[k2 * M + s + 1] * transRate_[(s + 1) * 3 + 1];
                sumTmp2(k2) += beta_mult_emit(k12) * PI_[k1 * M + s + 1] * transRate_[(s + 1) * 3 + 1];
                constTmp += beta_mult_emit(k12) * PI_[k1 * M + s + 1] * PI_[k2 * M + s + 1]
                            * transRate_[(s + 1) * 3 + 2];
            }
        }
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeBackwardInd(k12, s) =
                    (beta_mult_emit(k12) * transRate_[(s + 1) * 3 + 0] + sumTmp1(k1) + sumTmp2(k2) + constTmp)
                    * cs(s);
            }
        }
    }
    MyArr2D icluster = LikeForwardInd * LikeBackwardInd; // C x C x M
    return icluster;
}

#endif // COMMON_H_
