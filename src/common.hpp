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

using IntVec1D = std::vector<int>;
using IntVec2D = std::vector<IntVec1D>;
using FloatVec1D = std::vector<float>;
using FloatVec2D = std::vector<FloatVec1D>;
using DoubleVec1D = std::vector<double>;
using DoubleVec2D = std::vector<DoubleVec1D>;
using StringVec1D = std::vector<std::string>;
using MapStringInt1D = std::map<std::string, IntVec1D>;

using MyFloat = float; // use float if no accuracy drops
using MyFloat1D = std::vector<MyFloat>;
using MyFloat2D = std::vector<MyFloat1D>;
using MyMat2D = Eigen::Matrix<MyFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using MyMat1D = Eigen::Matrix<MyFloat, Eigen::Dynamic, 1, Eigen::ColMajor>;
using MyArr2D = Eigen::Array<MyFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using MyArr1D = Eigen::Array<MyFloat, Eigen::Dynamic, 1, Eigen::ColMajor>;

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
    int chunksize, nsamples, nsnps, nchunks, C; // number of clusters
    MyFloat2D PI, F, transRate; // M x C, 3 x M, fastphase pars
    IntVec1D ends; // chunk index where each chromo ends
    StringVec1D sampleids, chrs;
    IntVec2D pos; // store position of markers of each chunk
    MyFloat2D gls; // store gl(M, 3) of each chunk
};

inline auto calc_distRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    MyArr1D distRate(markers.size());
    // int nGen = 4 * Ne / C;
    // distRate(i) = (markers[i] - markers[i - 1]) * nGen * expRate / 1e8;
    distRate(0) = exp(-1e20);
    for(size_t i = 1; i < markers.size(); i++) distRate(i) = exp(-(markers[i] - markers[i - 1]) / 1e5);
    return distRate;
}

// check initialize_sigmaCurrent_m in STITCH
inline auto calc_transRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    MyArr1D distRate = calc_distRate(markers, C, Ne, expRate);
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
                    emitDip.col(k1 * C + k2) += gli.col(g1 + g2)
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
        // emitDip = emitDip.colwise() / emitDip.rowwise().maxCoeff(); // normalize it
        const double minEmission = 1e-10;
        emitDip = (emitDip < minEmission).select(minEmission, emitDip);
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
    MyArr2D emitDip = emissionCurIterInd(gli, F, false).transpose();
    MyArr2D LikeForwardInd(C * C, M); // likelihood of forward recursion for ind i, not log
    MyArr2D LikeBackwardInd(C * C, M); // likelihood of backward recursion for ind i, not log
    MyArr1D sumTmp1(C); // store sum over internal loop
    MyArr1D cs(M);
    double constTmp;

    // ======== forward recursion ===========
    int s{0};
    LikeForwardInd.col(s) =
        emitDip.col(s) * (PI.col(s).matrix() * PI.col(s).transpose().matrix()).reshaped().array();
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1 = LikeForwardInd.col(s - 1).reshaped(C, C).rowwise().sum() * transRate(1, s);
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeForwardInd(k12, s) =
                    emitDip(k12, s)
                    * (LikeForwardInd(k12, s - 1) * transRate(0, s) + PI(k1, s) * sumTmp1(k2)
                       + PI(k2, s) * sumTmp1(k1) + PI(k1, s) * PI(k2, s) * constTmp);
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
        constTmp = 0;
        auto beta_mult_emit = emitDip.row(s + 1).transpose() * LikeBackwardInd.col(s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI(k2, s + 1) * transRate(1, s + 1);
                constTmp += beta_mult_emit(k12) * PI(k1, s + 1) * PI(k2, s + 1) * transRate(2, s + 1);
            }
        }
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                // apply scaling
                LikeBackwardInd(k12, s) =
                    (beta_mult_emit(k12) * transRate(0, s + 1) + sumTmp1(k1) + sumTmp1(k2) + constTmp)
                    * cs(s);
            }
        }
    }
    MyArr2D icluster = LikeForwardInd * LikeBackwardInd; // C x C x M
    return icluster;
}

inline auto getClusterLikelihoods(int ind,
                                  MyArr2D & LikeForwardInd,
                                  MyArr2D & LikeBackwardInd,
                                  const MyFloat1D & GL,
                                  const MyFloat1D & transRate_,
                                  const MyFloat1D & PI_,
                                  const MyFloat1D & F_,
                                  bool gamma = false)
{
    const int C2 = LikeForwardInd.rows();
    const int M = LikeForwardInd.cols();
    int C = F_.size() / M;
    int g1, g2, k1, k2, k12;
    int igs = ind * M * 3;
    // ======== forward and backward recursion ===========
    MyArr2D emitDip(C2, M);
    MyArr1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    MyArr1D cs = MyArr1D::Zero(M);
    double constTmp;
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
                    emitDip(k12, s) += GL[igs + (g1 + g2) * M + s]
                                       * (g1 * F_[k1 * M + s] + (1 - g1) * (1 - F_[k1 * M + s]))
                                       * (g2 * F_[k2 * M + s] + (1 - g2) * (1 - F_[k2 * M + s]));
                }
            }
            LikeForwardInd(k12, s) = emitDip(k12, s) * PI_[s * C + k1] * PI_[s * C + k2];
            cs(s) += LikeForwardInd(k12, s);
        }
    }
    cs(s) = 1 / cs(s);
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1 = LikeForwardInd.col(s - 1).reshaped(C, C).rowwise().sum() * transRate_[s * 3 + 1];
        sumTmp2 = LikeForwardInd.col(s - 1).reshaped(C, C).colwise().sum() * transRate_[s * 3 + 1];
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
                        emitDip(k12, s) += GL[igs + (g1 + g2) * M + s]
                                           * (g1 * F_[k1 * M + s] + (1 - g1) * (1 - F_[k1 * M + s]))
                                           * (g2 * F_[k2 * M + s] + (1 - g2) * (1 - F_[k2 * M + s]));
                    }
                }
                LikeForwardInd(k12, s) =
                    emitDip(k12, s)
                    * (LikeForwardInd(k12, s - 1) * transRate_[s * 3 + 0] + PI_[s * C + k1] * sumTmp1(k2)
                       + PI_[s * C + k2] * sumTmp2(k1) + PI_[s * C + k1] * PI_[s * C + k2] * constTmp);
                cs(s) += LikeForwardInd(k12, s);
            }
        }
        cs(s) = 1 / cs(s);
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    // double indLike = LikeForwardInd.col(M - 1).sum(); // just 1
    // ======== backward recursion ===========
    s = M - 1;
    LikeBackwardInd.col(s).setConstant(cs(s)); // not log scale
    for(s = M - 2; s >= 0; s--)
    {
        auto beta_mult_emit = emitDip.col(s + 1) * LikeBackwardInd.col(s + 1);
        sumTmp1.setZero();
        sumTmp2.setZero();
        for(constTmp = 0, k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI_[(s + 1) * C + k2] * transRate_[(s + 1) * 3 + 1];
                sumTmp2(k2) += beta_mult_emit(k12) * PI_[(s + 1) * C + k1] * transRate_[(s + 1) * 3 + 1];
                constTmp += beta_mult_emit(k12) * PI_[(s + 1) * C + k1] * PI_[(s + 1) * C + k2]
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
    if(gamma)
    {
        // for(s = 0; s < M; s++) LikeForwardInd.col(s) /= cs(s);
        LikeForwardInd.rowwise() /= cs.transpose();
    }
}

#endif // COMMON_H_
