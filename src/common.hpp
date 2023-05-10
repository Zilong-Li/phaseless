/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/common.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#ifndef COMMON_H_
#define COMMON_H_

#include <Eigen/Dense>
#include <climits>
#include <clocale>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

using IntVec1D = std::vector<int>;
using IntVec2D = std::vector<IntVec1D>;
using FloatVec1D = std::vector<float>;
using FloatVec2D = std::vector<FloatVec1D>;
using DoubleVec1D = std::vector<double>;
using DoubleVec2D = std::vector<DoubleVec1D>;
using StringVec1D = std::vector<std::string>;
using MapStringInt1D = std::map<std::string, IntVec1D>;
using UMapStringInt = std::unordered_map<std::string, int>;

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

struct Options
{
    int ichunk{0}, chunksize{10000}, K{2}, C{10}, nadmix{1000}, nimpute{40}, nthreads{1}, seed{999};
    double ltol{1e-1}, qtol{1e-6}, info{0};
    bool noaccel{0}, noscreen{0}, single_chunk{0};
    std::filesystem::path out, in_beagle, in_vcf, in_bin;
    std::string samples{""}, region{""}, in_plink{""};
    std::string opts_in_effect{"Options in effect:\n   "};
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
** @param gli  genotype likelihoods of current individual i, nsnps x 3
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

inline auto getClusterLikelihoods(MyArr2D & LikeForwardInd,
                                  MyArr2D & LikeBackwardInd,
                                  const MyArr2D & gli,
                                  const MyArr2D & transRate,
                                  const MyArr2D & PI,
                                  const MyArr2D & F,
                                  const bool gamma = true)
{
    const int M = LikeForwardInd.cols();
    const int C = F.cols();
    MyArr2D emitDip = emissionCurIterInd(gli, F, false).transpose();
    MyArr1D sumTmp1(C), cs(M); // store sum over internal loop
    double constTmp;
    // ======== forward recursion ===========
    int z1, z2, z12;
    int s{0};
    LikeForwardInd.col(s) =
        emitDip.col(s) * (PI.col(s).matrix() * PI.col(s).transpose().matrix()).reshaped().array();
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1 = LikeForwardInd.col(s - 1).reshaped(C, C).rowwise().sum() * transRate(1, s);
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                LikeForwardInd(z12, s) =
                    emitDip(z12, s)
                    * (LikeForwardInd(z12, s - 1) * transRate(0, s) + PI(z1, s) * sumTmp1(z2)
                       + PI(z2, s) * sumTmp1(z1) + PI(z1, s) * PI(z2, s) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }

    // ======== backward recursion ===========
    s = M - 1; // set last site
    LikeBackwardInd.col(s).setConstant(cs(s)); // not log scale
    for(s = M - 2; s >= 0; s--)
    {
        auto beta_mult_emit = emitDip.col(s + 1) * LikeBackwardInd.col(s + 1);
        sumTmp1.setZero();
        for(constTmp = 0, z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                sumTmp1(z1) += beta_mult_emit(z12) * PI(z2, s + 1) * transRate(1, s + 1);
                constTmp += beta_mult_emit(z12) * PI(z1, s + 1) * PI(z2, s + 1) * transRate(2, s + 1);
            }
        }
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                // apply scaling
                LikeBackwardInd(z12, s) =
                    (beta_mult_emit(z12) * transRate(0, s + 1) + sumTmp1(z1) + sumTmp1(z2) + constTmp)
                    * cs(s);
            }
        }
    }

    if(gamma) LikeForwardInd.rowwise() /= cs.transpose();
    double indLogLikeForwardAll = log((1 / cs).sum()); // log likelhoods of the individual
    return indLogLikeForwardAll;
}

inline auto getClusterLikelihoods(int ind,
                                  MyArr2D & LikeForwardInd,
                                  MyArr2D & LikeBackwardInd,
                                  const MyFloat1D & GL,
                                  const MyFloat1D & transRate_,
                                  const MyFloat1D & PI,
                                  const MyFloat1D & F,
                                  bool gamma = true)
{
    const int C2 = LikeForwardInd.rows();
    const int M = LikeForwardInd.cols();
    int C = F.size() / M;
    int g1, g2, k1, k2, k12;
    int igs = ind * M * 3;
    // ======== forward and backward recursion ===========
    MyArr2D emitDip(C2, M);
    MyArr1D sumTmp1(C); // store sum over internal loop
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
                                       * (g1 * F[k1 * M + s] + (1 - g1) * (1 - F[k1 * M + s]))
                                       * (g2 * F[k2 * M + s] + (1 - g2) * (1 - F[k2 * M + s]));
                }
            }
            LikeForwardInd(k12, s) = emitDip(k12, s) * PI[s * C + k1] * PI[s * C + k2];
            cs(s) += LikeForwardInd(k12, s);
        }
    }
    cs(s) = 1 / cs(s);
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1 = LikeForwardInd.col(s - 1).reshaped(C, C).rowwise().sum() * transRate_[s * 3 + 1];
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
                                           * (g1 * F[k1 * M + s] + (1 - g1) * (1 - F[k1 * M + s]))
                                           * (g2 * F[k2 * M + s] + (1 - g2) * (1 - F[k2 * M + s]));
                    }
                }
                LikeForwardInd(k12, s) =
                    emitDip(k12, s)
                    * (LikeForwardInd(k12, s - 1) * transRate_[s * 3 + 0] + PI[s * C + k1] * sumTmp1(k2)
                       + PI[s * C + k2] * sumTmp1(k1) + PI[s * C + k1] * PI[s * C + k2] * constTmp);
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
        for(constTmp = 0, k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI[(s + 1) * C + k2] * transRate_[(s + 1) * 3 + 1];
                constTmp += beta_mult_emit(k12) * PI[(s + 1) * C + k1] * PI[(s + 1) * C + k2]
                            * transRate_[(s + 1) * 3 + 2];
            }
        }
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeBackwardInd(k12, s) =
                    (beta_mult_emit(k12) * transRate_[(s + 1) * 3 + 0] + sumTmp1(k1) + sumTmp1(k2) + constTmp)
                    * cs(s);
            }
        }
    }
    // for(s = 0; s < M; s++) LikeForwardInd.col(s) /= cs(s);
    if(gamma) LikeForwardInd.rowwise() /= cs.transpose();
}

inline auto calc_cluster_info(const int N, const MyArr2D & GZP1, const MyArr2D & GZP2)
{
    auto eij = GZP1 + GZP2 * 2;
    auto fij = GZP1 + GZP2 * 4;
    MyArr2D Info = 1 - (fij - eij.square()) / (eij * (1 - eij / (2 * N)));
    Info = (Info < 0).select(0, Info);
    Info = (Info > 1).select(1, Info);
    // Info = Info.isNaN().select(1, Info);
    return Info;
}

#endif // COMMON_H_
