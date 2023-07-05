/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/common.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#ifndef COMMON_H_
#define COMMON_H_

#include "log.hpp"
#include "timer.hpp"
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

// MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
Logger cao; // logger
Timer tim; // Timer
#else
extern Timer tim;
extern Logger cao;
#endif

// STD TYPES
using Int1D = std::vector<int>;
using Int2D = std::vector<Int1D>;
using Float1D = std::vector<float>;
using Float2D = std::vector<Float1D>;
using Double1D = std::vector<double>;
using Double2D = std::vector<Double1D>;
using String1D = std::vector<std::string>;
using MapStringInt1D = std::map<std::string, Int1D>;
using UMapStringInt = std::unordered_map<std::string, int>;

// EIGEN TYPES
using Mat2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Mat1D = Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>;
using Arr2D = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Arr1D = Eigen::Array<double, Eigen::Dynamic, 1, Eigen::ColMajor>;

// MY TYPES
using MyFloat = double; // use float if no accuracy drops
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
    int gridsize{0};
    double ltol{1e-1}, qtol{1e-6}, info{0};
    bool noaccel{0}, noscreen{0}, single_chunk{0}, debug{0};
    std::filesystem::path out, in_beagle, in_vcf, in_bin;
    std::string samples{""}, region{""}, in_plink{""};
    std::string opts_in_effect{"Options in effect:\n   "};
};

// all the genome info I need from fastphase
struct BigAss
{
    int chunksize, nsamples, nsnps, nchunks, C; // number of clusters
    MyFloat2D PI, F, transRate; // M x C, 3 x M, fastphase pars
    Int1D ends; // chunk index where each chromo ends
    String1D sampleids, chrs;
    Int2D pos; // store position of markers of each chunk
    MyFloat2D gls; // store gl(N, M*3) of each chunk
};

inline auto calc_position_distance(const Int1D & markers)
{
    Int1D dl(markers.size());
    dl[0] = 0;
    for(size_t i = 1; i < markers.size(); i++) dl[i] = markers[i] - markers[i - 1];
    return dl;
}

inline auto calc_distRate(const Int1D & markers, int C, double Ne = 20000, double expRate = 0.5)
{
    MyArr1D distRate(markers.size());
    double nGen = 4 * Ne / C;
    // distRate(i) = (markers[i] - markers[i - 1]) * nGen * expRate / 1e8;
    distRate(0) = 1; //  act as sentinel. so dim aligns with M
    for(size_t i = 1; i < markers.size(); i++)
        distRate(i) = std::exp(-(markers[i] - markers[i - 1]) * expRate * nGen / 1e5);
    return distRate;
}

// check initialize_sigmaCurrent_m in STITCH
// double nGen = 4 * Ne / C;
inline auto calc_transRate_diploid(const Int1D & dl, double nGen, double expRate = 0.5)
{
    MyArr2D transRate(3, dl.size());
    MyArr1D distRate(dl.size());
    distRate(0) = 1; //  act as sentinel. so dim aligns with M
    for(size_t i = 1; i < dl.size(); i++) distRate(i) = std::exp(-dl[i] * expRate * nGen / 100 / 1e6);
    transRate.row(0) = distRate.square();
    transRate.row(1) = distRate * (1 - distRate);
    transRate.row(2) = (1 - distRate).square();
    return transRate;
}

/*
** @param gli  genotype likelihoods of current individual i, (M, 3)
** @param F    cluster-specific allele frequence (M, C)
** @return emission probability (M, C2)
*/
inline auto get_emission_by_gl(const MyArr2D & gli, const MyArr2D & F, double minEmission = 1e-6)
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
    // emitDip = emitDip.colwise() / emitDip.rowwise().maxCoeff(); // normalize it
    // emitDip = (emitDip < minEmission).select(minEmission, emitDip);
    return emitDip;
}

/*
** @param alpha    forward probability, (C2,M)
** @param beta     backwards probability (C2,M)
** @param E        emission probability with individual genotype likelihood,(C2,M)
** @param R        transition probability (3,M)
** @param PI       cluster frequency (C,M)
** @return individual log likelihood
*/
inline auto forward_backwards_diploid(MyArr2D & alpha,
                                      MyArr2D & beta,
                                      const MyArr2D & E,
                                      const MyArr2D & R,
                                      const MyArr2D & PI)
{
    const int M = alpha.cols();
    const int C = PI.rows();
    MyArr1D sumTmp1(C), cs(M); // store sum over internal loop
    double constTmp;
    // ======== forward recursion ===========
    int z1, z2, z12;
    int s{0};
    alpha.col(s) = E.col(s) * (PI.col(s).matrix() * PI.col(s).transpose().matrix()).reshaped().array();
    cs(s) = 1.0 / alpha.col(s).sum();
    alpha.col(s) *= cs(s); // normalize it
    // alpha_s = emit * (alpha_(s-1) * R + pi(z1) * tmp1(z2) + pi(z2) * tmp2(z1) + P(switch into z1) *
    // P(switch into z2) * constTmp)
    for(s = 1; s < M; s++)
    {
        sumTmp1 = alpha.col(s - 1).reshaped(C, C).rowwise().sum() * R(1, s);
        // constTmp = alpha.col(s - 1).sum() * R(2, s);
        constTmp = R(2, s); // since alpha.col(s - 1).sum() = 1
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                alpha(z12, s) = E(z12, s)
                                * (alpha(z12, s - 1) * R(0, s) + PI(z1, s) * sumTmp1(z2)
                                   + PI(z2, s) * sumTmp1(z1) + PI(z1, s) * PI(z2, s) * constTmp);
            }
        }
        cs(s) = 1.0 / alpha.col(s).sum();
        alpha.col(s) *= cs(s); // normalize it
    }

    // ======== backward recursion ===========
    s = M - 1; // set last site
    beta.col(s).setConstant(1.0);
    for(s = M - 2; s >= 0; s--)
    {
        auto beta_mult_emit = E.col(s + 1) * beta.col(s + 1);
        sumTmp1.setZero();
        for(constTmp = 0, z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                sumTmp1(z1) += beta_mult_emit(z12) * PI(z2, s + 1) * R(1, s + 1);
                constTmp += beta_mult_emit(z12) * PI(z1, s + 1) * PI(z2, s + 1) * R(2, s + 1);
            }
        }
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                // apply scaling
                beta(z12, s) =
                    (beta_mult_emit(z12) * R(0, s + 1) + sumTmp1(z1) + sumTmp1(z2) + constTmp) * cs(s + 1);
            }
        }
    }
    return cs;
}

inline auto getClusterLikelihoods(int ind,
                                  MyArr2D & alpha,
                                  MyArr2D & beta,
                                  const MyFloat1D & GL,
                                  const MyFloat1D & R,
                                  const MyFloat1D & PI,
                                  const MyFloat1D & F,
                                  const double minEmission = 1e-6)
{
    const int C2 = alpha.rows();
    const int M = alpha.cols();
    const int C = F.size() / M;
    int g1, g2, k1, k2, k12;
    int igs = ind * M * 3;
    // ======== forward and backward recursion ===========
    MyArr2D emit(C2, M);
    MyArr1D sumTmp1(C); // store sum over internal loop
    MyArr1D cs = MyArr1D::Zero(M);
    double constTmp;
    int s{0};
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            emit(k12, s) = 0;
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emit(k12, s) += GL[igs + (g1 + g2) * M + s]
                                    * (g1 * F[k1 * M + s] + (1 - g1) * (1 - F[k1 * M + s]))
                                    * (g2 * F[k2 * M + s] + (1 - g2) * (1 - F[k2 * M + s]));
                }
            }
            // emit(k12, s) = emit(k12, s) < minEmission ? minEmission : emit(k12, s);
            alpha(k12, s) = emit(k12, s) * PI[s * C + k1] * PI[s * C + k2];
            cs(s) += alpha(k12, s);
        }
    }
    cs(s) = 1 / cs(s);
    alpha.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1 = alpha.col(s - 1).reshaped(C, C).rowwise().sum() * R[s * 3 + 1];
        constTmp = alpha.col(s - 1).sum() * R[s * 3 + 2];
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                emit(k12, s) = 0;
                for(g1 = 0; g1 <= 1; g1++)
                {
                    for(g2 = 0; g2 <= 1; g2++)
                    {
                        emit(k12, s) += GL[igs + (g1 + g2) * M + s]
                                        * (g1 * F[k1 * M + s] + (1 - g1) * (1 - F[k1 * M + s]))
                                        * (g2 * F[k2 * M + s] + (1 - g2) * (1 - F[k2 * M + s]));
                    }
                }
                // emit(k12, s) = emit(k12, s) < minEmission ? minEmission : emit(k12, s);
                alpha(k12, s) =
                    emit(k12, s)
                    * (alpha(k12, s - 1) * R[s * 3 + 0] + PI[s * C + k1] * sumTmp1(k2)
                       + PI[s * C + k2] * sumTmp1(k1) + PI[s * C + k1] * PI[s * C + k2] * constTmp);
                cs(s) += alpha(k12, s);
            }
        }
        cs(s) = 1 / cs(s);
        alpha.col(s) *= cs(s); // normalize it
    }
    // double indLike = LikeForwardInd.col(M - 1).sum(); // just 1
    // ======== backward recursion ===========
    s = M - 1;
    beta.col(s).setConstant(1.0);
    for(s = M - 2; s >= 0; s--)
    {
        auto beta_mult_emit = emit.col(s + 1) * beta.col(s + 1);
        sumTmp1.setZero();
        for(constTmp = 0, k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI[(s + 1) * C + k2] * R[(s + 1) * 3 + 1];
                constTmp +=
                    beta_mult_emit(k12) * PI[(s + 1) * C + k1] * PI[(s + 1) * C + k2] * R[(s + 1) * 3 + 2];
            }
        }
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                beta(k12, s) =
                    (beta_mult_emit(k12) * R[(s + 1) * 3 + 0] + sumTmp1(k1) + sumTmp1(k2) + constTmp)
                    * cs(s + 1);
            }
        }
    }
    return cs;
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

// @params GL genotype likelihoods, N x M x 3
inline auto estimate_af_by_gl(const MyFloat1D & GL, int N, int M, int niter = 100, double tol = 1e-4)
{
    Arr1D af_est = Arr1D::Constant(M, 0.25);
    Arr1D af_tmp = Arr1D::Zero(M);
    for(int it = 0; it < niter; it++)
    {
        for(int j = 0; j < M; j++)
        {
            af_tmp(j) = af_est(j);
            double p0, p1, p2, pt = 0.0;
            for(int i = 0; i < N; i++)
            {
                p0 = GL[i * M * 3 + 0 * M + j] * (1.0 - af_est(j)) * (1.0 - af_est(j));
                p1 = GL[i * M * 3 + 1 * M + j] * 2 * af_est(j) * (1.0 - af_est(j));
                p2 = GL[i * M * 3 + 2 * M + j] * af_est(j) * af_est(j);
                pt += (p1 + 2 * p2) / (2 * (p0 + p1 + p2));
            }
            af_est(j) = pt / (double)M;
        }
        double diff = sqrt((af_est - af_tmp).array().square().sum() / M);
        if(diff < tol)
            break;
        else if(it == niter - 1)
            cao.warn("EM for estimating AF may not converge\n");
    }

    return af_est;
}

#endif // COMMON_H_
