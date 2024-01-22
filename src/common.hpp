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
#include <cassert>
#include <climits>
#include <clocale>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <sys/utsname.h>
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

inline int SIG_COND = 1; // if we catch signal then quit program nicely

inline void handler(int s)
{
    SIG_COND = 0;
    cao.warn("Caught SIGNAL: ", s, ". will try to exit nicely. wait for current threads to finish");
}

// STD TYPES
using Bool1D = std::vector<bool>;
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
#ifdef USE_FLOAT
using MyFloat = float; // use float if no accuracy drops
#else
using MyFloat = double; // use float if no accuracy drops
#endif

using MyFloat1D = std::vector<MyFloat>;
using MyFloat2D = std::vector<MyFloat1D>;
using MyFloat3D = std::vector<MyFloat2D>;
using MyMat2D = Eigen::Matrix<MyFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using MyMat1D = Eigen::Matrix<MyFloat, Eigen::Dynamic, 1, Eigen::ColMajor>;
using MyArr2D = Eigen::Array<MyFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using MyArr1D = Eigen::Array<MyFloat, Eigen::Dynamic, 1, Eigen::ColMajor>;

template<typename MatrixType, typename RandomEngineType>
inline MatrixType RandomUniform(const Eigen::Index numRows,
                                const Eigen::Index numCols,
                                RandomEngineType & engine,
                                typename MatrixType::Scalar a = 0.0,
                                typename MatrixType::Scalar b = 1.0)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{a, b};
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
};

struct Options
{
    int ichunk{0}, chunksize{10000}, K{2}, C{10}, nadmix{1000}, nimpute{40}, nthreads{1}, seed{999};
    int gridsize{1};
    double ltol{1e-1}, info{0}, tol_pi{0.99}, tol_r{1e-5};
    double ptol{1e-6}; // threshold for P
    double ftol{1e-6}; // threshold for F
    double qtol{1e-6}; // threshold for Q
    bool noaccel{0}, noscreen{0}, single_chunk{0}, debug{0}, collapse{0};
    bool nQ{0}, nP{0}, nF{0}, nR{0}, aQ{0}, oVCF{0};
    std::string out, in_beagle, in_vcf, in_bin, in_impute, in_joint;
    std::string samples{""}, region{""}, in_plink{""}, in_qfile{""}, in_pfile{""}, in_rfile{""};
    std::string opts_in_effect{"Options in effect:\n   "};
};

// all the genome info I need from fastphase
struct BigAss
{
    int chunksize, nsamples, nsnps, nchunks;
    int B, G, C; // B: snps in a grid; G: total number of grids in a genome
    MyFloat2D PI, F, R, GammaAE; // M x C, 3 x M, fastphase pars
    Int1D ends; // chunk index where each chromo ends
    String1D sampleids, chrs;
    Int2D pos; // store position of markers of each chunk
    MyFloat2D gls; // store gl(N, M*3) of each chunk
};

// joint parameters
struct Pars
{
    void init(int K_,
              int C_,
              int M_,
              int N_,
              const MyArr1D & ier,
              const MyArr2D & iP,
              const MyArr2D & iQ,
              const std::vector<MyArr2D> & iF)
    {
        K = K_;
        C = C_;
        M = M_;
        N = N_;
        P = MyFloat1D(iP.data(), iP.data() + iP.size());
        Q = MyFloat1D(iQ.data(), iQ.data() + iQ.size());
        er = MyFloat1D(ier.data(), ier.data() + ier.size());
        for(size_t k = 0; k < iF.size(); k++) F.emplace_back(MyFloat1D(iF[k].data(), iF[k].data() + iF[k].size()));
    }
    int K, C, M, N;
    MyFloat1D P, Q;
    MyFloat2D F; // K x C x M, ancestral cluster frequency
    MyFloat1D er; // M, jumping rate
    Int2D pos; // store position of markers of each chunk
    MyFloat2D gls; // store gl(N, M*3) of each chunk
    String1D sampleids;
};

//******************************************************************************
//                               STRING UTILS
//******************************************************************************

inline std::string get_machine()
{
    struct utsname unameData;
    if(uname(&unameData) != 0)
    {
        perror("uname");
        exit(EXIT_FAILURE);
    }
    std::string machine{unameData.machine};
    std::string node{unameData.nodename};
    std::string release{unameData.release};
    std::string version{unameData.version};
    std::string sysname{unameData.sysname};
    return "Machine name: " + machine + "\nNode name: " + node + "\nOperating system release: " + release
           + "\nOperating system version: " + version + "\nOperating system name: " + sysname + "\n";
}

inline std::vector<std::string> split_string(const std::string & s, const std::string & separators)
{
    std::vector<std::string> ret;
    bool is_seperator[256] = {false};
    for(auto & ch : separators)
    {
        is_seperator[(unsigned int)ch] = true;
    }
    int begin = 0;
    for(int i = 0; i <= (int)s.size(); i++)
    {
        if(is_seperator[(uint8_t)s[i]] || i == (int)s.size())
        {
            ret.push_back(std::string(s.begin() + begin, s.begin() + i));
            begin = i + 1;
        }
    }
    return ret;
}

inline std::string trim_string(const std::string & s)
{
    int begin = 0, end = (int)s.size();
    while(begin < end && s[begin] == ' ') begin++;
    while(begin < end && s[end - 1] == ' ') end--;
    return std::string(s.begin() + begin, s.begin() + end);
}

inline bool ends_with(std::string const & str, std::string const & ending)
{
    if(ending.size() > str.size())
        return false;
    else
        return std::equal(ending.begin(), ending.end(), str.end() - ending.size());
}

inline bool starts_with(std::string const & str, std::string const & ending)
{
    if(ending.size() > str.size())
        return false;
    else
        return std::equal(ending.begin(), ending.end(), str.begin());
}

//******************************************************************************
//                               RECOMBINATION
//******************************************************************************

inline auto calc_position_distance(const Int1D & markers)
{
    Int1D dl(markers.size());
    dl[0] = 0;
    for(size_t i = 1; i < markers.size(); i++) dl[i] = markers[i] - markers[i - 1];
    return dl;
}

inline MyArr2D er2R(const MyArr1D & er)
{
    MyArr2D R(3, er.size());
    R.row(0) = er.square();
    R.row(1) = er * (1 - er);
    R.row(2) = (1 - er).square();
    return R;
}

inline void protect_er(MyArr1D & er)
{
    const double maxer{std::exp(-1e-9)};
    er = (er > maxer).select(maxer, er);
    er = (er < 1.0 - maxer).select(1.0 - maxer, er);
}

// check initialize_sigmaCurrent_m in STITCH
// double nGen = 4 * Ne / C;
inline MyArr1D calc_er(const Int1D & dl, double nGen, double expRate = 0.5)
{
    MyArr1D er(dl.size());
    // for(size_t i = 1; i < dl.size(); i++) distRate(i) = std::exp(-dl[i] / 1e6);
    for(size_t i = 1; i < dl.size(); i++) er(i) = std::exp(-dl[i] * expRate * nGen / 1e8);
    protect_er(er);
    return er;
}

// check initialize_sigmaCurrent_m in STITCH
// double nGen = 4 * Ne / C;
inline MyArr2D calc_transRate_diploid(const Int1D & dl, double nGen, double expRate = 0.5)
{
    auto er = calc_er(dl, nGen, expRate);
    return er2R(er);
}

/*
** @param gli  genotype likelihoods of current individual i, (M, 3)
** @param F    cluster-specific allele frequence (M, C)
** @return emission probability (M, C2)
*/
inline auto get_emission_by_gl(const MyArr2D & gli, const MyArr2D & F, double minEmission = 1e-10)
{
    int k1, k2, g1, g2;
    const int M = F.rows();
    const int C = F.cols();
    MyArr2D emitDip(M, C * C); // emission probabilies, nsnps x (C x C)
    for(k1 = 0; k1 < C; k1++)
        for(k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.col(g1 + g2) * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                                * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                }
            }
        }
    // emitDip = emitDip.colwise() / emitDip.rowwise().maxCoeff(); // normalize
    emitDip = (emitDip < minEmission).select(minEmission, emitDip);
    return emitDip;
}

/*
** @param gli  genotype likelihoods of current individual i, (M, 3)
** @param F    cluster-specific allele frequence (M, C)
** @return emission probability (M, C2)
*/
inline auto get_emission_by_grid(const MyFloat1D & GL,
                                 const MyFloat1D & F,
                                 int ind,
                                 int M,
                                 int B,
                                 double minEmission = 1e-10)
{
    const int C = F.size() / M;
    const int C2 = C * C;
    const int nGrids = B > 1 ? (M + B - 1) / B : M;
    MyArr2D emitGrid = MyArr2D::Ones(C2, nGrids);
    int z1, z2, z12, i, s, e, g, g1, g2;
    int igs = ind * M * 3;
    for(g = 0; g < nGrids; g++)
    {
        s = g * B;
        e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                for(i = s; i <= e; i++)
                {
                    double emit = 0;
                    for(g1 = 0; g1 <= 1; g1++)
                    {
                        for(g2 = 0; g2 <= 1; g2++)
                        {
                            emit += GL[igs + (g1 + g2) * M + i] * (g1 * F[z1 * M + i] + (1 - g1) * (1 - F[z1 * M + i]))
                                    * (g2 * F[z2 * M + i] + (1 - g2) * (1 - F[z2 * M + i]));
                        }
                    }
                    emitGrid(z12, g) *= emit;
                }
            }
        }
        // apply bounding
        // emitGrid.col(g) /= emitGrid.col(g).maxCoeff();
        // emitGrid.col(g) = (emitGrid.col(g) < minEmission).select(minEmission, emitGrid.col(g));
    }
    return emitGrid;
}

/*
** @param alpha    forward probability, (C2,M)
** @param beta     backwards probability (C2,M)
** @param E        emission probability with individual genotype
*likelihood,(C2,M)
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
    // alpha_s = emit * (alpha_(s-1) * R + pi(z1) * tmp1(z2) + pi(z2) * tmp2(z1)
    // + P(switch into z1) * P(switch into z2) * constTmp)
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
                                * (alpha(z12, s - 1) * R(0, s) + PI(z1, s) * sumTmp1(z2) + PI(z2, s) * sumTmp1(z1)
                                   + PI(z1, s) * PI(z2, s) * constTmp);
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
                beta(z12, s) = (beta_mult_emit(z12) * R(0, s + 1) + sumTmp1(z1) + sumTmp1(z2) + constTmp) * cs(s + 1);
            }
        }
    }
    return cs;
}

inline auto get_cluster_probability(int ind,
                                    const int M,
                                    MyArr2D & alpha,
                                    MyArr2D & beta,
                                    const MyFloat1D & GL,
                                    const MyFloat1D & R,
                                    const MyFloat1D & PI,
                                    const MyFloat1D & F,
                                    const double minEmission = 1e-6)
{
    const int C = F.size() / M;
    const int C2 = alpha.rows();
    const int nGrids = alpha.cols();
    const int B = (M + nGrids - 1) / nGrids;
    int g1, g2, z1, z2, z12;
    int igs = ind * M * 3;
    MyArr1D sumTmp1(C); // store sum over internal loop
    MyArr1D cs = MyArr1D::Zero(nGrids);
    double constTmp;
    // ======== forward and backward recursion ===========
    if(nGrids < M)
    {
        MyArr2D emitGrid = MyArr2D::Ones(C2, nGrids);
        int i, s, e, g{0};
        s = g * B;
        e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                for(i = s; i <= e; i++)
                {
                    double emit = 0;
                    for(g1 = 0; g1 <= 1; g1++)
                    {
                        for(g2 = 0; g2 <= 1; g2++)
                        {
                            emit += GL[igs + (g1 + g2) * M + i] * (g1 * F[z1 * M + i] + (1 - g1) * (1 - F[z1 * M + i]))
                                    * (g2 * F[z2 * M + i] + (1 - g2) * (1 - F[z2 * M + i]));
                        }
                    }
                    emitGrid(z12, g) *= emit;
                }
            }
        }
        // emitGrid.col(g) /= emitGrid.col(g).maxCoeff();
        // emitGrid.col(g) = (emitGrid.col(g) < minEmission).select(minEmission, emitGrid.col(g)); // apply bounding
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                alpha(z12, g) = emitGrid(z12, g) * PI[g * C + z1] * PI[g * C + z2];
                cs(g) += alpha(z12, g);
            }
        }
        cs(g) = 1 / cs(g);
        alpha.col(g) *= cs(g); // normalize it
        // now get the rest
        for(g = 1; g < nGrids; g++)
        {
            sumTmp1 = alpha.col(g - 1).reshaped(C, C).rowwise().sum() * R[g * 3 + 1];
            constTmp = alpha.col(g - 1).sum() * R[g * 3 + 2];
            s = g * B;
            e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
            for(z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    for(i = s; i <= e; i++)
                    {
                        double emit = 0;
                        for(g1 = 0; g1 <= 1; g1++)
                        {
                            for(g2 = 0; g2 <= 1; g2++)
                            {
                                emit += GL[igs + (g1 + g2) * M + i]
                                        * (g1 * F[z1 * M + i] + (1 - g1) * (1 - F[z1 * M + i]))
                                        * (g2 * F[z2 * M + i] + (1 - g2) * (1 - F[z2 * M + i]));
                            }
                        }
                        emitGrid(z12, g) *= emit;
                    }
                }
            }
            // emitGrid.col(g) /= emitGrid.col(g).maxCoeff();
            // emitGrid.col(g) = (emitGrid.col(g) < minEmission).select(minEmission, emitGrid.col(g)); // apply bounding
            for(z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    alpha(z12, g) = emitGrid(z12, g)
                                    * (alpha(z12, g - 1) * R[g * 3 + 0] + PI[g * C + z1] * sumTmp1(z2)
                                       + PI[g * C + z2] * sumTmp1(z1) + PI[g * C + z1] * PI[g * C + z2] * constTmp);
                    cs(g) += alpha(z12, g);
                }
            }
            cs(g) = 1 / cs(g);
            alpha.col(g) *= cs(g); // normalize it
        }
        // next backwards
        g = nGrids - 1;
        beta.col(g).setConstant(1.0);
        for(g = nGrids - 2; g >= 0; g--)
        {
            auto beta_mult_emit = emitGrid.col(g + 1) * beta.col(g + 1);
            sumTmp1.setZero();
            for(constTmp = 0, z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    sumTmp1(z1) += beta_mult_emit(z12) * PI[(g + 1) * C + z2] * R[(g + 1) * 3 + 1];
                    constTmp += beta_mult_emit(z12) * PI[(g + 1) * C + z1] * PI[(g + 1) * C + z2] * R[(g + 1) * 3 + 2];
                }
            }
            for(z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    beta(z12, g) =
                        (beta_mult_emit(z12) * R[(g + 1) * 3 + 0] + sumTmp1(z1) + sumTmp1(z2) + constTmp) * cs(g + 1);
                }
            }
        }
    }
    else if(nGrids == M)
    {
        MyArr2D emitSnp(C2, M);
        int s{0};
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                emitSnp(z12, s) = 0;
                for(g1 = 0; g1 <= 1; g1++)
                {
                    for(g2 = 0; g2 <= 1; g2++)
                    {
                        emitSnp(z12, s) += GL[igs + (g1 + g2) * M + s]
                                           * (g1 * F[z1 * M + s] + (1 - g1) * (1 - F[z1 * M + s]))
                                           * (g2 * F[z2 * M + s] + (1 - g2) * (1 - F[z2 * M + s]));
                    }
                }
                // emit(k12, s) = emit(k12, s) < minEmission ? minEmission :
                // emit(k12, s);
                alpha(z12, s) = emitSnp(z12, s) * PI[s * C + z1] * PI[s * C + z2];
                cs(s) += alpha(z12, s);
            }
        }
        cs(s) = 1 / cs(s);
        alpha.col(s) *= cs(s); // normalize it

        for(s = 1; s < M; s++)
        {
            sumTmp1 = alpha.col(s - 1).reshaped(C, C).rowwise().sum() * R[s * 3 + 1];
            constTmp = alpha.col(s - 1).sum() * R[s * 3 + 2];
            for(z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    emitSnp(z12, s) = 0;
                    for(g1 = 0; g1 <= 1; g1++)
                    {
                        for(g2 = 0; g2 <= 1; g2++)
                        {
                            emitSnp(z12, s) += GL[igs + (g1 + g2) * M + s]
                                               * (g1 * F[z1 * M + s] + (1 - g1) * (1 - F[z1 * M + s]))
                                               * (g2 * F[z2 * M + s] + (1 - g2) * (1 - F[z2 * M + s]));
                        }
                    }
                    alpha(z12, s) = emitSnp(z12, s)
                                    * (alpha(z12, s - 1) * R[s * 3 + 0] + PI[s * C + z1] * sumTmp1(z2)
                                       + PI[s * C + z2] * sumTmp1(z1) + PI[s * C + z1] * PI[s * C + z2] * constTmp);
                    cs(s) += alpha(z12, s);
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
            auto beta_mult_emit = emitSnp.col(s + 1) * beta.col(s + 1);
            sumTmp1.setZero();
            for(constTmp = 0, z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    sumTmp1(z1) += beta_mult_emit(z12) * PI[(s + 1) * C + z2] * R[(s + 1) * 3 + 1];
                    constTmp += beta_mult_emit(z12) * PI[(s + 1) * C + z1] * PI[(s + 1) * C + z2] * R[(s + 1) * 3 + 2];
                }
            }
            for(z1 = 0; z1 < C; z1++)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    z12 = z1 * C + z2;
                    beta(z12, s) =
                        (beta_mult_emit(z12) * R[(s + 1) * 3 + 0] + sumTmp1(z1) + sumTmp1(z2) + constTmp) * cs(s + 1);
                }
            }
        }
    }
    return cs;
}

inline auto get_cluster_frequency(MyArr2D & ae, const MyFloat1D & R_, const MyFloat1D & PI_)
{
    const int C2 = ae.rows();
    const int M = ae.cols();
    const int C = std::sqrt(C2);
    Eigen::Map<const MyArr2D> PI(PI_.data(), C, M);
    Eigen::Map<const MyArr2D> R(R_.data(), 3, M);

    int s{0};
    ae.col(s) = (PI.col(s).matrix() * PI.col(s).transpose().matrix()).reshaped().array();
    // ======== forward recursion ===========
    MyArr1D sumTmp1(C); // store sum over internal loop
    double constTmp;
    int z1, z2, z12;
    for(s = 1; s < M; s++)
    {
        sumTmp1 = ae.col(s - 1).reshaped(C, C).rowwise().sum() * R(1, s);
        constTmp = R(2, s); // since alpha.col(s - 1).sum() = 1
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                ae(z12, s) = (ae(z12, s - 1) * R(0, s) + PI(z1, s) * sumTmp1(z2) + PI(z2, s) * sumTmp1(z1)
                              + PI(z1, s) * PI(z2, s) * constTmp);
            }
        }
    }
    // (1.0 - ae.colwise().sum()).abs() < 1e-2 is OK. this may be due to
    // rounding error if we want to colsum equals 1.0. then normlize it
    // TODO: cluster frequency can be zero for certain cluster.
    // const double tol = 1e-6;
    // ae = (ae < tol).select(tol, ae);
    // ae = (ae > 1 - tol).select(1 - tol, ae);
    // ae.rowwise() /= ae.colwise().sum();
}

inline auto get_cluster_likelihoods(int ind,
                                    const int M,
                                    const int B,
                                    const MyFloat1D & GL,
                                    const MyFloat1D & R,
                                    const MyFloat1D & PI,
                                    const MyFloat1D & F,
                                    const double minEmission = 1e-10)
{
    const int C = F.size() / M;
    const int C2 = C * C;
    MyArr2D emitGrid = get_emission_by_grid(GL, F, ind, M, B, minEmission);
    const int nGrids = emitGrid.cols();
    MyArr2D alpha(C2, nGrids), beta(C2, nGrids), ae(C2, nGrids);
    int z1, z2, z12;
    MyArr1D sumTmp1(C); // store sum over internal loop for alpha
    MyArr1D sumTmp2(C); // store sum over internal loop for ae
    MyArr1D cs = MyArr1D::Zero(nGrids);
    double constTmp;
    // ======== forward and backward recursion ===========
    int g = 0;
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            alpha(z12, g) = emitGrid(z12, g) * PI[g * C + z1] * PI[g * C + z2];
            ae(z12, g) = PI[g * C + z1] * PI[g * C + z2];
            cs(g) += alpha(z12, g);
        }
    }
    cs(g) = 1 / cs(g);
    alpha.col(g) *= cs(g); // normalize it
    // now get the rest
    for(g = 1; g < nGrids; g++)
    {
        sumTmp1 = alpha.col(g - 1).reshaped(C, C).rowwise().sum() * R[g * 3 + 1];
        sumTmp2 = ae.col(g - 1).reshaped(C, C).rowwise().sum() * R[g * 3 + 1];
        constTmp = R[g * 3 + 2]; // since alpha.col(g).sum()==1
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                alpha(z12, g) = emitGrid(z12, g)
                                * (alpha(z12, g - 1) * R[g * 3 + 0] + PI[g * C + z1] * sumTmp1(z2)
                                   + PI[g * C + z2] * sumTmp1(z1) + PI[g * C + z1] * PI[g * C + z2] * constTmp);
                cs(g) += alpha(z12, g);
                ae(z12, g) = (ae(z12, g - 1) * R[g * 3 + 0] + PI[g * C + z1] * sumTmp2(z2)
                              + PI[g * C + z2] * sumTmp2(z1) + PI[g * C + z1] * PI[g * C + z2] * constTmp);
            }
        }
        cs(g) = 1 / cs(g);
        alpha.col(g) *= cs(g); // normalize it
    }
    // const double tol = 1e-6;
    // ae = (ae < tol).select(tol, ae);
    // ae = (ae > 1 - tol).select(1 - tol, ae);
    ae.rowwise() /= ae.colwise().sum();
    // next backwards
    g = nGrids - 1;
    beta.col(g).setConstant(1.0);
    for(g = nGrids - 2; g >= 0; g--)
    {
        auto beta_mult_emit = emitGrid.col(g + 1) * beta.col(g + 1);
        sumTmp1.setZero();
        for(constTmp = 0, z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                sumTmp1(z1) += beta_mult_emit(z12) * PI[(g + 1) * C + z2] * R[(g + 1) * 3 + 1];
                constTmp += beta_mult_emit(z12) * PI[(g + 1) * C + z1] * PI[(g + 1) * C + z2] * R[(g + 1) * 3 + 2];
            }
        }
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                beta(z12, g) =
                    (beta_mult_emit(z12) * R[(g + 1) * 3 + 0] + sumTmp1(z1) + sumTmp1(z2) + constTmp) * cs(g + 1);
            }
        }
    }
    return alpha * beta / ae;
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

inline auto divide_pos_into_grid(const Int1D & pos, int B)
{
    int M = pos.size();
    int G = (M + B - 1) / B;
    Int2D grids(G);
    int g, s, e;
    for(g = 0; g < G; g++)
    {
        s = g * B;
        e = g == G - 1 ? M - 1 : B * (g + 1) - 1;
        grids[g] = Int1D(pos.begin() + s, pos.begin() + e + 1);
    }
    return grids;
}

inline auto divide_pos_into_grid(const Int1D & pos, const Bool1D & collapse)
{
    assert(pos.size() == collapse.size());
    Int2D grids;
    for(auto i = 0; i < collapse.size(); i++)
    {
        if(collapse[i])
        {
            auto j = i + 1;
            while(collapse[j]) j++;
            grids.push_back(Int1D(pos.begin() + i, pos.begin() + j));
            i = j - 1;
        }
        else
        {
            grids.push_back(Int1D(1, pos[i]));
        }
    }

    return grids;
}

inline auto find_chunk_to_collapse(const MyArr2D & R, double tol_r = 1e-6)
{
    Bool1D collapse(R.cols(), false); // M sites
    for(auto i = 0; i < R.cols(); i++)
    {
        if(std::sqrt(R(2, i)) < tol_r) collapse[i] = true;
    }
    return collapse;
}

/*
** @params pos     snp position, first dim is each grid, second dim is snps in
*that grid
*/
inline auto calc_grid_distance(const Int2D & pos)
{
    Int1D dl(pos.size());
    dl[0] = 0;
    for(auto i = 1; i < pos.size(); i++)
    {
        dl[i] = pos[i][pos[i].size() / 2] - pos[i - 1][pos[i - 1].size() / 2];
    }
    return dl;
}

/*
** @param E original size of emission, full SNPs x C2
*/
inline auto collapse_emission_by_grid(const MyArr2D & E, const Int2D & grids, double minEmission = 1e-6)
{
    const int C2 = E.rows();
    const int G = grids.size();
    MyArr2D EG = MyArr2D::Ones(C2, G);
    int g, s, e, c, snp = 0;
    for(g = 0; g < G; g++)
    {
        // c = g == nGrids - 1 ? M - (nGrids - 1) * B : B;
        s = snp;
        e = snp + grids[g].size() - 1;
        snp = e + 1;
        for(c = s; c <= e; c++) EG.col(g) *= E.col(c);
        EG.col(g) /= EG.col(g).maxCoeff(); // rescale by maximum
        EG.col(g) = (EG.col(g) < minEmission).select(minEmission, EG.col(g)); // apply bounding
    }
    assert(snp == E.cols());

    return EG;
}

inline auto cat_stdvec_of_eigen(const std::vector<MyArr2D> & arr3)
{
    int K = arr3.size();
    int C = arr3[0].rows();
    int M = arr3[0].cols();
    MyArr2D out(C * K, M);
    for(int k = 0; k < K; k++) out.middleRows(k * C, C) = arr3[k];
    return out;
}

inline void load_csv(MyArr2D & Q, const std::string & path, bool colmajor = true, char sep = ' ')
{
    std::ifstream fin(path);
    std::string line;
    int i{0}, j{0};
    double val;
    while(std::getline(fin, line))
    {
        std::stringstream lineStream(line);
        std::string tok;
        j = 0;
        while(std::getline(lineStream, tok, sep))
        {
            try
            {
                val = std::stod(tok);
            }
            catch(const std::out_of_range & e)
            {
                val = 0;
            }
            if(colmajor)
                Q(j, i) = val;
            else
                Q(i, j) = val;
            ++j;
        }
        if(colmajor)
            assert(j == Q.rows());
        else
            assert(j == Q.cols());
        ++i;
    }
}

#endif // COMMON_H_
