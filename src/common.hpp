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
using Int1D = std::vector<int>;
using Int2D = std::vector<Int1D>;
using Int3D = std::vector<Int2D>;
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
using Bool1D = Eigen::Array<bool, Eigen::Dynamic, 1, Eigen::ColMajor>;
using Bool2D = std::vector<Bool1D>;

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

inline Eigen::IOFormat fmt6(6, Eigen::DontAlignCols, " ", "\n");
inline Eigen::IOFormat fmt10(10, Eigen::DontAlignCols, " ", "\n");

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
    int gridsize{1}, refillHaps{0};
    double ltol{1e-1}, info{0}, tol_pi{0.99}, tol_r{1e-5};
    double ptol{1e-6}; // threshold for P
    double ftol{1e-6}; // threshold for F
    double qtol{1e-6}; // threshold for Q
    bool noaccel{0}, noscreen{0}, single_chunk{0}, debug{0}, collapse{0};
    bool nQ{0}, nP{0}, nF{0}, nR{0}, aQ{0}, oVCF{0}, eHap{0}, oF{0}, cF{0};
    std::string out, in_beagle, in_vcf, in_bin, in_impute, in_joint;
    std::string samples{""}, region{""}, in_plink{""}, in_qfile{""}, in_pfile{""}, in_rfile{""};
    std::string opts_in_effect{"Options in effect:\n   "};
};

// all the genome info I need from fastphase
struct BigAss
{
    int chunksize, nsamples, nsnps, nchunks, B, C, G;
    MyFloat2D PI, F, R, AE; // M x C, 3 x M, fastphase pars
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

inline Int1D calc_position_distance(const Int1D & markers)
{
    Int1D dl(markers.size());
    dl[0] = 0;
    for(size_t i = 1; i < markers.size(); i++) dl[i] = markers[i] - markers[i - 1];
    return dl;
}

inline Int1D calc_grid_distance(const Int1D & pos, const Bool1D & collapse)
{
    // B = 1
    if((collapse == true).count() == 0) return calc_position_distance(pos);
    // B > 1, split pos into grids
    const int G = ((collapse == true).count() + 1) / 2;
    Int2D gpos(G);
    int s, e, g, m = pos.size();
    for(s = 0, e = 1, g = 0; g < G; g++)
    {
        for(;; e++)
            if(collapse[e] == true) break;
        if(g == G - 1) e = m - 1;
        gpos[g] = Int1D(pos.begin() + s, pos.begin() + e + 1);
        // cao.cerr("size of g:", gpos[g].size(), ",s:",s, ",e:",e);
        // for(auto j : gpos[g]) cao.cerr(j);
        s = e + 1 >= m ? m - 1 : e + 1;
        e = s + 1 >= m ? m - 1 : s + 1;
    }
    Int1D dl(G);
    dl[0] = 0;
    for(g = 1; g < G; g++)
    {
        // cao.cerr(g, ":", gpos[g][gpos[g].size() / 2], ",", gpos[g - 1][gpos[g - 1].size() / 2]);
        dl[g] = gpos[g][gpos[g].size() / 2] - gpos[g - 1][gpos[g - 1].size() / 2];
    }
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
    for(size_t i = 1; i < dl.size(); i++) er(i) = std::exp(-dl[i] / 1e6);
    // for(size_t i = 1; i < dl.size(); i++) er(i) = std::exp(-dl[i] * expRate * nGen / 1e8);
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
** @param P    cluster-specific allele frequence (M, C)
** @return emission probability (M, C2)
*/
inline MyArr2D get_emission_by_gl(const MyArr2D & gli, const MyArr2D & P, double minEmission = 1e-10)
{
    int k1, k2, g1, g2;
    const int M = P.rows();
    const int C = P.cols();
    MyArr2D emitDip(M, C * C); // emission probabilies, nsnps x (C x C)
    for(k1 = 0; k1 < C; k1++)
        for(k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.col(g1 + g2) * (g1 * P.col(k1) + (1 - g1) * (1 - P.col(k1)))
                                                * (g2 * P.col(k2) + (1 - g2) * (1 - P.col(k2)));
                }
            }
        }
    // emitDip = emitDip.colwise() / emitDip.rowwise().maxCoeff(); // normalize
    // emitDip = (emitDip < minEmission).select(minEmission, emitDip);
    return emitDip;
}

/*
** @param gli  genotype likelihoods of current individual i, (M, 3)
** @param P    cluster-specific allele frequence (M, C)
** @return emission probability (M, C2)
*/
inline MyArr2D get_emission_by_grid(const MyArr2D & gli,
                                    const MyArr2D & P,
                                    const Int2D & grids,
                                    double minEmission = 1e-10)
{
    const int C = P.cols();
    const int nGrids = grids.size();
    MyArr2D emitGrid = MyArr2D::Ones(C * C, nGrids);
    int z1, z2, z12, i, s, e, g, g1, g2, m;
    for(g = 0; g < nGrids; g++)
    {
        s = grids[g][0];
        e = grids[g][grids[g].size()];
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                for(i = s; i <= e; i++)
                {
                    double emit = 0;
                    m = i + grids[g].size();
                    for(g1 = 0; g1 <= 1; g1++)
                    {
                        for(g2 = 0; g2 <= 1; g2++)
                        {
                            emit += gli(m, g1 + g2) * (g1 * P(m, z1) + (1 - g1) * (1 - P(m, z1)))
                                    * (g2 * P(m, z2) + (1 - g2) * (1 - P(m, z2)));
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
** @param emit     emission probability(C2,M)
** @param R        jump probability (3,M)
** @param PI       cluster frequency (C,M)
** @return alpha, beta and scaling vector
*/
inline auto forward_backwards_diploid(const MyArr2D & emit, const MyArr2D & R, const MyArr2D & PI)
{
    const int M = emit.cols();
    const int C2 = emit.rows();
    const int C = PI.rows();
    MyArr2D alpha(C2, M), beta(C2, M);
    MyArr1D sumTmp1(C), cs(M); // store sum over internal loop
    double constTmp;
    // ======== forward recursion ===========
    int z1, z2, z12;
    int s{0};
    alpha.col(s) = emit.col(s) * (PI.col(s).matrix() * PI.col(s).transpose().matrix()).reshaped().array();
    cs(s) = 1.0 / alpha.col(s).sum();
    alpha.col(s) *= cs(s); // normalize it
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
                alpha(z12, s) = emit(z12, s)
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
        auto beta_mult_emit = emit.col(s + 1) * beta.col(s + 1);
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

    return std::tuple(alpha, beta, cs);
}

/// R: 3 x M; PI: C x M
inline auto get_cluster_frequency(const MyArr2D & R, const MyArr2D & PI)
{
    const int C = PI.rows();
    const int M = R.cols();
    MyArr2D ae(C * C, M);

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
    return ae;
}

inline auto get_cluster_likelihoods(const MyArr2D & gli,
                                    const MyArr2D & P,
                                    const MyArr2D & R,
                                    const MyArr2D & PI,
                                    const MyArr2D & AE,
                                    const double minEmission = 1e-10)
{
    MyArr2D emit = get_emission_by_gl(gli, P).transpose(); // CC x S
    const auto [alpha, beta, cs] = forward_backwards_diploid(emit, R, PI);
    // reuse emit
    emit = (alpha * beta) / AE;
    emit.rowwise() /= emit.colwise().sum(); // norm it
    return emit;
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
inline Arr1D estimate_af_by_gl(const MyFloat1D & GL, int N, int M, int niter = 100, double tol = 1e-4)
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

// grid size must be >=3
inline Bool1D find_grid_to_collapse(const Int1D & pos, int B)
{
    int M = pos.size();
    Bool1D collapse = Bool1D::Constant(M, false);
    if(B == 1) return collapse;
    int G = (M + B - 1) / B; // get ceiling number
    int g, s, e;
    for(g = 0; g < G; g++)
    {
        s = g * B;
        e = g == G - 1 ? M - 1 : B * (g + 1) - 1;
        collapse(s) = true;
        collapse(e) = true;
    }
    return collapse;
}

inline Bool1D find_grid_to_collapse(const MyArr2D & R, double tol_r = 1e-6)
{
    Bool1D collapse = Bool1D::Constant(R.cols(), false);
    for(auto i = 0; i < R.cols(); i++)
    {
        if(std::sqrt(R(2, i)) < tol_r) collapse(i) = true;
    }
    return collapse;
}

inline Int2D divide_pos_into_grid(const Int1D & pos, const Bool1D & collapse)
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

/*
** @param E original size of emission, full SNPs x C2
*/
inline MyArr2D collapse_emission_by_grid(const MyArr2D & E, const Int2D & grids, double minEmission = 1e-10)
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
        // EG.col(g) /= EG.col(g).maxCoeff(); // rescale by maximum
        // EG.col(g) = (EG.col(g) < minEmission).select(minEmission, EG.col(g)); // apply bounding
    }
    assert(snp == E.cols());

    return EG;
}

inline MyArr2D cat_stdvec_of_eigen(const std::vector<MyArr2D> & arr3)
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
