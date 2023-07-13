/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/fastphase.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef FASTPHASE_H_
#define FASTPHASE_H_

#include "common.hpp"
#include <cmath>
#include <mutex>

using fbd_res1 = std::tuple<double, MyArr2D, MyArr2D, MyArr2D, MyArr1D>;
using fbd_res2 = std::tuple<MyFloat1D, MyArr2D, MyArr2D, MyArr2D, MyArr2D>;

class FastPhaseK2
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK2(int m, int n, int c, int seed) : M(m), N(n), C(c), C2(c * c)
    {
        auto rng = std::default_random_engine{};
        rng.seed(seed);
        F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, alleleEmitThreshold,
                                                               1 - alleleEmitThreshold);
        GP.setZero(M * 3, N);
    }
    ~FastPhaseK2() {}

    // BOUNDING
    double minRate{0.1}, maxRate{100}; // threshold for R
    double clusterFreqThreshold{1e-4}; // threshold for PI
    double alleleEmitThreshold{1e-4}; // threshold for F(P)

    // FLAGS
    bool debug{0};

    // SHARED VARIBALES
    const int M, N, C, C2; // C2 = C x C
    int G, B; // G: number of grids after collapsing block
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr1D pi; // C, PI in first SNP
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr2D Ezj; // C x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D Ezg1, Ezg2; // C x M
    Int1D dist; // physical position distance between two markers
    MyArr1D AF;
    Int2D grids;
    double nGen;

    void initRecombination(const Int1D & pos, int B_ = 1, double Ne = 20000);
    void collapse_and_resize(const Int1D & pos, double tol_r = 1e-6);
    void initIteration();
    void updateIteration();
    void callGenoLoopC(int, int, int, const MyArr2D &, const MyArr1D &);
    double runWithOneThread(int, const MyFloat1D &);
    double forwardAndBackwardsLowRam(int, const MyFloat1D &, bool);
    double forwardAndBackwardsLowRamNormal(int, const MyFloat1D &, bool);
    double forwardAndBackwardsLowRamCollapse(int, const MyFloat1D &, bool);
    fbd_res1 forwardAndBackwardsHighRam(int, const MyFloat1D &, bool);
    fbd_res1 forwardAndBackwardsHighRamNormal(int, const MyFloat1D &, bool);
    fbd_res1 forwardAndBackwardsHighRamCollapse(int, const MyFloat1D &, bool);
};

fbd_res2 make_input_per_chunk(const std::unique_ptr<BigAss> & genome,
                              const int ic,
                              const int niters,
                              const int seed,
                              const bool collapse = false,
                              const double tol_pi = 0.99,
                              const double tol_r = 1e-5);

int run_impute_main(Options & opts);

#endif // FASTPHASE_H_
