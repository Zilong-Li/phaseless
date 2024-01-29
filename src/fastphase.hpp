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

class FastPhaseK2
{
  private:
    std::mutex mutex_it; // in case of race condition
    // randon engine
    std::default_random_engine rng = std::default_random_engine{};
    // BOUNDING
    double minRate{0.1}, maxRate{100}; // threshold for R
    double alleleEmitThreshold{1e-6}; // threshold for P
    double clusterFreqThreshold{1e-6}; // threshold for F
    double admixtureThreshold{1e-6}; // threshold for Q

  public:
    FastPhaseK2(int n, int m, int c, int seed) : N(n), M(m), C(c), CC(c * c)
    {
        rng.seed(seed);
        F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, alleleEmitThreshold, 1 - alleleEmitThreshold);
        PI = MyArr2D::Ones(C, M);
        PI.rowwise() /= PI.colwise().sum(); // normalize it per site
        R = MyArr2D(3, M);
        GP.setZero(M * 3, N);
    }
    ~FastPhaseK2() {}

    // FLAGS
    bool debug{0}, local{0}, post{1}, NQ{0}, NF{0}, NP{1}, NR{1};

    // SHARED VARIBALES
    const int N, M, C, CC; // CC = C x C
    int G, B; // G: number of grids after collapsing block
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr1D er; // M, jumping rate
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr2D Ezg1, Ezg2; // C x M
    MyArr2D Ezj; // C x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D HapSum; // C x M, sum(gammaK) for all inds
    double nGen;
    Int1D dist; // physical position distance between two markers
    Int1D pos_chunk; // store the start pos of each chunk in the full scale
    MyArr1D AF;

    void initRecombination(const Int2D & pos, std::string rfile = "", int B = 1, double Ne = 20000);
    void setFlags(double, double, double, bool, bool, bool, bool, bool);
    void protectPars();
    void initIteration();
    void updateIteration();
    void callGenoLoopC(int, int, int, const MyArr2D &, const MyArr1D &);
    double hmmIterWithJumps(const MyFloat1D &, const int, const int, bool);
    double runAllChunks(const MyFloat2D &, const int, bool);
};

int run_impute_main(Options & opts);

#endif // FASTPHASE_H_
