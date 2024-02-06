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
    double alleleEmitThreshold{1e-4}; // threshold for P
    double clusterFreqThreshold{1e-4}; // threshold for F
    double admixtureThreshold{1e-6}; // threshold for Q
    double minHapfreq{0.01}; // min haplotype frequency, or min(1/(10*C), 1/100)
    // FLAGS
    bool debug{0}, local{0}, post{1}, NQ{0}, NF{0}, NP{1}, NR{1};
    int G{0}, B{1}; // G: number of grids after collapsing block
    double nGen;

  public:
    FastPhaseK2(int n, int m, int c, int seed) : N(n), M(m), C(c), CC(c * c)
    {
        rng.seed(seed);
        minHapfreq = std::min(1.0 / (10 * C), minHapfreq);
        F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, alleleEmitThreshold, 1 - alleleEmitThreshold);
        GP.setZero(M * 3, N);
        collapse = Bool1D::Constant(M, false);
    }
    ~FastPhaseK2() {}


    // SHARED VARIBALES
    const int N, M, C, CC; // CC = C x C
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr2D PI; // C x M, cluster frequency
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr1D er; // M, jumping rate
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr2D Ezg1, Ezg2; // C x M
    MyArr2D Ezj; // C x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D HapSum; // C x M, sum(gammaK) for all inds
    Int1D dist; // physical position distance between two markers
    Int1D pos_chunk; // store the start pos of each chunk in the full scale
    Int3D grids; // physical position of each grid for each chunk
    Bool1D collapse;
    MyArr1D AF;

    void initRecombination(const Int2D & pos, std::string rfile = "", int B_ = 1, double Ne = 20000);
    void setFlags(double, double, double, bool, bool, bool, bool, bool);
    void refillHaps(int); // re-sample F for sites with hapfreq < minHapfreq
    void protectPars();
    void initIteration();
    void updateIteration();
    void callGenoLoopC(int, int, int, const MyArr1D &, const MyArr1D &, const MyArr1D &);
    double hmmIterWithJumps(const MyFloat1D &, const int, const int, bool);
    double runAllChunks(const MyFloat2D &, const int, bool);
};

int run_impute_main(Options & opts);

#endif // FASTPHASE_H_
