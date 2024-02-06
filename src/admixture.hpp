/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/admixture.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef ADMIXTURE_H_
#define ADMIXTURE_H_

#include "common.hpp"
#include <mutex>

class Admixture
{
  private:
    std::mutex mutex_it; // in case of race condition
    // randon engine
    std::default_random_engine rng = std::default_random_engine{};
    // BOUNDING
    double clusterFreqThreshold{1e-6}; // threshold for F
    double admixtureThreshold{1e-6}; // threshold for Q
    bool debug = false;
    bool nonewQ = false;
    bool cF = false;

  public:
    Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
    {
        rng.seed(seed);
        Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, admixtureThreshold, 1 - admixtureThreshold);
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
        F = RandomUniform<MyArr2D, std::default_random_engine>(C * K, M, rng, clusterFreqThreshold,
                                                               1 - clusterFreqThreshold);
        for(int k = 0; k < K; k++) F.middleRows(k * C, C).rowwise() /= F.middleRows(k * C, C).colwise().sum();
    }

    ~Admixture() {}

    // SHARED VARIBALES
    const int N, M, C, K; // M: number of grids in total,  C2 = C x C
    MyArr2D F; // (C x K) x M
    MyArr2D P; // C x M, for each k, F <= P
    MyArr2D Q; // K x N
    MyArr2D Ekc; // (C * K) x M, expected number of alleles per c per k
    MyArr2D NormF; // K x M

    void initIteration();
    void updateIteration();
    void protectPars();
    void constrainF();
    void setFlags(bool, bool, bool);
    void setStartPoint(const std::unique_ptr<BigAss> & genome, std::string qfile);
    double runNativeWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
    double runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
};

int run_admix_main(Options & opts);

#endif // ADMIXTURE_H_
