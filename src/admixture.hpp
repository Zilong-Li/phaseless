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

  public:
    Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
    {
        auto rng = std::default_random_engine{};
        rng.seed(seed);
        Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, admixtureThreshold, 1 - admixtureThreshold);
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
        F = RandomUniform<MyArr2D, std::default_random_engine>(C * K, M, rng, clusterFreqThreshold,
                                                               1 - clusterFreqThreshold);
        normalizeF();
    }

    Admixture(int n, int m, int c, int k, int seed, const std::string & qfile) : N(n), M(m), C(c), K(k)
    {
        auto rng = std::default_random_engine{};
        rng.seed(seed);
        // read Q file
        Q = MyArr2D(K, N);
        load_csv(qfile, Q);
        F = RandomUniform<MyArr2D, std::default_random_engine>(C * K, M, rng, clusterFreqThreshold,
                                                               1 - clusterFreqThreshold);
        normalizeF();
    }

    ~Admixture() {}

    // BOUNDING
    double clusterFreqThreshold{1e-6}; // threshold for F
    double admixtureThreshold{1e-6}; // threshold for Q
    bool debug = false;

    const int N, M, C, K; // M: number of grids in total,  C2 = C x C
    MyArr2D F; // (C x K) x M
    MyArr2D Q; // K x N
    MyArr2D Ekc; // (C * K) x M, expected number of alleles per c per k
    MyArr2D NormF; // K x M

    void initIteration();
    void updateIteration();
    void protectPars();
    void normalizeF();
    void writeQ(std::string out);
    double runNativeWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
    double runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
};

int run_admix_main(Options & opts);

#endif // ADMIXTURE_H_
