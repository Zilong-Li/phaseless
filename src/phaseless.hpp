#ifndef PHASELESS_H_
#define PHASELESS_H_

#include "common.hpp"
#include <mutex>

class Phaseless
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    Phaseless(int k, int c, int n, int m, int seed) : K(k), C(c), N(n), M(m), KK(k * k), CC(c * c)
    {
        auto rng = std::default_random_engine{};
        rng.seed(seed);
        P = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, alleleEmitThreshold, 1 - alleleEmitThreshold);
        Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, admixtureThreshold, 1 - admixtureThreshold);
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
        F.resize(K);
        for(int i = 0; i < K; i++)
        {
            F[i] = RandomUniform<MyArr2D, std::default_random_engine>(C, M, rng, clusterFreqThreshold,
                                                                      1 - clusterFreqThreshold);
            F[i].rowwise() /= F[i].colwise().sum();
        }
    }
    ~Phaseless() {}

    // BOUNDING
    double minRate{0.1}, maxRate{100}; // threshold for R
    double alleleEmitThreshold{1e-6}; // threshold for P
    double clusterFreqThreshold{1e-6}; // threshold for F
    double admixtureThreshold{1e-6}; // threshold for Q

    // FLAGS
    bool debug{0};

    // SHARED VARIBALES
    const int K, C, N, M, KK, CC; // CC = C x C, KK = K x K
    double nGen;
    Int1D pos_chunk; // store the start pos of each chunk in the full scale
    Int1D dist; // physical position distance between two markers
    MyArr1D er, et; // M, jumping rate
    std::vector<MyArr2D> F; // K x C x M, ancestral cluster frequency
    MyArr2D P; // M x C, ancestral cluster-specific allele frequence
    MyArr2D Q; // K x N, admixture proportions for all individuals
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    MyArr2D EclusterA1, EclusterA2; // C x M
    MyArr2D Eancestry; // K x N
    MyArr2D PclusterK; // C x K x M

    void initRecombination(const Int1D & pos, double Ne, int B);
    void initRecombination(const Int2D & pos, double Ne, int B);
    void setBoundaries(double, double, double);
    void protectPars();
    void initIteration();
    void updateIteration();
    double runForwardBackwards(const int, const int, const MyFloat1D &, bool);
    double runBigass(int, const MyFloat2D &, bool);
    void getForwardPrevSums(const MyArr2D &, MyArr2D &, MyArr2D &, MyArr2D &, MyArr1D &);
    void getBackwardPrevSums(const MyArr2D &, MyArr2D &, MyArr2D &, MyArr2D &, MyArr1D &, double &, int, int);
    void moveForward(int,
                     int,
                     MyArr2D &,
                     const MyArr1D &,
                     const MyArr2D &,
                     const MyArr2D &,
                     const MyArr2D &,
                     const MyArr2D &,
                     const MyArr1D &);
    void moveBackward(int,
                      int,
                      MyFloat,
                      MyArr2D &,
                      const MyArr2D &,
                      const MyArr2D &,
                      const MyArr2D &,
                      const MyArr2D &,
                      const MyArr1D &,
                      const double);
    void getPosterios(const int,
                      const int,
                      const MyArr2D &,
                      const MyArr2D &,
                      const MyArr1D &,
                      const std::vector<MyArr2D> &,
                      const std::vector<MyArr2D> &);
};

int run_phaseless_main(Options & opts);

#endif // PHASELESS_H_
