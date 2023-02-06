#ifndef ADMIXTURE_H_
#define ADMIXTURE_H_

#include "common.hpp"
#include <iostream>
#include <mutex>

class Admixture
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    Admixture(int n, int m, int c, int k, int seed);
    ~Admixture();
    const int N, M, C, K; // C2 = C x C
    ArrDouble2D F; // (M x C) x K
    ArrDouble2D Q; // K x N
    ArrDouble2D Ekg; // K x N, expected number of alleles per k
    ArrDouble2D Ekc; // (M x C) x K, expected number of alleles per c per k
    ArrDouble1D NormF; // K

    double updateQ(int ind, ArrDouble2D clusters);
    void updateF();
    void initIteration();
};

inline Admixture::Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<ArrDouble2D, std::default_random_engine>(M * C, K, rng, 0.05, 0.95);
    for(int s = 0; s < M; s++) // normalize it per snp per K
        F.middleRows(s * C, C).rowwise() /= F.middleRows(s * C, C).colwise().sum();
    Q = RandomUniform<ArrDouble2D, std::default_random_engine>(K, N, rng, 0.05, 0.95);
    Q = Q.rowwise() / Q.colwise().sum(); // normalize it per N individual
}

inline Admixture::~Admixture() {}

inline void Admixture::initIteration()
{
    Ekg.setZero(K, N);
    Ekc.setZero(M * C, K);
    NormF.setZero(K);
}

inline void Admixture::updateF()
{
    F = Ekc.rowwise() / NormF.transpose();
}

/*
** @params icluster cluster likelihoods for one individual, C x C x M
*/
inline double Admixture::updateQ(int ind, ArrDouble2D icluster)
{
    double norm;
    int s, k1, k2, k12, c1, c2, c12;
    ArrDouble2D w(C * C, K * K);
    // std::cout << std::this_thread::get_id() << ": " << ind << '\n';
    for(s = 0; s < M; s++)
    {
        norm = 0;
        for(k1 = 0; k1 < K; k1++)
        {
            for(k2 = 0; k2 < K; k2++)
            {
                k12 = k1 * K + k2;
                for(c1 = 0; c1 < C; c1++)
                {
                    for(c2 = 0; c2 < C; c2++)
                    {
                        c12 = c1 * C + c2;
                        w(c12, k12) =
                            icluster(c12, s) * F(s * C + c1, k1) * Q(k1, ind) * F(s * C + c2, k2) * Q(k2, ind);
                        norm += w(c12, k12);
                    }
                }
            }
        }
        norm = (icluster.col(s)
                * ((F.middleRows(s * C, C).matrix() * Q.col(ind).matrix())
                   * (F.middleRows(s * C, C).matrix() * Q.col(ind).matrix()).transpose())
                      .array()
                      .reshaped())
                   .sum();
        for(c1 = 0; c1 < C; c1++)
        {
            for(c2 = 0; c2 < C; c2++)
            {
                c12 = c1 * C + c2;
                for(k1 = 0; k1 < K; k1++)
                {
                    for(k2 = 0; k2 < K; k2++)
                    {
                        k12 = k1 * K + k2;
                        std::lock_guard<std::mutex> lock(mutex_it);
                        Ekg(k1, ind) += w(c12, k12) / norm;
                        Ekg(k2, ind) += w(c12, k12) / norm;
                        Ekc(s * C + c1, k1) += w(c12, k12) / norm;
                        Ekc(s * C + c2, k2) += w(c12, k12) / norm;
                        NormF(k1) += w(c12, k12) / norm;
                        NormF(k2) += w(c12, k12) / norm;
                    }
                }
            }
        }
    }

    Q.col(ind) = Ekg.col(ind) / (2 * M);
    return norm;
}

#endif // ADMIXTURE_H_
