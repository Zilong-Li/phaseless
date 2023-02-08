#ifndef ADMIXTURE_H_
#define ADMIXTURE_H_

#include "common.hpp"
#include <fstream>
#include <mutex>

class Admixture
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    Admixture(int n, int m, int c, int k, int seed);
    ~Admixture();
    const int N, M, C, K; // C2 = C x C
    MyArr2D FI; // (K x C) x M
    MyArr2D Q; // K x N
    MyArr2D Ekg; // (M x K) x N, expected number of alleles per k per n
    MyArr2D Ekc; // (K x C) x M, expected number of alleles per c per k
    MyArr2D NormF; // K x M

    void initIteration();
    void updateF();
    void writeQ(std::string out);
    double runWithClusterLikelihoods(int ind,
                                     const MyFloat1D & GL,
                                     const MyArr2D & transRate,
                                     const MyArr2D & PI,
                                     const MyArr2D & F);
};

inline Admixture::Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    FI = RandomUniform<MyArr2D, std::default_random_engine>(K * C, M, rng, 0.05, 0.95);
    for(int k = 0; k < K; k++) // normalize it per snp per k
        FI.middleRows(k * C, C).rowwise() /= FI.middleRows(k * C, C).colwise().sum();
    // for(int k = 0; k < K; k++) std::cout << FI.middleRows(k * C, C).colwise().sum() << "\n";
    Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, 0.05, 0.95);
    Q.rowwise() /= Q.colwise().sum(); // normalize it per individual
}

inline Admixture::~Admixture() {}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @param PI        PI from fastphase model
** @param F        F from fastphase model
** @return individual total likelihood
*/
inline double Admixture::runWithClusterLikelihoods(int ind,
                                                   const MyFloat1D & GL,
                                                   const MyArr2D & transRate,
                                                   const MyArr2D & PI,
                                                   const MyArr2D & F)
{
    auto icluster = getClusterLikelihoods(ind, GL, transRate, PI, F);
    const int iM = icluster.cols();
    double norm = 0, llike = 0;
    int c1, c2, c12;
    int k1, k2, k12, s;
    MyArr2D w(C * C, K * K);
    MyArr2D Ekc_i = MyArr2D::Zero(K * C, iM);
    MyArr2D NormF_i = MyArr2D::Zero(K, iM);
    for(s = 0; s < iM; s++)
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
                        w(c12, k12) = icluster(c12, s) * FI(k1 * C + c1, s) * Q(k1, ind) * FI(k2 * C + c2, s)
                                      * Q(k2, ind);
                        norm += w(c12, k12);
                    }
                }
            }
        }
        llike += log(norm);
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
                        Ekg(ind * K + k1, s) += w(c12, k12) / norm;
                        Ekg(ind * K + k2, s) += w(c12, k12) / norm;
                        Ekc_i(k1 * C + c1, s) += w(c12, k12) / norm;
                        Ekc_i(k2 * C + c2, s) += w(c12, k12) / norm;
                        NormF_i(k1, s) += w(c12, k12) / norm;
                        NormF_i(k2, s) += w(c12, k12) / norm;
                    }
                }
            }
        }
    }
    // update Q, Q.colwise().sum() should be 1
    for(int k = 0; k < K; k++) Q(k, ind) = Ekg.row(ind * K + k).sum() / (2 * M);
    {
        std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
        Ekc += Ekc_i;
        NormF += NormF_i;
    }
    return llike;
}

inline void Admixture::initIteration()
{
    Ekg.setZero(N * K, M);
    Ekc.setZero(K * C, M);
    NormF.setZero(K, M);
}

inline void Admixture::updateF()
{
    for(int k = 0; k < K; k++) FI.middleRows(k * C, C) = Ekc.middleRows(k * C, C).rowwise() / NormF.row(k);
    // for(int k = 0; k < K; k++) std::cout << FI.middleRows(k * C, C).sum() << "\n";
    // for(int k = 0; k < K; k++) std::cout << FI.middleRows(k * C, C).colwise().sum() << "\n";
}

inline void Admixture::writeQ(std::string out)
{
    std::ofstream ofs(out);
    if(!ofs) throw std::runtime_error(out + ": " + strerror(errno));
    Q = (Q * 1e6).round() / 1e6;
    ofs << std::fixed << Q.transpose() << "\n";
    ofs.close();
}
#endif // ADMIXTURE_H_
