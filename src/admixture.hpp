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
    ArrDouble2D Ekg; // (M x K) x N, expected number of alleles per k per n
    ArrDouble2D Ekc; // (M x C) x K, expected number of alleles per c per k
    ArrDouble1D NormF; // K

    double updateQ(int ind, ArrDouble2D clusters);
    double runWithClusterLikelihoods(int ind,
                                     const DoubleVec1D & GL,
                                     const ArrDouble2D & transRate,
                                     const ArrDouble2D & PI,
                                     const ArrDouble2D & PF);
    void updateF();
    void initIteration();
};

inline Admixture::Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<ArrDouble2D, std::default_random_engine>(M * C, K, rng, 0.05, 0.95);
    // normalize it per snp per K
    for(int s = 0; s < M; s++) F.middleRows(s * C, C).rowwise() /= F.middleRows(s * C, C).colwise().sum();
    Q = RandomUniform<ArrDouble2D, std::default_random_engine>(K, N, rng, 0.05, 0.95);
    Q = Q.rowwise() / Q.colwise().sum(); // normalize it per N individual
}

inline Admixture::~Admixture() {}

inline void Admixture::initIteration()
{
    Ekg.setZero(N * K, M);
    Ekc.setZero(K * C, M);
    NormF.setZero(K);
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @param PI        PI from fastphase model
** @param PF        F from fastphase model
** @return individual total likelihood
*/
inline double Admixture::runWithClusterLikelihoods(int ind,
                                                   const DoubleVec1D & GL,
                                                   const ArrDouble2D & transRate,
                                                   const ArrDouble2D & PI,
                                                   const ArrDouble2D & PF)
{
    int k1, k2, k12;
    // ======== forward and backward recursion ===========
    Eigen::Map<const ArrDouble2D> gli(GL.data() + ind * M * 3, 3, M);
    auto emitDip = emissionCurIterInd(gli, PF, false);
    ArrDouble2D LikeForwardInd(C * C, M); // likelihood of forward recursion for ind i, not log
    ArrDouble2D LikeBackwardInd(C * C, M); // likelihood of backward recursion for ind i, not log
    ArrDouble1D sumTmp1(C), sumTmp2(C); // store sum over internal loop
    ArrDouble1D cs(M);
    double constTmp;
    // ======== forward recursion ===========
    int s{0};
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            k12 = k1 * C + k2;
            LikeForwardInd(k12, s) = emitDip(s, k12) * PI(s, k1) * PI(s, k2);
        }
    }
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1.setZero();
        sumTmp2.setZero();
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += LikeForwardInd(k12, s - 1);
                sumTmp2(k2) += LikeForwardInd(k12, s - 1);
            }
        }
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        sumTmp1 *= transRate(1, s);
        sumTmp2 *= transRate(1, s);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeForwardInd(k12, s) =
                    emitDip(s, k12)
                    * (LikeForwardInd(k12, s - 1) * transRate(0, s) + PI(s, k1) * sumTmp1(k2)
                       + PI(s, k2) * sumTmp2(k1) + PI(s, k1) * PI(s, k2) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    // ======== backward recursion ===========
    s = M - 1;
    LikeBackwardInd.col(s).setOnes(); // not log scale
    LikeBackwardInd.col(s) *= cs(s);
    for(s = M - 2; s >= 0; s--)
    {
        sumTmp1.setZero();
        sumTmp2.setZero();
        constTmp = 0;
        auto beta_mult_emit = emitDip.row(s + 1).transpose() * LikeBackwardInd.col(s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                sumTmp1(k1) += beta_mult_emit(k12) * PI(s + 1, k2);
                sumTmp2(k2) += beta_mult_emit(k12) * PI(s + 1, k1);
                constTmp += beta_mult_emit(k12) * PI(s + 1, k1) * PI(s + 1, k2);
            }
        }
        sumTmp1 *= transRate(1, s + 1);
        sumTmp2 *= transRate(1, s + 1);
        constTmp *= transRate(2, s + 1);
        for(k1 = 0; k1 < C; k1++)
        {
            for(k2 = 0; k2 < C; k2++)
            {
                k12 = k1 * C + k2;
                LikeBackwardInd(k12, s) =
                    beta_mult_emit(k12) * transRate(0, s + 1) + sumTmp1(k1) + sumTmp2(k2) + constTmp;
            }
        }
        // apply scaling
        LikeBackwardInd.col(s) *= cs(s);
    }
    auto icluster = LikeForwardInd * LikeBackwardInd;

    //============= ready to update Q ============
    double norm;
    int c1, c2, c12;
    ArrDouble2D w(C * C, K * K);
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
                        w(c12, k12) = icluster(c12, s) * F(s * C + c1, k1) * Q(k1, ind) * F(s * C + c2, k2)
                                      * Q(k2, ind);
                        norm += w(c12, k12);
                    }
                }
            }
        }
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
                        Ekg(ind * K + k1, s) += w(c12, k12) / norm;
                        Ekg(ind * K + k2, s) += w(c12, k12) / norm;
                        Ekc(k1 * C + c1, s) += w(c12, k12) / norm;
                        Ekc(k2 * C + c2, s) += w(c12, k12) / norm;
                        NormF(k1) += w(c12, k12) / norm;
                        NormF(k2) += w(c12, k12) / norm;
                    }
                }
            }
        }
    }
    // update Q
    for(k1 = 0; k1 < K; k1++) Q(k1, ind) = Ekg.row(ind * K + k1).sum() / (2 * M);

    return norm;
}

inline void Admixture::updateF()
{
    F = Ekc.reshaped(M * C, K).rowwise() / NormF.transpose();
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
                        w(c12, k12) = icluster(c12, s) * F(s * C + c1, k1) * Q(k1, ind) * F(s * C + c2, k2)
                                      * Q(k2, ind);
                        norm += w(c12, k12);
                    }
                }
            }
        }
        // norm = (icluster.col(s)
        //         * ((F.middleRows(s * C, C).matrix() * Q.col(ind).matrix())
        //            * (F.middleRows(s * C, C).matrix() * Q.col(ind).matrix()).transpose())
        //               .array()
        //               .reshaped())
        //            .sum();
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
                        // w(c12, k12) =
                        //     icluster(c12, s) * F(s * C + c1, k1) * Q(k1, ind) * F(s * C + c2, k2) * Q(k2,
                        //     ind);
                        std::lock_guard<std::mutex> lock(mutex_it);
                        Ekg(k1 * M + s, ind) += w(c12, k12) / norm;
                        Ekg(k2 * M + s, ind) += w(c12, k12) / norm;
                        Ekc(s * C + c1, k1) += w(c12, k12) / norm;
                        Ekc(s * C + c2, k2) += w(c12, k12) / norm;
                        NormF(k1) += w(c12, k12) / norm;
                        NormF(k2) += w(c12, k12) / norm;
                    }
                }
            }
        }
    }
    // update Q
    for(k1 = 0; k1 < K; k1++) Q(k1, ind) = Ekg.middleRows(k1 * M, M).col(ind).sum() / (2 * M);

    return norm;
}

#endif // ADMIXTURE_H_
