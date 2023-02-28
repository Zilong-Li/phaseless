/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/admixture.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
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
    MyArr2D FI; // (C x K) x M
    MyArr2D Q; // K x N
    MyArr2D Ekg; // (M x K) x N, expected number of alleles per k per n
    MyArr2D Ekc; // (C * K) x M, expected number of alleles per c per k
    MyArr2D NormF; // K x M

    void initIteration(double tol = 1e-6);
    void updateIteration();
    void writeQ(std::string out);
    double runNaiveWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
    double runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
    double runWithSingleChunk(int, const MyFloat1D &, const MyArr2D &, const MyArr2D &, const MyArr2D &);
};

inline Admixture::Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    FI = RandomUniform<MyArr2D, std::default_random_engine>(C * K, M, rng, 0.05, 0.95);
    Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, 0.05, 0.95);
}

inline Admixture::~Admixture() {}

inline double Admixture::runNaiveWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D w((C * C + C) / 2, K * K);
    MyArr2D iEkc, LikeForwardInd, LikeBackwardInd;
    double norm = 0, llike = 0;
    int c1, c2, c12, cc;
    int k1, k2, k12, s;
    const bool logscale = false;
    for(int ic = 0, m = 0; ic < genome->nchunks; ic++)
    {
        int iM = genome->pos[ic].size();
        LikeForwardInd.setZero(C * C, iM); // likelihood of forward recursion for ind i, not log
        LikeBackwardInd.setZero(C * C, iM); // likelihood of backward recursion for ind i, not log
        getClusterLikelihoods(ind, LikeForwardInd, LikeBackwardInd, genome->gls[ic], genome->transRate[ic],
                              genome->PI[ic], genome->F[ic], true);
        iEkc.setZero(C * K, iM);
        for(s = 0; s < iM; s++, m++)
        {
            for(norm = 0, cc = 0, c1 = 0; c1 < C; c1++)
            {
                auto H1 = (Q.col(ind) * FI(Eigen::seqN(c1, K, C), m)).sum();
                for(c2 = c1; c2 < C; c2++)
                {
                    c12 = c1 * C + c2;
                    auto H2 = (Q.col(ind) * FI(Eigen::seqN(c2, K, C), m)).sum();
                    for(k1 = 0; k1 < K; k1++)
                    {
                        for(k2 = 0; k2 < K; k2++)
                        {
                            k12 = k1 * K + k2;
                            if(logscale)
                            {
                                w(cc, k12) = log(LikeForwardInd(c12, s)) + log(LikeBackwardInd(c12, s))
                                             + log(FI(k1 * C + c1, m)) + log(Q(k1, ind))
                                             + log(FI(k2 * C + c2, m)) + log(Q(k2, ind));
                                if(c1 != c2) w(cc, k12) -= log(2);
                                norm += exp(w(cc, k12));
                            }
                            else
                            {
                                w(cc, k12) = LikeForwardInd(c12, s) * LikeBackwardInd(c12, s)
                                             * FI(k1 * C + c1, m) * Q(k1, ind) * FI(k2 * C + c2, m)
                                             * Q(k2, ind);
                                w(cc, k12) /= H1 * H2;
                                if(c1 != c2) w(cc, k12) *= 2;
                                norm += w(cc, k12);
                            }
                        }
                    }
                    ++cc;
                }
            }
            llike += log(norm);
            if(logscale) w = w.exp();
            w /= norm;
            for(cc = 0, c1 = 0; c1 < C; c1++)
            {
                for(c2 = c1; c2 < C; c2++)
                {
                    for(k1 = 0; k1 < K; k1++)
                    {
                        for(k2 = 0; k2 < K; k2++)
                        {
                            k12 = k1 * K + k2;
                            Ekg(ind * K + k1, m) += w(cc, k12) ;
                            Ekg(ind * K + k2, m) += w(cc, k12) ;
                            iEkc(k1 * C + c1, s) += w(cc, k12) ;
                            iEkc(k2 * C + c2, s) += w(cc, k12) ;
                        }
                    }
                    ++cc;
                }
            }
            // iEkc = (iEkc < 2e-6).select(2e-6, iEkc);
            // if(s == 10) std::cerr << ind << ": " << iEkc.col(s).transpose() << std::endl;
            // if(m == 10)
            //     std::cerr << ind << ": " << Ekg.middleRows(ind * K, K).col(m).transpose() << std::endl;
        }
        {
            std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(m - iM, iM) += iEkc;
            // NormF.middleCols(m - iM, iM) += iNormF;
            NormF.middleCols(m - iM, iM) += Ekg.middleRows(ind * K, K).middleCols(m - iM, iM);
        }
    }
    // auto deno = Ekg.middleRows(ind * K, K).sum(); // 2 * M
    // update Q, Q.colwise().sum() should be 1
    for(int k = 0; k < K; k++) Q(k, ind) = Ekg.row(ind * K + k).sum() / (2 * M);

    return llike;
}

inline double Admixture::runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D kapa;
    MyArr2D LikeForwardInd, LikeBackwardInd;
    MyArr1D Hsz(K);
    double llike = 0;
    int c1, k1, s;
    for(int ic = 0, m = 0; ic < genome->nchunks; ic++)
    {
        int iM = genome->pos[ic].size();
        LikeForwardInd.setZero(C * C, iM);
        LikeBackwardInd.setZero(C * C, iM);
        getClusterLikelihoods(ind, LikeForwardInd, LikeBackwardInd, genome->gls[ic], genome->transRate[ic],
                              genome->PI[ic], genome->F[ic], true); // return gamma
        kapa.setZero(C * K, iM);
        for(s = 0; s < iM; s++, m++)
        {
            for(c1 = 0; c1 < C; c1++)
            {
                // auto Hsz = Q.col(ind) * FI.middleRows(c1 * K, K).col(m);
                Hsz = Q.col(ind) * FI(Eigen::seqN(c1, K, C), m);
                auto rho = (LikeForwardInd.middleRows(c1 * C, C).col(s)
                            * LikeBackwardInd.middleRows(c1 * C, C).col(s))
                               .sum();
                // kapa.middleRows(c1 * K, K).col(s) = Hsz * (rho / Hsz.sum()); // K x C x M layout
                kapa(Eigen::seqN(c1, K, C), s) = 2 * Hsz * (rho / Hsz.sum()); // C x K x M layout
            }
            // kapa = (kapa < 1e-6).select(1e-6, kapa);
            // std::cout << kapa.col(s).sum() << std::endl;
            for(k1 = 0; k1 < K; k1++) Ekg(ind * K + k1, m) = kapa.middleRows(k1 * C, C).col(s).sum();
            // if(s == 10) std::cerr << ind << ": " << kapa.col(s).transpose() << std::endl;
            // if(m == 10)
            //     std::cerr << ind << ": " << Ekg.middleRows(ind * K, K).col(m).transpose() << std::endl;
        }
        {
            std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(m - iM, iM) += kapa;
            NormF.middleCols(m - iM, iM) += Ekg.middleRows(ind * K, K).middleCols(m - iM, iM);
        }
    }
    // auto deno = Ekg.middleRows(ind * K, K).sum(); // 2 * M
    // assert(abs(deno - 2 * M) < 1e-4);
    // update Q, Q.colwise().sum() should be 1
    for(int k = 0; k < K; k++) Q(k, ind) = Ekg.row(ind * K + k).sum() / (2 * M);
    // std::cout << Q.col(ind).sum() << std::endl;

    return llike;
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param transRate (x^2, x(1-x), (1-x)^2),M x 3
** @param PI        PI from fastphase model
** @param F        F from fastphase model
** @return individual total likelihood
*/
inline double Admixture::runWithSingleChunk(int ind,
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
    MyArr2D iEkc = MyArr2D::Zero(K * C, iM);
    MyArr2D iNormF = MyArr2D::Zero(K, iM);
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
        w /= norm;
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
                        Ekg(ind * K + k1, s) += w(c12, k12);
                        Ekg(ind * K + k2, s) += w(c12, k12);
                        iEkc(k1 * C + c1, s) += w(c12, k12);
                        iEkc(k2 * C + c2, s) += w(c12, k12);
                        iNormF(k1, s) += w(c12, k12);
                        iNormF(k2, s) += w(c12, k12);
                    }
                }
            }
        }
    }
    // update Q, Q.colwise().sum() should be 1
    for(int k = 0; k < K; k++) Q(k, ind) = Ekg.row(ind * K + k).sum() / (2 * M);
    {
        std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
        Ekc += iEkc;
        NormF += iNormF;
    }
    return llike;
}

inline void Admixture::initIteration(double tol)
{
    if(Q.isNaN().any()) throw std::runtime_error("NaN in Q\n");
    Q = (Q < tol).select(tol, Q); // lower bound
    Q = (Q > 1 - tol).select(1 - tol, Q); // upper bound
    Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual

    if(FI.isNaN().any()) throw std::runtime_error("NaN in FI\n");
    FI = (FI < tol).select(tol, FI); // lower bound
    FI = (FI > 1 - tol).select(1 - tol, FI); // upper bound
    for(int k = 0; k < K; k++) // normalize FI per snp per k
        FI.middleRows(k * C, C).rowwise() /= FI.middleRows(k * C, C).colwise().sum();

    Ekg.setZero(N * K, M);
    Ekc.setZero(C * K, M);
    NormF.setZero(K, M);
}

inline void Admixture::updateIteration()
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
