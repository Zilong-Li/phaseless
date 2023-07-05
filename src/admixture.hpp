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

    bool debug = false;
    const int N, M, C, K; // M: number of grids in total,  C2 = C x C
    MyArr2D F; // (C x K) x M
    MyArr2D Q; // K x N
    MyArr2D Ekc; // (C * K) x M, expected number of alleles per c per k
    MyArr2D NormF; // K x M

    void initIteration(double tol = 1e-6);
    void updateIteration();
    void writeQ(std::string out);
    double runNativeWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
    double runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome);
};

inline Admixture::Admixture(int n, int m, int c, int k, int seed) : N(n), M(m), C(c), K(k)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(C * K, M, rng, 0.05, 0.95);
    Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, 0.05, 0.95);
}

inline Admixture::~Admixture() {}

// the complexity of this version should be O(CC + 2CK)
inline double Admixture::runOptimalWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D kapa, Ekg;
    MyArr2D alpha, beta;
    MyArr1D iQ = MyArr1D::Zero(K);
    MyArr1D Hz(C), gammaC(C);
    double norm = 0, llike = 0, tmp = 0;
    int c1, k1, s, c2, c12;
    for(int ic = 0, m = 0; ic < genome->nchunks; ic++)
    {
        const int nsnps = genome->pos[ic].size();
        const int nGrids = genome->B > 1 ? (nsnps + genome->B - 1) / genome->B : nsnps;
        alpha.setZero(C * C, nGrids);
        beta.setZero(C * C, nGrids);
        get_cluster_likelihood(ind, nsnps, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                               genome->F[ic]); // return gamma
        kapa.setZero(C * K, nGrids); // C x K x M layout
        Ekg.setZero(K, nGrids);
        for(s = 0; s < nGrids; s++, m++)
        {
            for(c1 = 0; c1 < C; c1++) Hz(c1) = (Q.col(ind) * F(Eigen::seqN(c1, K, C), m)).sum();
            gammaC = (alpha.col(s) * beta.col(s)).reshaped(C, C).colwise().sum();
            for(norm = 0, c1 = 0; c1 < C; c1++)
            {
                for(tmp = 0, c2 = 0; c2 < C; c2++)
                {
                    c12 = c1 * C + c2;
                    double xz = alpha(c12, s) * beta(c12, s) / (gammaC(c1) * gammaC(c2));
                    double zy = Hz(c1) * Hz(c2);
                    tmp += xz * zy;
                }
                norm += tmp;
                kapa(Eigen::seqN(c1, K, C), s) = (Q.col(ind) * F(Eigen::seqN(c1, K, C), m)) * tmp / Hz(c1);
            }
            llike += log(norm);
            kapa.col(s) /= kapa.col(s).sum();
            for(k1 = 0; k1 < K; k1++) Ekg(k1, s) = 2 * kapa.middleRows(k1 * C, C).col(s).sum();
        }
        iQ += Ekg.rowwise().sum();
        { // for update F
            std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(m - nGrids, nGrids) += 2 * kapa;
            NormF.middleCols(m - nGrids, nGrids) += Ekg;
        }
    }

    Q.col(ind) = iQ / (2 * M); // update Q, iQ.sum() should be 2M

    return llike;
}

// the complexity of this version is O(CCKK)
inline double Admixture::runNativeWithBigAss(int ind, const std::unique_ptr<BigAss> & genome)
{
    MyArr2D w((C * C + C) / 2, K * K);
    MyArr2D Ekg, iEkc, alpha, beta;
    MyArr1D gammaC(C);
    double norm = 0, llike = 0;
    int c1, c2, c12, cc;
    int k1, k2, k12, s;
    MyArr1D iQ = MyArr1D::Zero(K);
    for(int ic = 0, m = 0; ic < genome->nchunks; ic++)
    {
        const int nsnps = genome->pos[ic].size();
        const int nGrids = genome->B > 1 ? (nsnps + genome->B - 1) / genome->B : nsnps;
        alpha.setZero(C * C, nGrids);
        beta.setZero(C * C, nGrids);
        get_cluster_likelihood(ind, nsnps, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                               genome->F[ic]);
        iEkc.setZero(C * K, nGrids);
        Ekg.setZero(K, nGrids);
        for(s = 0; s < nGrids; s++, m++)
        {
            gammaC = (alpha.col(s) * beta.col(s)).reshaped(C, C).colwise().sum();
            for(norm = 0, cc = 0, c1 = 0; c1 < C; c1++)
            {
                for(c2 = c1; c2 < C; c2++)
                {
                    c12 = c1 * C + c2;
                    double xz = (double)alpha(c12, s) * beta(c12, s) / (gammaC(c1) * gammaC(c2));
                    for(k1 = 0; k1 < K; k1++)
                    {
                        for(k2 = 0; k2 < K; k2++)
                        {
                            k12 = k1 * K + k2;
                            w(cc, k12) = xz * F(k1 * C + c1, m) * Q(k1, ind) * F(k2 * C + c2, m) * Q(k2, ind);
                            if(c1 != c2) w(cc, k12) *= 2;
                            norm += w(cc, k12);
                        }
                    }
                    ++cc;
                }
            }
            llike += log(norm);
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
                            Ekg(k1, s) += w(cc, k12);
                            Ekg(k2, s) += w(cc, k12);
                            iEkc(k1 * C + c1, s) += w(cc, k12);
                            iEkc(k2 * C + c2, s) += w(cc, k12);
                        }
                    }
                    ++cc;
                }
            }
        }
        iQ += Ekg.rowwise().sum();
        {
            std::lock_guard<std::mutex> lock(mutex_it); // sum over all samples
            Ekc.middleCols(m - nGrids, nGrids) += iEkc;
            NormF.middleCols(m - nGrids, nGrids) += Ekg;
        }
    }
    // update Q, iQ.sum() should be 2M
    Q.col(ind) = iQ / (2 * M);

    return llike;
}

inline void Admixture::initIteration(double tol)
{
    if(Q.isNaN().any()) cao.error("NaN in Q\n");
    Q = (Q < tol).select(tol, Q); // lower bound
    Q = (Q > 1 - tol).select(1 - tol, Q); // upper bound
    Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual

    if(F.isNaN().any()) cao.error("NaN in F\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
    for(int k = 0; k < K; k++) // normalize F per snp per k
        F.middleRows(k * C, C).rowwise() /= F.middleRows(k * C, C).colwise().sum();

    Ekc.setZero(C * K, M);
    NormF.setZero(K, M);
}

inline void Admixture::updateIteration()
{
    for(int k = 0; k < K; k++) F.middleRows(k * C, C) = Ekc.middleRows(k * C, C).rowwise() / NormF.row(k);
}

inline void Admixture::writeQ(std::string out)
{
    std::ofstream ofs(out);
    if(!ofs) cao.error(out, strerror(errno));
    Q = (Q * 1e6).round() / 1e6;
    ofs << std::fixed << Q.transpose() << "\n";
    ofs.close();
}

#endif // ADMIXTURE_H_
