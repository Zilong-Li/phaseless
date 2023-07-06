#ifndef FASTPHASEK4_H_
#define FASTPHASEK4_H_

#include "../src/common.hpp"
#include <mutex>

class FastPhaseK4
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseK4(int n, int m, int c, int seed);
    ~FastPhaseK4();

    // SHARED VARIBALES
    std::ofstream ofs;
    const int N, M, C, C2; // C2 = C x C
    MyArr2D GP; // genotype probabilies for all individuals, N x (M x 3)
    MyArr2D PI; // nsnps x C
    MyArr2D F; // nsnps x C
    MyArr2D transHap; // nsnps x (C x C)
    MyArr2D transDip; // nsnps x (C x C x C x C)

    void updateClusterFreqPI(const MyArr2D &, double);
    void updateAlleleFreqWithinCluster(const MyArr2D &, double);
    void transitionCurIter(const MyArr1D &);
    auto callGenotypeInd(const MyArr2D &);
    void openClusterFile(std::string out);
    double forwardAndBackwards(int, const MyFloat1D &, MyArr2D &, MyArr2D &, bool);
};

inline FastPhaseK4::FastPhaseK4(int n, int m, int c, int seed)
: N(n), M(m), C(c), C2(c * c), GP(M * 3, N), transHap(M, C2), transDip(M, C2 * C2)
{
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, 0.0, 1.0);
    PI = PI.colwise() / PI.rowwise().sum(); // normalize it
}

inline FastPhaseK4::~FastPhaseK4()
{
    if(ofs.is_open()) ofs.close();
}

inline void FastPhaseK4::openClusterFile(std::string out)
{
    ofs.open(out, std::ios::binary);
    if(!ofs) throw std::runtime_error(out + ": " + strerror(errno));
    ofs.write((char *)&N, 4);
    ofs.write((char *)&M, 4);
    ofs.write((char *)&C, 4);
}

inline auto FastPhaseK4::callGenotypeInd(const MyArr2D & indPostProbsZandG)
{
    MyArr1D geno = MyArr1D::Zero(M * 3);
    int z1, z2, z12, g1, g2, g12, g3;
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    g3 = g1 + g2;
                    geno(Eigen::seqN(g3, M, 3)) += indPostProbsZandG.col(z12 * 4 + g12);
                }
            }
        }
    }
    return geno;
}

/*
** @param ind current individual i
** @param GL  genotype likelihood of all individuals in snp major form
** @return individual total likelihood
*/
inline double FastPhaseK4::forwardAndBackwards(int ind,
                                               const MyFloat1D & GL,
                                               MyArr2D & postProbsZ,
                                               MyArr2D & postProbsZandG,
                                               bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D emitDip(M, C2); // emission probabilies, nsnps x (C x C)
    int k1, k2, g1, g2, g12;
    for(k1 = 0; k1 < C; k1++)
    {
        for(k2 = 0; k2 < C; k2++)
        {
            emitDip.col(k1 * C + k2).setZero();
            for(g1 = 0; g1 <= 1; g1++)
            {
                for(g2 = 0; g2 <= 1; g2++)
                {
                    emitDip.col(k1 * C + k2) += gli.col(g1 + g2)
                                                * (g1 * F.col(k1) + (1 - g1) * (1 - F.col(k1)))
                                                * (g2 * F.col(k2) + (1 - g2) * (1 - F.col(k2)));
                }
            }
        }
    }
    emitDip = emitDip.log();
    MyArr2D logLikeForwardInd(C2, M); // log likelihood of forward recursion for ind i
    MyArr2D logLikeBackwardInd(C2, M); // log likelihood of backward recursion for ind i

    // ======== forward recursion ===========
    int z1, z2, z12, z12s, z12e;
    double protect_me = 0; // for underflow protection
    // first site
    int s{0};
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            logLikeForwardInd(z12, s) = emitDip(s, z12) + PI(s, z1) * PI(s, z2);
        }
    }

    // do forward recursion
    for(s = 1; s < M; s++)
    {
        // MyArr1D likeForwardTmp;
        auto likeForwardTmp = Eigen::exp(logLikeForwardInd.col(s - 1) - protect_me);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                z12s = (z1 * C + z2) * C2; // fist index
                z12e = z12s + C2 - 1; // last index
                logLikeForwardInd(z12, s) =
                    protect_me + emitDip(s, z12)
                    + log((likeForwardTmp * transDip(s, Eigen::seq(z12s, z12e)).transpose()).sum());
            }
        }
        // protect_me = logLikeForwardInd.row(s).maxCoeff();
        protect_me = logLikeForwardInd(C2 - 1, s);
    }
    // total likelhoods of the individual
    double indLogLikeForwardAll = protect_me + log((logLikeForwardInd.col(M - 1) - protect_me).exp().sum());

    // ======== backward recursion ===========
    // set last site
    logLikeBackwardInd.col(M - 1).setZero();
    // site M-2 to 0
    protect_me = 0;
    for(s = M - 2; s >= 0; s--)
    {
        auto likeBackwardTmp =
            Eigen::exp(emitDip.row(s + 1).transpose() + logLikeBackwardInd.col(s + 1) - protect_me);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                logLikeBackwardInd(z12, s) =
                    protect_me
                    + log((likeBackwardTmp * transDip(s + 1, Eigen::seqN(z12, C2, C2)).transpose()).sum());
            }
        }
        // protect_me = logLikeBackwardInd.row(s).maxCoeff();
        protect_me = logLikeBackwardInd(C2 - 1, s);
    }

    std::lock_guard<std::mutex> lock(mutex_it); // lock here if RAM cost really matters

    // ======== post decoding get p(Z|X, G) ===========
    // MyArr2D indPostProbsZ; // post probabilities for ind i, M x (C x C)
    MyArr2D indPostProbsZ = (logLikeBackwardInd + logLikeForwardInd - indLogLikeForwardAll).exp().transpose();
    // ======== post decoding get p(Z,G|X, theta) ===========
    // MyArr2D indPostProbsZandG, M x (C x C x 2 x 2)
    MyArr2D indPostProbsZandG(M, C2 * 4);
    MyArr1D tmpSum(M);
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            tmpSum.setZero();
            z12 = z1 * C + z2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    indPostProbsZandG.col(z12 * 4 + g12) = gli.col(g1 + g2)
                                                           * (g1 * F.col(z1) + (1 - g1) * (1 - F.col(z1)))
                                                           * (g2 * F.col(z2) + (1 - g2) * (1 - F.col(z2)));
                    tmpSum += indPostProbsZandG.col(z12 * 4 + g12);
                }
            }
            indPostProbsZandG.middleCols(z12 * 4, 4).colwise() *= indPostProbsZ.col(z12);
            indPostProbsZandG.middleCols(z12 * 4, 4).colwise() /= tmpSum;
        }
    }

    if(call_geno)
    {
        // std::lock_guard<std::mutex> lock(mutex_it);
        GP.col(ind) = callGenotypeInd(indPostProbsZandG);
        // output likelihood of each cluster
        MyArr2D likeCluster = (logLikeBackwardInd + logLikeForwardInd).exp().transpose();
        ofs.write((char *)likeCluster.data(), M * C2 * 8);
    }
    postProbsZ += indPostProbsZ;
    postProbsZandG += indPostProbsZandG;
    // std::cout << std::this_thread::get_id() << ": " << ind << '\n';

    return indLogLikeForwardAll;
}

/*
** @param distRate distance or recombination rate between two markers
*/
inline void FastPhaseK4::transitionCurIter(const MyArr1D & distRate)
{
    int to1, to2, from1, from2;
    int k4d, z2d1, z2d2;
    for(from1 = 0; from1 < C; from1++)
    {
        for(to1 = 0; to1 < C; to1++)
        {
            z2d1 = from1 * C + to1;
            if(from1 == to1)
                transHap.col(z2d1) = distRate + (1 - distRate) * PI.col(to1);
            else
                transHap.col(z2d1) = (1 - distRate) * PI.col(to1);
        }
    }

    for(from1 = 0; from1 < C; from1++)
    {
        for(from2 = 0; from2 < C; from2++)
        {
            for(to1 = 0; to1 < C; to1++)
            {
                for(to2 = 0; to2 < C; to2++)
                {
                    k4d = to2 + (to1 + (from2 + from1 * C) * C) * C;
                    z2d1 = from1 * C + to1;
                    z2d2 = from2 * C + to2;
                    transDip.col(k4d) = transHap.col(z2d1) * transHap.col(z2d2);
                }
            }
        }
    }
}

inline void FastPhaseK4::updateClusterFreqPI(const MyArr2D & postProbsZ, double tol)
{
    int z1, z2, z12;
    PI.setZero();
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            PI.col(z1) += postProbsZ.col(z12);
            PI.col(z2) += postProbsZ.col(z12);
        }
    }
    PI /= 2 * N;
    // map to domain
    if(PI.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    PI = (PI < tol).select(tol, PI); // lower bound
    PI = (PI > 1 - tol).select(1 - tol, PI); // upper bound
    // normalize it now
    PI = PI.colwise() / PI.rowwise().sum();
}

inline void FastPhaseK4::updateAlleleFreqWithinCluster(const MyArr2D & postProbsZandG, double tol)
{
    int z1, z2, z12, g1, g2, g12;
    MyArr2D Ekg = MyArr2D::Zero(M, C * 2);
    for(z1 = 0; z1 < C; z1++)
    {
        for(z2 = 0; z2 < C; z2++)
        {
            z12 = z1 * C + z2;
            for(g1 = 0; g1 < 2; g1++)
            {
                for(g2 = 0; g2 < 2; g2++)
                {
                    g12 = g1 * 2 + g2;
                    Ekg.col(z1 * 2 + g1) += postProbsZandG.col(z12 * 4 + g12);
                    Ekg.col(z2 * 2 + g2) += postProbsZandG.col(z12 * 4 + g12);
                }
            }
        }
    }
    F = Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2))
        / (Ekg(Eigen::all, Eigen::seq(1, Eigen::last, 2)) + Ekg(Eigen::all, Eigen::seq(0, Eigen::last, 2)));
    // map to domain but normalization
    if(F.isNaN().any()) throw std::runtime_error("NaN in PI\n");
    F = (F < tol).select(tol, F); // lower bound
    F = (F > 1 - tol).select(1 - tol, F); // upper bound
}

#endif // FASTPHASEK4_H_
