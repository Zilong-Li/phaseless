#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

using pars = std::tuple<double, MyArr2D, MyArr2D>;

inline auto getClusterLikelihoods2(MyArr2D & LikeForwardInd,
                                   MyArr2D & LikeBackwardInd,
                                   const MyArr2D & gli,
                                   const MyArr2D & transRate,
                                   const MyArr2D & PI,
                                   const MyArr2D & F,
                                   const bool gamma = true)
{
    const int M = LikeForwardInd.cols();
    const int C = F.cols();
    MyArr2D emitDip = emissionCurIterInd(gli, F, false).transpose();
    MyArr1D sumTmp1(C), sumTmp2(C), cs(M); // store sum over internal loop
    double constTmp;
    // ======== forward recursion ===========
    int z1, z2, z12;
    int s{0};
    LikeForwardInd.col(s) =
        emitDip.col(s) * (PI.col(s).matrix() * PI.col(s).transpose().matrix()).reshaped().array();
    cs(s) = 1 / LikeForwardInd.col(s).sum();
    LikeForwardInd.col(s) *= cs(s); // normalize it
    for(s = 1; s < M; s++)
    {
        sumTmp1 = LikeForwardInd.col(s - 1).reshaped(C, C).rowwise().sum() * transRate(1, s);
        sumTmp2 = LikeForwardInd.col(s - 1).reshaped(C, C).colwise().sum() * transRate(1, s);
        constTmp = LikeForwardInd.col(s - 1).sum() * transRate(2, s);
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                LikeForwardInd(z12, s) =
                    emitDip(z12, s)
                    * (LikeForwardInd(z12, s - 1) * transRate(0, s) + PI(z1, s) * sumTmp1(z2)
                       + PI(z2, s) * sumTmp2(z1) + PI(z1, s) * PI(z2, s) * constTmp);
            }
        }
        cs(s) = 1 / LikeForwardInd.col(s).sum();
        LikeForwardInd.col(s) *= cs(s); // normalize it
    }
    double indLike = LikeForwardInd.col(M - 1).sum(); // indLike should be 1 since we normalize them
    double indLogLikeForwardAll = log((indLike / cs).sum()); // log likelhoods of the individual

    // ======== backward recursion ===========
    s = M - 1; // set last site
    LikeBackwardInd.col(s).setConstant(cs(s)); // not log scale
    for(s = M - 2; s >= 0; s--)
    {
        auto beta_mult_emit = emitDip.col(s + 1) * LikeBackwardInd.col(s + 1);
        sumTmp1.setZero();
        sumTmp2.setZero();
        for(constTmp = 0, z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                sumTmp1(z1) += beta_mult_emit(z12) * PI(z2, s + 1) * transRate(1, s + 1);
                sumTmp2(z2) += beta_mult_emit(z12) * PI(z1, s + 1) * transRate(1, s + 1);
                constTmp += beta_mult_emit(z12) * PI(z1, s + 1) * PI(z2, s + 1) * transRate(2, s + 1);
            }
        }
        for(z1 = 0; z1 < C; z1++)
        {
            for(z2 = 0; z2 < C; z2++)
            {
                z12 = z1 * C + z2;
                // apply scaling
                LikeBackwardInd(z12, s) =
                    (beta_mult_emit(z12) * transRate(0, s + 1) + sumTmp1(z1) + sumTmp2(z2) + constTmp)
                    * cs(s);
            }
        }
    }

    if(gamma) LikeForwardInd.rowwise() /= cs.transpose();
    return indLogLikeForwardAll;
}

TEST_CASE("reconstruct alpha and beta from saved pars.bin", "[test-forward-backward]")
{
    int C{10}, seed{1}, chunksize{10000}, nimpute{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    ThreadPool poolit(1);
    vector<future<pars>> res;
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 faith(genome->pos[ic], genome->nsamples, C, seed);
        double prevlike{std::numeric_limits<double>::lowest()}, loglike;
        for(int it = 0; it <= nimpute; it++)
        {
            faith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == nimpute)
                    res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                    std::ref(genome->gls[ic]), true));
                else
                    res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                    std::ref(genome->gls[ic]), false));
            }
            loglike = 0;
            for(auto && ll : res)
            {
                const auto [l, iEk, iEkg] = ll.get();
                loglike += l;
                faith.Ek += iEk;
                faith.Ekg += iEkg;
            }
            res.clear(); // clear future and renew
            faith.updateIteration();
            // REQUIRE(loglike > prevlike); // works for double not float precision
            prevlike = loglike;
        }
        MyArr2D alpha, beta;
        const int iM = genome->pos[ic].size();
        for(int ind = 0; ind < genome->nsamples; ind++)
        {
            alpha.setZero(genome->C * genome->C, iM);
            beta.setZero(genome->C * genome->C, iM);
            Eigen::Map<const MyArr2D> gli(genome->gls[ic].data() + ind * iM * 3, iM, 3);
            getClusterLikelihoods(alpha, beta, gli, faith.J, faith.PI, faith.F);
            alpha *= beta;
            REQUIRE(((alpha.colwise().sum() - 1.0).abs() < 1e-5).all());
        }
    }
}

TEST_CASE("compare optmized fbd to native fbd", "[test-forward-backward]")
{
    int C{5}, seed{1}, chunksize{10000}, nphase{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    ThreadPool poolit(4);
    vector<future<pars>> res;
    int ic = 0;
    FastPhaseK2 faith(genome->pos[ic], genome->nsamples, C, seed);
    for(int it = 0; it <= nphase; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == nphase)
                res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                std::ref(genome->gls[ic]), true));
            else
                res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                std::ref(genome->gls[ic]), false));
        }
        for(auto && ll : res)
        {
            const auto [l, iEk, iEkg] = ll.get();
            faith.Ek += iEk;
            faith.Ekg += iEkg;
        }
        res.clear(); // clear future and renew
        faith.updateIteration();
    }
    MyArr2D alpha1, beta1, alpha2, beta2;
    const int iM = genome->pos[ic].size();
    for(int ind = 0; ind < genome->nsamples; ind++)
    {
        alpha1.setZero(genome->C * genome->C, iM);
        beta1.setZero(genome->C * genome->C, iM);
        Eigen::Map<const MyArr2D> gli(genome->gls[ic].data() + ind * iM * 3, iM, 3);
        getClusterLikelihoods(alpha1, beta1, gli, faith.J, faith.PI, faith.F);
        alpha1 *= beta1;
        alpha1.rowwise() /= alpha1.colwise().sum();
        alpha2.setZero(genome->C * genome->C, iM);
        beta2.setZero(genome->C * genome->C, iM);
        getClusterLikelihoods2(alpha2, beta2, gli, faith.J, faith.PI, faith.F);
        alpha2 *= beta2;
        alpha2.rowwise() /= alpha2.colwise().sum();
        auto delta = alpha1 - alpha2;
        cout << delta.abs().maxCoeff() << endl;
        REQUIRE((delta.abs() < 1e-5).all());
    }
}
