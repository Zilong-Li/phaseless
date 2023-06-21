#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

double make_input_per_chunk(int niters, int ic, const std::unique_ptr<BigAss> & genome, int seed)
{
    FastPhaseK2 nofaith(genome->pos[ic], genome->nsamples, genome->C, seed);
    auto transRate = calc_transRate(genome->pos[ic], genome->C);
    return nofaith.runWithOneThread(niters, genome->gls[ic]);
}

TEST_CASE("fastphasek2 forwardAndBackwardsLowRam", "[test-fastphasek2]")
{
    int N, M, C{5}, seed{1}, niters{5};
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    String1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_pos, N, M);
    auto ichr = chrs_pos.begin()->first;
    FastPhaseK2 nofaith(chrs_pos.begin()->second, N, C, seed);
    ThreadPool poolit(4);
    vector<future<double>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it < niters + 1; it++)
    {
        nofaith.initIteration();
        for(int i = 0; i < N; i++)
        {
            if(it == niters)
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &nofaith, i,
                                                  std::ref(genolikes), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &nofaith, i,
                                                  std::ref(genolikes), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        REQUIRE(loglike > prevlike);
        nofaith.updateIteration();
        prevlike = loglike;
    }
}

TEST_CASE("fastphasek2 forwardAndBackwardsHighRam", "[test-fastphasek2]")
{
    int N, M, C{5}, seed{1}, niters{5};
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    String1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_pos, N, M);
    auto ichr = chrs_pos.begin()->first;
    FastPhaseK2 nofaith(chrs_pos.begin()->second, N, C, seed);
    ThreadPool poolit(4);
    using pars = std::tuple<double, MyArr2D, MyArr2D>;
    vector<future<pars>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it < niters + 1; it++)
    {
        nofaith.initIteration();
        for(int i = 0; i < N; i++)
        {
            if(it == niters)
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &nofaith, i,
                                                  std::ref(genolikes), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &nofaith, i,
                                                  std::ref(genolikes), false));
        }
        loglike = 0;
        for(auto && ll : llike)
        {
            const auto [l, iEk, iEkg] = ll.get();
            loglike += l;
            nofaith.Ek += iEk;
            nofaith.Ekg += iEkg;
        }
        llike.clear(); // clear future and renew
        REQUIRE(loglike > prevlike);
        nofaith.updateIteration();
        prevlike = loglike;
    }
}

TEST_CASE("fastphasek2 runWithOneThread", "[test-fastphasek2]")
{
    int C{5}, seed{1}, niters{10}, chunksize{10000};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 nofaith(genome->pos[ic], genome->nsamples, C, seed);
        auto transRate = calc_transRate(genome->pos[ic], genome->C);
        cout << "diff: " << nofaith.runWithOneThread(niters, genome->gls[ic]) << endl;
    }
}

TEST_CASE("fastphasek2 runWithOneThread in threadpool", "[test-fastphasek2]")
{
    int C{5}, seed{1}, niters{10}, chunksize{10000};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    ThreadPool poolit(genome->nchunks);
    vector<future<double>> diff;
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        diff.emplace_back(poolit.enqueue(make_input_per_chunk, niters, ic, std::ref(genome), seed));
    }
    for(auto && ll : diff) cout << "diff: " << ll.get() << endl;
}

TEST_CASE("fastphasek4", "[test-fastphasek4]")
{
    double tol{1e-6};
    int N, M, C{5}, seed{1}, niters{5};
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    String1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_pos, N, M);
    auto ichr = chrs_pos.begin()->first;
    auto distRate = calc_distRate(chrs_pos[ichr], C);
    MyArr2D postProbsZ(M, C * C);
    MyArr2D postProbsZandG(M, C * C * 4);

    FastPhaseK4 nofaith(N, M, C, seed);
    ThreadPool poolit(4);
    vector<future<double>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it < niters + 1; it++)
    {
        nofaith.transitionCurIter(distRate);
        postProbsZ.setZero();
        postProbsZandG.setZero();
        for(int i = 0; i < N; i++)
        {
            if(it == niters)
                llike.emplace_back(poolit.enqueue(&FastPhaseK4::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(postProbsZ),
                                                  std::ref(postProbsZandG), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK4::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(postProbsZ),
                                                  std::ref(postProbsZandG), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        nofaith.updateClusterFreqPI(postProbsZ, tol);
        nofaith.updateAlleleFreqWithinCluster(postProbsZandG, tol);
        REQUIRE(loglike > prevlike);
        prevlike = loglike;
    }
}
