#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("fastphasek2 runWithOneThread", "[test-fastphasek2]")
{
    int C{5}, seed{1}, niters{10}, chunksize{10000};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), C, seed);
        nofaith.runWithOneThread(niters, genome->gls[ic], genome->pos[ic]);
    }
}

// TEST_CASE("fastphasek2 runWithOneThread in threadpool", "[test-fastphasek2]")
// {
//     int C{5}, seed{1}, niters{10}, chunksize{10000};
//     std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
//     genome->chunksize = chunksize, genome->C = C;
//     chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
//     ThreadPool poolit(genome->nchunks);
//     vector<future<double>> diff;
//     for(int ic = 0; ic < genome->nchunks; ic++)
//     {
//         diff.emplace_back(poolit.enqueue(&FastPhaseK2::runWithOneThread,
//                                          FastPhaseK2(genome->nsamples, genome->pos[ic].size(), C, seed),
//                                          niters, std::ref(genome->gls[ic]), std::ref(genome->pos[ic])));
//     }
//     for(auto && ll : diff) ll.get();
// }

TEST_CASE("fastphasek2 forwardAndBackwards", "[test-fastphasek2]")
{
    int N, M, C{5}, seed{1}, niters{5};
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    StringVec1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_pos, N, M);
    auto ichr = chrs_pos.begin()->first;
    auto transRate = calc_transRate(chrs_pos.begin()->second, C);

    FastPhaseK2 nofaith(N, M, C, seed);
    ThreadPool poolit(4);
    vector<future<double>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it < niters + 1; it++)
    {
        nofaith.initIteration();
        for(int i = 0; i < N; i++)
        {
            if(it == niters)
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        REQUIRE(loglike > prevlike);
        nofaith.updateIteration();
        prevlike = loglike;
    }
}
