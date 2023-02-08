#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("fastphasek4", "[test-fastphasek4]")
{
    double tol{1e-6};
    int N, M, C{5}, seed{1}, niters{5};
    MyFloat1D genolikes;
    uMapStringInt1D chrs_map;
    uMapStringUint chrs_starts;
    StringVec1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_map, chrs_starts, N, M);
    auto ichr = chrs_map.begin()->first;
    auto distRate = calc_distRate(chrs_map[ichr], C);
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
