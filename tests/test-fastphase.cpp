#define _DECLARE_TOOLBOX_HERE

#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"
#include "fastphasek4.hpp"

using namespace std;
using namespace Eigen;


TEST_CASE("fastphasek4", "[test-fastphasek4]")
{
    int C{10}, seed{1}, chunksize{10000}, niters{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    auto distRate = calc_er(genome->pos[ic], 4.0 * 20000 / C);
    MyArr2D postProbsZ(genome->nsnps, C * C);
    MyArr2D postProbsZandG(genome->nsnps, C * C * 4);
    FastPhaseK4 faith4(genome->nsamples, genome->nsnps, C, seed);
    ThreadPool pool(4);
    vector<future<double>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike, tol{1e-4};
    for(int it = 0; it <= niters; it++)
    {
        faith4.transitionCurIter(distRate);
        postProbsZ.setZero();
        postProbsZandG.setZero();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK4::forwardAndBackwards, &faith4, i,
                                                std::ref(genome->gls[ic]), std::ref(postProbsZ),
                                                std::ref(postProbsZandG), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK4::forwardAndBackwards, &faith4, i,
                                                std::ref(genome->gls[ic]), std::ref(postProbsZ),
                                                std::ref(postProbsZandG), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        faith4.updateClusterFreqPI(postProbsZ, tol);
        faith4.updateAlleleFreqWithinCluster(postProbsZandG, tol);
        // REQUIRE(loglike > prevlike);
        prevlike = loglike;
    }
}

