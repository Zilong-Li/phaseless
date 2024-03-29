#define _DECLARE_TOOLBOX_HERE

#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"
#include "fastphasek4.hpp"

using namespace std;
using namespace Eigen;

TEST_CASE("fastphasek2 forwardAndBackwardsHighRamCollapse", "[test-fastphasek2]")
{
    int C{10}, seed{1}, chunksize{10000}, niters{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    FastPhaseK2 faith(genome->pos[ic].size(), genome->nsamples, genome->C, seed);
    faith.initRecombination(genome->pos[ic]);
    faith.AF = estimate_af_by_gl(genome->gls[ic], genome->nsamples, genome->pos[ic].size()).cast<MyFloat>();
    ThreadPool pool(4);
    vector<future<fbd_res1>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                  std::ref(genome->gls[ic]), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                  std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : llike)
        {
            const auto [l, zj, zg1, zg2, gamma1] = ll.get();
            loglike += l;
            faith.Ezj += zj;
            faith.Ezg1 += zg1;
            faith.Ezg2 += zg2;
            faith.pi += gamma1;
        }
        llike.clear(); // clear future and renew
        // REQUIRE(loglike > prevlike);
        if(it != niters) faith.updateIteration();
        prevlike = loglike;
    }
    // start collapsing
    faith.collapse_and_resize(genome->pos[ic]);
    prevlike = std::numeric_limits<double>::lowest();
    niters = 2;
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                  std::ref(genome->gls[ic]), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                  std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : llike)
        {
            const auto [l, zj, zg1, zg2, gamma1] = ll.get();
            loglike += l;
            faith.Ezj += zj;
            faith.Ezg1 += zg1;
            faith.Ezg2 += zg2;
            faith.pi += gamma1;
        }
        llike.clear(); // clear future and renew
        // if(it > 1) REQUIRE(loglike > prevlike);
        prevlike = loglike;
        if(it != niters) faith.updateIteration();
    }
}

TEST_CASE("fastphasek2 forwardAndBackwardsLowRamNormal", "[test-fastphasek2]")
{
    int C{10}, seed{1}, chunksize{10000}, niters{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    FastPhaseK2 faith(genome->pos[ic].size(), genome->nsamples, genome->C, seed);
    faith.initRecombination(genome->pos[ic]);
    faith.AF = estimate_af_by_gl(genome->gls[ic], genome->nsamples, genome->pos[ic].size()).cast<MyFloat>();
    ThreadPool pool(4);
    vector<future<double>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &faith, i,
                                                  std::ref(genome->gls[ic]), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &faith, i,
                                                  std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        // REQUIRE(loglike > prevlike);
        if(it != niters) faith.updateIteration();
        prevlike = loglike;
    }
}

TEST_CASE("fastphasek2 forwardAndBackwardsLowRamCollapse", "[test-fastphasek2]")
{
    int C{10}, seed{999}, chunksize{10000}, niters{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    FastPhaseK2 faith(genome->pos[ic].size(), genome->nsamples, genome->C, seed);
    faith.initRecombination(genome->pos[ic]);
    faith.AF = estimate_af_by_gl(genome->gls[ic], genome->nsamples, genome->pos[ic].size()).cast<MyFloat>();
    ThreadPool pool(4);
    vector<future<double>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &faith, i,
                                                  std::ref(genome->gls[ic]), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &faith, i,
                                                  std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        // REQUIRE(loglike > prevlike);
        prevlike = loglike;
        if(it != niters) faith.updateIteration();
    }
    // start collapsing
    faith.collapse_and_resize(genome->pos[ic]);
    prevlike = std::numeric_limits<double>::lowest();
    niters = 2;
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &faith, i,
                                                  std::ref(genome->gls[ic]), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsLowRam, &faith, i,
                                                  std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        // if(it > 1) REQUIRE(loglike > prevlike);
        prevlike = loglike;
        if(it != niters) faith.updateIteration();
    }
}

TEST_CASE("fastphasek2 forwardAndBackwardsHighRamNormal", "[test-fastphasek2]")
{
    int C{10}, seed{1}, chunksize{10000}, niters{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    FastPhaseK2 faith(genome->pos[ic].size(), genome->nsamples, genome->C, seed);
    faith.initRecombination(genome->pos[ic]);
    faith.AF = estimate_af_by_gl(genome->gls[ic], genome->nsamples, genome->pos[ic].size()).cast<MyFloat>();
    ThreadPool pool(4);
    vector<future<fbd_res1>> llike;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                  std::ref(genome->gls[ic]), true));
            else
                llike.emplace_back(pool.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                  std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : llike)
        {
            const auto [l, zj, zg1, zg2, gamma1] = ll.get();
            loglike += l;
            faith.Ezj += zj;
            faith.Ezg1 += zg1;
            faith.Ezg2 += zg2;
            faith.pi += gamma1;
        }
        llike.clear(); // clear future and renew
        // REQUIRE(loglike > prevlike);
        if(it != niters) faith.updateIteration();
        prevlike = loglike;
    }
}

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

