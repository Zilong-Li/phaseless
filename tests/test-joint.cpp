#define _DECLARE_TOOLBOX_HERE

#include "../src/io.hpp"
#include "../src/phaseless.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("phaseless joint single chunk", "[test-joint]")
{
    cao.cerr("TEST: phaseless joint single chunk");
    int K{3}, C{5}, seed{1}, chunksize{INT_MAX}, nimpute{10};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    Phaseless faith(K, C, genome->nsamples, genome->nsnps, seed);
    faith.initRecombination(genome->pos, 20000, 1);
    ThreadPool pool(4);
    vector<future<double>> res;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; it <= nimpute; it++)
    {
        tim.clock();
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
            res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
        loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        diff = it ? loglike - prevlike : 0;
        res.clear(); // clear future and renew
        cao.cerr(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff = ", diff, ", time",
                 tim.reltime(), " sec");
        faith.updateIteration();
        prevlike = loglike;
    }
}

TEST_CASE("phaseless joint multiple chunks", "[test-joint]")
{
    cao.cerr("TEST: phaseless joint multiple chunks");
    int K{3}, C{5}, seed{1}, chunksize{INT_MAX}, nimpute{10};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    Phaseless faith(K, C, genome->nsamples, genome->nsnps, seed);
    faith.initRecombination(genome->pos, 20000, 1);
    ThreadPool pool(4);
    vector<future<double>> res;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; it <= nimpute; it++)
    {
        tim.clock();
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
            res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
        loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        diff = it ? loglike - prevlike : 0;
        cao.cerr(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff = ", diff, ", time",
                 tim.reltime(), " sec");
        faith.updateIteration();
        prevlike = loglike;
    }
}
