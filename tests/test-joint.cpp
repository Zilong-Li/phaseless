#include <climits>
#define _DECLARE_TOOLBOX_HERE

#include "../src/io.hpp"
#include "../src/phaseless.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("phaseless joint single chunk", "[test-joint]")
{
    int K{3}, C{5}, seed{1}, nadmix{10}, chunksize{INT_MAX}, nimpute{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    Phaseless faith(K, C, genome->nsamples, genome->pos[ic].size(), seed);
    faith.initRecombination(genome->pos[ic], 20000, 1);
    faith.debug = 1;
    for(int it = 0; it <= nimpute; it++)
    {
        tim.clock();
        faith.initIteration();
        double loglike = 0;
        for(int i = 0; i < genome->nsamples; i++)
            loglike += faith.runForwardBackwards(i, genome->gls[ic], false);
        cao.print(tim.date(), "run single chunk", ic, ", iteration", it,
                  ", likelihoods =", loglike, ",", tim.reltime(), " sec");
        if(it != nimpute) faith.updateIteration();
    }
}
