#define _DECLARE_TOOLBOX_HERE

#include "../src/collapse.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("FastPhaseGrid forwardAndBackwardsHighRam", "[test-collapse]")
{
    int C{10}, seed{1}, niters{20}, chunksize{10000}, B{5};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    auto [pos, central] = turn_pos_into_grid(genome->pos[ic], B);
    FastPhaseGrid faith(pos, central, genome->pos[ic].size(), genome->nsamples, C, seed);
    faith.debug = true;
    using pars1 = std::tuple<double, MyArr2D, MyArr2D, MyArr2D, MyArr1D>;
    vector<future<pars1>> res;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    ThreadPool poolit(4);
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                res.emplace_back(poolit.enqueue(&FastPhaseGrid::forwardAndBackwardsHighRam, &faith, i,
                                                std::ref(genome->gls[ic]), true));
            else
                res.emplace_back(poolit.enqueue(&FastPhaseGrid::forwardAndBackwardsHighRam, &faith, i,
                                                std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : res)
        {
            const auto [l, zj, zg1, zg2, gamma1] = ll.get();
            loglike += l;
            faith.Ezj += zj;
            faith.Ezg1 += zg1;
            faith.Ezg2 += zg2;
            faith.pi += gamma1;
        }
        res.clear(); // clear future and renew
        cao.cerr(tim.date(), "run single chunk", ic, ", iteration", it, ", likelihoods =", loglike);
        REQUIRE(loglike > prevlike);
        prevlike = loglike;
        if(it != niters) faith.updateIteration();
    }
}
