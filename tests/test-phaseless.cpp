#define _DECLARE_TOOLBOX_HERE

#include "../src/admixture.hpp"
#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"
#include <filesystem>

using namespace std;
using namespace Eigen;

TEST_CASE("phaseless naive vs dump dataset 1", "[test-phaseless]")
{
    int K{3}, C{10}, seed{1}, nadmix{10}, chunksize{10000}, nimpute{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    ThreadPool pool(genome->nchunks);
    vector<future<double>> res;
    filesystem::path outdir{"test.dir"};
    filesystem::create_directories(outdir);
    FastPhaseK2 faith(genome->nsamples, genome->nsnps, C, seed);
    for(int it = 0; it <= nimpute; it++)
    {
        faith.initIteration();
        for(int i = 0; i < faith.N; i++)
            res.emplace_back(pool.enqueue(&FastPhaseK2::runAllChunks, &faith, std::ref(genome->gls), i, false));
        double loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        faith.updateIteration();
    }
    double llike1, llike2;
    Admixture admixer1(genome->nsamples, genome->nsnps, genome->C, K, seed);
    Admixture admixer2(genome->nsamples, genome->nsnps, genome->C, K, seed);
    for(int it = 0; it < nadmix; it++)
    {
        admixer1.initIteration();
        llike1 = 0;
        for(int i = 0; i < genome->nsamples; i++) llike1 += admixer1.runOptimalWithBigAss(i, genome);
        admixer1.updateIteration();
        admixer2.initIteration();
        llike2 = 0;
        for(int i = 0; i < genome->nsamples; i++) llike2 += admixer2.runNativeWithBigAss(i, genome);
        admixer2.updateIteration();
        cerr << std::setprecision(6) << "llike2: " << llike2 << "\tllike1: " << llike1 << "\n";
        REQUIRE(abs(llike1 - llike2) < 1e-10);
    }
    REQUIRE(((admixer1.Q - admixer2.Q).abs() < 1e-6).all());
}
