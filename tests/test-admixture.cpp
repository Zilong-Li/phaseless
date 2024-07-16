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
    Options opts;
    opts.C = C, opts.chunksize = chunksize, opts.gridsize = 1, opts.seed = seed;
    opts.in_beagle = "../data/bgl.gz";
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    init_bigass(genome, opts);
    ThreadPool pool(4);
    vector<future<double>> res;
    FastPhaseK2 faith(genome->nsamples, genome->nsnps, C, seed);
    faith.initRecombination(genome->pos, opts.in_rfile);
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
    // reuse Ezj for AE
    faith.Ezj = get_cluster_frequency(faith.R, faith.PI);
    genome->collapse = Char1D(faith.collapse.data(), faith.collapse.data() + faith.collapse.size());
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        const int S = faith.pos_chunk[ic + 1] - faith.pos_chunk[ic];
        const int G = faith.grid_chunk[ic + 1] - faith.grid_chunk[ic];
        MyArr2D out = faith.Ezj.middleCols(faith.pos_chunk[ic], G);
        genome->AE.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.R.middleCols(faith.pos_chunk[ic], G);
        genome->R.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.PI.middleCols(faith.pos_chunk[ic], G);
        genome->PI.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
        out = faith.P.middleRows(faith.pos_chunk[ic], S);
        genome->P.emplace_back(MyFloat1D(out.data(), out.data() + out.size()));
    }
    double llike1, llike2;
    Admixture admixer1(genome->nsamples, genome->nsnps, genome->C, K, seed);
    admixer1.setStartPoint(genome, opts.in_qfile);
    Admixture admixer2(genome->nsamples, genome->nsnps, genome->C, K, seed);
    admixer2.setStartPoint(genome, opts.in_qfile);
    for(int it = 0; it < nadmix; it++)
    {
        llike1 = 0;
        admixer1.initIteration();
        for(int i = 0; i < genome->nsamples; i++) llike1 += admixer1.runOptimalWithBigAss(i, genome);
        admixer1.updateIteration();
        llike2 = 0;
        admixer2.initIteration();
        for(int i = 0; i < genome->nsamples; i++) llike2 += admixer2.runNativeWithBigAss(i, genome);
        admixer2.updateIteration();
        cerr << std::setprecision(6) << "llike2: " << llike2 << "\tllike1: " << llike1 << "\n";
        REQUIRE(abs(llike1 - llike2) < 1e-6);
    }
    REQUIRE(((admixer1.Q - admixer2.Q).abs() < 1e-6).all());
}
