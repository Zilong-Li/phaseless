#include "../src/admixture.hpp"
#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"
#include <filesystem>

using namespace std;
using namespace Eigen;

using pars = std::tuple<MyFloat1D, MyFloat1D, MyFloat1D>;

auto make_input_per_chunk(filesystem::path outdir,
                          const std::unique_ptr<BigAss> & genome,
                          const int ic,
                          const int niters,
                          const int seed)
{
    FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
    auto transRate = calc_transRate(genome->pos[ic], genome->C);
    nofaith.runWithOneThread(niters, genome->gls[ic], transRate);
    write_bcf_genotype_probability(nofaith.GP.data(), genome->chrs[ic], genome->pos[ic], genome->sampleids,
                                   outdir / string("chunk." + to_string(ic) + ".vcf.gz"));
    return std::tuple(MyFloat1D(nofaith.PI.data(), nofaith.PI.data() + nofaith.PI.size()),
                      MyFloat1D(nofaith.F.data(), nofaith.F.data() + nofaith.F.size()),
                      MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
}

TEST_CASE("phaseless naive vs optimal", "[test-phaseless]")
{
    int K{3}, C{5}, seed{1}, nadmix{10}, chunksize{10000}, nphase{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    ThreadPool poolit(genome->nchunks);
    vector<future<pars>> res;
    filesystem::path outdir{"t.phaseless"};
    filesystem::create_directories(outdir);
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), C, seed);
        auto transRate = calc_transRate(genome->pos[ic], genome->C);
        res.emplace_back(poolit.enqueue(make_input_per_chunk, outdir, std::ref(genome), ic, nphase, seed));
    }
    for(auto && ll : res)
    {
        const auto [PI, F, transRate] = ll.get();
        genome->PI.emplace_back(PI);
        genome->F.emplace_back(F);
        genome->transRate.emplace_back(transRate);
    }
    res.clear(); // clear future and renew
    Admixture admixer1(genome->nsamples, genome->nsnps, genome->C, K, seed);
    for(int it = 0; it < nadmix; it++)
    {
        admixer1.initIteration();
        for(int i = 0; i < genome->nsamples; i++) admixer1.runNaiveWithBigAss(i, genome);
        admixer1.updateIteration();
    }
    Admixture admixer2(genome->nsamples, genome->nsnps, genome->C, K, seed);
    for(int it = 0; it < nadmix; it++)
    {
        admixer2.initIteration();
        for(int i = 0; i < genome->nsamples; i++) admixer2.runOptimalWithBigAss(i, genome);
        admixer2.updateIteration();
    }
    REQUIRE(((admixer1.Q - admixer2.Q).abs() < 1e-1).all());
    // REQUIRE(((admixer1.FI - admixer2.FI).abs() > 1e-1).count() == 0);
}

TEST_CASE("phaseless normal iteration with make_input_per_chunk", "[test-phaseless]")
{
    int K{3}, C{5}, seed{1}, nadmix{10}, chunksize{10000}, nphase{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    ThreadPool poolit(genome->nchunks);
    vector<future<pars>> res;
    filesystem::path outdir{"t.phaseless"};
    filesystem::create_directories(outdir);
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), C, seed);
        auto transRate = calc_transRate(genome->pos[ic], genome->C);
        res.emplace_back(poolit.enqueue(make_input_per_chunk, outdir, std::ref(genome), ic, nphase, seed));
    }
    for(auto && ll : res)
    {
        const auto [PI, F, transRate] = ll.get();
        genome->PI.emplace_back(PI);
        genome->F.emplace_back(F);
        genome->transRate.emplace_back(transRate);
    }
    res.clear(); // clear future and renew
    Admixture admixer(genome->nsamples, genome->nsnps, genome->C, K, seed);
    vector<future<double>> llike;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; it < nadmix; it++)
    {
        admixer.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
            llike.emplace_back(poolit.enqueue(&Admixture::runNaiveWithBigAss, &admixer, i, std::ref(genome)));
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        diff = it ? loglike - prevlike : 0;
        REQUIRE(loglike > prevlike);
        if(diff > 0 && diff < 0.1) break;
        prevlike = loglike;
        admixer.updateIteration();
    }
}
