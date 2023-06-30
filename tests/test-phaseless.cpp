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
    FastPhaseK2 faith(genome->pos[ic], genome->nsamples, genome->C, seed);
    faith.runWithOneThread(niters, genome->gls[ic]);
    write_bcf_genotype_probability(faith.GP.data(), genome->chrs[ic], genome->pos[ic], genome->sampleids,
                                   outdir / string("chunk." + to_string(ic) + ".vcf.gz"));
    std::cerr << faith.R << "\n";
    return std::tuple(MyFloat1D(faith.PI.data(), faith.PI.data() + faith.PI.size()),
                      MyFloat1D(faith.F.data(), faith.F.data() + faith.F.size()),
                      MyFloat1D(faith.R.data(), faith.R.data() + faith.R.size()));
}

TEST_CASE("phaseless naive vs dump", "[test-phaseless]")
{
    int K{3}, C{10}, seed{1}, nadmix{10}, chunksize{10000}, nphase{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    ThreadPool poolit(genome->nchunks);
    vector<future<pars>> res;
    filesystem::path outdir{"t.phaseless"};
    filesystem::create_directories(outdir);
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 faith(genome->pos[ic], genome->nsamples, C, seed);
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
        cerr << std::setprecision(6) << llike2 << "\t" << llike1 << "\n";
        REQUIRE(abs(llike1 - llike2) < 1e-2);
    }
    REQUIRE(((admixer1.Q - admixer2.Q).abs() < 1e-6).all());
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
        FastPhaseK2 faith(genome->pos[ic], genome->nsamples, C, seed);
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
            llike.emplace_back(poolit.enqueue(&Admixture::runNativeWithBigAss, &admixer, i, std::ref(genome)));
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
