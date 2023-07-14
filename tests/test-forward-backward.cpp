#define _DECLARE_TOOLBOX_HERE

#include "../src/fastphase.hpp"
#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("reconstruct alpha and beta from saved pars.bin", "[test-forward-backward]")
{
    int C{10}, seed{1}, chunksize{10000}, niters{40};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    ThreadPool poolit(4);
    vector<future<fbd_res1>> res;
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 faith(genome->pos[ic].size(), genome->nsamples, C, seed);
        faith.initRecombination(genome->pos[ic]);
        faith.AF =
            estimate_af_by_gl(genome->gls[ic], genome->nsamples, genome->pos[ic].size()).cast<MyFloat>();
        faith.debug = true;
        double prevlike{std::numeric_limits<double>::lowest()}, loglike;
        for(int it = 0; it <= niters; it++)
        {
            faith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == niters)
                    res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                    std::ref(genome->gls[ic]), true));
                else
                    res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
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
            REQUIRE(loglike > prevlike); // might only work for double not float precision
            prevlike = loglike;
            if(it != niters) faith.updateIteration();
        }
        MyArr2D alpha, beta;
        const int iM = genome->pos[ic].size();
        for(int ind = 0; ind < genome->nsamples; ind++)
        {
            alpha.setZero(genome->C * genome->C, iM);
            beta.setZero(genome->C * genome->C, iM);
            Eigen::Map<const MyArr2D> gli(genome->gls[ic].data() + ind * iM * 3, iM, 3);
            MyArr2D emit = get_emission_by_gl(gli, faith.F).transpose(); // C2 x M
            forward_backwards_diploid(alpha, beta, emit, faith.R, faith.PI);
            alpha *= beta;
            REQUIRE(((alpha.colwise().sum() - 1.0).abs() < 1e-6).all());
        }
    }
}
