#include "fastphase.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"
#include <alpaca/alpaca.h>

using namespace std;
using namespace vcfpp;

using pars1 = std::tuple<double, MyArr2D, MyArr2D>;
using pars2 = std::tuple<MyFloat1D, MyArr2D, MyArr2D, MyArr2D, MyArr2D>;

inline auto make_input_per_chunk(const std::unique_ptr<BigAss> & genome,
                                 const int ic,
                                 const int niters,
                                 const int seed)
{
    FastPhaseK2 faith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
    auto transRate = calc_transRate(genome->pos[ic], genome->C);
    faith.runWithOneThread(niters, genome->gls[ic], transRate);
    auto Info = calc_cluster_info(faith.N, faith.GZP1, faith.GZP2);
    return std::tuple(MyFloat1D(faith.GP.data(), faith.GP.data() + faith.GP.size()), Info, transRate,
                      faith.PI, faith.F);
}

// inline void filter_input_per_chunk(filesystem::path out,
//                                    const std::unique_ptr<BigAss> & genome,
//                                    const int ic,
//                                    const int niters,
//                                    const int seed,
//                                    const double info)
// {
//     FastPhaseK2 faith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
//     auto transRate = calc_transRate(genome->pos[ic], genome->C);
//     faith.runWithOneThread(niters, genome->gls[ic], transRate);
//     auto idx2rm = filter_sites_per_chunk(faith.GP.data(), info, genome->nsamples, genome->pos[ic].size());
//     thin_bigass_per_chunk(ic, idx2rm, genome);
// }

inline int run_impute_main(Options & opts)
{
    Logger cao(opts.out.string() + "log", !opts.noscreen);
    cao.print(opts.opts_in_effect);
    Timer tm;
    cao.warn(tm.date(), "-> running fastphase");
    int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < allthreads ? opts.nthreads : allthreads;
    cao.print(tm.date(), allthreads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool poolit(opts.nthreads);

    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = opts.chunksize, genome->C = opts.C;
    tm.clock();
    chunk_beagle_genotype_likelihoods(genome, opts.in_beagle);
    cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks, ", seed =", opts.seed);
    cao.done(tm.date(), "elapsed time for parsing beagle file", std::fixed, tm.reltime(), " secs");
    if(opts.info > 0)
    {
        // TODO implement this at some point
        return 1;
        // cao.warn(tm.date(), "--info", opts.info, " is applied, which will do two rounds imputation");
        // vector<future<void>> res;
        // for(int ic = 0; ic < genome->nchunks; ic++)
        //     res.emplace_back(poolit.enqueue(filter_input_per_chunk, opts.out, std::ref(genome), ic,
        //                                     opts.nimpute, opts.seed, opts.info));
        // for(auto && ll : res) ll.get();
        // update_bigass_inplace(genome);
    }
    auto bw = make_bcfwriter(opts.out.string() + "all.vcf.gz", genome->chrs, genome->sampleids);
    std::ofstream oinfo(opts.out.string() + "all.info");
    std::ofstream opi(opts.out.string() + "all.pi");
    Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
    if(opts.single_chunk)
    {
        vector<future<pars1>> res;
        for(int ic = 0; ic < genome->nchunks; ic++)
        {
            FastPhaseK2 faith(genome->nsamples, genome->pos[ic].size(), opts.C, opts.seed);
            auto transRate = calc_transRate(genome->pos[ic], opts.C);
            for(int it = 0; it <= opts.nimpute; it++)
            {
                tm.clock();
                faith.initIteration();
                for(int i = 0; i < genome->nsamples; i++)
                {
                    if(it == opts.nimpute)
                        res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                        std::ref(genome->gls[ic]), std::ref(transRate),
                                                        true));
                    else
                        res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                        std::ref(genome->gls[ic]), std::ref(transRate),
                                                        false));
                }
                double loglike = 0;
                for(auto && ll : res)
                {
                    const auto [l, iEk, iEkg] = ll.get();
                    loglike += l;
                    faith.Ek += iEk;
                    faith.Ekg += iEkg;
                }
                res.clear(); // clear future and renew
                cao.print(tm.date(), "run single chunk", ic, ", iteration", it, ", likelihoods =", loglike,
                          ",", tm.reltime(), " sec");
                faith.updateIteration();
            }
            tm.clock();
            write_bigass_to_bcf(bw, faith.GP.data(), genome->chrs[ic], genome->pos[ic]);
            genome->transRate.emplace_back(MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
            genome->PI.emplace_back(MyFloat1D(faith.PI.data(), faith.PI.data() + faith.PI.size()));
            genome->F.emplace_back(MyFloat1D(faith.F.data(), faith.F.data() + faith.F.size()));
            auto ClusterInfo = calc_cluster_info(faith.N, faith.GZP1, faith.GZP2);
            oinfo << ClusterInfo.format(fmt) << "\n";
            opi << faith.PI.transpose().format(fmt) << "\n";
            cao.done(tm.date(), "chunk", ic, " done. outputting elapsed", tm.reltime(), " secs");
        }
    }
    else
    {
        if(genome->nchunks < opts.nthreads)
            cao.warn(tm.date(), "nchunks < nthreads. only", genome->nchunks, " threads will be working");
        vector<future<pars2>> res;
        for(int ic = 0; ic < genome->nchunks; ic++)
            res.emplace_back(
                poolit.enqueue(make_input_per_chunk, std::ref(genome), ic, opts.nimpute, opts.seed));
        int ic = 0;
        for(auto && ll : res)
        {
            const auto [GP, ClusterInfo, transRate, PI, F] = ll.get();
            write_bigass_to_bcf(bw, GP.data(), genome->chrs[ic], genome->pos[ic]);
            genome->transRate.emplace_back(MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
            genome->PI.emplace_back(MyFloat1D(PI.data(), PI.data() + PI.size()));
            genome->F.emplace_back(MyFloat1D(F.data(), F.data() + F.size()));
            oinfo << ClusterInfo.format(fmt) << "\n";
            opi << PI.transpose().format(fmt) << "\n";
            cao.print(tm.date(), "chunk", ic++, " imputation done and outputting");
        }
    }
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::ofstream ofs(opts.out.string() + "pars.bin", std::ios::out | std::ios::binary);
    auto bytes_written = alpaca::serialize<OPTIONS, BigAss>(*genome, ofs);
    ofs.close();
    assert(std::filesystem::file_size(opts.out.string() + "pars.bin") == bytes_written);
    cao.done(tm.date(), "imputation done and outputting.", bytes_written, " bytes written to file");
    return 0;
}
