/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/phaseless.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "admixture.hpp"
#include "cli.hpp"
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

int run_bootstrap(const std::unique_ptr<BigAss> & genome, ThreadPool & poolit, Logger & cao, int K, int nseeds)
{
    std::vector<std::future<double>> res;
    std::vector<double> llikes;
    for(int seed = 0; seed < nseeds; seed++)
    {
        Admixture admixer(genome->nsamples, genome->nsnps, genome->C, K, seed);
        for(int i = 0; i < genome->nsamples; i++)
            res.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
        double loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        llikes.push_back(loglike);
    }
    int seed = 0;
    for(auto ll : llikes) cao.print("seed =", seed++, ", log likelihoods =", ll);
    auto it = std::max_element(llikes.begin(), llikes.end());
    return std::distance(llikes.begin(), it);
}

auto make_input_per_chunk(const std::unique_ptr<BigAss> & genome,
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

void filter_input_per_chunk(filesystem::path out,
                            const std::unique_ptr<BigAss> & genome,
                            const int ic,
                            const int niters,
                            const int seed,
                            const double info)
{
    FastPhaseK2 faith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
    auto transRate = calc_transRate(genome->pos[ic], genome->C);
    faith.runWithOneThread(niters, genome->gls[ic], transRate);
    auto idx2rm = filter_sites_per_chunk(faith.GP.data(), info, genome->nsamples, genome->pos[ic].size());
    thin_bigass_per_chunk(ic, idx2rm, genome);
}

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ===========================
    //
    auto opts = parsecli(argc, argv);
    Logger cao(opts.out.string() + "log", !opts.noscreen);
    cao.print(opts.opts_in_effect);
    Timer tm;
    cao.warn(tm.date(), "-> running fastphase");
    int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < allthreads ? opts.nthreads : allthreads;
    cao.print(tm.date(), allthreads, " concurrent threads are supported. use", opts.nthreads, " threads");

    // ======================== core calculation part ===========================
    ThreadPool poolit(opts.nthreads);
    std::unique_ptr<BigAss> genome;
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;

    if(opts.run_impute)
    {
        genome = std::make_unique<BigAss>();
        genome->chunksize = opts.chunksize, genome->C = opts.C;
        tm.clock();
        chunk_beagle_genotype_likelihoods(genome, opts.in_beagle);
        cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
                  ", M =", genome->nsnps, ", nchunks =", genome->nchunks);
        cao.done(tm.date(), "elapsed time for parsing beagle file", std::fixed, tm.reltime(), " secs");
        if(opts.info > 0)
        {
            // TODO implement this at some point
            return 1;
            cao.warn(tm.date(), "--info", opts.info, " is applied, which will do two rounds imputation");
            vector<future<void>> res;
            for(int ic = 0; ic < genome->nchunks; ic++)
                res.emplace_back(poolit.enqueue(filter_input_per_chunk, opts.out, std::ref(genome), ic,
                                                opts.nimpute, opts.seed, opts.info));
            for(auto && ll : res) ll.get();
            update_bigass_inplace(genome);
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
                            res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith,
                                                            i, std::ref(genome->gls[ic]), std::ref(transRate),
                                                            true));
                        else
                            res.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith,
                                                            i, std::ref(genome->gls[ic]), std::ref(transRate),
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
                    cao.print(tm.date(), "run single chunk", ic, ", iteration", it,
                              ", log likelihoods =", loglike, ",", tm.reltime(), " sec");
                    faith.updateIteration();
                }
                tm.clock();
                write_bigass_to_bcf(bw, faith.GP.data(), genome->chrs[ic], genome->pos[ic]);
                genome->transRate.emplace_back(
                    MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
                genome->PI.emplace_back(MyFloat1D(faith.PI.data(), faith.PI.data() + faith.PI.size()));
                genome->F.emplace_back(MyFloat1D(faith.F.data(), faith.F.data() + faith.F.size()));
                auto ClusterInfo = calc_cluster_info(faith.N, faith.GZP1, faith.GZP2);
                oinfo << ClusterInfo.format(fmt) << "\n";
                opi << faith.PI.transpose().format(fmt) << "\n";
                cao.done(tm.date(), "chunk", ic, "done. outputting ", tm.reltime(), " secs");
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
                genome->transRate.emplace_back(
                    MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
                genome->PI.emplace_back(MyFloat1D(PI.data(), PI.data() + PI.size()));
                genome->F.emplace_back(MyFloat1D(F.data(), F.data() + F.size()));
                oinfo << ClusterInfo.format(fmt) << "\n";
                opi << PI.transpose().format(fmt) << "\n";
                cao.print(tm.date(), "chunk", ic++, " imputation done and outputting");
            }
        }
        std::ofstream ofs(opts.out.string() + "pars.bin", std::ios::out | std::ios::binary);
        auto bytes_written = alpaca::serialize<OPTIONS, BigAss>(*genome, ofs);
        ofs.close();
        assert(std::filesystem::file_size(opts.out.string() + "pars.bin") == bytes_written);
        cao.done(tm.date(), "imputation done and outputting.", bytes_written, " bytes written to file");
        return 0;
    }

    if(opts.run_admix)
    { // Deserialize from file
        auto filesize = std::filesystem::file_size(opts.in_bin);
        std::error_code ec;
        std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
        genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
        ifs.close();
        assert((bool)ec == false);
        cao.done(tm.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
        cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
                  ", M =", genome->nsnps, ", nchunks =", genome->nchunks);
        assert(opts.K < genome->C);

        cao.warn(tm.date(), "-> running admixture with seed =", opts.seed);
        Admixture admixer(genome->nsamples, genome->nsnps, genome->C, opts.K, opts.seed);
        vector<future<double>> llike;
        if(!opts.noaccel)
        {
            MyArr2D F0, Q0, F1, Q1;
            const int istep{4};
            double alpha, diff, stepMax{4}, alphaMax{1280};
            for(int it = 0; it < opts.nadmix / 3; it++)
            {
                // first accel iteration
                admixer.initIteration();
                F0 = admixer.F;
                Q0 = admixer.Q;
                for(int i = 0; i < genome->nsamples; i++)
                    llike.emplace_back(
                        poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
                for(auto && ll : llike) ll.get();
                llike.clear(); // clear future and renew
                admixer.updateIteration();
                // second accel iteration
                admixer.initIteration();
                F1 = admixer.F;
                Q1 = admixer.Q;
                diff = (Q1 - Q0).square().sum();
                tm.clock();
                for(int i = 0; i < genome->nsamples; i++)
                    llike.emplace_back(
                        poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
                double loglike = 0;
                for(auto && ll : llike) loglike += ll.get();
                llike.clear(); // clear future and renew
                admixer.updateIteration();
                cao.print(tm.date(), "SqS3 iteration", it * 3 + 1, ", diff(Q) =", std::scientific, diff,
                          ", likelihoods =", std::fixed, loglike, ",", tm.reltime(), " sec");
                if(diff < opts.qtol)
                {
                    cao.print(tm.date(), "hit stopping criteria, diff(Q) =", std::scientific, diff, " <",
                              opts.qtol);
                    break;
                }
                // accel iteration with steplen
                admixer.initIteration();
                alpha =
                    ((F1 - F0).square().sum() + (Q1 - Q0).square().sum())
                    / ((admixer.F - 2 * F1 + F0).square().sum() + (admixer.Q - 2 * Q1 + Q0).square().sum());
                alpha = max(1.0, sqrt(alpha));
                if(alpha >= stepMax)
                {
                    alpha = min(stepMax, alphaMax);
                    stepMax = min(stepMax * istep, alphaMax);
                }
                admixer.F = F0 + 2 * alpha * (F1 - F0) + alpha * alpha * (admixer.F - 2 * F1 + F0);
                admixer.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (admixer.Q - 2 * Q1 + Q0);
                admixer.initIteration();
                for(int i = 0; i < genome->nsamples; i++)
                    llike.emplace_back(
                        poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
                for(auto && ll : llike) ll.get();
                llike.clear(); // clear future and renew
                admixer.updateIteration();
            }
        }
        else
        {
            MyArr2D Q0;
            double diff, loglike;
            for(int it = 0; it < opts.nadmix; it++)
            {
                tm.clock();
                admixer.initIteration();
                Q0 = admixer.Q;
                for(int i = 0; i < genome->nsamples; i++)
                    llike.emplace_back(
                        poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));

                loglike = 0;
                for(auto && ll : llike) loglike += ll.get();
                llike.clear(); // clear future and renew
                admixer.updateIteration();
                diff = (admixer.Q - Q0).square().sum();
                cao.print(tm.date(), "normal iteration", it, ", diff(Q) =", std::scientific, diff,
                          ", likelihoods =", std::fixed, loglike, ",", tm.reltime(), " sec");
                if(diff < opts.qtol) break;
            }
        }
        cao.done(tm.date(), "admixture done and outputting");
        admixer.writeQ(opts.out.string() + "Q");
        // if(admixer.debug) admixer.writeBin(out.string() + "qf.bin", genome);
        cao.done(tm.date(), "-> good job. have a nice day, bye!");

        return 0;
    }

    if(opts.run_pars)
    {
        auto filesize = std::filesystem::file_size(opts.in_bin);
        std::error_code ec;
        std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
        genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
        ifs.close();
        assert((bool)ec == false);
        const int ic = opts.ichunk;
        if(ic >= genome->nchunks)
        {
            cao.error("the chunk to extract", ic, "(0-baded) is not less than total chunks", genome->nchunks);
            return 1;
        }
        const int iM = genome->pos[ic].size();
        std::ofstream ofs(opts.out.string() + "haplike.bin", std::ios::binary);
        ofs.write((char *)&genome->C, 4);
        ofs.write((char *)&genome->nsamples, 4);
        ofs.write((char *)&iM, 4);
        // haplike is p(X_is ,Z_is | \theta) = alpha * beta
        MyArr2D alpha, beta;
        for(int ind = 0; ind < genome->nsamples; ind++)
        {
            alpha.setZero(genome->C * genome->C, iM);
            beta.setZero(genome->C * genome->C, iM);
            getClusterLikelihoods(ind, alpha, beta, genome->gls[ic], genome->transRate[ic], genome->PI[ic],
                                  genome->F[ic]);
            alpha *= beta;
            ofs.write((char *)alpha.data(), genome->C * genome->C * iM * 4);
        }
        return 0;
    }
}
