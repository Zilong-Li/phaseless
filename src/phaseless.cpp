/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/phaseless.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "admixture.hpp"
#include "alpaca/alpaca.h"
#include "fastphase.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"
#include <filesystem>

using namespace std;

using pars = std::tuple<MyFloat1D, MyFloat1D, MyFloat1D>;
using pars2 = std::tuple<int, MyFloat1D, MyFloat1D, MyFloat1D>;

auto filter_input_per_chunk(filesystem::path out,
                            const std::unique_ptr<BigAss> & genome,
                            const int ic,
                            const int niters,
                            const int seed,
                            const double info)
{
    int nsnp2rm{0};
    {
        FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
        auto transRate = calc_transRate(genome->pos[ic], genome->C);
        nofaith.runWithOneThread(niters, genome->gls[ic], transRate);
        auto idx2rm = write_bcf_genotype_probability(
            nofaith.GP.data(), genome->chrs[ic], genome->pos[ic], genome->sampleids,
            out.string() + string("chunk." + to_string(ic) + ".vcf.gz"), info);
        nsnp2rm = thin_bigass_per_chunk(ic, idx2rm, genome);
    }
    FastPhaseK2 faith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
    auto transRate = calc_transRate(genome->pos[ic], genome->C);
    faith.runWithOneThread(niters, genome->gls[ic], transRate);
    return std::tuple(nsnp2rm, MyFloat1D(faith.PI.data(), faith.PI.data() + faith.PI.size()),
                      MyFloat1D(faith.F.data(), faith.F.data() + faith.F.size()),
                      MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
}

auto make_input_per_chunk(filesystem::path out,
                          const std::unique_ptr<BigAss> & genome,
                          const int ic,
                          const int niters,
                          const int seed)
{
    FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), genome->C, seed);
    auto transRate = calc_transRate(genome->pos[ic], genome->C);
    nofaith.runWithOneThread(niters, genome->gls[ic], transRate);
    write_bcf_genotype_probability(nofaith.GP.data(), genome->chrs[ic], genome->pos[ic], genome->sampleids,
                                   out.string() + string("chunk." + to_string(ic) + ".vcf.gz"));
    return std::tuple(MyFloat1D(nofaith.PI.data(), nofaith.PI.data() + nofaith.PI.size()),
                      MyFloat1D(nofaith.F.data(), nofaith.F.data() + nofaith.F.size()),
                      MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
}

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -g beagle.gz -o dir -c 10 -k 3 -n 20\n"
                  << "\nOptions:\n"
                  << "     -a      accelerated EM with SqS3 scheme [1]\n"
                  << "     -b      input binary file with all parameters\n"
                  << "     -c      number of ancestral haplotype clusters\n"
                  << "     -f      input vcf/bcf format\n"
                  << "     -g      gziped beagle format\n"
                  << "     -i      maximum iterations of admixture [1000]\n"
                  << "     -I      maximum iterations of imputation [40]\n"
                  << "     -k      number of ancestry in admixture model\n"
                  << "     -n      number of threads\n"
                  << "     -o      output prefix\n"
                  << "     -p      print out log to screen [0]\n"
                  << "     -P      run phasing/imputation only [0]\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      number of sites of each chunk [100000]\n"
                  << "     -info   filter and re-impute sites with low info [0]\n"
                  << "     -qtol   tolerance of stopping criteria [1e-6]\n"
                  << "     -seed   for reproducing results [1]\n"
                  << std::endl;
        return 1;
    }

    filesystem::path out, in_beagle, in_vcf, in_bin;
    string samples = "-", region = "";
    int chunksize{100000}, accel{1}, phase_only{0}, K{1}, C{0}, nadmix{1000}, nphase{40}, nthreads{4};
    int seed{1}, isscreen{0};
    double qtol{1e-6}, diff, info{0};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-a") accel = stoi(args[++i]);
        if(args[i] == "-b") in_bin = args[++i];
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-k") K = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") out.assign(args[++i]);
        if(args[i] == "-g") in_beagle.assign(args[++i]);
        if(args[i] == "-i") nadmix = stoi(args[++i]);
        if(args[i] == "-I") nphase = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-p") isscreen = stoi(args[++i]);
        if(args[i] == "-P") phase_only = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") chunksize = stoi(args[++i]);
        if(args[i] == "-info") info = stod(args[++i]);
        if(args[i] == "-qtol") qtol = stod(args[++i]);
        if(args[i] == "-seed") seed = stoi(args[++i]);
    }

    Logger cao(out.string() + "phaseless.log");
    cao.is_screen = isscreen;
    cao.cao << "Options in effect:\n";
    for(size_t i = 0; i < args.size(); i++) // print out options in effect
        i % 2 ? cao.cao << args[i] + "\n" : cao.cao << "  " + args[i] + " ";
    Timer tm;
    cao.warn(tm.date(), "-> running fastphase");
    int allthreads = std::thread::hardware_concurrency();
    nthreads = nthreads < allthreads ? nthreads : allthreads;
    cao.print(tm.date(), allthreads, " concurrent threads are supported. use", nthreads, " threads");

    // ========= core calculation part ===========================================
    ThreadPool poolit(nthreads);
    std::unique_ptr<BigAss> genome;
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    if(in_bin.empty())
    {
        genome = std::make_unique<BigAss>();
        genome->chunksize = chunksize, genome->C = C;
        tm.clock();
        chunk_beagle_genotype_likelihoods(genome, in_beagle, true);
        cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
                  ", M =", genome->nsnps, ", nchunks =", genome->nchunks);
        cao.done(tm.date(), "elapsed time for parsing beagle file", std::fixed, tm.reltime(), " secs");
        if(info > 0)
        {
            cao.warn(tm.date(), "-info", info, " is applied, which will do two rounds imputation");
            vector<future<pars2>> res;
            for(int ic = 0; ic < genome->nchunks; ic++)
                res.emplace_back(
                    poolit.enqueue(filter_input_per_chunk, out, std::ref(genome), ic, nphase, seed, info));
            int ic = 0;
            for(auto && ll : res)
            {
                const auto [n, PI, F, transRate] = ll.get();
                genome->nsnps -= n;
                genome->PI.emplace_back(PI);
                genome->F.emplace_back(F);
                genome->transRate.emplace_back(transRate);
                cao.print(tm.date(), "chunk", ic++, " imputation done and outputting");
            }
        }
        else
        {
            vector<future<pars>> res;
            for(int ic = 0; ic < genome->nchunks; ic++)
                res.emplace_back(
                    poolit.enqueue(make_input_per_chunk, out, std::ref(genome), ic, nphase, seed));
            int ic = 0;
            for(auto && ll : res)
            {
                const auto [PI, F, transRate] = ll.get();
                genome->PI.emplace_back(PI);
                genome->F.emplace_back(F);
                genome->transRate.emplace_back(transRate);
                cao.print(tm.date(), "chunk", ic++, " imputation done and outputting");
            }
        }
        std::ofstream ofs(out.string() + "pars.bin", std::ios::out | std::ios::binary);
        auto bytes_written = alpaca::serialize<OPTIONS, BigAss>(*genome, ofs);
        ofs.close();
        assert(std::filesystem::file_size(out.string() + "pars.bin") == bytes_written);
        cao.done(tm.date(), "imputation done and outputting.", bytes_written, " bytes written to file");
    }
    else
    { // Deserialize from file
        auto filesize = std::filesystem::file_size(in_bin);
        std::error_code ec;
        std::ifstream ifs(in_bin, std::ios::in | std::ios::binary);
        genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
        ifs.close();
        assert((bool)ec == false);
        cao.done(tm.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
        cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
                  ", M =", genome->nsnps, ", nchunks =", genome->nchunks);
        assert(K < genome->C);
    }
    if(phase_only) return 0;

    cao.warn(tm.date(), "-> running admixture");
    Admixture admixer(genome->nsamples, genome->nsnps, genome->C, K, seed);
    vector<future<double>> llike;
    if(accel)
    {
        MyArr2D F0, Q0, F1, Q1;
        const int istep{4};
        double alpha, stepMax{4}, alphaMax{1280};
        for(int it = 0; it < nadmix / 3; it++)
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
            if(diff < qtol) break;
            tm.clock();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(
                    poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            for(auto && ll : llike) ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            cao.print(tm.date(), "SqS3 iteration", it * 3 + 1, ", diff(Q) =", std::scientific, diff, ",",
                      std::fixed, tm.reltime(), " sec");
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
        for(int it = 0; it < nadmix; it++)
        {
            tm.clock();
            admixer.initIteration();
            Q0 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(
                    poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            for(auto && ll : llike) ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            diff = (admixer.Q - Q0).square().sum();
            cao.print(tm.date(), "normal iteration", it, ", diff(Q) =", std::scientific, diff, ",",
                      std::fixed, tm.reltime(), " sec");
            if(diff < qtol) break;
        }
    }
    cao.done(tm.date(), "admixture done and outputting");
    admixer.writeQ(out.string() + "admixture.Q");
    cao.done(tm.date(), "-> good job. have a nice day, bye!");

    return 0;
}
