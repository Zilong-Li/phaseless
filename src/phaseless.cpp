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
                  << "     -o      output directory\n"
                  << "     -P      run phasing/imputation only [0]\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      number of sites of each chunk [100000]\n"
                  << "     -seed   for reproducing results [1]\n"
                  << std::endl;
        return 1;
    }

    filesystem::path outdir, in_beagle, in_vcf, in_bin;
    string samples = "-", region = "";
    int chunksize{100000}, accel{1}, phase_only{0}, K{1}, C{0}, nadmix{1000}, nphase{40}, nthreads{4},
        seed{1};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-a") accel = stoi(args[++i]);
        if(args[i] == "-b") in_bin = args[++i];
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-k") K = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") outdir.assign(args[++i]);
        if(args[i] == "-g") in_beagle.assign(args[++i]);
        if(args[i] == "-i") nadmix = stoi(args[++i]);
        if(args[i] == "-I") nphase = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-P") phase_only = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") chunksize = stoi(args[++i]);
        if(args[i] == "-seed") seed = stoi(args[++i]);
    }

    if(!outdir.empty())
    {
        assert((K > 0) && (C > 0) && (C > K));
        filesystem::create_directories(outdir);
    }
    else if(!in_bin.empty())
    {
        assert(K > 0);
        outdir = in_bin.parent_path();
    }
    else
        throw invalid_argument("please check the input and output options\n");

    Logger cao(outdir / string("phaseless.k" + to_string(K) + ".log"));
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
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    std::unique_ptr<BigAss> genome;
    if(in_bin.empty())
    {
        genome = std::make_unique<BigAss>();
        genome->chunksize = chunksize, genome->C = C;
        tm.clock();
        chunk_beagle_genotype_likelihoods(genome, in_beagle);
        cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
                  ", M =", genome->nsnps, ", nchunks =", genome->nchunks);
        cao.done(tm.date(), "elapsed time for parsing beagle file", tm.reltime(), "secs");
        vector<future<pars>> res;
        int ic;
        for(ic = 0; ic < genome->nchunks; ic++)
            res.emplace_back(
                poolit.enqueue(make_input_per_chunk, outdir, std::ref(genome), ic, nphase, seed));
        ic = 0;
        for(auto && ll : res)
        {
            const auto [PI, F, transRate] = ll.get();
            genome->PI.emplace_back(PI);
            genome->F.emplace_back(F);
            genome->transRate.emplace_back(transRate);
            cao.print(tm.date(), "chunk", ic++, " imputation done and outputting.");
        }
        res.clear(); // clear future and renew
        std::ofstream ofs(outdir / "pars.bin", std::ios::out | std::ios::binary);
        auto bytes_written = alpaca::serialize<BigAss>(*genome, ofs);
        ofs.close();
        cao.done(tm.date(), "imputation done and outputting.", bytes_written, "bytes written to file");
    }
    else
    {
        // Deserialize from file
        auto filesize = std::filesystem::file_size(in_bin);
        std::error_code ec;
        std::ifstream ifs(in_bin, std::ios::in | std::ios::binary);
        genome = std::make_unique<BigAss>(alpaca::deserialize<BigAss>(ifs, filesize, ec));
        ifs.close();
        cao.done(tm.date(), filesize, "bytes deserialized from file. skip imputation");
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
        MyArr2D FI0, Q0, FI1, Q1;
        const int istep{4};
        double alpha, stepMax{4}, alphaMax{1280};
        for(int it = 0; it < nadmix; it++)
        {
            tm.clock();
            // first accel iteration
            admixer.initIteration();
            FI0 = admixer.FI;
            Q0 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runWithBigAss, &admixer, i, std::ref(genome)));
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            diff = it ? loglike - prevlike : 0;
            prevlike = loglike;
            cao.print(tm.date(), "SqS3 iteration", it++, ", log likelihoods =", loglike, ", diff =", diff,
                      ",", tm.reltime(), " sec");
            if(diff > 0 && diff < 0.1) break;
            // second accel iteration
            tm.clock();
            admixer.initIteration();
            FI1 = admixer.FI;
            Q1 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runWithBigAss, &admixer, i, std::ref(genome)));
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            diff = it ? loglike - prevlike : 0;
            prevlike = loglike;
            cao.print(tm.date(), "SqS3 iteration", it++, ", log likelihoods =", loglike, ", diff =", diff,
                      ",", tm.reltime(), " sec");
            if(diff > 0 && diff < 0.1) break;
            // accel iteration with steplen
            tm.clock();
            admixer.initIteration();
            alpha =
                ((FI1 - FI0).square().sum() + (Q1 - Q0).square().sum())
                / ((admixer.FI - 2 * FI1 + FI0).square().sum() + (admixer.Q - 2 * Q1 + Q0).square().sum());
            alpha = max(1.0, sqrt(alpha));
            if(alpha >= stepMax)
            {
                alpha = min(stepMax, alphaMax);
                stepMax = min(stepMax * istep, alphaMax);
            }
            admixer.FI = FI0 + 2 * alpha * (FI1 - FI0) + alpha * alpha * (admixer.FI - 2 * FI1 + FI0);
            admixer.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (admixer.Q - 2 * Q1 + Q0);
            admixer.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runWithBigAss, &admixer, i, std::ref(genome)));
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            diff = it ? loglike - prevlike : 0;
            prevlike = loglike;
            cao.print(tm.date(), "SqS3 iteration", it, ", log likelihoods =", loglike, ", diff =", diff, ",",
                      tm.reltime(), " sec");
            if(diff > 0 && diff < 0.1) break;
        }
    }
    else
    {
        for(int it = 0; it < nadmix; it++)
        {
            tm.clock();
            admixer.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(poolit.enqueue(&Admixture::runWithBigAss, &admixer, i, std::ref(genome)));
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            diff = it ? loglike - prevlike : 0;
            cao.print(tm.date(), "iteration", it, ", log likelihoods =", loglike, ", diff =", diff, ",",
                      tm.reltime(), " sec");
            if(diff > 0 && diff < 0.1) break;
            prevlike = loglike;
            admixer.updateIteration();
        }
    }
    cao.done(tm.date(), "admixture done and outputting.");
    admixer.writeQ(outdir / string("admixture.k" + to_string(K) + ".Q"));
    cao.done(tm.date(), "-> good job. have a nice day, bye!");

    return 0;
}
