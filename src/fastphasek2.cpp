/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/fastphasek2.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "alpaca/alpaca.h"
#include "fastphase.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"
#include <filesystem>

using namespace std;

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -g beagle.gz -o dir -c 10 -n 20\n"
                  << "\nOptions:\n"
                  << "     -b      save all parameters as binary file [1]\n"
                  << "     -c      number of ancestral haplotype clusters\n"
                  << "     -f      input vcf/bcf format\n"
                  << "     -g      gziped beagle format\n"
                  << "     -i      number of iterations of imputation [40]\n"
                  << "     -n      number of threads\n"
                  << "     -o      output directory\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      size of each chunk in sites unit [100000]\n"
                  << "     -seed   for reproducing results [1]\n"
                  << std::endl;
        return 1;
    }

    filesystem::path outdir, in_beagle, in_vcf, in_bin;
    std::string samples = "-", region = "";
    int C{0}, niters{40}, nthreads{4}, seed{1}, chunksize{100000};
    double info{0};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") outdir.assign(args[++i]);
        if(args[i] == "-g") in_beagle.assign(args[++i]);
        if(args[i] == "-i") niters = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") chunksize = stoi(args[++i]);
        if(args[i] == "-seed") seed = stoi(args[++i]);
    }
    assert(C > 0);
    filesystem::create_directories(outdir);

    Logger cao(outdir / "imputation.log");
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
    vector<future<double>> llike;
    double loglike, diff, prevlike;
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    tm.clock();
    chunk_beagle_genotype_likelihoods(genome, in_beagle);
    cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks);
    cao.done(tm.date(), "elapsed time for parsing beagle file", tm.reltime(), " secs");
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), C, seed);
        auto transRate = calc_transRate(genome->pos[ic], C);
        for(int it = 0; it <= niters; it++)
        {
            tm.clock();
            nofaith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == niters)
                    llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                      std::ref(genome->gls[ic]), std::ref(transRate), true));
                else
                    llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                      std::ref(genome->gls[ic]), std::ref(transRate), false));
            }
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            diff = it ? loglike - prevlike : 0;
            cao.print(tm.date(), "run chunk", ic, ", iteration", it, ", log likelihoods =", loglike,
                      ", diff =", diff, ",", tm.reltime(), " sec");
            prevlike = loglike;
            // if(diff > 0 && diff < 0.1) break;
            nofaith.updateIteration();
        }
        auto idx2rm = write_bcf_genotype_probability(
            nofaith.GP.data(), genome->chrs[ic], genome->pos[ic], genome->sampleids,
            outdir / string("chunk." + to_string(ic) + ".vcf.gz"), info);
        thin_bigass(ic, idx2rm, genome, nofaith.PI, nofaith.F, transRate);
    }
    std::ofstream ofs(outdir / "pars.bin", std::ios::out | std::ios::binary);
    auto bytes_written = alpaca::serialize<BigAss>(*genome, ofs);
    ofs.close();
    cao.done(tm.date(), "imputation done and outputting.", bytes_written, "bytes written to file");

    return 0;
}
