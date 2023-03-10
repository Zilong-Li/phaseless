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

using pars = std::tuple<double, MyArr2D, MyArr2D>;

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -g beagle.gz -o out -c 10 -n 20\n"
                  << "\nOptions:\n"
                  << "     -b      save all parameters as binary file [1]\n"
                  << "     -c      number of ancestral haplotype clusters\n"
                  << "     -f      input vcf/bcf format\n"
                  << "     -g      gziped beagle format\n"
                  << "     -i      number of iterations of imputation [40]\n"
                  << "     -n      number of threads\n"
                  << "     -o      output prefix\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      size of each chunk in sites unit [100000]\n"
                  << "     -seed   for reproducing results [1]\n"
                  << std::endl;
        return 1;
    }

    filesystem::path out, in_beagle, in_vcf, in_bin;
    std::string samples = "-", region = "";
    int C{0}, niters{40}, nthreads{4}, seed{1}, chunksize{100000};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") out.assign(args[++i]);
        if(args[i] == "-g") in_beagle.assign(args[++i]);
        if(args[i] == "-i") niters = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") chunksize = stoi(args[++i]);
        if(args[i] == "-seed") seed = stoi(args[++i]);
    }
    assert(C > 0);
    // filesystem::create_directories(outdir);

    Logger cao(out.string() + "imputation.log");
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
    vector<future<pars>> llike;
    double loglike, diff, prevlike;
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    tm.clock();
    chunk_beagle_genotype_likelihoods(genome, in_beagle);
    cao.print(tm.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks);
    cao.done(tm.date(), "elapsed time for parsing beagle file", tm.reltime(), " secs");
    auto bw = make_bcfwriter(out.string() + "all.vcf.gz", genome->chrs, genome->sampleids);
    // std::ofstream ofs(out.string() + "all.info");
    // Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 faith(genome->nsamples, genome->pos[ic].size(), C, seed);
        auto transRate = calc_transRate(genome->pos[ic], C);
        for(int it = 0; it <= niters; it++)
        {
            tm.clock();
            faith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == niters)
                    llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                      std::ref(genome->gls[ic]), std::ref(transRate), true));
                else
                    llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwardsHighRam, &faith, i,
                                                      std::ref(genome->gls[ic]), std::ref(transRate), false));
            }
            loglike = 0;
            for(auto && ll : llike)
            {
                const auto [l, iEk, iEkg] = ll.get();
                loglike += l;
                faith.Ek += iEk;
                faith.Ekg += iEkg;
            }
            llike.clear(); // clear future and renew
            diff = it ? loglike - prevlike : 0;
            cao.print(tm.date(), "run chunk", ic, ", iteration", it, ", log likelihoods =", loglike,
                      ", diff =", diff, ",", tm.reltime(), " sec");
            prevlike = loglike;
            // if(diff > 0 && diff < 0.1) break;
            faith.updateIteration();
        }
        tm.clock();
        write_bigass_to_bcf(bw, faith.GP.data(), genome->chrs[ic], genome->pos[ic]);
        genome->transRate.emplace_back(MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
        genome->PI.emplace_back(MyFloat1D(faith.PI.data(), faith.PI.data() + faith.PI.size()));
        genome->F.emplace_back(MyFloat1D(faith.F.data(), faith.F.data() + faith.F.size()));
        cao.done(tm.date(), "chunk", ic, "done. outputting ", tm.reltime(), " secs");
    }
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::ofstream ofbin(out.string() + "pars.bin", std::ios::out | std::ios::binary);
    auto bytes_written = alpaca::serialize<OPTIONS, BigAss>(*genome, ofbin);
    ofbin.close();
    assert(std::filesystem::file_size(out.string() + "pars.bin") == bytes_written);
    cao.done(tm.date(), "imputation done and outputting.", bytes_written, "bytes written to file");

    return 0;
}
