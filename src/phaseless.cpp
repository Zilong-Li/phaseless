#include "admixture.hpp"
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
                  << "     " + (std::string)argv[0] + " -g beagle.gz -o dir -c 10 -k 3 -i 100 -n 10\n"
                  << "\nOptions:\n"
                  << "     -b      input binary file with all parameters\n"
                  << "     -c      number of ancestral clusters\n"
                  << "     -f      input vcf/bcf format\n"
                  << "     -g      gziped beagle format\n"
                  << "     -i      number of iterations of admixture [100]\n"
                  << "     -I      number of iterations of imputation [40]\n"
                  << "     -k      number of ancestry\n"
                  << "     -n      number of threads\n"
                  << "     -o      output directory\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      size of each chunk in sites unit\n"
                  << "     -seed   for reproducing results [1]\n"
                  << std::endl;
        return 1;
    }

    filesystem::path outdir, in_beagle, in_vcf, in_bin;
    string samples = "-", region = "";
    int K{0}, C{0}, niters_admix{100}, niters_impute{40}, nthreads{4}, seed{1};
    int chunksize{100000};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-b") in_bin = args[++i];
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-k") K = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") outdir.assign(args[++i]);
        if(args[i] == "-g") in_beagle.assign(args[++i]);
        if(args[i] == "-i") niters_admix = stoi(args[++i]);
        if(args[i] == "-I") niters_impute = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") chunksize = stoi(args[++i]);
        if(args[i] == "-seed") seed = stoi(args[++i]);
    }
    assert((K > 0) && (C > 0) && (C > K));

    if(!outdir.empty())
        filesystem::create_directories(outdir);
    else if(!in_bin.empty())
        outdir = in_bin.parent_path();
    else
        throw invalid_argument("please check the input and output options\n");

    Logger cao(outdir / "phaseless.log");
    cao.cao << "Options in effect:\n";
    for(size_t i = 0; i < args.size(); i++) // print out options in effect
    {
        if(i % 2)
            cao.cao << args[i] + "\n";
        else
            cao.cao << "  " + args[i] + " ";
    }
    Timer tm;
    cao.warn(tm.date(), "-> running fastphase");
    cout.flags(std::ios::fixed | std::ios::right);

    // ========= core calculation part ===========================================
    int allthreads = std::thread::hardware_concurrency();
    nthreads = nthreads < allthreads ? nthreads : allthreads;
    cao.print(tm.date(), allthreads, " concurrent threads are supported. use", nthreads, " threads");
    ThreadPool poolit(nthreads);
    vector<future<double>> llike;
    double loglike{0}, diff, prevlike;

    std::unique_ptr<BigAss> genome;
    if(in_bin.empty())
    {
        genome = std::make_unique<BigAss>();
        genome->chunksize = chunksize;
        tm.clock();
        chunk_beagle_genotype_likelihoods(genome, in_beagle);
        cao.print(tm.date(), "parsing input -> C =", C, ", N =", genome->nsamples, ", M =", genome->nsnps,
                  ", nchunks =", genome->nchunks);
        cao.done(tm.date(), "elapsed time for parsing beagle file", tm.reltime(), " secs");
        for(int ic = 0; ic < genome->nchunks; ic++)
        {
            FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), C, seed);
            auto transRate = calc_transRate(genome->pos[ic], C);
            genome->transRate.emplace_back(MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
            for(int it = 0; it <= niters_impute; it++)
            {
                tm.clock();
                nofaith.initIteration();
                for(int i = 0; i < genome->nsamples; i++)
                {
                    if(it == niters_impute)
                        llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                          std::ref(genome->gls[ic]), std::ref(transRate),
                                                          true));
                    else
                        llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                          std::ref(genome->gls[ic]), std::ref(transRate),
                                                          false));
                }
                loglike = 0;
                for(auto && ll : llike) loglike += ll.get();
                llike.clear(); // clear future and renew
                if(it)
                    diff = loglike - prevlike;
                else
                    diff = 0;
                cao.print(tm.date(), "run chunk", ic, ", iteration", it, ", log likelihoods =", loglike,
                          ", diff =", diff, ",", tm.reltime(), " sec");
                prevlike = loglike;
                // if(diff > 0 && diff < 0.1) break;
                nofaith.updateIteration();
            }
            write_bcf_genotype_probability(nofaith.GP.data(), genome->chrs[ic], genome->pos[ic],
                                           genome->sampleids,
                                           outdir / string("chunk." + to_string(ic) + ".vcf.gz"));
            genome->PI.emplace_back(MyFloat1D(nofaith.PI.data(), nofaith.PI.data() + nofaith.PI.size()));
            genome->F.emplace_back(MyFloat1D(nofaith.F.data(), nofaith.F.data() + nofaith.F.size()));
        }
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
        cao.done(tm.date(), filesize, "bytes deserialized from file. skip imputation");
        genome = std::make_unique<BigAss>(alpaca::deserialize<BigAss>(ifs, filesize, ec));
        ifs.close();
        cao.print(tm.date(), "parsing input -> C =", C, ", N =", genome->nsamples, ", M =", genome->nsnps,
                  ", nchunks =", genome->nchunks);
    }

    cao.warn(tm.date(), "-> running admixture");
    Admixture admixer(genome->nsamples, genome->nsnps, C, K, seed);
    for(int it = 0, prevlike = 0; it <= niters_admix; it++)
    {
        tm.clock();
        admixer.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
            llike.emplace_back(poolit.enqueue(&Admixture::runWithBigAss, &admixer, i, std::ref(genome)));
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        if(it)
            diff = loglike - prevlike;
        else
            diff = 0;
        cao.print(tm.date(), "iteration", it, ", log likelihoods =", loglike, ", diff =", diff, ",",
                  tm.reltime(), " sec");
        if(diff > 0 && diff < 0.1) break;
        prevlike = loglike;
        admixer.updateIteration();
    }
    cao.done(tm.date(), "admixture done and outputting.");
    admixer.writeQ(outdir / "admixture.q");
    cao.done(tm.date(), "-> good job. have a nice day, bye!");

    return 0;
}
