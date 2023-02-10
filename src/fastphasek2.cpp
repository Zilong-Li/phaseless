#include "fastphase.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0]
                         + " -g beagle.gz -c 10 -i 20 -n 10 -o out.vcf.gz -b cluster.bin\n"
                  << "\nOptions:\n"
                  << "     -b      binary file with clusters likelihood\n"
                  << "     -c      number of ancestral clusters\n"
                  << "     -f      input vcf/bcf format\n"
                  << "     -g      gziped beagle format\n"
                  << "     -i      number of iterations of imputation\n"
                  << "     -n      number of threads\n"
                  << "     -o      output compressed/uncompressed vcf/bcf\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      samples in vcf/bcf to subset\n"
                  << "     -seed   for reproducing results [1]\n"
                  << std::endl;
        return 1;
    }

    std::string in_beagle, in_vcf, out_vcf, out_cluster;
    std::string samples = "-", region = "";
    int C{0}, niters{40}, nthreads{4}, seed{1};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-b") out_cluster = args[++i];
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") out_vcf = args[++i];
        if(args[i] == "-g") in_beagle = args[++i];
        if(args[i] == "-i") niters = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") samples = args[++i];
        if(args[i] == "-seed") seed = stoi(args[++i]);
    }
    assert(C > 0);

    Logger cao(out_vcf + ".log");
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
    int allthreads = std::thread::hardware_concurrency();
    // nthreads = nthreads < genome->nsamples ? nthreads : genome->nsamples;
    nthreads = nthreads < allthreads ? nthreads : allthreads;
    cao.print(tm.date(), allthreads, " concurrent threads are supported. use", nthreads, " threads");

    // ========= core calculation part ===========================================
    int N, M;
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    StringVec1D sampleids;
    tm.clock();
    read_beagle_genotype_likelihoods(in_beagle, genolikes, sampleids, chrs_pos, N, M);
    cao.print(tm.date(), "parsing input -> C =", C, ", N =", N, ", M =", M);
    assert(chrs_pos.size() == 1);
    auto ichr = chrs_pos.begin()->first;
    auto transRate = calc_transRate(chrs_pos.begin()->second, C);

    double loglike{0}, prevlike, diff;
    FastPhaseK2 nofaith(N, M, C, seed);
    ThreadPool poolit(nthreads);
    vector<future<double>> llike;
    for(int it = 0; it < niters + 1; it++)
    {
        tm.clock();
        nofaith.initIteration();
        for(int i = 0; i < N; i++)
        {
            if(it == niters)
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate), false));
        }
        loglike = 0;
        for(auto && l : llike) loglike += l.get();
        llike.clear(); // clear future and renew
        if(it)
            diff = loglike - prevlike;
        else
            diff = 0;
        cao.print(tm.date(), "iteration", it, ", log likelihoods =", loglike, ", diff =", diff, ",",
                  tm.reltime(), " sec");
        prevlike = loglike;
        // if(diff > 0 && diff < 0.1) break;
        nofaith.updateIteration();
    }
    cao.done(tm.date(), "imputation done and outputting.");
    write_bcf_genotype_probability(nofaith.GP.data(), ichr, chrs_pos.begin()->second, sampleids, out_vcf);
    cao.done(tm.date(), "-> good job. have a nice day, bye!");

    return 0;
}
