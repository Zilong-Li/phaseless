#include "admixture.hpp"
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
                  << "     -b      save parameters in binary file\n"
                  << "     -c      number of ancestral clusters\n"
                  << "     -f      input vcf/bcf format\n"
                  << "     -g      gziped beagle format\n"
                  << "     -i      number of iterations of admixture [100]\n"
                  << "     -I      number of iterations of imputation [40]\n"
                  << "     -k      number of ancestry\n"
                  << "     -n      number of threads\n"
                  << "     -o      output compressed/uncompressed vcf/bcf\n"
                  << "     -q      output estimated ancestry proportion\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      size of each chunk in sites unit\n"
                  << "     -seed   for reproducing results [1]\n"
                  << "     -tol    tolerance value [1e-6]\n"
                  << std::endl;
        return 1;
    }

    std::string in_beagle, in_vcf, out_vcf, out_admixture, out_bin;
    std::string samples = "-", region = "";
    int K{0}, C{0}, niters_admix{100}, niters_impute{40}, nthreads{4}, seed{1};
    int chunksize{100};
    double tol{1e-6};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-b") out_bin = args[++i];
        if(args[i] == "-c") C = stoi(args[++i]);
        if(args[i] == "-k") K = stoi(args[++i]);
        if(args[i] == "-f") in_vcf = args[++i];
        if(args[i] == "-o") out_vcf = args[++i];
        if(args[i] == "-q") out_admixture = args[++i];
        if(args[i] == "-g") in_beagle = args[++i];
        if(args[i] == "-i") niters_admix = stoi(args[++i]);
        if(args[i] == "-I") niters_impute = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") chunksize = stoi(args[++i]);
        if(args[i] == "-seed") seed = stoi(args[++i]);
        if(args[i] == "-tol") tol = stod(args[++i]);
    }
    assert((K > 0) && (C > 0) && (C > K));

    Logger log(out_admixture + ".log");
    log.cao << "Options in effect:\n";
    for(size_t i = 0; i < args.size(); i++) // print out options in effect
    {
        if(i % 2)
            log.cao << args[i] + "\n";
        else
            log.cao << "  " + args[i] + " ";
    }
    Timer tm;
    log.warn(tm.date() + "-> running fastphase\n");
    // log.cao.precision(4);

    // ========= core calculation part ===========================================
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>(chunksize);
    tm.clock();
    chunk_beagle_genotype_likelihoods(genome, in_beagle);
    log.done(tm.date()) << "parsing input -> C:" << C << ", N:" << genome->nsamples << ", M:" << genome->nsnps
                        << ", nchunks:" << genome->nchunks << "; " << tm.reltime() << " ms" << endl;
    auto allthreads = std::thread::hardware_concurrency();
    nthreads = nthreads < genome->nsamples ? nthreads : genome->nsamples;
    nthreads = nthreads < allthreads ? nthreads : allthreads;
    log.done(tm.date()) << allthreads << " concurrent threads are supported. use " << nthreads
                        << " threads\n";
    ThreadPool poolit(nthreads);
    vector<future<double>> llike;

    double loglike{0};
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        FastPhaseK2 nofaith(genome->nsamples, genome->pos[ic].size(), C, seed);
        auto transRate = calc_transRate(genome->pos[ic], C);
        for(int it = 0; it <= niters_impute; it++)
        {
            tm.clock();
            nofaith.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
            {
                if(it == niters_impute)
                    llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                      std::ref(genome->gls[ic]), std::ref(transRate), true));
                else
                    llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                      std::ref(genome->gls[ic]), std::ref(transRate), false));
            }
            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            nofaith.updateClusterFreqPI(tol);
            nofaith.updateAlleleFreqWithinCluster(tol);
            log.done(tm.date()) << setw(2) << "run chunk " << ic << ", iteration " << setw(2) << it
                                << ", log likelihoods: " << std::fixed << loglike << "; " << tm.reltime()
                                << " ms" << endl;
        }
    }

    // write_bcf_genotype_probability(nofaith.GP.data(), out_vcf, in_vcf, sampleids, chrs_pos[ichr], ichr, N,
    // M);
    // log.done(tm.date()) << "imputation done and outputting.\n";

    // log.warn(tm.date() + "-> running admixture\n");
    // Admixture admixer(N, M, C, K, seed);
    // double loglike_prev = 0, diff;
    // for(int it = 0; it <= niters_admix; it++)
    // {
    //     tm.clock();
    //     admixer.initIteration();
    //     for(int i = 0; i < N; i++)
    //         llike.emplace_back(poolit.enqueue(&Admixture::runWithClusterLikelihoods, &admixer, i,
    //                                           std::ref(genolikes), std::ref(transRate),
    //                                           std::ref(nofaith.PI), std::ref(nofaith.F)));
    //     loglike = 0;
    //     for(auto && ll : llike) loglike += ll.get();
    //     llike.clear(); // clear future and renew
    //     diff = loglike - loglike_prev;
    //     log.done(tm.date()) << "iteration " << setw(3) << it << ", log likelihoods: " << std::fixed <<
    //     loglike
    //                         << ", diff=" << diff << ". " << tm.reltime() << " ms" << endl;
    //     if(diff > 0 && diff < 0.1) break;
    //     loglike_prev = loglike;
    //     admixer.updateF();
    // }
    // admixer.writeQ(out_admixture);
    // log.done(tm.date()) << "admixture done and outputting.\n";
    // log.warn(tm.date() + "-> have a nice day, bye!\n");

    return 0;
}
