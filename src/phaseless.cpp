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
                  << "     -b      binary file with clusters likelihood\n"
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
                  << "     -s      samples in vcf/bcf to subset\n"
                  << "     -seed   for reproducing results [1]\n"
                  << "     -tol    tolerance value [1e-6]\n"
                  << std::endl;
        return 1;
    }

    std::string in_beagle, in_vcf, out_vcf, out_admixture, out_cluster;
    std::string samples = "-", region = "";
    int K{0}, C{0}, niters_admix{100}, niters_impute{40}, nthreads{4}, seed{1};
    double tol{1e-6};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-b") out_cluster = args[++i];
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
        if(args[i] == "-s") samples = args[++i];
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
    int N, M;
    DoubleVec1D genolikes;
    StringIntVecMapU chrs_map;
    StringIntMapU chrs_starts;
    StringVec1D sampleids;
    std::string ichr;
    tm.clock();
    read_beagle_genotype_likelihoods(in_beagle, genolikes, sampleids, chrs_map, chrs_starts, N, M);
    log.done(tm.date()) << "parsing input -> N:" << N << ", M:" << M << ", C:" << C << "; " << tm.reltime()
                        << " ms" << endl;
    assert(chrs_map.size() == 1);
    ichr = chrs_map.begin()->first;
    auto transRate = calc_transRate(chrs_map[ichr], C);
    nthreads = nthreads < N ? nthreads : N;

    double loglike{0};
    FastPhaseK2 nofaith(N, M, C, seed);
    if(!out_cluster.empty()) nofaith.openClusterFile(out_cluster);
    ThreadPool poolit(nthreads);
    vector<future<double>> llike;
    for(int it = 0; it <= niters_impute; it++)
    {
        tm.clock();
        nofaith.initIteration();
        for(int i = 0; i < N; i++)
        {
            if(it == niters_impute)
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate), false));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        nofaith.updateClusterFreqPI(tol);
        nofaith.updateAlleleFreqWithinCluster(tol);
        log.done(tm.date()) << "iteration " << setw(2) << it << ", log likelihoods: " << std::fixed << loglike
                            << "; " << tm.reltime() << " ms" << endl;
    }
    write_bcf_genotype_probability(nofaith.GP.data(), out_vcf, in_vcf, sampleids, chrs_map[ichr], ichr, N, M);
    log.done(tm.date()) << "imputation done and outputting.\n";

    log.warn(tm.date() + "-> running admixture\n");
    Admixture admixer(N, M, C, K, seed);
    double loglike_prev = 0, diff;
    for(int it = 0; it <= niters_admix; it++)
    {
        tm.clock();
        admixer.initIteration();
        for(int i = 0; i < N; i++)
            llike.emplace_back(poolit.enqueue(&Admixture::runWithClusterLikelihoods, &admixer, i,
                                              std::ref(genolikes), std::ref(transRate), std::ref(nofaith.PI),
                                              std::ref(nofaith.F)));
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        diff = loglike - loglike_prev;
        log.done(tm.date()) << "iteration " << setw(3) << it << ", log likelihoods: " << std::fixed << loglike
                            << ", diff=" << diff << ". " << tm.reltime() << " ms" << endl;
        if(diff > 0 && diff < 0.1) break;
        loglike_prev = loglike;
        admixer.updateF();
    }
    admixer.writeQ(out_admixture);
    log.done(tm.date()) << "admixture done and outputting.\n";
    log.warn(tm.date() + "-> have a nice day, bye!\n");

    return 0;
}
