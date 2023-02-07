#include "admixture.hpp"
#include "fastphase.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"

using namespace std;

// check initialize_sigmaCurrent_m in STITCH
ArrDouble2D calc_transRate(const IntVec1D & markers, int C, int Ne = 20000, double expRate = 0.5)
{
    ArrDouble1D distRate(markers.size());
    distRate(0) = exp(-1e20);
    // int nGen = 4 * Ne / C;
    for(size_t i = 1; i < markers.size(); i++) distRate(i) = exp(-(markers[i] - markers[i - 1]) / 1e6);
    ArrDouble2D transRate(3, markers.size());
    transRate.row(0) = distRate.square();
    transRate.row(1) = distRate * (1 - distRate);
    transRate.row(2) = (1 - distRate).square();
    return transRate;
}

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
                  << "     -i      number of iterations\n"
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
    int K{0}, C{0}, niters{40}, nthreads{4}, seed{1};
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
        if(args[i] == "-i") niters = stoi(args[++i]);
        if(args[i] == "-n") nthreads = stoi(args[++i]);
        if(args[i] == "-r") region = args[++i];
        if(args[i] == "-s") samples = args[++i];
        if(args[i] == "-seed") seed = stoi(args[++i]);
        if(args[i] == "-tol") tol = stod(args[++i]);
    }
    assert((K > 0) && (C > 0) && (C > K));

    Logger log(out_vcf + ".log");
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
    log.cao.precision(4);

    // ========= core calculation part ===========================================
    int N, M;
    DoubleVec1D genolikes;
    IntVec1D markers;
    StringVec1D sampleids, chrs;

    tm.clock();
    read_beagle_genotype_likelihoods(in_beagle, genolikes, sampleids, chrs, markers, N, M);
    log.done(tm.date()) << "parsing input -> N:" << N << ", M:" << M << ", C:" << C << "; " << tm.reltime()
                        << " ms" << endl;
    auto transRate = calc_transRate(markers, C);
    nthreads = nthreads < N ? nthreads : N;

    double loglike{0};
    FastPhaseK2 nofaith(N, M, C, seed);
    if(!out_cluster.empty()) nofaith.openClusterFile(out_cluster);
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
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        nofaith.updateClusterFreqPI(tol);
        nofaith.updateAlleleFreqWithinCluster(tol);
        log.done(tm.date()) << "iteration " << setw(2) << it << ", log likelihoods: " << std::fixed << loglike
                            << "; " << tm.reltime() << " ms" << endl;
    }
    write_bcf_genotype_probability(nofaith.GP.data(), out_vcf, in_vcf, sampleids, markers, chrs[0], N, M);
    log.done(tm.date()) << "imputation done and outputting.\n";
    log.warn(tm.date() + "-> running admixture\n");
    Admixture admixer(N, M, C, K, seed);
    for(int it = 0; it < niters + 1; it++)
    {
        admixer.initIteration();
        for(int i = 0; i < N; i++)
        {
            llike.emplace_back(poolit.enqueue(&Admixture::runWithClusterLikelihoods, &admixer, i,
                                              std::ref(genolikes), std::ref(transRate), std::ref(nofaith.PI),
                                              std::ref(nofaith.F)));
        }
        loglike = 0;
        for(auto && ll : llike) loglike += ll.get();
        llike.clear(); // clear future and renew
        admixer.updateF();
        log.done(tm.date()) << "iteration " << setw(2) << it << ", log likelihoods: " << std::fixed << loglike
                            << "; " << tm.reltime() << " ms" << endl;
    }
    admixer.writeQ(out_admixture);
    log.done(tm.date()) << "admixture done and outputting.\n";
    log.warn(tm.date() + "-> have a nice day, bye!\n");

    return 0;
}
