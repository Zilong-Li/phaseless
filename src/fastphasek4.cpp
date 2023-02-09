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
                  << "     -i      number of iterations\n"
                  << "     -n      number of threads\n"
                  << "     -o      output compressed/uncompressed vcf/bcf\n"
                  << "     -r      region in vcf/bcf to subset\n"
                  << "     -s      samples in vcf/bcf to subset\n"
                  << "     -seed   for reproducing results [1]\n"
                  << "     -tol    tolerance value [1e-6]\n"
                  << std::endl;
        return 1;
    }

    std::string in_beagle, in_vcf, out_vcf, out_cluster;
    std::string samples = "-", region = "";
    int C{0}, niters{40}, nthreads{4}, seed{1};
    double tol{1e-6};
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
        if(args[i] == "-tol") tol = stod(args[++i]);
    }
    assert(C > 0);

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
    log.warn(tm.date() + "-> program start\n");
    // log.cao.precision(4);

    // ========= core calculation part ===========================================
    int N, M;
    MyFloat1D genolikes;
    StringVec1D sampleids;
    MapStringInt1D chrs_pos;
    tm.clock();
    read_beagle_genotype_likelihoods(in_beagle, genolikes, sampleids, chrs_pos, N, M);
    log.done(tm.date()) << "parsing input -> N:" << N << ", M:" << M << ", C:" << C << "; " << tm.reltime()
                        << " ms" << endl;
    assert(chrs_pos.size() == 1);
    auto ichr = chrs_pos.begin()->first;
    auto distRate = calc_distRate(chrs_pos.begin()->second, C);

    double loglike{0};
    MyArr2D postProbsZ(M, C * C);
    MyArr2D postProbsZandG(M, C * C * 4);

    FastPhaseK4 nofaith(N, M, C, seed);
    if(!out_cluster.empty()) nofaith.openClusterFile(out_cluster);
    nthreads = nthreads < N ? nthreads : N;
    ThreadPool poolit(nthreads);
    vector<future<double>> llike;
    for(int it = 0; it < niters + 1; it++)
    {
        tm.clock();
        nofaith.transitionCurIter(distRate);
        postProbsZ.setZero();
        postProbsZandG.setZero();
        for(int i = 0; i < N; i++)
        {
            // loglike += nofaith.forwardAndBackwards(i, genolikes, postProbsZ, postProbsZandG);
            if(it == niters)
                llike.emplace_back(poolit.enqueue(&FastPhaseK4::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(postProbsZ),
                                                  std::ref(postProbsZandG), true));
            else
                llike.emplace_back(poolit.enqueue(&FastPhaseK4::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(postProbsZ),
                                                  std::ref(postProbsZandG), false));
        }
        loglike = 0;
        for(auto && l : llike) loglike += l.get();
        llike.clear(); // clear future and renew
        nofaith.updateClusterFreqPI(postProbsZ, tol);
        nofaith.updateAlleleFreqWithinCluster(postProbsZandG, tol);
        log.done(tm.date()) << "iteration " << setw(2) << it << ", log likelihoods: " << std::fixed << loglike
                            << "; " << tm.reltime() << " ms" << endl;
    }
    write_bcf_genotype_probability(nofaith.GP.data(), ichr, chrs_pos.begin()->second, sampleids, out_vcf);
    log.warn(tm.date() + "-> have a nice day, bye!\n");

    return 0;
}
