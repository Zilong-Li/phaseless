// -*- compile-command: "g++ fastphasek4.cpp -o fastphasek4 -std=c++17 -g -O3 -Wall -lz && ./fastphasek4 -f ../data/gl.vcf.gz"; -*-
#include "fastphase.hpp"
#include "io.hpp"
#include "threadpool.hpp"

using namespace std;


// check initialize_sigmaCurrent_m in STITCH
ArrDouble2D calc_transRate(const IntVec1D& markers, int C, int Ne = 20000, double expRate = 0.5)
{
    ArrDouble1D distRate(markers.size());
    distRate(0) = exp(-1e20);
    // int nGen = 4 * Ne / C;
    for (size_t i = 1; i < markers.size(); i++)
        distRate(i) = exp(-(markers[i] - markers[i - 1]) / 1e6);
    ArrDouble2D transRate(3, markers.size());
    transRate.row(0) = distRate.square();
    transRate.row(1) = distRate * (1 - distRate);
    transRate.row(2) = (1 - distRate).square();
    return transRate;
}

int main(int argc, char* argv[])
{
    // ========= helper message and parameters parsing ============================
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] +
                         " -g beagle.gz -c 10 -i 20 -n 10 -o out.vcf.gz -b cluster.bin\n"
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

    std::string cluster_out, beagle_in = "", vcf_in = "", vcf_out = "", samples = "-", region = "";
    int C{3}, niters{1}, nthreads{4}, seed{1};
    double tol{1e-6};
    for (size_t i = 0; i < args.size(); i++)
    {
        if (args[i] == "-b")
            cluster_out = args[++i];
        if (args[i] == "-c")
            C = stoi(args[++i]);
        if (args[i] == "-f")
            vcf_in = args[++i];
        if (args[i] == "-o")
            vcf_out = args[++i];
        if (args[i] == "-g")
            beagle_in = args[++i];
        if (args[i] == "-i")
            niters = stoi(args[++i]);
        if (args[i] == "-n")
            nthreads = stoi(args[++i]);
        if (args[i] == "-r")
            region = args[++i];
        if (args[i] == "-s")
            samples = args[++i];
        if (args[i] == "-seed")
            seed = stoi(args[++i]);
        if (args[i] == "-tol")
            tol = stod(args[++i]);
    }

    // ========= core calculation part ===========================================
    int N, M;
    DoubleVec1D genolikes;
    IntVec1D markers;
    StringVec1D sampleids, chrs;

    // read_bcf_genotype_likelihoods(genolikes, markers, N, M, vcffile, samples, region);
    // cout << N << endl;
    // cout << M << endl;
    // cout << Eigen::Map<ArrDouble2D>(genolikes.data(), N * 3, M) << endl;

    read_beagle_genotype_likelihoods(beagle_in, genolikes, sampleids, chrs, markers, N, M);
    cout << "N:" << N << ",M:" << M << ",C:" << C << endl;
    auto transRate = calc_transRate(markers, C);

    double loglike{0};
    ArrDouble2D postProbsZ(M, C * C);
    ArrDouble2D postProbsZandG(M, C * C * 4);

    fastPhaseK2 nofaith(N, M, C, seed, cluster_out);
    nthreads = nthreads < N ? nthreads : N;
    ThreadPool poolit(nthreads);
    vector<future<double>> llike;
    for (int it = 0; it < niters + 1; it++)
    {
        postProbsZ.setZero();
        postProbsZandG.setZero();
        for (int i = 0; i < N; i++)
        {
            if (it == niters)
                llike.emplace_back(poolit.enqueue(&fastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate),
                                                  std::ref(postProbsZ), std::ref(postProbsZandG), true));
            else
                llike.emplace_back(poolit.enqueue(&fastPhaseK2::forwardAndBackwards, &nofaith, i,
                                                  std::ref(genolikes), std::ref(transRate),
                                                  std::ref(postProbsZ), std::ref(postProbsZandG), false));
        }
        loglike = 0;
        for (auto&& l : llike)
            loglike += l.get();
        cout << "iteration " << it << ", log likelihoods: " << loglike << endl;
        nofaith.updateClusterFreqPI(postProbsZ, tol);
        nofaith.updateAlleleFreqWithinCluster(postProbsZandG, tol);
        llike.clear(); // clear future and renew
    }

    write_bcf_genotype_probability(nofaith.GP.data(), vcf_out, vcf_in, sampleids, markers, chrs[0], N, M);

    return 0;
}
