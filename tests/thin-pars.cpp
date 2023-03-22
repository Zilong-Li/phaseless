#include "alpaca/alpaca.h"
#include "common.hpp"
#include "log.hpp"
#include "timer.hpp"
#include <algorithm>
#include <filesystem>

using namespace std;

// pick random number of chunks from the bigass genome
template<class ForwardIterator>
ForwardIterator pick_chunk_by_step(ForwardIterator first, ForwardIterator last, int step)
{
    ForwardIterator current = first;
    ForwardIterator result = first;
    while(current != last)
    {
        if(std::distance(first, current) % step != 0)
        {
            if(result != current) *result = std::move(*current);
            ++result;
        }
        ++current;
    }
    return result;
}

int main(int argc, char * argv[])
{
    // ===================== helper message and parameters parsing ===============
    std::vector<std::string> args(argv + 1, argv + argc);
    if(argc <= 1 || args[0] == "-h" || args[0] == "-help")
    {
        std::cout << "Author: Zilong-Li (zilong.dk@gmail.com)\n"
                  << "Usage example:\n"
                  << "     " + (std::string)argv[0] + " -b par.bin -o out -f 0.5\n"
                  << "\nOptions:\n"
                  << "     -b      input binary file with all parameters\n"
                  << "     -f      fraction to be subset\n"
                  << "     -o      output prefix\n"
                  << std::endl;
        return 1;
    }

    filesystem::path out, in_beagle, in_vcf, in_bin;
    double frac{1};
    for(size_t i = 0; i < args.size(); i++)
    {
        if(args[i] == "-b") in_bin = args[++i];
        if(args[i] == "-f") frac = stod(args[++i]);
        if(args[i] == "-o") out.assign(args[++i]);
    }

    // ===================== ===  core calculation part   ========================
    Timer tm;
    Logger cao(out.string() + "pars.log");
    cao.warn(tm.date(), "-> running thin-pars");
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    auto filesize = std::filesystem::file_size(in_bin);
    std::error_code ec;
    std::ifstream ifs(in_bin, std::ios::in | std::ios::binary);
    auto genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    cao.done(tm.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
    cao.print("parsing input:\nC =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks, ", chunksize =", genome->chunksize);
    int nchunks_ = genome->nchunks * frac;
    int step = genome->nchunks / nchunks_;
    genome->chrs.erase(pick_chunk_by_step(genome->chrs.begin(), genome->chrs.end(), step),
                       genome->chrs.end());
    genome->pos.erase(pick_chunk_by_step(genome->pos.begin(), genome->pos.end(), step), genome->pos.end());
    genome->gls.erase(pick_chunk_by_step(genome->gls.begin(), genome->gls.end(), step), genome->gls.end());
    genome->PI.erase(pick_chunk_by_step(genome->PI.begin(), genome->PI.end(), step), genome->PI.end());
    genome->F.erase(pick_chunk_by_step(genome->F.begin(), genome->F.end(), step), genome->F.end());
    genome->transRate.erase(pick_chunk_by_step(genome->transRate.begin(), genome->transRate.end(), step),
                            genome->transRate.end());
    genome->nsnps = 0;
    for(auto & p : genome->pos) genome->nsnps += p.size();
    genome->nchunks = genome->pos.size();
    std::ofstream ofs(out.string() + "pars.bin", std::ios::out | std::ios::binary);
    auto bytes_written = alpaca::serialize<OPTIONS, BigAss>(*genome, ofs);
    ofs.close();
    assert(std::filesystem::file_size(out.string() + "pars.bin") == bytes_written);
    cao.done(tm.date(), "thinnig par.bin done.", bytes_written, " bytes written to file");

    return 0;
}
