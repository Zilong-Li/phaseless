#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>

inline int run_parse_main(Options & opts)
{
    Logger cao(opts.out.string() + "log", !opts.noscreen);
    cao.print(opts.opts_in_effect);
    auto filesize = std::filesystem::file_size(opts.in_bin);
    std::error_code ec;
    std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome =
        std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    const int ic = opts.ichunk;
    if(ic >= genome->nchunks)
    {
        cao.error("the chunk to extract", ic, "(0-baded) is not less than total chunks", genome->nchunks);
        return 1;
    }
    const int iM = genome->pos[ic].size();
    std::ofstream ofs(opts.out.string() + "haplike.bin", std::ios::binary);
    ofs.write((char *)&genome->C, 4);
    ofs.write((char *)&genome->nsamples, 4);
    ofs.write((char *)&iM, 4);
    // haplike is p(Z_is |X_is , theta) = alpha * beta / p(X|theta) = gamma
    MyArr2D alpha, beta;
    for(int ind = 0; ind < genome->nsamples; ind++)
    {
        alpha.setZero(genome->C * genome->C, iM);
        beta.setZero(genome->C * genome->C, iM);
        getClusterLikelihoods(ind, alpha, beta, genome->gls[ic], genome->transRate[ic], genome->PI[ic],
                              genome->F[ic]);
        alpha *= beta;
        ofs.write((char *)alpha.data(), genome->C * genome->C * iM * 4);
    }
    return 0;
}
