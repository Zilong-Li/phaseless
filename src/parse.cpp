#include "common.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>
#include <stdexcept>

inline int run_parse_main(Options & opts)
{
    cao.cao.open(opts.out.string() + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    auto filesize = std::filesystem::file_size(opts.in_bin);
    std::error_code ec;
    std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome =
        std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    int ic = opts.ichunk;
    if(ic >= genome->nchunks)
    {
        cao.error("the chunk ", ic, " (0-based) to be extracted is not less than total chunks",
                  genome->nchunks);
        return 1;
    }
    else if(ic < 0)
    {
        cao.warn("the chunk ", ic, " to be extracted is less than 0. all chunks will be extracted!");
    }
    Int1D ids;
    if(!opts.samples.empty())
    {
        UMapStringInt ids_m;
        int i = 0;
        for(auto & id : genome->sampleids) ids_m[id] = i++;
        std::ifstream ifs(opts.samples);
        if(!ifs.is_open()) cao.error("can not open the file: ", opts.samples);
        std::string line;
        while(getline(ifs, line)) ids.push_back(ids_m[line]);
    }
    else
    {
        for(int ind = 0; ind < genome->nsamples; ind++) ids.push_back(ind);
    }
    // haplike is p(Z_is |X_is , theta) = alpha * beta / p(X|theta) = gamma
    std::ofstream ofs(opts.out.string() + ".haplike.bin", std::ios::binary);
    int nsamples = ids.size();
    ofs.write((char *)&genome->C, 4);
    ofs.write((char *)&nsamples, 4);
    MyArr2D alpha, beta;
    if(ic < 0)
    {
        ofs.write((char *)&genome->G, 4);
        for(auto ind : ids)
        {
            for(ic = 0; ic < genome->nchunks; ic++)
            {
                const int iM = genome->pos[ic].size();
                const int nGrids = genome->B > 1 ? (iM + genome->B - 1) / genome->B : iM;
                alpha.setZero(genome->C * genome->C, nGrids);
                beta.setZero(genome->C * genome->C, nGrids);
                get_cluster_likelihood(ind, iM, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                                       genome->F[ic]);
                alpha *= beta;
                if(!((1 - alpha.colwise().sum()).abs() < 1e-6).all()) cao.error("gamma sum is not 1.0!\n");
                ofs.write((char *)alpha.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
            }
        }
    }
    else
    {
        const int iM = genome->pos[ic].size();
        const int nGrids = genome->B > 1 ? (iM + genome->B - 1) / genome->B : iM;
        ofs.write((char *)&nGrids, 4);
        for(auto ind : ids)
        {
            alpha.setZero(genome->C * genome->C, nGrids);
            beta.setZero(genome->C * genome->C, nGrids);
            get_cluster_likelihood(ind, iM, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                                   genome->F[ic]);
            alpha *= beta;
            if(!((1 - alpha.colwise().sum()).abs() < 1e-6).all()) cao.error("gamma sum is not 1.0!\n");
            ofs.write((char *)alpha.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
        }
    }
    return 0;
}
