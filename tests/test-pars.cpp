#include "alpaca/alpaca.h"
#include "common.hpp"
#include "log.hpp"
#include "timer.hpp"
#include <filesystem>

using namespace std;

int main(int argc, char * argv[])
{
    Timer tm;
    Logger cao("test-pars.log");
    cao.warn(tm.date(), "-> running test-pars");
    filesystem::path in_bin{argv[1]};
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    auto filesize = std::filesystem::file_size(in_bin);
    std::error_code ec;
    std::ifstream ifs(in_bin, std::ios::in | std::ios::binary);
    auto genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    cao.done(tm.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
    cao.print("parsing input:\nC =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks, ", chunksize =", genome->chunksize);
    for(auto i : genome->ends)
    {
        cao.print("end:", i);
    }
    for(int i = 0; i < genome->nchunks; i++)
    {
        for(auto pos : genome->pos[i])
        {
            cao.print("chunk:", i, ",chr:", genome->chrs[i], ",pos:", pos);
        }
    }
    return 0;
}
