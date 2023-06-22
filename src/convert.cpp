#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"

using namespace std;

inline auto make_beagle_header(std::string fam)
{
    std::ifstream ifs(fam);
    std::string line;
    std::string hdr{"marker\tallele1\tallele2"};
    while(getline(ifs, line))
    {
        std::stringstream ss(line);
        std::vector<std::string> token(std::istream_iterator<std::string>{ss}, std::istream_iterator<std::string>{});
        hdr += "\t" + token[1] + "\t" + token[1] + "\t" + token[1];
    }
    hdr += "\n";
    return hdr;
}

inline int run_convert_main(Options & opts)
{
    Logger cao(opts.out.string() + "log", !opts.noscreen);
    cao.print(opts.opts_in_effect);
    std::setlocale(LC_ALL, "C"); // use minial locale
    std::ios_base::sync_with_stdio(false); // don't sync
    uint64_t nsamples = count_lines(opts.in_plink + ".fam");
    uint64_t nsnps = count_lines(opts.in_plink + ".bim");
    cao.print("nsamples :", nsamples, ", nsnps :", nsnps);
    std::string hdr = make_beagle_header(opts.in_plink + ".fam");
    int nchunks = (nsnps + opts.chunksize - 1) / opts.chunksize;
    std::ifstream ifs_bed(opts.in_plink + ".bed", std::ios::in | std::ios::binary);
    uint8_t header[3];
    ifs_bed.read(reinterpret_cast<char *>(&header[0]), 3);
    if((header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
        throw std::invalid_argument("Incorrect magic number in plink bed file.\n");
    ThreadPool poolit(opts.nthreads);
    vector<future<string>> res;
    std::ifstream ifs_bim(opts.in_plink + ".bim", std::ios::in);
    for(int ic = 0; ic < nchunks; ic++)
    {
        const auto [bed, marker] = read_plink_bed(ifs_bed, ifs_bim, nsamples, opts.chunksize);
        res.emplace_back(poolit.enqueue(convert_geno2like, bed, marker, nsamples));
    }
    string out{opts.out.string() + ".gz"};
    gzFile gzfp = gzopen(out.c_str(), "wb");
    gzwrite(gzfp, hdr.c_str(), hdr.size());
    for(auto && ll : res)
    {
        auto gl = ll.get();
        gzwrite(gzfp, gl.c_str(), gl.size());
    }
    gzclose(gzfp);
    return 0;
}
