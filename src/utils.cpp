#include "utils.hpp"

#include "io.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>

using namespace std;

std::string make_beagle_header(std::string fam)
{
    std::ifstream ifs(fam);
    std::string line;
    std::string hdr{"marker\tallele1\tallele2"};
    while(getline(ifs, line))
    {
        std::stringstream ss(line);
        std::vector<std::string> token(std::istream_iterator<std::string>{ss},
                                       std::istream_iterator<std::string>{});
        hdr += "\t" + token[1] + "\t" + token[1] + "\t" + token[1];
    }
    hdr += "\n";
    return hdr;
}

int run_parse_main(Options & opts)
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
    std::ofstream ofs_alpha(opts.out.string() + ".alpha.bin", std::ios::binary);
    std::ofstream ofs_beta(opts.out.string() + ".beta.bin", std::ios::binary);
    int nsamples = ids.size();
    ofs_alpha.write((char *)&genome->C, 4);
    ofs_alpha.write((char *)&nsamples, 4);
    ofs_beta.write((char *)&genome->C, 4);
    ofs_beta.write((char *)&nsamples, 4);
    MyArr2D alpha, beta;
    if(ic < 0)
    {
        ofs_alpha.write((char *)&genome->G, 4);
        ofs_beta.write((char *)&genome->G, 4);
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
                if(!((1 - (alpha * beta).colwise().sum()).abs() < 1e-6).all())
                    cao.error("gamma sum is not 1.0!\n");
                ofs_alpha.write((char *)alpha.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
                ofs_beta.write((char *)beta.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
            }
        }
    }
    else
    {
        const int iM = genome->pos[ic].size();
        const int nGrids = genome->B > 1 ? (iM + genome->B - 1) / genome->B : iM;
        ofs_alpha.write((char *)&nGrids, 4);
        ofs_beta.write((char *)&nGrids, 4);
        std::ofstream ofs_ae(opts.out.string() + ".cluster.freq");
        Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
        MyArr2D ae = MyArr2D::Zero(genome->C, nGrids);
        for(auto ind : ids)
        {
            alpha.setZero(genome->C * genome->C, nGrids);
            beta.setZero(genome->C * genome->C, nGrids);
            get_cluster_likelihood(ind, iM, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                                   genome->F[ic]);
            if(!((1 - (alpha * beta).colwise().sum()).abs() < 1e-6).all())
                cao.error("gamma sum is not 1.0!\n");
            ofs_alpha.write((char *)alpha.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
            ofs_beta.write((char *)beta.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
            for(int i = 0; i < nGrids; i++)
                ae.col(i) += (alpha.col(i) * beta.col(i)).reshaped(genome->C, genome->C).colwise().sum();
        }
        ae.rowwise() /= ae.colwise().sum();
        ofs_ae << ae.transpose().format(fmt) << "\n";
    }
    return 0;
}

int run_convert_main(Options & opts)
{
    cao.cao.open(opts.out.string() + ".log");
    cao.is_screen = !opts.noscreen;
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
        cao.error("Incorrect magic number in plink bed file.\n");
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
