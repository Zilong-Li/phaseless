#ifndef PHASELESS_IO_H_
#define PHASELESS_IO_H_

#include "common.hpp"
#include "vcfpp.h"
#include <cmath>
#include <zlib.h>

inline auto make_bcfwriter(std::string vcfout, const String1D & chrs, const String1D & sampleids)
{
    vcfpp::BcfWriter bw(vcfout, "VCF4.1");
    // add contig
    for(auto const & chr : chrs) bw.header.addContig(chr);
    // add GT,GP,DS tag into the header
    bw.header.addFORMAT("GT", "1", "String", "Unphased genotype");
    bw.header.addFORMAT("GP", "3", "Float", "Posterior genotype probability of 0/0, 0/1, and 1/1");
    bw.header.addFORMAT("DS", "1", "Float", "Diploid dosage");
    bw.header.addINFO("INFO", "1", "Float", "INFO score");
    bw.header.addINFO("EAF", "1", "Float", "Estimated allele frequency");
    // add all samples in the header
    for(auto const & id : sampleids) bw.header.addSample(id);
    bw.writeHeader();
    return bw;
}

inline void write_bigass_to_bcf(vcfpp::BcfWriter & bw,
                                const MyFloat * GP,
                                std::string chr,
                                const Int1D & markers)
{
    int M = markers.size();
    int N = bw.header.nSamples();
    vcfpp::BcfRecord var(bw.header); // construct a variant record from the header
    var.setCHR(chr.c_str());
    double thetaHat, info, eaf, eij, fij, a0, a1;
    int i, m;
    Float1D gp(N * 3), ds(N);
    Int1D gt(N * 2);
    for(m = 0; m < M; m++)
    {
        for(eij = 0, fij = 0, i = 0; i < N; i++)
        {
            gp[i * 3 + 0] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 0]) / 1e3;
            gp[i * 3 + 1] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 1]) / 1e3;
            gp[i * 3 + 2] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 2]) / 1e3;
            ds[i] = gp[i * 3 + 1] + gp[i * 3 + 2] * 2;
            gt[i * 2 + 0] = !(gp[i * 3 + 0] > gp[i * 3 + 1] && gp[i * 3 + 0] > gp[i * 3 + 2]);
            gt[i * 2 + 1] = (gp[i * 3 + 2] > gp[i * 3 + 1] && gt[i * 2 + 0]);
            a0 = GP[i * M * 3 + m * 3 + 1] + GP[i * M * 3 + m * 3 + 2] * 2;
            a1 = GP[i * M * 3 + m * 3 + 1] + GP[i * M * 3 + m * 3 + 2] * 4;
            eij += a0;
            fij += a1 - a0 * a0;
        }
        eaf = eij / 2 / N;
        thetaHat = std::lround(1e2 * eaf) / 1e2;
        if(thetaHat == 0 || thetaHat == 1)
            info = 1;
        else if(info < 0)
            info = 0;
        else
        {
            info = 1 - fij / (2 * N * eaf * (1 - eaf));
            info = std::lround(1e3 * info) / 1e3;
            info = info > 1 ? 1 : info;
        }
        eaf = std::lround(1e6 * eaf) / 1e6;
        var.setPOS(markers[m]);
        var.setGenotypes(gt);
        var.setFORMAT("GP", gp);
        var.setFORMAT("DS", ds);
        var.setINFO("INFO", info);
        var.setINFO("EAF", eaf);
        bw.writeRecord(var);
    }
}

inline Int1D filter_sites_per_chunk(const MyFloat * GP, double tol, int N, int M)
{
    Float1D gp(N * 3);
    double thetaHat, info, eaf, eij, fij, a0, a1;
    int i, m;
    Int1D idx2rm;
    for(m = 0; m < M; m++)
    {
        for(eij = 0, fij = 0, i = 0; i < N; i++)
        {
            gp[i * 3 + 0] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 0]) / 1e3;
            gp[i * 3 + 1] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 1]) / 1e3;
            gp[i * 3 + 2] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 2]) / 1e3;
            a0 = gp[i * 3 + 1] + gp[i * 3 + 2] * 2;
            a1 = gp[i * 3 + 1] + gp[i * 3 + 2] * 4;
            eij += a0;
            fij += a1 - a0 * a0;
        }
        eaf = eij / 2 / N;
        info = 1 - fij / (eij * (1 - eaf));
        thetaHat = std::lround(1e2 * eaf) / 1e2;
        if(thetaHat == 0 || thetaHat == 1)
            info = 1;
        else if(info < 0)
            info = 0;
        else
            info = std::lround(1e3 * info) / 1e3;
        eaf = std::lround(1e6 * eaf) / 1e6;
        if(info < tol || info == 1) idx2rm.push_back(m);
    }
    return idx2rm;
}

/*
** @GP maps Eigen matrix layout, (3 x nsnps) x nsamples
*/
inline Int1D write_bcf_genotype_probability(const MyFloat * GP,
                                            std::string chr,
                                            const Int1D & markers,
                                            const String1D & sampleids,
                                            std::string vcfout,
                                            double infotol = 0)
{
    int N = sampleids.size();
    int M = markers.size();
    Float1D gp(N * 3), ds(N);
    Int1D gt(N * 2);
    vcfpp::BcfWriter bw(vcfout, "VCF4.1");
    bw.header.addContig(chr);
    // add GT,GP,DS tag into the header
    bw.header.addFORMAT("GT", "1", "String", "Unphased genotype");
    bw.header.addFORMAT("GP", "3", "Float", "Posterior genotype probability of 0/0, 0/1, and 1/1");
    bw.header.addFORMAT("DS", "1", "Float", "Diploid dosage");
    bw.header.addINFO("INFO", "1", "Float", "INFO score");
    bw.header.addINFO("EAF", "1", "Float", "Estimated allele frequency");
    // add all samples in the header
    for(int i = 0; i < N; i++)
    {
        if(sampleids.empty())
            bw.header.addSample(std::string("id_" + std::to_string(i)));
        else
            bw.header.addSample(sampleids[i]);
    }
    bw.writeHeader();
    vcfpp::BcfRecord var(bw.header); // construct a variant record from the header
    double thetaHat, info, eaf, eij, fij, a0, a1;
    int i, m;
    Int1D idx2rm;
    for(m = 0; m < M; m++)
    {
        var.setPOS(markers[m]);
        for(eij = 0, fij = 0, i = 0; i < N; i++)
        {
            gp[i * 3 + 0] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 0]) / 1e3;
            gp[i * 3 + 1] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 1]) / 1e3;
            gp[i * 3 + 2] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 2]) / 1e3;
            ds[i] = gp[i * 3 + 1] + gp[i * 3 + 2] * 2;
            gt[i * 2 + 0] = !(gp[i * 3 + 0] > gp[i * 3 + 1] && gp[i * 3 + 0] > gp[i * 3 + 2]);
            gt[i * 2 + 1] = (gp[i * 3 + 2] > gp[i * 3 + 1] && gt[i * 2 + 0]);
            a0 = gp[i * 3 + 1] + gp[i * 3 + 2] * 2;
            a1 = gp[i * 3 + 1] + gp[i * 3 + 2] * 4;
            eij += a0;
            fij += a1 - a0 * a0;
        }
        eaf = eij / 2 / N;
        info = 1 - fij / (eij * (1 - eaf));
        thetaHat = std::lround(1e2 * eaf) / 1e2;
        if(thetaHat == 0 || thetaHat == 1)
            info = 1;
        else if(info < 0)
            info = 0;
        else
            info = std::lround(1e3 * info) / 1e3;
        eaf = std::lround(1e6 * eaf) / 1e6;
        var.setGenotypes(gt);
        var.setFORMAT("GP", gp);
        var.setFORMAT("DS", ds);
        var.setINFO("INFO", info);
        var.setINFO("EAF", eaf);
        bw.writeRecord(var);
        if(infotol > 0 && info < infotol) idx2rm.push_back(m);
    }
    return idx2rm;
}

inline void read_bcf_genotype_likelihoods(const std::string & vcffile,
                                          const std::string & samples,
                                          const std::string & region,
                                          MyFloat1D & GL,
                                          Int1D & markers,
                                          int & nsamples,
                                          int & nsnps,
                                          bool snp_major = true)
{
    vcfpp::BcfReader vcf(vcffile, samples, region);
    vcfpp::BcfRecord var(vcf.header);
    nsamples = vcf.nsamples;
    nsnps = 0;
    Int1D PL;
    Int2D PL2;
    GL.clear();
    while(vcf.getNextVariant(var))
    {
        nsnps++;
        markers.push_back(var.POS());
        var.getFORMAT("PL", PL);
        if(snp_major)
            PL2.push_back(PL);
        else
            GL.insert(GL.end(), PL.begin(), PL.end());
    }
    // now do transpose if snp major is wanted
    if(snp_major)
    {
        GL.resize(nsamples * nsnps * 3);
        for(int i = 0; i < nsamples; i++)
        {
            for(int j = 0; j < nsnps; j++)
            {
                GL[i * nsnps * 3 + j * 3 + 0] = PL2[j][i * 3 + 0];
                GL[i * nsnps * 3 + j * 3 + 1] = PL2[j][i * 3 + 1];
                GL[i * nsnps * 3 + j * 3 + 2] = PL2[j][i * 3 + 2];
            }
        }
    }
}

/*
** @brief gz line gets
** @param gz   file hander returned by gzopen
** @param buf  buffer used for storing data
** @param size buffer size for realloc buffer
** @return length and buffer of current line
*/
inline int zlgets(gzFile gz, char ** buf, uint64_t * size)
{
    int rlen = 0;
    char * tok = gzgets(gz, *buf + rlen, *size - rlen); // return buf or NULL
    if(!tok) return rlen;
    int tmp = tok ? strlen(tok) : 0;
    if(tok[tmp - 1] != '\n')
    {
        // expand buf size if no end-of-line found
        rlen += tmp;
        *size *= 2;
        *buf = (char *)realloc(*buf, *size);
    }
    rlen += tmp;
    return rlen;
}

/*
** @return GL genotype likelihoods, nsnps x (nsamples x 3)
*/
inline void read_beagle_genotype_likelihoods(const std::string & beagle,
                                             MyFloat1D & GL,
                                             String1D & sampleids,
                                             MapStringInt1D & chrs,
                                             int & nsamples,
                                             int & nsnps,
                                             bool snp_major = true)
{
    // VARIBLES
    gzFile fp = nullptr;
    char *original, *buffer, *tok, *chr, *pos;
    uint64_t bufsize = (uint64_t)128 * 1024 * 1024;
    int i, j;
    const char * delims = "\t \n";
    Float2D GL2; // return 2D either sample-major or snp-major

    // PARSE BEAGLE FILE
    fp = gzopen(beagle.c_str(), "r");
    original = buffer = (char *)calloc(bufsize, sizeof(char));
    zlgets(fp, &buffer, &bufsize);
    if(buffer != original) original = buffer;
    int nCol = 1;
    strtok_r(buffer, delims, &buffer);
    while((tok = strtok_r(NULL, delims, &buffer)))
        if(nCol++ % 3 == 0) sampleids.push_back(std::string(tok));
    if(nCol % 3) throw std::runtime_error("Number of columns should be a multiple of 3.\n");
    nsamples = nCol / 3 - 1;
    Float1D gli(nsamples * 3); // current line
    // now fill in GL and get nsnps
    nsnps = 0;
    GL.clear();
    buffer = original;
    while(zlgets(fp, &buffer, &bufsize))
    {
        if(buffer != original) original = buffer;
        tok = strtok_r(buffer, delims, &buffer); // id: chr_pos
        chr = strtok(tok, "_");
        pos = strtok(NULL, "_");
        chrs[chr].push_back(std::stoi(pos));
        tok = strtok_r(NULL, delims, &buffer); // ref
        tok = strtok_r(NULL, delims, &buffer); // alt
        for(i = 0; i < nsamples; i++)
        {
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 0] = strtod(tok, NULL);
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 1] = strtod(tok, NULL);
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 2] = strtod(tok, NULL);
        }
        buffer = original;
        nsnps++;
        if(snp_major)
            GL2.push_back(gli);
        else
            GL.insert(GL.end(), gli.begin(), gli.end());
    }
    gzclose(fp);
    // now do transpose if snp major is wanted
    if(snp_major)
    {
        GL.resize(nsamples * nsnps * 3);
        for(i = 0; i < nsamples; i++)
        {
            for(j = 0; j < nsnps; j++)
            {
                GL[i * nsnps * 3 + 0 * nsnps + j] = GL2[j][i * 3 + 0];
                GL[i * nsnps * 3 + 1 * nsnps + j] = GL2[j][i * 3 + 1];
                GL[i * nsnps * 3 + 2 * nsnps + j] = GL2[j][i * 3 + 2];
            }
        }
    }
}

inline void chunk_beagle_genotype_likelihoods(const std::unique_ptr<BigAss> & genome,
                                              const std::string & beagle)
{
    // VARIBLES
    gzFile fp = nullptr;
    char *original, *buffer, *tok, *pos;
    std::string chr0 = "", chr1 = "";
    uint64_t bufsize = (uint64_t)128 * 1024 * 1024;
    const char * delims = "\t \n";

    // PARSE BEAGLE FILE
    fp = gzopen(beagle.c_str(), "r");
    original = buffer = (char *)calloc(bufsize, sizeof(char));
    zlgets(fp, &buffer, &bufsize);
    if(buffer != original) original = buffer;
    int nCol = 1;
    strtok_r(buffer, delims, &buffer);
    while((tok = strtok_r(NULL, delims, &buffer)))
        if(nCol++ % 3 == 0) genome->sampleids.push_back(std::string(tok));
    if(nCol % 3) throw std::runtime_error("Number of columns should be a multiple of 3.\n");
    genome->nsamples = nCol / 3 - 1;
    MyFloat2D glchunk; // gl for one chunk
    MyFloat1D gli(genome->nsamples * 3); // current line
    genome->nsnps = 0;
    buffer = original;
    Int1D markers;
    genome->nchunks = 0;
    bool samechr;
    int i, j, im, isnp{0};
    while(zlgets(fp, &buffer, &bufsize))
    {
        if(buffer != original) original = buffer;
        tok = strtok_r(buffer, delims, &buffer); // id: chr_pos
        chr1 = (std::string)strtok(tok, "_");
        if(chr0.empty() || chr0 == chr1)
            samechr = true;
        else
            samechr = false;
        pos = strtok(NULL, "_");
        tok = strtok_r(NULL, delims, &buffer); // ref
        tok = strtok_r(NULL, delims, &buffer); // alt
        for(i = 0; i < genome->nsamples; i++)
        {
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 0] = strtod(tok, NULL);
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 1] = strtod(tok, NULL);
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 2] = strtod(tok, NULL);
        }
        if(samechr || isnp == 0)
        {
            markers.push_back(std::stoi(pos));
            glchunk.push_back(gli);
        }
        if(((++isnp % genome->chunksize == 0) || (samechr == false)) && isnp != 1)
        {
            genome->nchunks++;
            genome->chrs.push_back(chr0);
            im = markers.size();
            genome->pos.push_back(markers);
            markers.clear();
            // transpose glchunk into genome->gl then clear it
            MyFloat1D gl(genome->nsamples * im * 3);
            for(i = 0; i < genome->nsamples; i++)
            {
                for(j = 0; j < im; j++)
                {
                    gl[i * im * 3 + 0 * im + j] = glchunk[j][i * 3 + 0];
                    gl[i * im * 3 + 1 * im + j] = glchunk[j][i * 3 + 1];
                    gl[i * im * 3 + 2 * im + j] = glchunk[j][i * 3 + 2];
                }
            }
            glchunk.clear();
            genome->gls.push_back(gl);
            if((isnp % genome->chunksize != 0) && samechr == false)
            {
                genome->ends.push_back(genome->nchunks - 1);
                markers.push_back(std::stoi(pos));
                glchunk.push_back(gli);
                isnp = 1;
            }
            else
            {
                isnp = 0;
            }
        }
        genome->nsnps++;
        chr0 = chr1;
        buffer = original;
    }
    gzclose(fp);
    if(!markers.empty())
    { // add the rest into genome
        im = markers.size();
        genome->pos.push_back(markers);
        markers.clear();
        genome->nchunks++;
        genome->chrs.push_back(chr0);
        // transpose glchunk into genome->gl then clear it
        MyFloat1D gl(genome->nsamples * im * 3);
        for(i = 0; i < genome->nsamples; i++)
        {
            for(j = 0; j < im; j++)
            {
                gl[i * im * 3 + 0 * im + j] = glchunk[j][i * 3 + 0];
                gl[i * im * 3 + 1 * im + j] = glchunk[j][i * 3 + 1];
                gl[i * im * 3 + 2 * im + j] = glchunk[j][i * 3 + 2];
            }
        }
        glchunk.clear();
        genome->gls.push_back(gl);
        genome->ends.push_back(genome->nchunks - 1);
    }
    // now evenly split last two chunks of each chromosome
}

inline void update_bigass_inplace(const std::unique_ptr<BigAss> & genome)
{
    int ic, ndiff;
    for(ic = 0, genome->nsnps = 0; ic < genome->nchunks; ic++) genome->nsnps += genome->pos[ic].size();
    for(ic = 0; ic < genome->nchunks; ic++)
    {
        if(ic == 0) continue; // assume the first chunksize is not greater than the defined
        if((int)genome->pos[ic - 1].size() < genome->chunksize)
        {
            ndiff = genome->chunksize - genome->pos[ic - 1].size();
            if((int) genome->pos[ic].size() >= ndiff)
            {
                genome->pos[ic - 1].insert(genome->pos[ic - 1].end(), genome->pos[ic].begin(),
                                           genome->pos[ic].begin() + ndiff);
                Int1D(genome->pos[ic].begin() + ndiff, genome->pos[ic].end()).swap(genome->pos[ic]);
            }
            else
                break;
        }
    }
}

inline size_t count_lines(std::string fpath)
{
    std::ifstream ifs(fpath);
    size_t count = 0;
    std::string line;
    while(getline(ifs, line)) count++;
    return count;
}

inline auto read_plink_bed(std::string plink)
{
    std::setlocale(LC_ALL, "C"); // use minial locale
    std::ios_base::sync_with_stdio(false); // don't sync
    uint64_t nsamples = count_lines(plink + ".fam");
    uint64_t nsnps = count_lines(plink + ".bim");
    uint64_t bed_bytes_per_snp = (nsamples + 3) >> 2; // get ceiling(nsamples/4)
    std::ifstream ifs(plink + ".bed", std::ios::in | std::ios::binary);
    uint8_t header[3];
    ifs.read(reinterpret_cast<char *>(&header[0]), 3);
    if((header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
        throw std::invalid_argument("Incorrect magic number in plink bed file.\n");
    // read all into inbed
    std::vector<uint8_t> bed(bed_bytes_per_snp * nsnps);
    ifs.read(reinterpret_cast<char *>(&bed[0]), bed_bytes_per_snp * nsnps);
    return bed;
}

inline auto read_plink_bed(std::ifstream & ifs_bed,
                           std::ifstream & ifs_bim,
                           const uint64_t nsamples,
                           const uint64_t nsnps)
{
    std::string line;
    std::vector<std::string> marker;
    uint64_t i = 0;
    while(i < nsnps && getline(ifs_bim, line))
    {
        std::stringstream ss(line);
        std::vector<std::string> token(std::istream_iterator<std::string>{ss},
                                       std::istream_iterator<std::string>{});
        marker.push_back(token[0] + "_" + token[3] + "\t" + token[4] + "\t" + token[5]);
        i++;
    }
    uint64_t bed_bytes_per_snp = (nsamples + 3) >> 2; // get ceiling(nsamples/4)
    uint64_t bed_size = bed_bytes_per_snp * i;
    std::vector<uint8_t> bed(bed_size);
    ifs_bed.read(reinterpret_cast<char *>(bed.data()), bed_size);
    return std::tuple(bed, marker);
}

inline auto convert_geno2like(std::vector<uint8_t> bed,
                              std::vector<std::string> marker,
                              const uint64_t nsamples)
{
    uint64_t bed_bytes_per_snp = (nsamples + 3) >> 2; // get ceiling(nsamples/4)
    std::string res;
    uint8_t buf, g;
    uint64_t i, b, j, k;
    /**
     * 00 ->  0 Homozygous for first allele in .bim file
     * 01 ->  1 Missing genotype
     * 10 ->  2 Heterozygous (copy of A1)
     * 11 ->  3 Homozygous for second allele in .bim file
     */
    for(i = 0; i < marker.size(); i++)
    {
        std::string gl = marker[i];
        for(j = 0, b = 0; b < bed_bytes_per_snp; b++)
        {
            buf = bed[i * bed_bytes_per_snp + b];
            for(k = 0; k < 4; ++k, ++j)
            {
                if(j < nsamples)
                {
                    g = buf & 3;
                    if(g == 0)
                        gl += "\t1\t0\t0";
                    else if(g == 1)
                        gl += "\t0\t0\t0";
                    else if(g == 2)
                        gl += "\t0\t1\t0";
                    else if(g == 3)
                        gl += "\t0\t0\t1";
                    buf >>= 2; // shift packed data and throw away genotype just processed.
                }
            }
        }
        gl += "\n";
        res += gl;
    }
    return res;
}

inline void init_bigass(const std::unique_ptr<BigAss> & genome, const Options & opts)
{
    genome->chunksize = opts.chunksize, genome->C = opts.C, genome->B = opts.gridsize;
    tim.clock();
    chunk_beagle_genotype_likelihoods(genome, opts.in_beagle);
    int G{0};
    for(int ic = 0; ic < genome->nchunks; ic++)
    {
        int nsnps = genome->pos[ic].size();
        int nGrids = genome->B > 1 ? (nsnps + genome->B - 1) / genome->B : nsnps;
        G += nGrids;
    }
    genome->G = G;
    if(genome->B == 1 && genome->G != genome->nsnps)
        cao.error("number of grids should be same as snps if B=1");
    cao.print(tim.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
              ", M =", genome->nsnps, ", nchunks =", genome->nchunks, ", B =", opts.gridsize,
              ", seed =", opts.seed);
    cao.done(tim.date(), "elapsed time for parsing beagle file", std::fixed, tim.reltime(), " secs");
}

#endif // PHASELESS_IO_H_
