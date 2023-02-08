#ifndef PHASELESS_IO_H_
#define PHASELESS_IO_H_

#include "vcfpp.h"
#include <cmath>
#include <map>
#include <unordered_map>
#include <zlib.h>

using IntVec1D = std::vector<int>;
using IntVec2D = std::vector<IntVec1D>;
using FloatVec1D = std::vector<float>;
using FloatVec2D = std::vector<FloatVec1D>;
using DoubleVec1D = std::vector<double>;
using DoubleVec2D = std::vector<DoubleVec1D>;
using StringVec1D = std::vector<std::string>;
using StringIntVecMapU = std::unordered_map<std::string, IntVec1D>;
using StringIntMapU = std::unordered_map<std::string, uint64_t>;
using StringIntPairMapU = std::unordered_map<std::string, std::pair<int, int>>;

using MyFloat1D = FloatVec1D;
using MyFloat2D = FloatVec2D;

/*
** @GP maps Eigen matrix layout, (3 x nsnps) x nsamples
*/
inline void write_bcf_genotype_probability(float * GP,
                                           const std::string & vcfout,
                                           const std::string & vcfin,
                                           const StringVec1D & sampleids,
                                           const IntVec1D & markers,
                                           const std::string & chr,
                                           int N,
                                           int M)
{
    FloatVec1D gp(N * 3), ds(N);
    IntVec1D gt(N * 2);
    if(vcfin.empty())
    {
        vcfpp::BcfWriter bw(vcfout, "VCF4.1");
        bw.header.addContig(chr);
        // add GT,GP,DS tag into the header
        bw.header.addFORMAT("GT", "1", "String", "Unphased genotype");
        bw.header.addFORMAT("GP", "3", "Float", "Posterior genotype probability of 0/0, 0/1, and 1/1");
        bw.header.addFORMAT("DS", "1", "Float", "Diploid dosage");
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
        for(int m = 0; m < M; m++)
        {
            var.setPOS(markers[m]);
            for(int i = 0; i < N; i++)
            {
                gp[i * 3 + 0] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 0]) / 1e3;
                gp[i * 3 + 1] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 1]) / 1e3;
                gp[i * 3 + 2] = std::lround(1e3 * GP[i * M * 3 + m * 3 + 2]) / 1e3;
                ds[i] = gp[i * 3 + 1] + gp[i * 3 + 2] * 2;
                gt[i * 2 + 0] = !(gp[i * 3 + 0] > gp[i * 3 + 1] && gp[i * 3 + 0] > gp[i * 3 + 2]);
                gt[i * 2 + 1] = (gp[i * 3 + 2] > gp[i * 3 + 1] && gt[i * 2 + 0]);
            }
            var.setGenotypes(gt);
            var.setFORMAT("GP", gp);
            var.setFORMAT("DS", ds);
            bw.writeRecord(var);
        }
    }
}

inline void read_bcf_genotype_likelihoods(const std::string & vcffile,
                                          const std::string & samples,
                                          const std::string & region,
                                          MyFloat1D & GL,
                                          IntVec1D & markers,
                                          int & nsamples,
                                          int & nsnps,
                                          bool snp_major = true)
{
    vcfpp::BcfReader vcf(vcffile, samples, region);
    vcfpp::BcfRecord var(vcf.header);
    nsamples = vcf.nsamples;
    nsnps = 0;
    IntVec1D PL;
    IntVec2D PL2;
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
                                             StringVec1D & sampleids,
                                             StringIntVecMapU & chrs,
                                             StringIntMapU & starts,
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
    FloatVec2D GL2; // return 2D either sample-major or snp-major

    // PARSE BEAGLE FILE
    fp = gzopen(beagle.c_str(), "r");
    original = buffer = (char *)calloc(bufsize, sizeof(char));
    zlgets(fp, &buffer, &bufsize);
    if(buffer != original) original = buffer;
    int nCol = 1;
    strtok_r(buffer, delims, &buffer);
    while((tok = strtok_r(NULL, delims, &buffer)))
    {
        if(nCol % 3 == 0) sampleids.push_back(std::string(tok));
        nCol++;
    }
    if(nCol % 3) throw std::runtime_error("Number of columns should be a multiple of 3.\n");
    nsamples = nCol / 3 - 1;
    FloatVec1D gli(nsamples * 3); // current line
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
        if(starts.count(chr) == 0) starts[chr] = nsnps; // assume chr ordered
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
                GL[i * nsnps * 3 + j * 3 + 0] = GL2[j][i * 3 + 0];
                GL[i * nsnps * 3 + j * 3 + 1] = GL2[j][i * 3 + 1];
                GL[i * nsnps * 3 + j * 3 + 2] = GL2[j][i * 3 + 2];
            }
        }
    }
}

/*
** @params GL genotype likelihoods, sample-majored
** @return GL chunked genotype likelihoods, snp-majored
*/
inline auto subset_genotype_likelihoods(const std::string & ichr,
                                        const MyFloat1D & GL,
                                        const StringIntVecMapU & chrs,
                                        const StringIntMapU & starts,
                                        int N)
{
    int M = chrs.at(ichr).size();
    MyFloat1D chunkGL(M * N * 3);
    int i, j;
    uint64_t step = (uint64_t)starts.at(ichr) * N * 3;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < M; j++)
        {
            chunkGL[i * M * 3 + j * 3 + 0] = GL[step + j * N * 3 + i * 3 + 0];
            chunkGL[i * M * 3 + j * 3 + 1] = GL[step + j * N * 3 + i * 3 + 1];
            chunkGL[i * M * 3 + j * 3 + 2] = GL[step + j * N * 3 + i * 3 + 2];
        }
    }

    return chunkGL;
}

#endif // PHASELESS_IO_H_
