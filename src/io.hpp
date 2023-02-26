#ifndef PHASELESS_IO_H_
#define PHASELESS_IO_H_

#include "common.hpp"
#include "vcfpp.h"
#include <cmath>
#include <zlib.h>

/*
** @GP maps Eigen matrix layout, (3 x nsnps) x nsamples
*/
inline IntVec1D write_bcf_genotype_probability(MyFloat * GP,
                                               std::string chr,
                                               const IntVec1D & markers,
                                               const StringVec1D & sampleids,
                                               std::string vcfout,
                                               double infotol = 0)
{
    int N = sampleids.size();
    int M = markers.size();
    FloatVec1D gp(N * 3), ds(N);
    IntVec1D gt(N * 2);
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
    IntVec1D idx2rm;
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
        if(infotol > 0 && info > infotol) idx2rm.push_back(m);
    }
    return idx2rm;
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
    FloatVec2D GL2; // return 2D either sample-major or snp-major

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
                                              const std::string & beagle,
                                              bool discard = false)
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
    IntVec1D markers;
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
    if((!markers.empty() && discard == false) || genome->nchunks == 0)
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
    }
    else
    { // discard the rest small piece
        genome->nsnps -= markers.size();
    }
}

inline void thin_bigass(int ichunk,
                        const IntVec1D & idx2rm,
                        const std::unique_ptr<BigAss> & genome,
                        MyArr2D & PI,
                        MyArr2D & F,
                        MyArr2D & transRate)
{
    if(!idx2rm.empty())
    {
        IntVec1D idx2keep;
        int M = genome->pos[ichunk].size();
        for(int i = 0, j = 0; i < M; i++)
        {
            if(idx2rm[j] == i)
                j++;
            else
                idx2keep.push_back(i);
        }
        PI = PI(idx2keep, Eigen::all);
        F = F(idx2keep, Eigen::all);
        transRate = transRate(Eigen::all, idx2keep);
        int im = idx2keep.size();
        MyFloat1D gls(genome->nsamples * im * 3);
        for(int i = 0; i < genome->nsamples; i++)
        {
            for(int j = 0; j < im; j++)
            {
                gls[i * im * 3 + 0 * im + j] = genome->gls[ichunk][i * M * 3 + 0 * M + idx2keep[j]];
                gls[i * im * 3 + 1 * im + j] = genome->gls[ichunk][i * M * 3 + 1 * M + idx2keep[j]];
                gls[i * im * 3 + 2 * im + j] = genome->gls[ichunk][i * M * 3 + 2 * M + idx2keep[j]];
            }
        }
        genome->gls[ichunk] = gls;
        IntVec1D pos(im);
        for(int j = 0; j < im; j++) pos[j] = genome->pos[ichunk][idx2keep[j]];
        genome->pos[ichunk] = pos;
        genome->nsnps -= idx2rm.size();
    }
    genome->transRate.emplace_back(MyFloat1D(transRate.data(), transRate.data() + transRate.size()));
    genome->PI.emplace_back(MyFloat1D(PI.data(), PI.data() + PI.size()));
    genome->F.emplace_back(MyFloat1D(F.data(), F.data() + F.size()));
}

#endif // PHASELESS_IO_H_
