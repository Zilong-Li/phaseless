#ifndef PHASELESS_IO_H_
#define PHASELESS_IO_H_

#include "common.hpp"
#include "vcfpp.h"
#include <cmath>
#include <zlib.h>

struct VariantMetadata
{
    std::vector<String1D> ids;
    std::vector<String1D> refs;
    std::vector<String1D> alts;
};

vcfpp::BcfWriter make_bcfwriter(std::string vcfout, const String1D & chrs, const String1D & sampleids);
MyArr2D extract_gp_chunk(const MyArr2D & GP, int start, int count);
void write_bigass_to_bcf(vcfpp::BcfWriter & bw,
                          const MyFloat * GP,
                          std::string chr,
                          const Int1D & markers,
                          const String1D & ids,
                          const String1D & refs,
                          const String1D & alts);
Int1D filter_sites_per_chunk(const MyFloat * GP, double tol, int N, int M);
Int1D write_bcf_genotype_probability(const MyFloat * GP,
                                     std::string chr,
                                     const Int1D & markers,
                                     const String1D & sampleids,
                                     std::string vcfout,
                                     double infotol);
void chunk_bcf_genotype_likelihoods(const std::unique_ptr<BigAss> & genome,
                                    const std::string & vcffile,
                                    const std::string & region,
                                    const std::string & samples,
                                    VariantMetadata * metadata = nullptr);
int zlgets(gzFile gz, char ** buf, uint64_t * size);
void read_beagle_genotype_likelihoods(const std::string & beagle,
                                      MyFloat1D & GL,
                                      String1D & sampleids,
                                      MapStringInt1D & chrs,
                                      int & nsamples,
                                      int & nsnps,
                                      bool snp_major);
void chunk_beagle_genotype_likelihoods(const std::unique_ptr<BigAss> & genome,
                                       const std::string & beagle,
                                       VariantMetadata * metadata = nullptr);
void update_bigass_inplace(const std::unique_ptr<BigAss> & genome);
size_t count_lines(std::string fpath);
std::vector<uint8_t> read_plink_bed(std::string plink);
std::tuple<std::vector<uint8_t>, std::vector<std::string>> read_plink_bed(std::ifstream & ifs_bed,
                                                                          std::ifstream & ifs_bim,
                                                                          const uint64_t nsamples,
                                                                          const uint64_t nsnps);
std::string convert_geno2like(std::vector<uint8_t> bed,
                              std::vector<std::string> marker,
                              const uint64_t nsamples);
void init_bigass(const std::unique_ptr<BigAss> & genome, const Options & opts, VariantMetadata * metadata = nullptr);
void update_bigass(const std::unique_ptr<BigAss> & genome, const Options & opts);

#endif // PHASELESS_IO_H_
