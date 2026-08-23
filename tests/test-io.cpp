#define _DECLARE_TOOLBOX_HERE

#include "catch.hh"
#include "../src/io.hpp"

using namespace std;

TEST_CASE("test beagle parser for single contig", "[test-io]")
{
    int N, M;
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    String1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_pos, N, M, false);
    REQUIRE(N == 60);
    REQUIRE(M == 1000);
    REQUIRE(chrs_pos.size() == 1);
    auto ichr = chrs_pos.begin()->first;
    REQUIRE(ichr == "5");
}

TEST_CASE("test beagle parser for multiple contigs", "[test-io]")
{
    int N, M;
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    String1D sampleids;
    read_beagle_genotype_likelihoods("../data/all.bgl.gz", genolikes, sampleids, chrs_pos, N, M, false);
    REQUIRE(N == 60);
    REQUIRE(M == 1000);
    REQUIRE(chrs_pos.size() == 2);
    REQUIRE(chrs_pos.begin()->first == "4");
    REQUIRE(next(chrs_pos.begin())->first == "5");
}

TEST_CASE("test beagle parser for bigass genome with balanced chunksize", "[test-io]")
{
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = 100;
    VariantMetadata metadata;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz", &metadata);
    REQUIRE(genome->nsamples == 60);
    REQUIRE(genome->nsnps == 1000);
    REQUIRE(genome->nchunks == genome->pos.size());
    REQUIRE(metadata.ids.size() == genome->pos.size());
    REQUIRE(metadata.refs.size() == genome->pos.size());
    REQUIRE(metadata.alts.size() == genome->pos.size());
    for(int i = 0; i < genome->nchunks; i++) REQUIRE(genome->pos[i].size() == 100);
    REQUIRE(metadata.ids[0][0] == "4_10001");
    REQUIRE(metadata.refs[0][0] == "3");
    REQUIRE(metadata.alts[0][0] == "1");
}

TEST_CASE("chunked GP export preserves site, sample, and allele metadata", "[test-io]")
{
    MyArr2D GP(6, 2);
    GP.col(0) << 1, 0, 0, 0, 0, 1;
    GP.col(1) << 0, 1, 0, 1, 0, 0;
    const MyArr2D second_site = extract_gp_chunk(GP, 1, 1);

    const std::string path = (std::filesystem::temp_directory_path() / "phaseless-test-output.vcf").string();
    std::filesystem::remove(path);
    {
        auto bw = make_bcfwriter(path, {"7"}, {"sample0", "sample1"});
        write_bigass_to_bcf(bw, second_site.data(), "7", {202}, {"rs202"}, {"A"}, {"G"});
    }

    vcfpp::BcfReader reader(path);
    vcfpp::BcfRecord record(reader.header);
    REQUIRE(reader.getNextVariant(record));
    REQUIRE(record.CHROM() == "7");
    REQUIRE(record.POS() == 202);
    REQUIRE(record.ID() == "rs202");
    REQUIRE(record.REF() == "A");
    REQUIRE(record.ALT() == "G");
    Float1D gp;
    record.getFORMAT("GP", gp);
    const Float1D expected{0, 0, 1, 1, 0, 0};
    REQUIRE(gp == expected);
    REQUIRE_FALSE(reader.getNextVariant(record));
    std::filesystem::remove(path);
}

TEST_CASE("test beagle parser for bigass genome with unbalanced chunksize - small", "[test-io]")
{
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = 87;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    REQUIRE(genome->nsamples == 60);
    REQUIRE(genome->nsnps == 1000);
    REQUIRE(genome->nchunks == genome->pos.size());
    // for(int i = 0; i < genome->nchunks; i++) cout << genome->chrs[i] << ": " << genome->pos[i].size() <<
    // endl;
}

TEST_CASE("test beagle parser for bigass genome with unbalanced chunksize - big", "[test-io]")
{
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = 107;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    REQUIRE(genome->nsamples == 60);
    REQUIRE(genome->nsnps == 1000);
    REQUIRE(genome->nchunks == genome->pos.size());
    for(int i = 0; i < genome->nchunks; i++) cout << genome->chrs[i] << ": " << genome->pos[i].size() << endl;
}
