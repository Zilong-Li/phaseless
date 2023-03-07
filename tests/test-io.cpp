#include "../src/io.hpp"
#include "catch.hh"

using namespace std;

TEST_CASE("test beagle parser for single contig", "[test-io]")
{
    int N, M;
    MyFloat1D genolikes;
    MapStringInt1D chrs_pos;
    StringVec1D sampleids;
    read_beagle_genotype_likelihoods("../data/bgl.gz", genolikes, sampleids, chrs_pos, N, M);
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
    StringVec1D sampleids;
    read_beagle_genotype_likelihoods("../data/all.bgl.gz", genolikes, sampleids, chrs_pos, N, M);
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
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    REQUIRE(genome->nsamples == 60);
    REQUIRE(genome->nsnps == 1000);
    REQUIRE(genome->nchunks == genome->pos.size());
    for(int i = 0; i < genome->nchunks; i++) REQUIRE(genome->pos[i].size() == 100);
}

TEST_CASE("test beagle parser for bigass genome with unbalanced chunksize - small", "[test-io]")
{
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = 87;
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    REQUIRE(genome->nsamples == 60);
    REQUIRE(genome->nsnps == 1000);
    REQUIRE(genome->nchunks == genome->pos.size());
    // for(int i = 0; i < genome->nchunks; i++) cout << genome->chrs[i] << ": " << genome->pos[i].size() << endl;
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
