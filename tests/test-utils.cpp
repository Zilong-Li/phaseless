#define _DECLARE_TOOLBOX_HERE

#include "../src/common.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("calc_grid_distance", "[test-utils]")
{
    // Int1D pos{1,2,3,4,5,6,7,8,9,10,11};
    Int1D pos{1, 2, 3, 4, 5, 6, 7, 8, 9};
    const auto b = find_grid_to_collapse(pos, 4);
    
    Char1D cb(b.size());
    for(int i = 0; i < b.size();i++) cb[i] = (b(i)==true);
    for(auto i : cb) cao.cerr(i==1);
    auto d = calc_grid_distance(pos, b);
    // cao.cerr("size:", d.size());
}

TEST_CASE("runtime safeguards", "[test-utils]")
{
    REQUIRE(resolve_thread_count(8, 4) == 4);
    REQUIRE(resolve_thread_count(8, 0) == 1);
    REQUIRE_THROWS_AS(resolve_thread_count(0, 4), std::invalid_argument);

    REQUIRE(likelihood_converged(1e-5, 1e-4));
    REQUIRE_FALSE(likelihood_converged(-1e-5, 1e-4));
    REQUIRE_FALSE(likelihood_converged(NAN, 1e-4));
}

TEST_CASE("allele frequency EM normalizes over samples", "[test-utils]")
{
    const int N = 2, M = 2;
    // Sample-major blocks, with each genotype block containing M sites.
    const MyFloat1D GL{1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0};
    const Arr1D af = estimate_af_by_gl(GL, N, M);
    REQUIRE(af(0) == Approx(0.5));
    REQUIRE(af(1) == Approx(0.5));
}
