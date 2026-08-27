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
    REQUIRE(resolve_thread_count(-1, 4) == 4);
    REQUIRE(resolve_thread_count(-1, 0) == 1);
    REQUIRE_THROWS_AS(resolve_thread_count(0, 4), std::invalid_argument);
    REQUIRE_THROWS_AS(resolve_thread_count(-2, 4), std::invalid_argument);

    REQUIRE(likelihood_converged(1e-5, 1e-4));
    REQUIRE_FALSE(likelihood_converged(-1e-5, 1e-4));
    REQUIRE_FALSE(likelihood_converged(NAN, 1e-4));
}

TEST_CASE("whole-genome likelihood convergence metrics", "[test-utils]")
{
    const auto metrics = assess_likelihood_convergence(-105, -110, -120, 100);
    REQUIRE(metrics.delta == Approx(5));
    REQUIRE(metrics.relative_change == Approx(5.0 / 105));
    REQUIRE(metrics.aitken_rate == Approx(0.5));
    REQUIRE(metrics.improvement_per_observation == Approx(0.05));
    REQUIRE(metrics.aitken_gap_per_observation == Approx(0.05));
    REQUIRE(metrics.monotone);
    REQUIRE(metrics.aitken_valid);

    const auto decreased = assess_likelihood_convergence(-111, -110, -120, 100);
    REQUIRE_FALSE(decreased.monotone);
    REQUIRE_FALSE(decreased.aitken_valid);
    REQUIRE(decreased.improvement_per_observation == Approx(0.01));
    REQUIRE(std::isnan(decreased.aitken_gap_per_observation));

    const auto noncontractive = assess_likelihood_convergence(-100, -105, -109, 100);
    REQUIRE(noncontractive.aitken_rate == Approx(1.25));
    REQUIRE_FALSE(noncontractive.aitken_valid);
    REQUIRE(noncontractive.improvement_per_observation == Approx(0.05));
    REQUIRE(std::isnan(noncontractive.aitken_gap_per_observation));
}

TEST_CASE("adaptive shared-haplotype convergence uses likelihood and posterior profiles", "[test-utils]")
{
    REQUIRE(initialization_converged(5e-5, 1e-3, 1e-4, 2e-3));
    REQUIRE_FALSE(initialization_converged(2e-4, 1e-3, 1e-4, 2e-3));
    REQUIRE_FALSE(initialization_converged(5e-5, 3e-3, 1e-4, 2e-3));
    REQUIRE_FALSE(initialization_converged(NAN, 1e-3, 1e-4, 2e-3));
}

TEST_CASE("SqS3 handoff uses long rejection and efficiency windows", "[test-utils]")
{
    const auto too_short = assess_sqs3_handoff(8, 4, 8, 0.0, 1e-5, 1e-6);
    REQUIRE_FALSE(too_short.persistent_rejection);
    REQUIRE_FALSE(too_short.low_efficiency);

    const auto rejected = assess_sqs3_handoff(24, 5, 16, 4.0, 5e-5, 1e-6);
    REQUIRE(rejected.rejection_rate == Approx(5.0 / 24));
    REQUIRE(rejected.persistent_rejection);
    REQUIRE_FALSE(rejected.low_efficiency);

    const auto inefficient = assess_sqs3_handoff(24, 2, 16, 0.4, 2e-5, 1e-6);
    REQUIRE_FALSE(inefficient.persistent_rejection);
    REQUIRE(inefficient.mean_realized_gain == Approx(0.025));
    REQUIRE(inefficient.low_efficiency);
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
