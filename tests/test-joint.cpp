#define _DECLARE_TOOLBOX_HERE

#include "../src/io.hpp"
#include "../src/phaseless.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("joint ancestry responsibility uses ancestry-specific cluster frequency", "[test-joint]")
{
    const int K{2}, C{2}, N{1}, M{1};
    Phaseless faith(K, C, N, M, 1);
    faith.pos_chunk = {0, M};
    faith.Q << 0.5, 0.5;
    faith.F[0] << 0.9, 0.1;
    faith.F[1] << 0.1, 0.9;
    faith.P.setConstant(0.5);
    faith.initIteration();

    MyArr2D gli = MyArr2D::Ones(M, 3);
    MyArr2D emit = MyArr2D::Ones(diploid_unordered_state_count(C), M);
    MyArr2D H(C, M);
    H << 0.5, 0.5;
    MyArr1D cs = MyArr1D::Ones(M);
    MyArr2D alpha = MyArr2D::Zero(diploid_unordered_state_count(C), M);
    alpha(0, 0) = 0.8;
    alpha(diploid_unordered_state_index(1, 1, C), 0) = 0.2;
    MyArr2D beta = MyArr2D::Ones(diploid_unordered_state_count(C), M);

    faith.getPosterios(0, 0, gli, emit, H, cs, alpha, beta, false);

    REQUIRE(faith.Eancestry(0, 0) == Approx(0.74));
    REQUIRE(faith.Eancestry(1, 0) == Approx(0.26));
    REQUIRE(faith.EclusterK.col(0).sum() == Approx(1.0));
    REQUIRE(faith.EclusterUsage.col(0).sum() == Approx(2.0));

    faith.NF = faith.NP = faith.NR = true;
    faith.updateIteration();
    REQUIRE(faith.Q(0, 0) == Approx(0.74));
    REQUIRE(faith.Q(1, 0) == Approx(0.26));
}

TEST_CASE("joint ancestry counts only cluster-refresh events", "[test-joint]")
{
    const int K{2}, C{2}, N{1}, M{2};
    Phaseless faith(K, C, N, M, 1);
    faith.pos_chunk = {0, M};
    faith.Q << 0.6, 0.4;
    faith.F[0].setConstant(0.5);
    faith.F[1].setConstant(0.5);
    faith.P.setConstant(0.5);
    faith.R.setZero();
    faith.R.row(0).setOnes(); // both copies remain in their previous clusters
    faith.initIteration();

    MyArr2D gli = MyArr2D::Ones(M, 3);
    MyArr2D emit = MyArr2D::Ones(diploid_unordered_state_count(C), M);
    MyArr2D H = MyArr2D::Constant(C, M, 0.5);
    MyArr1D cs = MyArr1D::Ones(M);
    MyArr2D alpha = MyArr2D::Zero(diploid_unordered_state_count(C), M);
    alpha.row(0).setOnes();
    MyArr2D beta = MyArr2D::Ones(diploid_unordered_state_count(C), M);

    faith.getPosterios(0, 0, gli, emit, H, cs, alpha, beta, false);

    REQUIRE(faith.Eancestry.col(0).sum() == Approx(1.0));
    REQUIRE(faith.EclusterK.col(0).sum() == Approx(1.0));
    REQUIRE(faith.EclusterK.col(1).sum() == Approx(0.0).margin(1e-12));
    REQUIRE(faith.EclusterUsage.col(0).sum() == Approx(2.0));
    REQUIRE(faith.EclusterUsage.col(1).sum() == Approx(2.0));
}

TEST_CASE("joint Q update applies a symmetric ancestry pseudocount", "[test-joint]")
{
    Phaseless faith(3, 2, 1, 1, 5);
    faith.NP = faith.NF = faith.NR = true;
    faith.NQ = false;
    faith.Eancestry.resize(3, 1);
    faith.Eancestry << 9.0, 1.0, 0.0;
    faith.setAdmixturePseudocount(0.5);

    faith.updateIteration();

    REQUIRE(faith.Q(0, 0) == Approx(9.5 / 11.5));
    REQUIRE(faith.Q(1, 0) == Approx(1.5 / 11.5));
    REQUIRE(faith.Q(2, 0) == Approx(0.5 / 11.5));
    REQUIRE(faith.Q.col(0).sum() == Approx(1.0));
    REQUIRE_THROWS_AS(faith.setAdmixturePseudocount(-0.1), std::invalid_argument);
}

TEST_CASE("joint P update shrinks low-occupancy cells toward pooled site frequency", "[test-joint]")
{
    Phaseless faith(1, 2, 1, 2, 7);
    faith.initIteration();
    faith.NQ = faith.NF = faith.NR = true;
    faith.NP = false;
    faith.EclusterA1 << 9.0, 0.0,
                        0.0, 10.0;
    faith.EclusterA2 << 1.0, 10.0,
                        0.0, 0.0;
    faith.setEmissionShrinkage(1.0);

    faith.updateIteration();

    REQUIRE(faith.P(0, 0) == Approx(0.1));
    REQUIRE(faith.P(0, 1) == Approx(0.1));
    REQUIRE(faith.P(1, 0) == Approx(10.5 / 11.0));
    REQUIRE(faith.P(1, 1) == Approx(0.5 / 11.0));
    REQUIRE_THROWS_AS(faith.setEmissionShrinkage(-0.1), std::invalid_argument);
}

TEST_CASE("joint STITCH heuristic aligns swapped cluster labels", "[test-joint]")
{
    constexpr int K{2}, C{2}, N{2}, M{4};
    Phaseless faith(K, C, N, M, 7);
    faith.pos_chunk = {0, M};
    faith.Q << 1.0, 0.0, 0.0, 1.0;
    faith.F[0] << 0.9, 0.9, 0.1, 0.1, 0.1, 0.1, 0.9, 0.9;
    faith.F[1] << 0.1, 0.1, 0.9, 0.9, 0.9, 0.9, 0.1, 0.1;
    faith.P << 0.1, 0.8, 0.2, 0.7, 0.7, 0.2, 0.8, 0.1;

    const JointHeuristicReport report = faith.alignClusterLabels(2, 0);

    REQUIRE(report.relabelled_boundaries == 1);
    REQUIRE(report.reset_sites == 1);
    REQUIRE(faith.F[0](0, 2) == Approx(0.9));
    REQUIRE(faith.F[1](0, 2) == Approx(0.1));
    REQUIRE(faith.P(3, 0) == Approx(0.1));
    REQUIRE(faith.P(3, 1) == Approx(0.8));
}

TEST_CASE("joint initialization derives Q and F from posterior cluster profiles", "[test-joint]")
{
    constexpr int K{2}, C{2}, N{4}, M{3};
    Phaseless faith(K, C, N, M, 23);
    faith.initializeSharedHaplotypeStart();
    for(int k = 1; k < K; ++k) REQUIRE((faith.F[k] - faith.F[0]).abs().maxCoeff() == Approx(0.0));
    faith.EindividualClusterUsage.resize(C, N);
    faith.EindividualClusterUsage << 95.0, 90.0, 10.0, 5.0,
                                      5.0, 10.0, 90.0, 95.0;

    REQUIRE(faith.initializeAncestryFromPosterior(0.05, 0));
    REQUIRE((faith.F[1] - faith.F[0]).abs().maxCoeff() > 0.0);
    REQUIRE(std::abs(faith.Q(0, 0) - faith.Q(0, 1)) < 0.1);
    REQUIRE(std::abs(faith.Q(0, 2) - faith.Q(0, 3)) < 0.1);
    REQUIRE(std::abs(faith.Q(0, 0) - faith.Q(0, 3)) > 0.8);
    for(int i = 0; i < N; ++i) REQUIRE(faith.Q.col(i).sum() == Approx(1.0));
    for(int k = 0; k < K; ++k)
        for(int m = 0; m < M; ++m) REQUIRE(faith.F[k].col(m).sum() == Approx(1.0));
}

TEST_CASE("low-cost initialization profile pruning retains cluster-informative SNPs", "[test-joint]")
{
    constexpr int K{1}, C{2}, N{10}, M{6};
    Phaseless faith(K, C, N, M, 31);
    faith.pos_chunk = {0, M};
    faith.EclusterUsage.setConstant(C, M, 10.0);
    faith.P << 0.50, 0.50,
               0.10, 0.90,
               0.20, 0.80,
               0.50, 0.50,
               0.01, 0.99,
               0.45, 0.55;

    const auto report = faith.configureInformativeProfileSites(3, 0.80, 1);

    REQUIRE_FALSE(report.used_all_sites_fallback);
    REQUIRE(report.blocks == 2);
    REQUIRE(report.informative_blocks == 2);
    REQUIRE(report.informative_sites == 4);
    REQUIRE(report.retained_sites == 3);
    REQUIRE(faith.initializationProfileWeights(0) == Approx(0.0));
    REQUIRE(faith.initializationProfileWeights(1) == Approx(1.5));
    REQUIRE(faith.initializationProfileWeights(2) == Approx(1.5));
    REQUIRE(faith.initializationProfileWeights(3) == Approx(0.0));
    REQUIRE(faith.initializationProfileWeights(4) == Approx(3.0));
    REQUIRE(faith.initializationProfileWeights(5) == Approx(0.0));

    faith.P.setConstant(0.5);
    const auto fallback = faith.configureInformativeProfileSites(3, 0.95, 1);
    REQUIRE(fallback.used_all_sites_fallback);
    REQUIRE((faith.initializationProfileWeights == 1.0).all());
}

TEST_CASE("initialization profile weights affect only per-individual profile aggregation", "[test-joint]")
{
    constexpr int K{1}, C{2}, N{1}, M{3};
    Phaseless faith(K, C, N, M, 37);
    faith.pos_chunk = {0, M};
    faith.initializationProfileWeights.resize(M);
    faith.initializationProfileWeights << 0.0, 3.0, 0.0;
    faith.initIteration();

    const MyArr2D gli = MyArr2D::Ones(M, 3);
    const MyArr2D emit = MyArr2D::Ones(diploid_unordered_state_count(C), M);
    const MyArr2D H = MyArr2D::Constant(C, M, 0.5);
    const MyArr1D cs = MyArr1D::Ones(M);
    MyArr2D alpha = MyArr2D::Zero(diploid_unordered_state_count(C), M);
    alpha.row(diploid_unordered_state_index(0, 0, C)).setOnes();
    const MyArr2D beta = MyArr2D::Ones(diploid_unordered_state_count(C), M);

    faith.getPosterios(0, 0, gli, emit, H, cs, alpha, beta, false);

    REQUIRE(faith.EindividualClusterUsage(0, 0) == Approx(6.0));
    REQUIRE(faith.EindividualClusterUsage(1, 0) == Approx(0.0));
    REQUIRE((faith.EclusterUsage.row(0) == 2.0).all());
    REQUIRE((faith.EclusterUsage.row(1) == 0.0).all());
}

TEST_CASE("posterior initialization restarts perturb the clustering view", "[test-joint]")
{
    constexpr int K{3}, C{4}, N{8}, M{2};
    Phaseless first(K, C, N, M, 41);
    Phaseless second(K, C, N, M, 41);
    MyArr2D profiles(C, N);
    profiles << 70, 60, 45, 35, 25, 15, 10, 5,
                10, 20, 35, 45, 50, 55, 20, 10,
                15, 15, 10, 10, 15, 20, 55, 70,
                 5,  5, 10, 10, 10, 10, 15, 15;
    first.initializeSharedHaplotypeStart();
    second.initializeSharedHaplotypeStart();
    first.EindividualClusterUsage = profiles;
    second.EindividualClusterUsage = profiles;

    REQUIRE(first.initializeAncestryFromPosterior(0.05, 1));
    REQUIRE(second.initializeAncestryFromPosterior(0.05, 2));
    REQUIRE((first.Q - second.Q).abs().maxCoeff() > 1e-6);
    for(int i = 0; i < N; ++i)
    {
        REQUIRE(first.Q.col(i).sum() == Approx(1.0));
        REQUIRE(second.Q.col(i).sum() == Approx(1.0));
    }
}

TEST_CASE("phase-stage alignment uses centered posterior occupancy profiles", "[test-joint]")
{
    constexpr int K{1}, C{2}, N{4}, M{4};
    Phaseless faith(K, C, N, M, 29);
    faith.pos_chunk = {0, M};
    faith.configurePhaseAlignment(2);
    REQUIRE(faith.phaseAlignmentBoundaries == Int1D{2});

    faith.P << 0.1, 0.8, 0.2, 0.7, 0.7, 0.2, 0.8, 0.1;
    faith.F[0] << 0.8, 0.7, 0.2, 0.1, 0.2, 0.3, 0.8, 0.9;
    faith.EclusterUsage.setZero(C, M);
    faith.EclusterUsage.col(1) << 4.0, 4.0;
    faith.EclusterUsage.col(2) << 4.0, 4.0;
    faith.EphaseAlignmentCross[0] << 2.0, 6.0, 6.0, 2.0;
    faith.EphaseAlignmentSquares.col(0).setConstant(8.0);

    const JointHeuristicReport report = faith.alignPhaseClusterLabels(0);

    REQUIRE(report.relabelled_boundaries == 1);
    REQUIRE(report.reset_sites == 1);
    REQUIRE(faith.P(3, 0) == Approx(0.1));
    REQUIRE(faith.P(3, 1) == Approx(0.8));
    REQUIRE(faith.F[0](0, 3) == Approx(0.9));
    REQUIRE(faith.F[0](1, 3) == Approx(0.1));
}

TEST_CASE("joint STITCH heuristic revives a low-usage cluster with donor plus noise", "[test-joint]")
{
    constexpr int K{1}, C{2}, N{10}, M{6};
    Phaseless faith(K, C, N, M, 11);
    faith.pos_chunk = {0, M};
    faith.P.col(0).setConstant(0.25);
    faith.P.col(1).setConstant(0.95);
    faith.EclusterUsage.resize(C, M);
    faith.EclusterUsage.row(0).setConstant(19.98);
    faith.EclusterUsage.row(1).setConstant(0.02);

    const JointHeuristicReport report = faith.reviveUnusedClusters(0.005, 2, 0.8);

    REQUIRE(report.revived_intervals == 1);
    REQUIRE(report.revived_sites == M);
    for(int m = 0; m < M; ++m)
    {
        REQUIRE(faith.P(m, 1) >= 0.2);
        REQUIRE(faith.P(m, 1) <= 0.4);
    }
}

TEST_CASE("joint sufficient statistics reduce individuals deterministically", "[test-joint]")
{
    constexpr int K{1}, C{1}, N{3}, M{1};
    Phaseless ordered(K, C, N, M, 19);
    Phaseless reversed(K, C, N, M, 19);
    ordered.pos_chunk = reversed.pos_chunk = {0, M};
    ordered.Q.setOnes();
    reversed.Q.setOnes();
    ordered.F[0].setOnes();
    reversed.F[0].setOnes();
    ordered.P.setConstant(0.5);
    reversed.P.setConstant(0.5);
    ordered.initIteration();
    reversed.initIteration();

    const MyArr2D gli = MyArr2D::Ones(M, 3);
    const MyArr2D emit = MyArr2D::Ones(1, M);
    const MyArr2D H = MyArr2D::Ones(C, M);
    const MyArr1D cs = MyArr1D::Ones(M);
    const MyArr2D beta = MyArr2D::Ones(1, M);
    std::array<MyArr2D, N> alpha{MyArr2D::Constant(1, M, 1e16), MyArr2D::Ones(1, M),
                                 MyArr2D::Ones(1, M)};

    for(int ind = 0; ind < N; ++ind)
        ordered.getPosterios(ind, 0, gli, emit, H, cs, alpha[ind], beta, false);

    std::array<std::future<void>, N> futures;
    for(int launch = 0; launch < N; ++launch)
    {
        const int ind = N - 1 - launch;
        futures[launch] = std::async(std::launch::async, [&, ind]
        { reversed.getPosterios(ind, 0, gli, emit, H, cs, alpha[ind], beta, false); });
    }
    for(auto & result : futures) result.get();

    REQUIRE((reversed.EclusterUsage == ordered.EclusterUsage).all());
    REQUIRE((reversed.EclusterA1 == ordered.EclusterA1).all());
    REQUIRE((reversed.EclusterA2 == ordered.EclusterA2).all());
    REQUIRE((reversed.EclusterK == ordered.EclusterK).all());
}

TEST_CASE("symmetric diploid recursion matches ordered recursion", "[test-joint]")
{
    constexpr int C{3}, M{4};
    MyArr2D gli(M, 3), P(M, C), PI(C, M), R(3, M);
    gli << 0.75, 0.20, 0.05, 0.10, 0.65, 0.25, 0.30, 0.45, 0.25, 0.05, 0.25, 0.70;
    P << 0.10, 0.45, 0.80, 0.20, 0.55, 0.75, 0.15, 0.60, 0.85, 0.25, 0.50, 0.90;
    PI << 0.50, 0.20, 0.35, 0.45, 0.30, 0.50, 0.25, 0.20, 0.20, 0.30, 0.40, 0.35;
    R << 1.0, 0.64, 0.49, 0.81, 0.0, 0.16, 0.21, 0.09, 0.0, 0.04, 0.09, 0.01;

    const MyArr2D ordered_emit = get_emission_by_gl(gli, P);
    const MyArr2D symmetric_emit = get_emission_by_gl_symmetric(gli, P);
    const auto [ordered_alpha, ordered_beta, ordered_cs] = forward_backwards_diploid(ordered_emit, R, PI);
    const auto [symmetric_alpha, symmetric_beta, symmetric_cs] =
        forward_backwards_diploid_symmetric(symmetric_emit, R, PI);

    REQUIRE((symmetric_cs - ordered_cs).abs().maxCoeff() < 1e-12);
    for(int z1 = 0; z1 < C; ++z1)
        for(int z2 = z1; z2 < C; ++z2)
        {
            const int unordered = diploid_unordered_state_index(z1, z2, C);
            const int ordered = z1 * C + z2;
            REQUIRE((symmetric_emit.row(unordered) - ordered_emit.row(ordered)).abs().maxCoeff() < 1e-12);
            REQUIRE((symmetric_alpha.row(unordered) - ordered_alpha.row(ordered)).abs().maxCoeff() < 1e-12);
            REQUIRE((symmetric_beta.row(unordered) - ordered_beta.row(ordered)).abs().maxCoeff() < 1e-12);
        }

    for(int s = 0; s < M; ++s)
    {
        double posterior_sum = 0;
        for(int z1 = 0; z1 < C; ++z1)
            for(int z2 = z1; z2 < C; ++z2)
            {
                const int state = diploid_unordered_state_index(z1, z2, C);
                posterior_sum +=
                    (z1 == z2 ? 1.0 : 2.0) * symmetric_alpha(state, s) * symmetric_beta(state, s);
            }
        REQUIRE(posterior_sum == Approx(1.0).margin(1e-12));
    }
}

TEST_CASE("symmetric off-diagonal state uses context-specific multiplicity", "[test-joint]")
{
    constexpr int K{1}, C{2}, N{1}, M{1};
    Phaseless faith(K, C, N, M, 1);
    faith.pos_chunk = {0, M};
    faith.Q.setOnes();
    faith.F[0].setConstant(0.5);
    faith.P << 0.2, 0.8;
    faith.initIteration();
    faith.GP.setZero(M * 3, N);

    const MyArr2D gli = MyArr2D::Ones(M, 3);
    const MyArr2D emit = MyArr2D::Ones(diploid_unordered_state_count(C), M);
    const MyArr2D H = MyArr2D::Constant(C, M, 0.5);
    const MyArr1D cs = MyArr1D::Ones(M);
    MyArr2D alpha = MyArr2D::Zero(diploid_unordered_state_count(C), M);
    alpha(diploid_unordered_state_index(0, 1, C), 0) = 0.5;
    const MyArr2D beta = MyArr2D::Ones(diploid_unordered_state_count(C), M);

    faith.getPosterios(0, 0, gli, emit, H, cs, alpha, beta, true);

    REQUIRE(faith.EclusterA1(0, 0) == Approx(0.4));
    REQUIRE(faith.EclusterA1(1, 0) == Approx(0.1));
    REQUIRE(faith.EclusterA2(0, 0) == Approx(0.1));
    REQUIRE(faith.EclusterA2(1, 0) == Approx(0.4));
    REQUIRE(faith.EclusterK.col(0).sum() == Approx(1.0));
    REQUIRE(faith.GP(0, 0) == Approx(0.16));
    REQUIRE(faith.GP(1, 0) == Approx(0.68));
    REQUIRE(faith.GP(2, 0) == Approx(0.16));
}

TEST_CASE("phaseless joint single chunk", "[test-joint]")
{
    cao.cerr("TEST: phaseless joint single chunk");
    int K{3}, C{5}, seed{1}, chunksize{INT_MAX}, nimpute{10};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    Phaseless faith(K, C, genome->nsamples, genome->nsnps, seed);
    faith.initRecombination(genome->pos, "", 20000, 1);
    ThreadPool pool(4);
    vector<future<double>> res;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; it <= nimpute; it++)
    {
        tim.clock();
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
            res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
        loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        diff = it ? loglike - prevlike : 0;
        res.clear(); // clear future and renew
        cao.cerr(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff = ", diff, ", time",
                 tim.reltime(), " sec");
        faith.updateIteration();
        prevlike = loglike;
    }
}

TEST_CASE("phaseless joint multiple chunks", "[test-joint]")
{
    cao.cerr("TEST: phaseless joint multiple chunks");
    int K{3}, C{5}, seed{1}, chunksize{INT_MAX}, nimpute{10};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C, genome->B = 1;
    tim.clock();
    chunk_beagle_genotype_likelihoods(genome, "../data/all.bgl.gz");
    cao.print(tim.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples, ", M =", genome->nsnps,
              ", nchunks =", genome->nchunks);
    Phaseless faith(K, C, genome->nsamples, genome->nsnps, seed);
    faith.initRecombination(genome->pos, "");
    ThreadPool pool(4);
    vector<future<double>> res;
    double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
    for(int it = 0; it <= nimpute; it++)
    {
        tim.clock();
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
            res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), false));
        loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        diff = it ? loglike - prevlike : 0;
        cao.cerr(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff = ", diff, ", time",
                 tim.reltime(), " sec");
        faith.updateIteration();
        prevlike = loglike;
    }
}
