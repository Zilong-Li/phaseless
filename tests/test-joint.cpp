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
