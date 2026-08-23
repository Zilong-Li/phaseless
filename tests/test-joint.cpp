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
    MyArr2D emit = MyArr2D::Ones(C * C, M);
    MyArr2D H(C, M);
    H << 0.5, 0.5;
    MyArr1D cs = MyArr1D::Ones(M);
    MyArr2D alpha = MyArr2D::Zero(C * C, M);
    alpha(0, 0) = 0.8;
    alpha(3, 0) = 0.2;
    MyArr2D beta = MyArr2D::Ones(C * C, M);

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
    MyArr2D emit = MyArr2D::Ones(C * C, M);
    MyArr2D H = MyArr2D::Constant(C, M, 0.5);
    MyArr1D cs = MyArr1D::Ones(M);
    MyArr2D alpha = MyArr2D::Zero(C * C, M);
    alpha.row(0).setOnes();
    MyArr2D beta = MyArr2D::Ones(C * C, M);

    faith.getPosterios(0, 0, gli, emit, H, cs, alpha, beta, false);

    REQUIRE(faith.Eancestry.col(0).sum() == Approx(1.0));
    REQUIRE(faith.EclusterK.col(0).sum() == Approx(1.0));
    REQUIRE(faith.EclusterK.col(1).sum() == Approx(0.0).margin(1e-12));
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
