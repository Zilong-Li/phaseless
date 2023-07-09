#define _DECLARE_TOOLBOX_HERE

#include "../src/io.hpp"
#include "../src/threadpool.hpp"
#include "catch.hh"
#include <mutex>

using namespace std;
using namespace Eigen;

/*
** this class works on Grids level
*/
class FastPhaseGrid
{
  private:
    std::mutex mutex_it; // in case of race condition

  public:
    FastPhaseGrid(const Int2D & pos, int m, int n, int c, int seed);
    ~FastPhaseGrid();

    // BOUNDING
    double minRate{0.1}, maxRate{100}; // threshold for R
    double clusterFreqThreshold{1e-4}; // threshold for PI
    double alleleEmitThreshold{1e-4}; // threshold for F
    double nGen, Ne;

    // FLAGS
    bool debug{0};

    // SHARED VARIBALES
    const int nGrids; // number of grids in total
    const int B; // number of snps in each grid
    const int M; // number of snps in total
    const int N, C, C2; // N: samples, C: cluster C2: C*C
    MyArr2D F; // M x C, cluster-specific allele frequence
    MyArr2D R; // 3 x M, jumping / recombination rate
    MyArr1D pi; // C, PI in first SNP
    MyArr2D PI; // C x nGrids, cluster frequency
    MyArr2D Ezj; // nGrids x M, E(Z=z,J=1|X,par), expectation of switch into state k
    MyArr2D Ezg1, Ezg2; // C x M
    MyArr2D GP; // N x (M x 3), genotype probabilies for all individuals
    Int1D dist; // physical position distance between central markers in two grids

    void initIteration();
    void updateIteration();
    auto forwardAndBackwardsHighRam(int, const MyFloat1D &, bool);
};

inline FastPhaseGrid::FastPhaseGrid(const Int2D & pos, int m, int n, int c, int seed)
: nGrids(pos.size()), B(pos[0].size()), M(m), N(n), C(c), C2(c * c)
{
    GP.setZero(M * 3, N);
    auto rng = std::default_random_engine{};
    rng.seed(seed);
    F = RandomUniform<MyArr2D, std::default_random_engine>(M, C, rng, alleleEmitThreshold,
                                                           1 - alleleEmitThreshold);
    PI = MyArr2D::Ones(C, nGrids);
    PI.rowwise() /= PI.colwise().sum(); // normalize it per grid
    Ne = 20000; // for human
    nGen = 4 * Ne / C;
    dist = calc_grid_distance(pos);
    R = calc_transRate_diploid(dist, nGen);
}

inline FastPhaseGrid::~FastPhaseGrid() {}

inline void FastPhaseGrid::initIteration()
{
    if(debug) cao.cerr(R);
    // initial temp variables
    pi.setZero(C); // reset pi at first SNP
    Ezj.setZero(C, nGrids); // reset post(Z,j)
    Ezg1.setZero(C, M); // reset pos(Z,g)
    Ezg2.setZero(C, M); // reset pos(Z,g)
}

inline void FastPhaseGrid::updateIteration()
{
    MyArr1D er = 1.0 - Ezj.colwise().sum() / N;
    // morgans per SNP assuming T=100, 0.5 cM/Mb
    for(int i = 0; i < er.size(); i++)
    {
        double miner = std::exp(-nGen * maxRate * dist[i] / 1e8);
        double maxer = std::exp(-nGen * minRate * dist[i] / 1e8);
        er(i) = er(i) < miner ? miner : er(i);
        er(i) = er(i) > maxer ? maxer : er(i);
    }
    R.row(0) = er.square();
    R.row(1) = (1 - er) * er;
    R.row(2) = (1 - er).square();

    // update F
    F = (Ezg2 / (Ezg1 + Ezg2)).transpose();
    if(F.isNaN().any()) cao.cerr("NaN in F in FastPhaseGrid model. will fill it with AF");
    // map F to domain but no normalization
    F = (F < alleleEmitThreshold).select(alleleEmitThreshold, F); // lower bound
    F = (F > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, F); // upper bound

    // update PI(C, M) except the first snp
    // first we normalize Ezj so that each col sum to 1
    Ezj.col(0) = pi / pi.sum(); // now update the first SNP
    Ezj.rowwise() /= Ezj.colwise().sum();
    if(Ezj.isNaN().any() || (Ezj < clusterFreqThreshold).any())
    {
        // std::cerr << "reset values below threshold\n";
        Ezj = (Ezj < clusterFreqThreshold).select(0, Ezj); // reset to 0 first
        for(int i = 0; i < nGrids; i++)
        {
            // for columns with an entry below 0
            // each 0 entry becomes threshold
            // then rest re-scaled so whole thing has sum 1
            if(auto c = (Ezj.col(i) == 0).count() > 0)
            {
                double xsum = 1 - c * clusterFreqThreshold;
                double csum = Ezj.col(i).sum();
                Ezj.col(i) = (Ezj.col(i) > 0).select(Ezj.col(i) * xsum / csum, clusterFreqThreshold);
            }
        }
    }
    PI = Ezj;

    if(Ezj.isNaN().any()) cao.cerr(Ezj, "NaN in PI from FastPhaseGrid\n");
    if(debug && !((1 - PI.colwise().sum()).abs() < 1e-3).all())
        cao.cerr(PI.colwise().sum(), "\ncolsum of PI is not 1.0!\n");
}

/*
** @param ind       current individual i
** @param GL        genotype likelihood of all individuals in snp major form
** @param call_geno boolean, call genotype of not
** @return individual total likelihood
*/
inline auto FastPhaseGrid::forwardAndBackwardsHighRam(int ind, const MyFloat1D & GL, bool call_geno)
{
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * M * 3, M, 3);
    MyArr2D alpha(C2, nGrids), beta(C2, nGrids);
    MyArr2D emit = get_emission_by_gl(gli, F).transpose(); // C2 x M
    MyArr2D emitGrids = collapse_emission_by_grid(emit, B, nGrids);
    auto cs = forward_backwards_diploid(alpha, beta, emitGrids, R, PI);
    if(debug && !((1 - ((alpha * beta).colwise().sum())).abs() < 1e-4).all())
        cao.cerr((alpha * beta).colwise().sum() / cs.transpose(), "\ngamma sum is not 1.0!\n");

    MyArr1D gamma = alpha.col(0) * beta.col(0); //  gamma in first snp
    MyArr1D gamma1 = gamma.reshaped(C, C).colwise().sum();
    MyArr2D ind_post_zj = MyArr2D::Zero(C, nGrids);
    MyArr2D ind_post_zg1 = MyArr2D::Zero(C, M);
    MyArr2D ind_post_zg2 = MyArr2D::Zero(C, M);
    MyArr1D tmp_zg(4);

    int z1, z2, s, e, i, g = 0;
    // now get expectation of post(Z,J)
    s = g * B;
    e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
    for(z1 = 0; z1 < C; z1++)
    {
        for(i = s; i <= e; i++)
        {
            auto gamma_div_emit = gamma / emit.col(i);
            ind_post_zg1(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(i, z1))
                                   * (gli(i, 0) * (1 - F.row(i)) + gli(i, 1) * F.row(i)).transpose())
                                      .sum();
            ind_post_zg2(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(i, z1))
                                   * (gli(i, 1) * (1 - F.row(i)) + gli(i, 2) * F.row(i)).transpose())
                                      .sum();
            if(call_geno)
            {
                for(z2 = 0; z2 < C; z2++)
                {
                    int z12 = z1 * C + z2;
                    tmp_zg(0) = gli(i, 0) * (1 - F(i, z1)) * (1 - F(i, z2));
                    tmp_zg(1) = gli(i, 1) * (1 - F(i, z1)) * F(i, z2);
                    tmp_zg(2) = gli(i, 1) * F(i, z1) * (1 - F(i, z2));
                    tmp_zg(3) = gli(i, 2) * F(i, z1) * F(i, z2);
                    GP(3 * i + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                    GP(3 * i + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                    GP(3 * i + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
                }
            }
        }
    }
    for(g = 1; g < nGrids; g++)
    {
        const MyArr1D beta_mult_emit = emitGrids.col(g) * beta.col(g); // C2
        gamma = (alpha.col(g) * beta.col(g)); // C2
        MyArr1D alphatmp(C);
        s = g * B;
        e = g == nGrids - 1 ? M - 1 : B * (g + 1) - 1;
        for(z1 = 0; z1 < C; z1++)
        {
            alphatmp(z1) = alpha(Eigen::seqN(z1, C, C), g - 1).sum() * R(1, g);
            for(i = s; i <= e; i++)
            {
                auto gamma_div_emit = gamma / emit.col(i);
                ind_post_zg1(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - F(i, z1))
                                       * (gli(i, 0) * (1 - F.row(i)) + gli(i, 1) * F.row(i)).transpose())
                                          .sum();
                ind_post_zg2(z1, i) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (F(i, z1))
                                       * (gli(i, 1) * (1 - F.row(i)) + gli(i, 2) * F.row(i)).transpose())
                                          .sum();
                if(call_geno)
                {
                    for(z2 = 0; z2 < C; z2++)
                    {
                        int z12 = z1 * C + z2;
                        tmp_zg(0) = gli(i, 0) * (1 - F(i, z1)) * (1 - F(i, z2));
                        tmp_zg(1) = gli(i, 1) * (1 - F(i, z1)) * F(i, z2);
                        tmp_zg(2) = gli(i, 1) * F(i, z1) * (1 - F(i, z2));
                        tmp_zg(3) = gli(i, 2) * F(i, z1) * F(i, z2);
                        GP(3 * i + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
                        GP(3 * i + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
                        GP(3 * i + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
                    }
                }
            }
        }
        alphatmp += PI.col(g) * R(2, g) * 1.0;
        for(z1 = 0; z1 < C; z1++)
            ind_post_zj(z1, g) = cs(g) * (PI(z1, g) * alphatmp * beta_mult_emit(Eigen::seqN(z1, C, C))).sum();
    }

    double indLogLike = (1 / cs).log().sum();
    return std::tuple(indLogLike, ind_post_zj, ind_post_zg1, ind_post_zg2, gamma1);
}

TEST_CASE("FastPhaseGrid forwardAndBackwardsHighRam", "[test-collapse]")
{
    int C{10}, seed{1}, niters{20}, chunksize{10000}, B{5};
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    genome->chunksize = chunksize, genome->C = C;
    chunk_beagle_genotype_likelihoods(genome, "../data/bgl.gz");
    int ic = 0;
    auto pos = divide_pos_into_grid(genome->pos[ic], B);
    FastPhaseGrid faith(pos, genome->pos[ic].size(), genome->nsamples, C, seed);
    faith.debug = true;
    using pars1 = std::tuple<double, MyArr2D, MyArr2D, MyArr2D, MyArr1D>;
    vector<future<pars1>> res;
    double prevlike{std::numeric_limits<double>::lowest()}, loglike;
    ThreadPool poolit(4);
    for(int it = 0; it <= niters; it++)
    {
        faith.initIteration();
        for(int i = 0; i < genome->nsamples; i++)
        {
            if(it == niters)
                res.emplace_back(poolit.enqueue(&FastPhaseGrid::forwardAndBackwardsHighRam, &faith, i,
                                                std::ref(genome->gls[ic]), true));
            else
                res.emplace_back(poolit.enqueue(&FastPhaseGrid::forwardAndBackwardsHighRam, &faith, i,
                                                std::ref(genome->gls[ic]), false));
        }
        loglike = 0;
        for(auto && ll : res)
        {
            const auto [l, zj, zg1, zg2, gamma1] = ll.get();
            loglike += l;
            faith.Ezj += zj;
            faith.Ezg1 += zg1;
            faith.Ezg2 += zg2;
            faith.pi += gamma1;
        }
        res.clear(); // clear future and renew
        cao.cerr(tim.date(), "run single chunk", ic, ", iteration", it, ", likelihoods =", loglike);
        REQUIRE(loglike > prevlike);
        prevlike = loglike;
        if(it != niters) faith.updateIteration();
    }
}
