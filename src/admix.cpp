#include "admixture.hpp"
#include "io.hpp"
#include "log.hpp"
#include "threadpool.hpp"
#include "timer.hpp"
#include <alpaca/alpaca.h>

using namespace std;

inline int run_bootstrap(const std::unique_ptr<BigAss> & genome,
                         ThreadPool & poolit,
                         Logger & cao,
                         int K,
                         int nseeds)
{
    std::vector<std::future<double>> res;
    std::vector<double> llikes;
    for(int seed = 0; seed < nseeds; seed++)
    {
        Admixture admixer(genome->nsamples, genome->nsnps, genome->C, K, seed);
        for(int i = 0; i < genome->nsamples; i++)
            res.emplace_back(poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
        double loglike = 0;
        for(auto && ll : res) loglike += ll.get();
        res.clear(); // clear future and renew
        llikes.push_back(loglike);
    }
    int seed = 0;
    for(auto ll : llikes) cao.print("seed =", seed++, ", likelihoods =", ll);
    auto it = std::max_element(llikes.begin(), llikes.end());
    return std::distance(llikes.begin(), it);
}

inline int run_admix_main(Options & opts)
{
    cao.cao.open(opts.out.string() + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running admixture");
    int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = opts.nthreads < allthreads ? opts.nthreads : allthreads;
    cao.print(tim.date(), allthreads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool poolit(opts.nthreads);
    // Deserialize from file
    auto filesize = std::filesystem::file_size(opts.in_bin);
    std::error_code ec;
    std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome =
        std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    cao.done(tim.date(), filesize, " bytes deserialized from file. skip imputation, ec", ec);
    cao.print(tim.date(), "parsing input -> C =", genome->C, ", N =", genome->nsamples,
              ", M =", genome->nsnps, ", nchunks =", genome->nchunks);
    assert(opts.K < genome->C);

    cao.warn(tim.date(), "-> running admixture with seed =", opts.seed);
    Admixture admixer(genome->nsamples, genome->nsnps, genome->C, opts.K, opts.seed);
    vector<future<double>> llike;
    if(!opts.noaccel)
    {
        MyArr2D F0, Q0, F1, Q1;
        const int istep{4};
        double alpha, qdiff, ldiff, stepMax{4}, alphaMax{1280};
        double prevlike{std::numeric_limits<double>::lowest()};
        for(int it = 0; it < opts.nadmix / 3; it++)
        {
            // first accel iteration
            admixer.initIteration();
            F0 = admixer.F;
            Q0 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(
                    poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            for(auto && ll : llike) ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            // second accel iteration
            admixer.initIteration();
            F1 = admixer.F;
            Q1 = admixer.Q;
            qdiff = (Q1 - Q0).square().sum();
            tim.clock();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(
                    poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            double loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            ldiff = loglike - prevlike;
            prevlike = loglike;
            cao.print(tim.date(), "SqS3 iteration", it * 3 + 1, ", diff(Q) =", std::scientific, qdiff,
                      ", likelihoods =", std::fixed, loglike, ",", tim.reltime(), " sec");
            if(qdiff < opts.qtol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff(Q) =", std::scientific, qdiff, " <",
                          opts.qtol);
                break;
            }
            if(ldiff < opts.ltol)
            {
                cao.print(tim.date(), "hit stopping criteria, diff(loglikelihood) =", std::scientific, ldiff,
                          " <", opts.ltol);
                break;
            }
            // accel iteration with steplen
            admixer.initIteration();
            alpha = ((F1 - F0).square().sum() + (Q1 - Q0).square().sum())
                    / ((admixer.F - 2 * F1 + F0).square().sum() + (admixer.Q - 2 * Q1 + Q0).square().sum());
            alpha = max(1.0, sqrt(alpha));
            if(alpha >= stepMax)
            {
                alpha = min(stepMax, alphaMax);
                stepMax = min(stepMax * istep, alphaMax);
            }
            admixer.F = F0 + 2 * alpha * (F1 - F0) + alpha * alpha * (admixer.F - 2 * F1 + F0);
            admixer.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (admixer.Q - 2 * Q1 + Q0);
            admixer.initIteration();
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(
                    poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));
            for(auto && ll : llike) ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
        }
    }
    else
    {
        MyArr2D Q0;
        double qdiff, ldiff, loglike;
        double prevlike{std::numeric_limits<double>::lowest()};
        for(int it = 0; it < opts.nadmix; it++)
        {
            tim.clock();
            admixer.initIteration();
            Q0 = admixer.Q;
            for(int i = 0; i < genome->nsamples; i++)
                llike.emplace_back(
                    poolit.enqueue(&Admixture::runOptimalWithBigAss, &admixer, i, std::ref(genome)));

            loglike = 0;
            for(auto && ll : llike) loglike += ll.get();
            llike.clear(); // clear future and renew
            admixer.updateIteration();
            ldiff = loglike - prevlike;
            prevlike = loglike;
            qdiff = (admixer.Q - Q0).square().sum();
            cao.print(tim.date(), "normal iteration", it, ", diff(Q) =", std::scientific, qdiff,
                      ", likelihoods =", std::fixed, loglike, ",", tim.reltime(), " sec");
            if(qdiff < opts.qtol) break;
            if(ldiff < opts.ltol) break;
        }
    }
    cao.done(tim.date(), "admixture done and outputting");
    admixer.writeQ(opts.out.string() + ".Q");
    cao.done(tim.date(), "-> good job. have a nice day, bye!");

    return 0;
}
