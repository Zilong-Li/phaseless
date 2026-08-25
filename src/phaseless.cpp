/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/phaseless.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "phaseless.hpp"

#include "common.hpp"
#include "io.hpp"
#include "joint_cuda.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>

using namespace std;

namespace
{
struct JointParameterSnapshot
{
    MyArr2D Q;
    MyArr2D P;
    MyArr2D F;
    MyArr1D er;
};

struct JointParameterChange
{
    double q_max{0};
    double p_rms{0};
    double f_rms{0};
    double r_rms{0};
};

JointParameterSnapshot snapshot_parameters(const Phaseless & model)
{
    return {model.Q, model.P, cat_stdvec_of_eigen(model.F), model.er};
}

double rms_change(const MyArr2D & current, const MyArr2D & previous)
{
    return std::sqrt((current - previous).square().mean());
}

double rms_change(const MyArr1D & current, const MyArr1D & previous)
{
    return std::sqrt((current - previous).square().mean());
}

JointParameterChange parameter_change(const Phaseless & model, const JointParameterSnapshot & previous)
{
    JointParameterChange out;
    if(!model.NQ) out.q_max = (model.Q - previous.Q).abs().maxCoeff();
    if(!model.NP) out.p_rms = rms_change(model.P, previous.P);
    if(!model.NF) out.f_rms = rms_change(cat_stdvec_of_eigen(model.F), previous.F);
    if(!model.NR) out.r_rms = rms_change(model.er, previous.er);
    return out;
}
} // namespace

void Phaseless::initRecombination(const Int1D & pos, std::string rfile, int B, double Ne)
{
    nGen = 4 * Ne / C;
    dist = calc_position_distance(pos);
    if(rfile.empty())
        R = calc_transRate_diploid(dist, nGen);
    else
        load_csv(R, rfile, true);
    er = R.row(0).sqrt();
    protect_er(er);
    R = er2R(er);
}

void Phaseless::initRecombination(const Int2D & pos, std::string rfile, int B, double Ne)
{
    nGen = 4 * Ne / C;
    int nchunks = pos.size();
    pos_chunk.resize(nchunks + 1);
    int i{0}, ss{0};
    dist.reserve(M);
    for(i = 0; i < nchunks; i++)
    {
        pos_chunk[i] = ss;
        auto tmp = calc_position_distance(pos[i]);
        dist.insert(dist.end(), tmp.begin(), tmp.end());
        R.middleCols(ss, pos[i].size()) = calc_transRate_diploid(tmp, nGen);
        ss += pos[i].size();
    }
    pos_chunk[nchunks] = ss; // add sentinel
    if(!rfile.empty()) load_csv(R, rfile, true);
    er = R.row(0).sqrt();
    protect_er(er);
    R = er2R(er);
}

void Phaseless::setFlags(double tol_p, double tol_f, double tol_q, bool debug_, bool nQ, bool nP, bool nF, bool nR)
{
    alleleEmitThreshold = tol_p;
    clusterFreqThreshold = tol_f;
    admixtureThreshold = tol_q;
    debug = debug_;
    NQ = nQ;
    NP = nP;
    NF = nF;
    NR = nR;
}

void Phaseless::setStartPoint(std::string qfile, std::string pfile)
{
    if(!qfile.empty()) load_csv(Q, qfile, true);
    if(!pfile.empty()) load_csv(P, pfile, false);
}

void Phaseless::setStartPoint(const std::unique_ptr<Pars> & par)
{
    er = Eigen::Map<MyArr1D>(par->er.data(), M);
    R = er2R(er);
    Q = Eigen::Map<MyArr2D>(par->Q.data(), K, N);
    P = Eigen::Map<MyArr2D>(par->P.data(), M, C);
    for(int i = 0; i < K; i++) F[i] = Eigen::Map<MyArr2D>(par->F[i].data(), C, M);
}

void Phaseless::protectPars()
{
    // if we accelerate pars, protect them!
    if(!NQ)
    { // protect Q
        if(debug && Q.isNaN().any()) cao.error("NaN in Q in Phaseless model. reset it to the threshold");
        Q = (Q < admixtureThreshold).select(admixtureThreshold, Q); // lower bound
        Q = (Q > 1 - admixtureThreshold).select(1 - admixtureThreshold, Q); // upper bound
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
        if(debug && !(1.0 - Q.colwise().sum() == 0).any()) cao.warn("Q colsum is not 1.0");
    }
    if(!NP)
    { // protect P
        if(P.isNaN().any()) cao.warn("NaN in P in Phaseless model. will fill it with AF");
        P = (P < alleleEmitThreshold).select(alleleEmitThreshold, P); // lower bound
        P = (P > 1 - alleleEmitThreshold).select(1 - alleleEmitThreshold, P); // upper bound
    }
    if(!NF)
    { // protect F
        for(int k = 0; k < K; k++)
        {
            // could cluster jump be zero?
            if(F[k].isNaN().any()) cao.error("NaN in F in Phaseless model. reset it to the threshold. k =", k);
            F[k] = (F[k] < clusterFreqThreshold).select(clusterFreqThreshold, F[k]);
            F[k] = (F[k] > 1 - clusterFreqThreshold).select(1 - clusterFreqThreshold, F[k]);
            // re-normalize F per site. hope should work well. otherwise do the complicated.
            F[k].rowwise() /= F[k].colwise().sum();
        }
    }
    if(!NR)
    {
        for(int i = 0; i < er.size(); i++)
        {
            const double miner = std::exp(-nGen * maxRate * dist[i] / 100 / 1e6);
            const double maxer = std::exp(-nGen * minRate * dist[i] / 100 / 1e6);
            er(i) = er(i) < miner ? miner : er(i);
            er(i) = er(i) > maxer ? maxer : er(i);
        }
        protect_er(er);
        R = er2R(er);
    }
}

void Phaseless::initIteration()
{
    EclusterK.setZero(C * K, M);
    EclusterA1.setZero(C, M);
    EclusterA2.setZero(C, M);
    Eancestry.setZero(K, N);
}

void Phaseless::updateIteration()
{
    // update P
    if(!NP) P = (EclusterA2 / (EclusterA1 + EclusterA2)).transpose();
    if(!NQ)
    { // update Q
        Q = Eancestry;
        Q.rowwise() /= Q.colwise().sum(); // normalize Q per individual
    }
    if(!NF)
    { // update F
        for(int k = 0; k < K; k++)
        {
            F[k] = EclusterK.middleRows(k * C, C); // C x M
            F[k].rowwise() /= F[k].colwise().sum(); // normalize F per site per K
        }
    }
    if(!NR) er = 1.0 - EclusterK.colwise().sum() / N;
    protectPars();
}

void Phaseless::callGenoLoopC(int ind, int s, int z1, const MyArr2D & gli, const MyArr1D & gamma_div_emit)
{
    MyArr1D tmp_zg(4);
    for(int z2 = 0; z2 < C; z2++)
    {
        const int z12 = diploid_unordered_state_index(z1, z2, C);
        tmp_zg(0) = gli(s, 0) * (1 - P(s, z1)) * (1 - P(s, z2));
        tmp_zg(1) = gli(s, 1) * (1 - P(s, z1)) * P(s, z2);
        tmp_zg(2) = gli(s, 1) * P(s, z1) * (1 - P(s, z2));
        tmp_zg(3) = gli(s, 2) * P(s, z1) * P(s, z2);
        GP(3 * s + 0, ind) += gamma_div_emit(z12) * tmp_zg(0);
        GP(3 * s + 1, ind) += gamma_div_emit(z12) * (tmp_zg(1) + tmp_zg(2));
        GP(3 * s + 2, ind) += gamma_div_emit(z12) * tmp_zg(3);
    }
}

void Phaseless::getPosterios(const int ind,
                             const int ic,
                             const MyArr2D & gli,
                             const MyArr2D & emit,
                             const MyArr2D & H,
                             const MyArr1D & cs,
                             const MyArr2D & alpha,
                             const MyArr2D & beta,
                             bool finalIter)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    int m{0}, s{0}, z1{0}, z2{0}, y1{0}, zz{0};
    MyArr2D ind_post_zg1(C, S), ind_post_zg2(C, S);
    MyArr2D ind_post_zy(C * K, S);
    MyArr1D gamma_div_emit(diploid_unordered_state_count(C));
    ind_post_zy.setZero();
    for(s = 0; s < S; s++)
    {
        m = s + pos_chunk[ic];
        gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // what if emit is 0
        for(z1 = 0; z1 < C; z1++)
        {
            if(finalIter) callGenoLoopC(ind, m, z1, gli, gamma_div_emit);
            double post_zg1 = 0, post_zg2 = 0;
            for(z2 = 0; z2 < C; ++z2)
            {
                const double weight = gamma_div_emit(diploid_unordered_state_index(z1, z2, C));
                post_zg1 += weight * (1 - P(m, z1)) * (gli(s, 0) * (1 - P(m, z2)) + gli(s, 1) * P(m, z2));
                post_zg2 += weight * P(m, z1) * (gli(s, 1) * (1 - P(m, z2)) + gli(s, 2) * P(m, z2));
            }
            ind_post_zg1(z1, s) = post_zg1;
            ind_post_zg2(z1, s) = post_zg2;
            if(s == 0)
            {
                double tmp = 0;
                for(z2 = 0; z2 < C; ++z2)
                {
                    const int state = diploid_unordered_state_index(z1, z2, C);
                    tmp += alpha(state, 0) * beta(state, 0);
                }
                for(y1 = 0; y1 < K; y1++)
                {
                    // The first site is a compulsory cluster refresh.  Given
                    // Z=c, its ancestry responsibility is Q_k F_ck / H_c.
                    ind_post_zy(y1 * C + z1, 0) = tmp * Q(y1, ind) * F[y1](z1, m) / H(z1, 0);
                }
            }
        }
        if(s == 0) continue;
        MyArr1D alphaprev(C); // previous alpha colsums
        for(z1 = 0; z1 < C; z1++)
        {
            alphaprev(z1) = 0;
            for(z2 = 0; z2 < C; ++z2) alphaprev(z1) += alpha(diploid_unordered_state_index(z1, z2, C), s - 1);
        }
        for(z1 = 0; z1 < C; z1++)
        {
            double tmp{0};
            for(z2 = 0; z2 < C; z2++)
            {
                zz = diploid_unordered_state_index(z1, z2, C);
                double eb = emit(zz, s) * beta(zz, s);
                tmp += eb * (R(1, m) * alphaprev(z2) + R(2, m) * H(z2, s));
            }
            for(y1 = 0; y1 < K; y1++) ind_post_zy(y1 * C + z1, s) = tmp * Q(y1, ind) * F[y1](z1, m) * cs(s);
        }
    }
    // Q is the ancestry distribution at cluster-refresh events.  No-refresh
    // transitions contain no ancestry draw and therefore contribute no count.
    for(y1 = 0; y1 < K; y1++) Eancestry(y1, ind) += ind_post_zy.middleRows(y1 * C, C).sum();
    { // sum over all samples for updates
        std::scoped_lock<std::mutex> lock(mutex_it);
        EclusterA1.middleCols(pos_chunk[ic], S) += ind_post_zg1;
        EclusterA2.middleCols(pos_chunk[ic], S) += ind_post_zg2;
        EclusterK.middleCols(pos_chunk[ic], S) += ind_post_zy;
    }
}

double Phaseless::runForwardBackwards(const int ind, const int ic, const MyFloat1D & GL, bool finalIter)
{
    const int S = pos_chunk[ic + 1] - pos_chunk[ic];
    Eigen::Map<const MyArr2D> gli(GL.data() + ind * S * 3, S, 3);
    MyArr2D emit = get_emission_by_gl_symmetric(gli, P.middleRows(pos_chunk[ic], S));
    // first get H ie old PI in fastphase
    MyArr2D H = MyArr2D::Zero(C, S);
    int z1, y1, s;
    for(s = 0; s < S; s++)
        for(z1 = 0; z1 < C; z1++)
            for(y1 = 0; y1 < K; y1++) H(z1, s) += Q(y1, ind) * F[y1](z1, s + pos_chunk[ic]);
    // cs is 1 / colsum(alpha)
    const auto [alpha, beta, cs] =
        forward_backwards_diploid_symmetric(emit, R.middleCols(pos_chunk[ic], S), H);
    // get posterios
    getPosterios(ind, ic, gli, emit, H, cs, alpha, beta, finalIter);
    return (1 / cs).log().sum();
}

double Phaseless::runBigass(int ind, const MyFloat2D & GL, bool finalIter)
{
    if(pos_chunk.size() == 0) cao.error("please run initRecombination first");
    int nchunks = GL.size();
    double loglike{0};
    for(int ic = 0; ic < nchunks; ic++)
    {
        loglike += runForwardBackwards(ind, ic, GL[ic], finalIter);
    }
    return loglike;
}

int run_phaseless_main(Options & opts)
{
    cao.cao.open(opts.out + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    cao.warn(tim.date(), "-> running phaseless");
    const unsigned int allthreads = std::thread::hardware_concurrency();
    opts.nthreads = resolve_thread_count(opts.nthreads, allthreads);
    cao.print(tim.date(), allthreads, " concurrent threads are available. use", opts.nthreads, " threads");
    ThreadPool pool(opts.nthreads);
    if(opts.gpu)
    {
        std::string reason;
        if(!joint_cuda_available(reason)) throw std::runtime_error("cannot use --gpu: " + reason);
        cao.print(tim.date(), "using NVIDIA CUDA for the joint-model E step");
    }

    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>();
    VariantMetadata metadata;
    init_bigass(genome, opts, opts.oVCF ? &metadata : nullptr);
    vector<future<double>> res;
    Phaseless faith(opts.K, opts.C, genome->nsamples, genome->nsnps, opts.seed);
    faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.nP, opts.nF, opts.nR);
    faith.setStartPoint(opts.in_qfile, opts.in_pfile);
    faith.initRecombination(genome->pos, opts.in_rfile);
    auto evaluate_e_step = [&](bool final_iteration)
    {
        if(opts.gpu) return joint_cuda_e_step(faith, genome->gls, final_iteration);
        double value = 0;
        for(int i = 0; i < faith.N; i++)
            res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), final_iteration));
        for(auto && ll : res) value += ll.get();
        res.clear();
        return value;
    };
    constexpr double monotonicity_tol{1e-10};
    const double gap_tol = opts.conv_gap_tol;
    const double relative_tol = opts.conv_relative_tol;
    const double parameter_tol = opts.conv_parameter_tol;
    const int stable_iterations_required = opts.conv_stable_iterations;
    const double observations = std::max(1.0, static_cast<double>(faith.N) * faith.M);
    double loglike{NAN}, previous_like{NAN}, previous_previous_like{NAN};
    JointParameterSnapshot previous_parameters;
    bool have_previous_parameters{false};
    bool did_converge{false};
    int stable_iterations{0};
    cao.print(tim.date(), std::scientific, "joint convergence: gap/observation < ", gap_tol,
              ", relative likelihood change < ", relative_tol, ", parameter change < ", parameter_tol, " for ",
              stable_iterations_required, " accepted iterations");

    auto joint_converged = [&](int iteration, double current_like)
    {
        if(!have_previous_parameters)
        {
            previous_like = current_like;
            previous_parameters = snapshot_parameters(faith);
            have_previous_parameters = true;
            return false;
        }

        const auto likelihood = assess_likelihood_convergence(current_like, previous_like, previous_previous_like,
                                                               observations, monotonicity_tol);
        const auto parameters = parameter_change(faith, previous_parameters);
        const bool likelihood_stable = std::isfinite(previous_previous_like)
                                    && likelihood.gap_per_observation < gap_tol
                                    && likelihood.relative_change < relative_tol;
        const bool parameters_stable = parameters.q_max < parameter_tol && parameters.p_rms < parameter_tol
                                    && parameters.f_rms < parameter_tol && parameters.r_rms < parameter_tol;
        const bool stable = likelihood.monotone && likelihood_stable && parameters_stable;
        stable_iterations = stable ? stable_iterations + 1 : 0;

        cao.print(tim.date(), "accepted iteration", iteration, ", likelihood =", current_like, ", delta =",
                  likelihood.delta, ", relative =", std::scientific, likelihood.relative_change,
                  ", gap/observation =", likelihood.gap_per_observation, ", Aitken =", likelihood.aitken_rate,
                  ", dQ(max) =", parameters.q_max, ", dP(rms) =", parameters.p_rms, ", dF(rms) =",
                  parameters.f_rms, ", dR(rms) =", parameters.r_rms, ", stable =", stable_iterations, "/",
                  stable_iterations_required);
        if(!likelihood.monotone)
            cao.warn(tim.date(), "accepted joint-model likelihood decreased by", likelihood.delta,
                     "; convergence counter reset");

        previous_previous_like = previous_like;
        previous_like = current_like;
        previous_parameters = snapshot_parameters(faith);
        return stable_iterations >= stable_iterations_required;
    };

    if(opts.noaccel)
    {
        for(int it = 0; SIG_COND && it <= opts.nimpute; it++)
        {
            tim.clock();
            faith.initIteration();
            loglike = evaluate_e_step(false);
            cao.print(tim.date(), "run whole genome, iteration", it, ", likelihood =", loglike, ", time",
                      tim.reltime(), " sec");
            if(joint_converged(it, loglike))
            {
                did_converge = true;
                cao.print(tim.date(), "joint model converged after", stable_iterations,
                          " consecutive stable accepted iterations");
                break;
            }
            if(it == opts.nimpute) break;
            faith.updateIteration();
        }
    }
    else
    {
        MyArr2D Q0, Q1, Q2, Qt;
        MyArr2D F0, F1, F2, Ft;
        const int istep{4};
        double alpha{0}, stepMax{4}, alphaMax{1280}, logcheck{0};
        const int max_outer_iterations = opts.nimpute / 4;
        for(int it = 0; SIG_COND && it <= max_outer_iterations; it++)
        {
            // Evaluate the current accepted state, then take the first normal EM step.
            tim.clock();
            faith.initIteration();
            Q0 = faith.Q;
            F0 = cat_stdvec_of_eigen(faith.F);
            loglike = evaluate_e_step(false);
            if(joint_converged(it, loglike))
            {
                did_converge = true;
                cao.print(tim.date(), "joint model converged after", stable_iterations,
                          " consecutive stable accepted iterations");
                break;
            }
            if(it == max_outer_iterations) break;
            faith.updateIteration();
            // second normal iter
            faith.initIteration();
            Q1 = faith.Q;
            F1 = cat_stdvec_of_eigen(faith.F);
            loglike = evaluate_e_step(false);
            faith.updateIteration();
            cao.print(tim.date(), "SqS3 outer iteration", it, ", accepted likelihood =", previous_like,
                      ", second EM likelihood =", loglike, ", time", tim.reltime(), " sec");
            // save for later comparison
            Qt = faith.Q;
            Ft = cat_stdvec_of_eigen(faith.F);
            // calculate alpha based on first two pars
            if(opts.aQ)
            {

                alpha = ((Q1 - Q0).square().sum()) / ((faith.Q - 2 * Q1 + Q0).square().sum());
            }
            else
            {
                alpha = ((F1 - F0).square().sum() + (Q1 - Q0).square().sum())
                        / ((cat_stdvec_of_eigen(faith.F) - 2 * F1 + F0).square().sum()
                           + (faith.Q - 2 * Q1 + Q0).square().sum());
            }
            alpha = max(1.0, sqrt(alpha));
            if(alpha >= stepMax)
            {
                alpha = min(stepMax, alphaMax);
                stepMax = min(stepMax * istep, alphaMax);
            }
            // third accel iter
            // update Q and F using the second em iter
            faith.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (faith.Q - 2 * Q1 + Q0);
            for(int k = 0; k < faith.K; k++)
                faith.F[k] = F0.middleRows(k * faith.C, faith.C)
                             + 2 * alpha * (F1.middleRows(k * faith.C, faith.C) - F0.middleRows(k * faith.C, faith.C))
                             + alpha * alpha
                                   * (faith.F[k] - 2 * F1.middleRows(k * faith.C, faith.C)
                                      + F0.middleRows(k * faith.C, faith.C));
            faith.protectPars();
            faith.initIteration();
            loglike = evaluate_e_step(false);
            faith.updateIteration();
            // save current pars
            Q2 = faith.Q;
            F2 = cat_stdvec_of_eigen(faith.F);
            // check if normal third iter is better
            faith.Q = Qt;
            for(int k = 0; k < faith.K; k++) faith.F[k] = Ft.middleRows(k * faith.C, faith.C);
            faith.initIteration();
            logcheck = evaluate_e_step(false);
            faith.updateIteration();
            if(loglike < logcheck - monotonicity_tol * observations)
            {
                stepMax = istep;
                cao.warn(tim.date(), "reset stepMax to 4, normal EM yields better likelihoods than the accelerated EM.",
                         logcheck, " -", loglike, ">", monotonicity_tol * observations);
            }
            else
            {
                faith.Q = Q2;
                for(int k = 0; k < faith.K; k++) faith.F[k] = F2.middleRows(k * faith.C, faith.C);
            }
        }
    }
    if(!did_converge)
        cao.warn(tim.date(), "joint model reached the iteration limit before satisfying the convergence criterion");
    std::ofstream oanc(opts.out + ".Q");
    oanc << std::fixed << faith.Q.transpose().format(fmt10) << "\n";
    oanc.close();
    std::ofstream op(opts.out + ".P");
    op << faith.P.format(fmt6) << "\n";
    if(opts.oF)
    {
        std::ofstream of(opts.out + ".F");
        for(size_t k = 0; k < faith.F.size(); k++) of << faith.F[k].format(fmt6) << "\n";
    }
    std::unique_ptr<Pars> par = std::make_unique<Pars>();
    par->init(faith.K, faith.C, faith.M, faith.N, faith.er, faith.P, faith.Q, faith.F);
    par->pos = genome->pos;
    par->gls = genome->gls;
    std::ofstream opar(opts.out + ".pars.bin", std::ios::out | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    auto bytes_written = alpaca::serialize<OPTIONS, Pars>(*par, opar);
    opar.close();
    assert(std::filesystem::file_size(opts.out + ".pars.bin") == bytes_written);
    cao.done(tim.date(), "joint model done and outputting.", bytes_written, " bytes written to file");
    if(opts.oVCF)
    {
        faith.initIteration();
        cao.done(tim.date(), "run one more iteration to output vcf.");
        faith.GP.setZero(faith.M * 3, faith.N);
        loglike = evaluate_e_step(true);
        auto bw = make_bcfwriter(opts.out + ".vcf.gz", genome->chrs, genome->sampleids);
        for(int ic = 0; ic < genome->nchunks; ic++)
        {
            const int S = faith.pos_chunk[ic + 1] - faith.pos_chunk[ic];
            MyArr2D out = extract_gp_chunk(faith.GP, faith.pos_chunk[ic], S);
            write_bigass_to_bcf(bw, out.data(), genome->chrs[ic], genome->pos[ic], metadata.ids[ic],
                                metadata.refs[ic], metadata.alts[ic]);
        }
    }

    return 0;
}
