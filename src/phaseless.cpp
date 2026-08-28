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
#include <deque>

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

void restore_parameters(Phaseless & model, const JointParameterSnapshot & snapshot)
{
    model.Q = snapshot.Q;
    model.P = snapshot.P;
    for(int k = 0; k < model.K; ++k) model.F[k] = snapshot.F.middleRows(k * model.C, model.C);
    model.er = snapshot.er;
    model.R = er2R(model.er);
}

std::vector<int> maximum_weight_assignment(const MyArr2D & score)
{
    const int n = score.rows();
    if(score.cols() != n) throw std::invalid_argument("cluster assignment score must be square");
    std::vector<double> u(n + 1), v(n + 1);
    std::vector<int> p(n + 1), way(n + 1);
    for(int i = 1; i <= n; ++i)
    {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(n + 1, std::numeric_limits<double>::infinity());
        std::vector<bool> used(n + 1, false);
        do
        {
            used[j0] = true;
            const int i0 = p[j0];
            double delta = std::numeric_limits<double>::infinity();
            int j1 = 0;
            for(int j = 1; j <= n; ++j)
            {
                if(used[j]) continue;
                const double cur = -score(i0 - 1, j - 1) - u[i0] - v[j];
                if(cur < minv[j])
                {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if(minv[j] < delta)
                {
                    delta = minv[j];
                    j1 = j;
                }
            }
            for(int j = 0; j <= n; ++j)
            {
                if(used[j])
                {
                    u[p[j]] += delta;
                    v[j] -= delta;
                }
                else
                    minv[j] -= delta;
            }
            j0 = j1;
        } while(p[j0] != 0);
        do
        {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while(j0 != 0);
    }
    std::vector<int> assignment(n);
    for(int j = 1; j <= n; ++j) assignment[p[j] - 1] = j - 1;
    return assignment;
}

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

void Phaseless::setAdmixturePseudocount(double value)
{
    if(!std::isfinite(value) || value < 0)
        throw std::invalid_argument("Q pseudocount must be finite and non-negative");
    admixturePseudocount = value;
}

void Phaseless::setEmissionShrinkage(double value)
{
    if(!std::isfinite(value) || value < 0)
        throw std::invalid_argument("P shrinkage must be finite and non-negative");
    emissionShrinkage = value;
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

void Phaseless::initializeSharedHaplotypeStart()
{
    MyArr2D common = MyArr2D::Zero(C, M);
    for(const auto & frequencies : F) common += frequencies;
    common /= static_cast<double>(K);
    common.rowwise() /= common.colwise().sum();
    for(auto & frequencies : F) frequencies = common;
    Q.setConstant(1.0 / static_cast<double>(K));
}

bool Phaseless::initializeAncestryFromPosterior(double noise, int restart)
{
    if(noise < 0 || noise > 1) throw std::invalid_argument("ancestry initialization noise must be in [0, 1]");
    if(restart < 0) throw std::invalid_argument("ancestry initialization restart must be non-negative");
    if(K == 1)
    {
        Q.setOnes();
        return true;
    }
    auto jittered_symmetric_fallback = [&]()
    {
        const MyArr2D common = F.front();
        Q = RandomUniform<MyArr2D, std::default_random_engine>(K, N, rng, admixtureThreshold,
                                                                1 - admixtureThreshold);
        Q.rowwise() /= Q.colwise().sum();
        for(auto & frequencies : F)
        {
            const MyArr2D jitter = RandomUniform<MyArr2D, std::default_random_engine>(
                C, M, rng, 1 - noise, 1 + noise);
            frequencies = common * jitter;
            frequencies.rowwise() /= frequencies.colwise().sum();
        }
        protectPars();
        return false;
    };
    if(EindividualClusterUsage.rows() != C || EindividualClusterUsage.cols() != N || N == 0)
        return jittered_symmetric_fallback();

    MyArr2D profiles = EindividualClusterUsage;
    for(int i = 0; i < N; ++i)
    {
        const double total = profiles.col(i).sum();
        if(!std::isfinite(total) || total <= 0) return jittered_symmetric_fallback();
        profiles.col(i) /= total;
    }
    const MyArr1D mean_profile = profiles.rowwise().mean();
    const double total_variation = (profiles.colwise() - mean_profile).square().mean();
    if(!std::isfinite(total_variation) || total_variation <= std::numeric_limits<double>::epsilon())
        return jittered_symmetric_fallback();

    // Perturb only the clustering view for later restarts.  The posterior
    // profiles themselves remain unchanged when constructing Q and F.  This
    // avoids spending nominal restarts on the same deterministic partition.
    MyArr2D clustering_profiles = profiles;
    if(restart > 0 && noise > 0)
    {
        const double amplitude = std::min(0.25, noise * (1.0 + 0.25 * restart));
        clustering_profiles *= RandomUniform<MyArr2D, std::default_random_engine>(
            C, N, rng, 1 - amplitude, 1 + amplitude);
        clustering_profiles.rowwise() /= clustering_profiles.colwise().sum();
    }

    MyArr2D centroids(C, K);
    Int1D selected(K, 0);
    if(restart == 0)
    {
        MyArr1D distance(N);
        for(int i = 0; i < N; ++i) distance(i) = (profiles.col(i) - mean_profile).square().sum();
        Eigen::Index first = 0;
        distance.maxCoeff(&first);
        selected[0] = static_cast<int>(first);
    }
    else
    {
        MyArr1D weights(N);
        for(int i = 0; i < N; ++i)
            weights(i) = (clustering_profiles.col(i) - mean_profile).square().sum()
                       + std::numeric_limits<double>::epsilon();
        MyFloat1D sampling_weights(weights.data(), weights.data() + weights.size());
        std::discrete_distribution<int> choose_first(sampling_weights.begin(), sampling_weights.end());
        selected[0] = choose_first(rng);
    }
    centroids.col(0) = clustering_profiles.col(selected[0]);
    MyArr1D nearest = MyArr1D::Constant(N, std::numeric_limits<double>::infinity());
    for(int k = 1; k < K; ++k)
    {
        for(int i = 0; i < N; ++i)
            nearest(i) = std::min<double>(nearest(i),
                                          (clustering_profiles.col(i) - centroids.col(k - 1)).square().sum());
        int choice = k % N;
        if(restart == 0)
        {
            Eigen::Index farthest = 0;
            nearest.maxCoeff(&farthest);
            choice = static_cast<int>(farthest);
        }
        else if(nearest.sum() > std::numeric_limits<double>::epsilon())
        {
            MyFloat1D weights(nearest.data(), nearest.data() + nearest.size());
            std::discrete_distribution<int> choose(weights.begin(), weights.end());
            choice = choose(rng);
        }
        selected[k] = choice;
        centroids.col(k) = clustering_profiles.col(choice);
    }

    Int1D labels(N, -1);
    for(int iteration = 0; iteration < 50; ++iteration)
    {
        bool changed = false;
        for(int i = 0; i < N; ++i)
        {
            int best = 0;
            double best_distance = std::numeric_limits<double>::infinity();
            for(int k = 0; k < K; ++k)
            {
                const double distance = (clustering_profiles.col(i) - centroids.col(k)).square().sum();
                if(distance < best_distance)
                {
                    best_distance = distance;
                    best = k;
                }
            }
            changed = changed || labels[i] != best;
            labels[i] = best;
        }
        MyArr2D next = MyArr2D::Zero(C, K);
        Int1D counts(K, 0);
        for(int i = 0; i < N; ++i)
        {
            next.col(labels[i]) += clustering_profiles.col(i);
            ++counts[labels[i]];
        }
        for(int k = 0; k < K; ++k)
        {
            if(counts[k])
                next.col(k) /= counts[k];
            else
                next.col(k) = clustering_profiles.col(selected[k]);
        }
        centroids = std::move(next);
        if(!changed) break;
    }

    // Convert the randomized partition back to centroids on the unperturbed
    // posterior profiles before deriving model parameters.
    centroids.setZero();
    Int1D final_counts(K, 0);
    for(int i = 0; i < N; ++i)
    {
        centroids.col(labels[i]) += profiles.col(i);
        ++final_counts[labels[i]];
    }
    for(int k = 0; k < K; ++k)
        if(final_counts[k]) centroids.col(k) /= final_counts[k];
        else centroids.col(k) = profiles.col(selected[k]);

    double within = 0;
    for(int i = 0; i < N; ++i) within += (profiles.col(i) - centroids.col(labels[i])).square().sum();
    const double base_temperature = std::max({within / std::max(1, N),
                                              (0.02 + noise) * total_variation * C,
                                              100 * std::numeric_limits<double>::epsilon()});
    const int temperature_level = (restart + 1) / 2;
    const double temperature_scale = restart == 0 ? 1.0
                                   : restart % 2 ? std::pow(1.5, temperature_level)
                                                 : std::pow(1.5, -temperature_level);
    const double temperature = base_temperature * temperature_scale;
    for(int i = 0; i < N; ++i)
    {
        for(int k = 0; k < K; ++k)
        {
            const double distance = (profiles.col(i) - centroids.col(k)).square().sum();
            Q(k, i) = std::exp(-0.5 * distance / temperature);
        }
        Q.col(i) += admixtureThreshold;
        Q.col(i) /= Q.col(i).sum();
    }

    const MyArr2D common = F.front();
    for(int k = 0; k < K; ++k)
    {
        const double weight = Q.row(k).sum();
        MyArr1D ancestry_profile = MyArr1D::Zero(C);
        for(int i = 0; i < N; ++i) ancestry_profile += Q(k, i) * profiles.col(i);
        ancestry_profile /= std::max(weight, std::numeric_limits<double>::epsilon());
        MyArr1D ratio = (ancestry_profile + clusterFreqThreshold)
                      / (mean_profile + clusterFreqThreshold);
        ratio = ratio.sqrt();
        const MyArr2D jitter = RandomUniform<MyArr2D, std::default_random_engine>(
            C, M, rng, 1 - noise, 1 + noise);
        F[k] = common.colwise() * ratio;
        F[k] *= jitter;
        F[k].rowwise() /= F[k].colwise().sum();
    }
    protectPars();
    return true;
}

void Phaseless::configurePhaseAlignment(int boundary_stride)
{
    phaseAlignmentBoundaries.clear();
    phaseAlignmentIndex.assign(M, -1);
    EphaseAlignmentCross.clear();
    EphaseAlignmentSquares.resize(0, 0);
    if(boundary_stride <= 0) return;
    if(pos_chunk.empty()) throw std::logic_error("phase alignment requires initialized recombination chunks");
    for(size_t ic = 0; ic + 1 < pos_chunk.size(); ++ic)
        for(int boundary = pos_chunk[ic] + boundary_stride; boundary < pos_chunk[ic + 1]; boundary += boundary_stride)
        {
            phaseAlignmentIndex[boundary] = static_cast<int>(phaseAlignmentBoundaries.size());
            phaseAlignmentBoundaries.push_back(boundary);
            EphaseAlignmentCross.emplace_back(MyArr2D::Zero(C, C));
        }
    EphaseAlignmentSquares.setZero(2 * C, phaseAlignmentBoundaries.size());
}

JointHeuristicReport Phaseless::alignPhaseClusterLabels(int reset_radius)
{
    if(reset_radius < 0) throw std::invalid_argument("heuristic reset radius cannot be negative");
    JointHeuristicReport report;
    if(C < 2 || phaseAlignmentBoundaries.empty()) return report;
    Int1D current_to_original(C);
    std::iota(current_to_original.begin(), current_to_original.end(), 0);
    size_t current_chunk = 0;

    for(size_t b = 0; b < phaseAlignmentBoundaries.size(); ++b)
    {
        const int boundary = phaseAlignmentBoundaries[b];
        while(current_chunk + 1 < pos_chunk.size() && boundary >= pos_chunk[current_chunk + 1])
        {
            ++current_chunk;
            std::iota(current_to_original.begin(), current_to_original.end(), 0);
        }
        MyArr2D score = MyArr2D::Zero(C, C);
        for(int left = 0; left < C; ++left)
            for(int right = 0; right < C; ++right)
            {
                const int original_left = current_to_original[left];
                const int original_right = current_to_original[right];
                const double left_sum = EclusterUsage(original_left, boundary - 1);
                const double right_sum = EclusterUsage(original_right, boundary);
                const double covariance = EphaseAlignmentCross[b](original_left, original_right)
                                        - left_sum * right_sum / N;
                const double left_variance = EphaseAlignmentSquares(original_left, b)
                                           - left_sum * left_sum / N;
                const double right_variance = EphaseAlignmentSquares(C + original_right, b)
                                            - right_sum * right_sum / N;
                const double scale = std::sqrt(std::max(0.0, left_variance) * std::max(0.0, right_variance));
                if(scale > std::numeric_limits<double>::epsilon()) score(left, right) = covariance / scale;
            }
        const auto assignment = maximum_weight_assignment(score);
        double identity_score = 0, assigned_score = 0;
        for(int c = 0; c < C; ++c)
        {
            identity_score += score(c, c);
            assigned_score += score(c, assignment[c]);
        }
        const double improvement_tolerance = 1e-12 * std::max(1.0, std::abs(identity_score));
        if(assigned_score <= identity_score + improvement_tolerance) continue;

        const int chunk_end = pos_chunk[current_chunk + 1];
        const int tail_length = chunk_end - boundary;
        const MyArr2D old_p = P.middleRows(boundary, tail_length);
        for(int c = 0; c < C; ++c) P.middleRows(boundary, tail_length).col(c) = old_p.col(assignment[c]);
        for(int k = 0; k < K; ++k)
        {
            const MyArr2D old_f = F[k].middleCols(boundary, tail_length);
            for(int c = 0; c < C; ++c)
                F[k].middleCols(boundary, tail_length).row(c) = old_f.row(assignment[c]);
        }
        Int1D next_to_original(C);
        for(int c = 0; c < C; ++c) next_to_original[c] = current_to_original[assignment[c]];
        current_to_original = std::move(next_to_original);

        const int reset_start = std::max(pos_chunk[current_chunk], boundary - reset_radius);
        const int reset_end = std::min(chunk_end, boundary + reset_radius + 1);
        P.middleRows(reset_start, reset_end - reset_start) =
            RandomUniform<MyArr2D, std::default_random_engine>(reset_end - reset_start, C, rng,
                                                              alleleEmitThreshold, 1 - alleleEmitThreshold);
        report.reset_sites += reset_end - reset_start;
        ++report.relabelled_boundaries;
    }
    protectPars();
    return report;
}

void Phaseless::initIteration()
{
    EclusterK.setZero(C * K, M);
    EclusterA1.setZero(C, M);
    EclusterA2.setZero(C, M);
    Eancestry.setZero(K, N);
    EclusterUsage.setZero(C, M);
    EindividualClusterUsage.setZero(C, N);
    for(auto & cross : EphaseAlignmentCross) cross.setZero();
    if(EphaseAlignmentSquares.size()) EphaseAlignmentSquares.setZero();
    next_merge_ind.assign(pos_chunk.empty() ? 0 : pos_chunk.size() - 1, 0);
}

void Phaseless::updateIteration()
{
    // update P
    if(!NP)
    {
        if(emissionShrinkage == 0)
            P = (EclusterA2 / (EclusterA1 + EclusterA2)).transpose();
        else
        {
            // Treat the pooled site frequency as the prior mean and the
            // configured strength as an effective chromosome-copy count.
            // This stabilizes low-occupancy cluster/site cells without
            // pulling rare or common sites indiscriminately toward 0.5.
            for(int site = 0; site < M; ++site)
            {
                const double pooled_alt = EclusterA2.col(site).sum();
                const double pooled_total = pooled_alt + EclusterA1.col(site).sum();
                const double pooled_frequency = pooled_total > 0 ? pooled_alt / pooled_total : 0.5;
                for(int cluster = 0; cluster < C; ++cluster)
                {
                    const double cluster_total = EclusterA1(cluster, site) + EclusterA2(cluster, site);
                    P(site, cluster) = (EclusterA2(cluster, site)
                                      + emissionShrinkage * pooled_frequency)
                                     / (cluster_total + emissionShrinkage);
                }
            }
        }
    }
    if(!NQ)
    { // update Q
        Q = Eancestry + admixturePseudocount;
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

JointHeuristicReport Phaseless::alignClusterLabels(int boundary_stride, int reset_radius)
{
    if(boundary_stride < 1) throw std::invalid_argument("heuristic block size must be positive");
    if(reset_radius < 0) throw std::invalid_argument("heuristic reset radius cannot be negative");
    JointHeuristicReport report;
    if(C < 2 || M < 2) return report;
    const MyArr2D individual_ancestry_gram = (Q.matrix() * Q.matrix().transpose()).array();

    for(size_t ic = 0; ic + 1 < pos_chunk.size(); ++ic)
    {
        const int chunk_start = pos_chunk[ic];
        const int chunk_end = pos_chunk[ic + 1];
        for(int boundary = chunk_start + boundary_stride; boundary < chunk_end; boundary += boundary_stride)
        {
            MyArr2D left = MyArr2D::Zero(C, K);
            MyArr2D right = MyArr2D::Zero(C, K);
            for(int k = 0; k < K; ++k)
            {
                left.col(k) = F[k].col(boundary - 1);
                right.col(k) = F[k].col(boundary);
            }
            const MyArr2D score =
                (left.matrix() * individual_ancestry_gram.matrix() * right.matrix().transpose()).array();
            const auto assignment = maximum_weight_assignment(score);
            double identity_score = 0, assigned_score = 0;
            for(int c = 0; c < C; ++c)
            {
                identity_score += score(c, c);
                assigned_score += score(c, assignment[c]);
            }
            const double improvement_tolerance = 1e-12 * std::max(1.0, std::abs(identity_score));
            if(assigned_score <= identity_score + improvement_tolerance) continue;

            const int tail_length = chunk_end - boundary;
            const MyArr2D old_p = P.middleRows(boundary, tail_length);
            for(int c = 0; c < C; ++c) P.middleRows(boundary, tail_length).col(c) = old_p.col(assignment[c]);
            for(int k = 0; k < K; ++k)
            {
                const MyArr2D old_f = F[k].middleCols(boundary, tail_length);
                for(int c = 0; c < C; ++c)
                    F[k].middleCols(boundary, tail_length).row(c) = old_f.row(assignment[c]);
            }

            const int reset_start = std::max(chunk_start, boundary - reset_radius);
            const int reset_end = std::min(chunk_end, boundary + reset_radius + 1);
            P.middleRows(reset_start, reset_end - reset_start) =
                RandomUniform<MyArr2D, std::default_random_engine>(reset_end - reset_start, C, rng,
                                                                  alleleEmitThreshold, 1 - alleleEmitThreshold);
            report.reset_sites += reset_end - reset_start;
            ++report.relabelled_boundaries;
        }
    }
    protectPars();
    return report;
}

JointHeuristicReport Phaseless::reviveUnusedClusters(double min_usage, int bin_size, double donor_weight)
{
    if(min_usage < 0 || min_usage >= 1) throw std::invalid_argument("heuristic minimum usage must be in [0, 1)");
    if(bin_size < 1) throw std::invalid_argument("heuristic block size must be positive");
    if(donor_weight < 0 || donor_weight > 1)
        throw std::invalid_argument("heuristic donor weight must be in [0, 1]");
    JointHeuristicReport report;
    if(C < 2 || N < 1) return report;

    const double usage_scale = 2.0 * N;
    for(size_t ic = 0; ic + 1 < pos_chunk.size(); ++ic)
    {
        const int chunk_start = pos_chunk[ic];
        const int chunk_end = pos_chunk[ic + 1];
        const int bins = (chunk_end - chunk_start + bin_size - 1) / bin_size;
        for(int c = 0; c < C; ++c)
        {
            int b = 0;
            while(b < bins)
            {
                const int bin_start = chunk_start + b * bin_size;
                const int bin_end = std::min(chunk_end, bin_start + bin_size);
                const double average = EclusterUsage.row(c).segment(bin_start, bin_end - bin_start).mean()
                                             / usage_scale;
                if(average >= min_usage)
                {
                    ++b;
                    continue;
                }
                const int run_first_bin = b;
                do
                {
                    ++b;
                    if(b == bins) break;
                    const int next_start = chunk_start + b * bin_size;
                    const int next_end = std::min(chunk_end, next_start + bin_size);
                    const double next_average =
                        EclusterUsage.row(c).segment(next_start, next_end - next_start).mean() / usage_scale;
                    if(next_average >= min_usage) break;
                } while(true);
                const int run_start = chunk_start + run_first_bin * bin_size;
                const int run_end = std::min(chunk_end, chunk_start + b * bin_size);

                MyArr1D donor_usage = EclusterUsage.middleCols(run_start, run_end - run_start).rowwise().sum();
                donor_usage(c) = 0;
                int donor = 0;
                const double donor_total = donor_usage.sum();
                if(donor_total > 0)
                {
                    MyFloat1D weights(donor_usage.data(), donor_usage.data() + donor_usage.size());
                    std::discrete_distribution<int> distribution(weights.begin(), weights.end());
                    donor = distribution(rng);
                }
                else
                {
                    std::uniform_int_distribution<int> distribution(0, C - 2);
                    donor = distribution(rng);
                    if(donor >= c) ++donor;
                }
                const int length = run_end - run_start;
                const MyArr1D noise = RandomUniform<MyArr1D, std::default_random_engine>(
                    length, 1, rng, alleleEmitThreshold, 1 - alleleEmitThreshold);
                P.col(c).segment(run_start, length) =
                    donor_weight * P.col(donor).segment(run_start, length) + (1 - donor_weight) * noise;
                ++report.revived_intervals;
                report.revived_sites += length;
            }
        }
    }
    protectPars();
    return report;
}

InitializationProfilePruningReport Phaseless::configureInformativeProfileSites(int block_size,
                                                                                double information_fraction,
                                                                                int minimum_sites)
{
    if(block_size < 1) throw std::invalid_argument("initialization profile block size must be positive");
    if(!std::isfinite(information_fraction) || information_fraction <= 0 || information_fraction > 1)
        throw std::invalid_argument("initialization profile information fraction must be in (0, 1]");
    if(minimum_sites < 1 || minimum_sites > block_size)
        throw std::invalid_argument("initialization profile minimum sites must be in [1, block size]");
    if(P.rows() != M || P.cols() != C || EclusterUsage.rows() != C || EclusterUsage.cols() != M)
        throw std::logic_error("initialization profile pruning requires P and posterior cluster usage");

    InitializationProfilePruningReport report;
    MyArr1D information = MyArr1D::Zero(M);
    auto binary_entropy = [](double probability)
    {
        if(probability <= 0 || probability >= 1) return 0.0;
        return -probability * std::log(probability)
             - (1 - probability) * std::log1p(-probability);
    };
    for(int site = 0; site < M; ++site)
    {
        const double total_usage = EclusterUsage.col(site).sum();
        if(!std::isfinite(total_usage) || total_usage <= 0) continue;
        double pooled_frequency = 0;
        double conditional_entropy = 0;
        for(int cluster = 0; cluster < C; ++cluster)
        {
            const double weight = EclusterUsage(cluster, site) / total_usage;
            const double frequency = std::clamp<double>(P(site, cluster), 0, 1);
            pooled_frequency += weight * frequency;
            conditional_entropy += weight * binary_entropy(frequency);
        }
        const double score = binary_entropy(pooled_frequency) - conditional_entropy;
        if(std::isfinite(score) && score > 0) information(site) = score;
    }

    initializationProfileWeights = MyArr1D::Zero(M);
    std::vector<std::pair<int, int>> chunks;
    if(pos_chunk.size() >= 2)
    {
        for(size_t chunk = 0; chunk + 1 < pos_chunk.size(); ++chunk)
            chunks.emplace_back(pos_chunk[chunk], pos_chunk[chunk + 1]);
    }
    else
        chunks.emplace_back(0, M);

    constexpr double information_epsilon = 64 * std::numeric_limits<double>::epsilon();
    for(const auto & [chunk_start, chunk_end] : chunks)
        for(int block_start = chunk_start; block_start < chunk_end; block_start += block_size)
        {
            const int block_end = std::min(chunk_end, block_start + block_size);
            ++report.blocks;
            std::vector<int> ranked;
            double block_information = 0;
            for(int site = block_start; site < block_end; ++site)
                if(information(site) > information_epsilon)
                {
                    ranked.push_back(site);
                    block_information += information(site);
                }
            report.informative_sites += static_cast<int>(ranked.size());
            report.total_information += block_information;
            if(ranked.empty() || block_information <= 0) continue;
            ++report.informative_blocks;
            std::stable_sort(ranked.begin(), ranked.end(), [&](int left, int right)
            {
                return information(left) > information(right);
            });
            const double target = information_fraction * block_information;
            double retained = 0;
            int keep = 0;
            const int required = std::min<int>(minimum_sites, ranked.size());
            while(keep < static_cast<int>(ranked.size()) && (keep < required || retained < target))
            {
                retained += information(ranked[keep]);
                ++keep;
            }
            // Preserve the block's former total profile weight while using
            // only its informative representatives. Uninformative blocks
            // intentionally receive zero weight.
            const double site_weight = static_cast<double>(block_end - block_start) / keep;
            for(int rank = 0; rank < keep; ++rank) initializationProfileWeights(ranked[rank]) = site_weight;
            report.retained_sites += keep;
            report.retained_information += retained;
        }

    if(report.retained_sites == 0)
    {
        initializationProfileWeights.setOnes();
        report.retained_sites = M;
        report.used_all_sites_fallback = true;
    }
    return report;
}

void Phaseless::clearInformativeProfileSites()
{
    initializationProfileWeights.resize(0);
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
    MyArr2D ind_cluster_usage = MyArr2D::Zero(C, S);
    MyArr2D ind_post_zy(C * K, S);
    MyArr1D gamma_div_emit(diploid_unordered_state_count(C));
    ind_post_zy.setZero();
    for(s = 0; s < S; s++)
    {
        m = s + pos_chunk[ic];
        gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // what if emit is 0
        for(z1 = 0; z1 < C; ++z1)
            for(z2 = z1; z2 < C; ++z2)
            {
                const double posterior = 2 * alpha(diploid_unordered_state_index(z1, z2, C), s) * beta(
                                               diploid_unordered_state_index(z1, z2, C), s);
                ind_cluster_usage(z1, s) += posterior;
                if(z2 != z1) ind_cluster_usage(z2, s) += posterior;
            }
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
    if(initializationProfileWeights.size() == M)
        EindividualClusterUsage.col(ind) +=
            (ind_cluster_usage.matrix()
             * initializationProfileWeights.segment(pos_chunk[ic], S).matrix()).array();
    else
        EindividualClusterUsage.col(ind) += ind_cluster_usage.rowwise().sum();
    { // Sum individuals in a fixed order so seeded heuristic runs are reproducible across thread schedules.
        std::unique_lock<std::mutex> lock(mutex_it);
        merge_cv.wait(lock, [&] { return ind == next_merge_ind[ic]; });
        EclusterA1.middleCols(pos_chunk[ic], S) += ind_post_zg1;
        EclusterA2.middleCols(pos_chunk[ic], S) += ind_post_zg2;
        EclusterK.middleCols(pos_chunk[ic], S) += ind_post_zy;
        EclusterUsage.middleCols(pos_chunk[ic], S) += ind_cluster_usage;
        if(!EphaseAlignmentCross.empty())
            for(int local_boundary = 1; local_boundary < S; ++local_boundary)
            {
                const int boundary = pos_chunk[ic] + local_boundary;
                const int boundary_index = phaseAlignmentIndex[boundary];
                if(boundary_index < 0) continue;
                const MyArr1D left = ind_cluster_usage.col(local_boundary - 1);
                const MyArr1D right = ind_cluster_usage.col(local_boundary);
                EphaseAlignmentCross[boundary_index] += (left.matrix() * right.matrix().transpose()).array();
                EphaseAlignmentSquares.col(boundary_index).head(C) += left.square();
                EphaseAlignmentSquares.col(boundary_index).tail(C) += right.square();
            }
        ++next_merge_ind[ic];
        lock.unlock();
        merge_cv.notify_all();
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
    faith.setAdmixturePseudocount(opts.q_pseudocount);
    faith.setEmissionShrinkage(opts.p_shrinkage);
    cao.print(tim.date(), "Q update pseudocount =", opts.q_pseudocount,
              ", P site-frequency shrinkage =", opts.p_shrinkage);
    faith.setStartPoint(opts.in_qfile, opts.in_pfile);
    faith.initRecombination(genome->pos, opts.in_rfile);
    bool initialization_cpu_e_step{false};
    auto evaluate_e_step = [&](bool final_iteration)
    {
        if(opts.gpu && !initialization_cpu_e_step) return joint_cuda_e_step(faith, genome->gls, final_iteration);
        double value = 0;
        for(int i = 0; i < faith.N; i++)
            res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(genome->gls), final_iteration));
        for(auto && ll : res) value += ll.get();
        res.clear();
        return value;
    };
    auto run_initialization_stage = [&](const char * name, int scans)
    {
        double stage_like = NAN;
        for(int it = 0; SIG_COND && it < scans; ++it)
        {
            tim.clock();
            faith.initIteration();
            stage_like = evaluate_e_step(false);
            faith.updateIteration();
            cao.print(tim.date(), name, "scan", it + 1, "/", scans, ", likelihood =", stage_like,
                      ", time", tim.reltime(), "sec");
        }
        return stage_like;
    };
    const bool posterior_initialization = !opts.random_init && opts.in_qfile.empty();
    JointParameterSnapshot initialization_checkpoint;
    double initialization_checkpoint_like{-std::numeric_limits<double>::infinity()};
    bool have_initialization_checkpoint{false};
    if(!posterior_initialization && !opts.in_qfile.empty())
        cao.print(tim.date(), "using explicit Q start; posterior-driven initialization is disabled");
    if(posterior_initialization)
    {
        cao.print(tim.date(), "posterior-driven initialization: shared haplotype scans =",
                  opts.init_haplotype_min_iterations, "-", opts.init_haplotype_iterations,
                  " (adaptive), ancestry refinement scans =", opts.init_ancestry_iterations,
                  ", ancestry starts =", opts.init_restarts, ", initialization noise =", opts.init_noise,
                  ", Q pseudocount =", opts.q_pseudocount, ", STITCH heuristics =",
                  opts.stitch_heuristics ? "integrated" : "disabled", ", no block warm-up");

        // Learn the haplotype map without asking a random ancestry split to
        // organize the cluster labels at the same time.
        faith.initializeSharedHaplotypeStart();
        faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, true, opts.nP, false, opts.nR);
        initialization_cpu_e_step = true;
        if(opts.gpu)
            cao.warn(tim.date(), "posterior-driven initialization uses CPU E-steps; GPU resumes afterward");
        if(opts.stitch_heuristics)
            faith.configurePhaseAlignment(opts.heuristic_block_size);
        // Heuristics belong exclusively to this shared-haplotype fit. End
        // them once the unperturbed profile first looks stable, or early
        // enough to guarantee a clean tail before the scan limit.
        constexpr int initialization_heuristic_cooldown{8};
        const int last_initialization_heuristic_scan =
            std::max(0, opts.init_haplotype_iterations - initialization_heuristic_cooldown);
        bool initialization_heuristics_active = opts.stitch_heuristics;
        int initialization_heuristics_end_scan{0};
        double previous_stage_like{NAN};
        MyArr2D previous_profile;
        bool have_previous_profile{false};
        int initialization_stable_scans{0};
        int completed_haplotype_scans{0};
        bool shared_haplotype_converged{false};
        for(int it = 0; SIG_COND && it < opts.init_haplotype_iterations; ++it)
        {
            tim.clock();
            faith.initIteration();
            const double stage_like = evaluate_e_step(false);
            MyArr2D current_profile = faith.EindividualClusterUsage;
            for(int individual = 0; individual < current_profile.cols(); ++individual)
            {
                const double total = current_profile.col(individual).sum();
                if(std::isfinite(total) && total > 0) current_profile.col(individual) /= total;
                else current_profile.col(individual).setZero();
            }
            faith.updateIteration();
            JointHeuristicReport heuristic_report;
            const int scan = it + 1;
            if(initialization_heuristics_active && scan <= last_initialization_heuristic_scan
               && it >= 4 && it % 4 == 0)
            {
                const auto report = faith.alignPhaseClusterLabels(opts.heuristic_reset_radius);
                heuristic_report.relabelled_boundaries += report.relabelled_boundaries;
                heuristic_report.reset_sites += report.reset_sites;
                cao.warn(tim.date(), "initialization posterior label alignment", scan, ": relabelled",
                         report.relabelled_boundaries, "boundaries and reset", report.reset_sites, "SNPs");
            }
            if(initialization_heuristics_active && scan <= last_initialization_heuristic_scan
               && it >= 6 && (it - 2) % 4 == 0)
            {
                const auto report = faith.reviveUnusedClusters(opts.heuristic_min_usage,
                                                                opts.heuristic_block_size,
                                                                opts.heuristic_donor_weight);
                heuristic_report.revived_intervals += report.revived_intervals;
                heuristic_report.revived_sites += report.revived_sites;
                cao.warn(tim.date(), "initialization cluster revival", scan, ": revived", report.revived_intervals,
                         "intervals covering", report.revived_sites, "SNPs");
            }
            const double relative_change = std::isfinite(previous_stage_like)
                ? std::abs(stage_like - previous_stage_like) / std::max(1.0, std::abs(stage_like))
                : NAN;
            const double profile_rms = have_previous_profile
                ? std::sqrt((current_profile - previous_profile).square().mean())
                : NAN;
            const bool profile_stable = !heuristic_report.changed()
                                     && initialization_converged(relative_change, profile_rms,
                                                                 opts.init_haplotype_relative_tol,
                                                                 opts.init_haplotype_profile_tol);
            if(initialization_heuristics_active && scan >= opts.init_haplotype_min_iterations
               && profile_stable)
            {
                initialization_heuristics_active = false;
                initialization_heuristics_end_scan = scan;
                initialization_stable_scans = 0;
                cao.print(tim.date(), "shared-haplotype profile stabilized; ending STITCH heuristics at scan",
                          scan, "and starting clean posterior-profile cooldown");
            }
            else if(initialization_heuristics_active && scan >= last_initialization_heuristic_scan)
            {
                initialization_heuristics_active = false;
                initialization_heuristics_end_scan = scan;
                initialization_stable_scans = 0;
                cao.print(tim.date(), "ending STITCH heuristics at scan", scan,
                          "to reserve the posterior-profile cooldown");
            }
            const bool eligible = scan >= opts.init_haplotype_min_iterations
                               && (!opts.stitch_heuristics
                                   || (!initialization_heuristics_active
                                       && scan > initialization_heuristics_end_scan));
            const bool stable = eligible && profile_stable;
            initialization_stable_scans = stable
                ? std::min(opts.init_haplotype_stable_iterations, initialization_stable_scans + 1)
                : 0;
            const bool heuristic_cooldown_complete = !opts.stitch_heuristics
                || (!initialization_heuristics_active
                    && scan >= initialization_heuristics_end_scan + initialization_heuristic_cooldown);
            completed_haplotype_scans = scan;
            cao.print(tim.date(), "shared-haplotype initialization scan", scan, "/",
                      opts.init_haplotype_iterations, ", likelihood =", stage_like, ", relative =",
                      relative_change, ", profile RMS =", profile_rms, ", stable =",
                      initialization_stable_scans, "/", opts.init_haplotype_stable_iterations,
                      ", time", tim.reltime(), "sec");
            previous_stage_like = stage_like;
            previous_profile = std::move(current_profile);
            have_previous_profile = true;
            if(initialization_stable_scans >= opts.init_haplotype_stable_iterations
               && heuristic_cooldown_complete)
            {
                shared_haplotype_converged = true;
                cao.print(tim.date(), "shared-haplotype initialization converged after", scan, "scans");
                break;
            }
        }
        faith.configurePhaseAlignment(0);
        if(completed_haplotype_scans == opts.init_haplotype_iterations && !shared_haplotype_converged)
            cao.warn(tim.date(), "shared-haplotype initialization reached its scan limit before adaptive convergence");
        if(opts.init_profile_pruning)
        {
            const auto report = faith.configureInformativeProfileSites(opts.init_profile_block_size,
                                                                        opts.init_profile_information_fraction,
                                                                        opts.init_profile_min_snps);
            if(report.used_all_sites_fallback)
                cao.warn(tim.date(), "shared-haplotype SNP information scores were all zero;",
                         "using every SNP in the ancestry initialization profile");
            else
                cao.print(tim.date(), "initialization-profile SNP pruning retained", report.retained_sites,
                          "/", faith.M, "sites across", report.informative_blocks, "/", report.blocks,
                          "informative blocks and", 100 * report.retained_information / report.total_information,
                          "% of cluster-allele information");
        }
        else
            faith.clearInformativeProfileSites();
        // Re-evaluate once after the last M-step. This makes the profiles used
        // to seed ancestry correspond exactly to the saved shared state. The
        // optional low-cost site weights affect only profile aggregation; all
        // SNPs remain in this HMM E-step and in final genotype imputation.
        faith.initIteration();
        const double final_shared_like = evaluate_e_step(false);
        const JointParameterSnapshot shared_parameters = snapshot_parameters(faith);
        const MyArr2D posterior_profiles = faith.EindividualClusterUsage;
        faith.clearInformativeProfileSites();
        cao.print(tim.date(), "shared-haplotype final profile likelihood =", final_shared_like);
        initialization_cpu_e_step = false;

        // Cluster genome-wide posterior haplotype occupancy, use the resulting
        // soft groups for Q, and tilt the shared F by each group's profile.
        faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, true, opts.nF, true);
        JointParameterSnapshot best_initialization;
        double best_initialization_like = -std::numeric_limits<double>::infinity();
        for(int restart = 0; SIG_COND && restart < opts.init_restarts; ++restart)
        {
            restore_parameters(faith, shared_parameters);
            faith.EindividualClusterUsage = posterior_profiles;
            const bool informative = faith.initializeAncestryFromPosterior(opts.init_noise, restart);
            if(!informative)
                cao.warn(tim.date(), "posterior cluster profiles are uninformative; using a jittered symmetric start");
            cao.print(tim.date(), "posterior-driven ancestry start", restart + 1, "/", opts.init_restarts);
            run_initialization_stage("fixed-haplotype ancestry refinement", opts.init_ancestry_iterations);
            faith.initIteration();
            const double candidate_like = evaluate_e_step(false);
            cao.print(tim.date(), "fixed-haplotype ancestry candidate", restart + 1,
                      ", observed likelihood =", candidate_like);
            if(candidate_like > best_initialization_like)
            {
                best_initialization_like = candidate_like;
                best_initialization = snapshot_parameters(faith);
            }
        }
        if(std::isfinite(best_initialization_like))
        {
            restore_parameters(faith, best_initialization);
            initialization_checkpoint = best_initialization;
            initialization_checkpoint_like = best_initialization_like;
            have_initialization_checkpoint = true;
        }
        else
            restore_parameters(faith, shared_parameters);

        faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.nP, opts.nF, opts.nR);
        cao.print(tim.date(), "posterior-driven initialization complete; selected candidate likelihood =",
                  best_initialization_like, ", releasing all requested joint parameter blocks");
    }
    constexpr double monotonicity_tol{1e-10};
    const double improvement_tol = opts.conv_gap_tol;
    const double relative_tol = opts.conv_relative_tol;
    const double parameter_tol = opts.conv_parameter_tol;
    const int stable_iterations_required = opts.conv_stable_iterations;
    const double observations = std::max(1.0, static_cast<double>(faith.N) * faith.M);
    double loglike{NAN}, previous_like{NAN}, previous_previous_like{NAN};
    JointParameterSnapshot previous_parameters;
    JointParameterSnapshot best_parameters;
    bool have_previous_parameters{false};
    bool have_best_parameters{false};
    bool did_converge{false};
    int stable_iterations{0};
    double best_like{-std::numeric_limits<double>::infinity()};
    if(have_initialization_checkpoint)
    {
        best_parameters = initialization_checkpoint;
        best_like = initialization_checkpoint_like;
        have_best_parameters = true;
    }
    cao.print(tim.date(), std::scientific, "joint convergence: improvement/observation < ", improvement_tol,
              ", relative likelihood change < ", relative_tol, ", parameter change < ", parameter_tol, " for ",
              stable_iterations_required, " accepted iterations");

    auto reset_convergence_history = [&]()
    {
        loglike = NAN;
        previous_like = NAN;
        previous_previous_like = NAN;
        have_previous_parameters = false;
        stable_iterations = 0;
    };

    auto preserve_best = [&](double current_like)
    {
        if(!std::isfinite(current_like) || current_like <= best_like) return;
        best_like = current_like;
        best_parameters = snapshot_parameters(faith);
        have_best_parameters = true;
    };

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
                                    && likelihood.improvement_per_observation < improvement_tol
                                    && likelihood.relative_change < relative_tol;
        const bool parameters_stable = parameters.q_max < parameter_tol && parameters.p_rms < parameter_tol
                                    && parameters.f_rms < parameter_tol && parameters.r_rms < parameter_tol;
        const bool stable = likelihood.monotone && likelihood_stable && parameters_stable;
        stable_iterations = stable ? stable_iterations + 1 : 0;

        cao.print(tim.date(), "accepted iteration", iteration, ", likelihood =", current_like, ", delta =",
                  likelihood.delta, ", relative =", std::scientific, likelihood.relative_change,
                  ", improvement/observation =", likelihood.improvement_per_observation,
                  ", Aitken rate =", likelihood.aitken_rate, ", Aitken gap/observation =",
                  likelihood.aitken_gap_per_observation,
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

    if(opts.stitch_heuristics)
    {
        if(posterior_initialization)
            cao.print(tim.date(), "STITCH heuristics completed inside shared-haplotype initialization;",
                      "joint optimization will not perturb clusters a second time");
        else
            cao.warn(tim.date(), "--stitch-heuristics requires posterior-driven initialization and was not applied");
    }

    if(opts.noaccel)
    {
        for(int it = 0; SIG_COND && it <= opts.nimpute; it++)
        {
            tim.clock();
            faith.initIteration();
            loglike = evaluate_e_step(false);
            preserve_best(loglike);
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
    const bool run_acceleration = !opts.noaccel;
    int ordinary_finish_start{-1};
    if(run_acceleration && SIG_COND && !did_converge)
    {
        const int istep{4};
        double alpha{0}, stepMax{4}, alphaMax{1280};
        const int acceleration_scans = opts.nimpute;
        const int max_outer_iterations = acceleration_scans / 4;
        const int iteration_offset = 0;
        const int finish_reserve = std::min(acceleration_scans,
                                            std::max(12, 4 * (stable_iterations_required + 2)));
        const int forced_finish_scan = opts.nimpute - finish_reserve;
        double accelerated_previous_like{NAN};
        JointParameterSnapshot accelerated_previous_parameters;
        bool have_accelerated_previous{false};
        std::deque<bool> recent_rejections;
        std::deque<double> recent_realized_gains;
        for(int it = 0; SIG_COND && it <= max_outer_iterations; it++)
        {
            // Evaluate the current accepted state, then take the first normal EM step.
            tim.clock();
            faith.initIteration();
            const JointParameterSnapshot sqs3_x0 = snapshot_parameters(faith);
            loglike = evaluate_e_step(false);
            preserve_best(loglike);
            const int scan = iteration_offset + 4 * it;
            double accelerated_relative{NAN};
            JointParameterChange accelerated_change;
            if(have_accelerated_previous)
            {
                accelerated_relative = std::abs(loglike - accelerated_previous_like)
                                     / std::max(1.0, std::abs(loglike));
                accelerated_change = parameter_change(faith, accelerated_previous_parameters);
                cao.print(tim.date(), "accelerated checkpoint", scan, ", likelihood =", loglike,
                          ", relative =", accelerated_relative, ", dQ(max) =", accelerated_change.q_max,
                          ", dP(rms) =", accelerated_change.p_rms, ", dF(rms) =",
                          accelerated_change.f_rms, ", dR(rms) =", accelerated_change.r_rms);
            }
            accelerated_previous_like = loglike;
            accelerated_previous_parameters = snapshot_parameters(faith);
            have_accelerated_previous = true;
            const double accepted_like = loglike;

            const int rejected = static_cast<int>(std::count(recent_rejections.begin(),
                                                              recent_rejections.end(), true));
            const double realized_gain_sum = std::accumulate(recent_realized_gains.begin(),
                                                              recent_realized_gains.end(), 0.0);
            const auto handoff = assess_sqs3_handoff(static_cast<int>(recent_rejections.size()), rejected,
                                                     static_cast<int>(recent_realized_gains.size()),
                                                     realized_gain_sum, accelerated_relative, relative_tol);
            const bool reserve_plain_em = scan >= forced_finish_scan;
            if(handoff.persistent_rejection || handoff.low_efficiency || reserve_plain_em
               || it == max_outer_iterations)
            {
                ordinary_finish_start = scan;
                const char * reason = handoff.persistent_rejection ? "persistent rejected SqS3 proposals"
                                     : handoff.low_efficiency ? "low realized SqS3 gain"
                                     : "reserved ordinary-EM convergence budget";
                cao.print(tim.date(), "switching from SqS3 to ordinary-EM finishing at scan", scan,
                          " because of ", reason, "; rejection rate =", handoff.rejection_rate,
                          ", mean realized gain =", handoff.mean_realized_gain, ", ",
                          opts.nimpute - scan, " scans remain");
                break;
            }
            faith.updateIteration();
            // second normal iter
            faith.initIteration();
            const JointParameterSnapshot sqs3_x1 = snapshot_parameters(faith);
            loglike = evaluate_e_step(false);
            faith.updateIteration();
            cao.print(tim.date(), "SqS3 outer iteration", it, ", accepted likelihood =", accepted_like,
                      ", second EM likelihood =", loglike, ", time", tim.reltime(), " sec");
            const JointParameterSnapshot normal_candidate = snapshot_parameters(faith);
            // Decide when to finish from the underlying ordinary EM map, not
            // from the deliberately enlarged distance between extrapolated
            // checkpoints.  The likelihood change is x0 -> x1 and the
            // parameter change is x1 -> x2, so both measure one ordinary map.
            const double ordinary_relative = std::abs(loglike - accepted_like)
                                           / std::max(1.0, std::abs(loglike));
            const JointParameterChange ordinary_change = parameter_change(faith, sqs3_x1);
            const bool ordinary_near_convergence = ordinary_em_near_convergence(
                loglike >= accepted_like - monotonicity_tol * observations,
                ordinary_relative, ordinary_change.q_max, ordinary_change.p_rms,
                ordinary_change.f_rms, ordinary_change.r_rms, relative_tol, parameter_tol);
            cao.print(tim.date(), "ordinary-EM map checkpoint", scan + 2, ", relative =", ordinary_relative,
                      ", dQ(max) =", ordinary_change.q_max, ", dP(rms) =", ordinary_change.p_rms,
                      ", dF(rms) =", ordinary_change.f_rms, ", dR(rms) =", ordinary_change.r_rms);
            if(ordinary_near_convergence)
            {
                ordinary_finish_start = scan + 2;
                cao.print(tim.date(), "switching from SqS3 to ordinary-EM finishing at scan",
                          ordinary_finish_start, " because the ordinary EM map reached the near-convergence",
                          " threshold; rejection rate =", handoff.rejection_rate,
                          ", mean realized gain =", handoff.mean_realized_gain, ", ",
                          opts.nimpute - ordinary_finish_start, " scans remain");
                break;
            }
            SqS3StepMoments step_moments;
            if(opts.aQ)
                add_sqs3_block_moments(step_moments, sqs3_x0.Q, sqs3_x1.Q, normal_candidate.Q);
            else
            {
                // Estimate one common step from every active parameter block.
                // Block means give Q, P, F, and r equal geometric weight even
                // though their arrays have very different sizes.
                if(!faith.NQ)
                    add_sqs3_block_moments(step_moments, sqs3_x0.Q, sqs3_x1.Q, normal_candidate.Q);
                if(!faith.NP)
                    add_sqs3_block_moments(step_moments, sqs3_x0.P, sqs3_x1.P, normal_candidate.P);
                if(!faith.NF)
                    add_sqs3_block_moments(step_moments, sqs3_x0.F, sqs3_x1.F, normal_candidate.F);
                if(!faith.NR)
                    add_sqs3_block_moments(step_moments, sqs3_x0.er, sqs3_x1.er, normal_candidate.er);
            }
            alpha = sqs3_step_length(step_moments);
            if(alpha >= stepMax)
            {
                alpha = min(stepMax, alphaMax);
                stepMax = min(stepMax * istep, alphaMax);
            }
            // Extrapolate the same complete active state used to estimate the
            // step length. Frozen blocks remain at their second-EM values.
            if(!faith.NQ)
                faith.Q = sqs3_x0.Q + 2 * alpha * (sqs3_x1.Q - sqs3_x0.Q)
                        + alpha * alpha * (normal_candidate.Q - 2 * sqs3_x1.Q + sqs3_x0.Q);
            if(!faith.NP)
                faith.P = sqs3_x0.P + 2 * alpha * (sqs3_x1.P - sqs3_x0.P)
                        + alpha * alpha * (normal_candidate.P - 2 * sqs3_x1.P + sqs3_x0.P);
            if(!faith.NF)
            {
                const MyArr2D extrapolated_F = sqs3_x0.F + 2 * alpha * (sqs3_x1.F - sqs3_x0.F)
                                             + alpha * alpha
                                                   * (normal_candidate.F - 2 * sqs3_x1.F + sqs3_x0.F);
                for(int k = 0; k < faith.K; ++k)
                    faith.F[k] = extrapolated_F.middleRows(k * faith.C, faith.C);
            }
            if(!faith.NR)
                faith.er = sqs3_x0.er + 2 * alpha * (sqs3_x1.er - sqs3_x0.er)
                         + alpha * alpha * (normal_candidate.er - 2 * sqs3_x1.er + sqs3_x0.er);
            faith.protectPars();
            const JointParameterSnapshot accelerated_candidate = snapshot_parameters(faith);
            faith.initIteration();
            const double accelerated_like = evaluate_e_step(false);

            // Compare complete, evaluated parameter states. Rejected SqS3
            // proposals cannot leak their P/r trajectory into the EM fallback.
            restore_parameters(faith, normal_candidate);
            faith.initIteration();
            const double normal_like = evaluate_e_step(false);
            const bool rejected_acceleration = !std::isfinite(accelerated_like)
                                            || accelerated_like < normal_like
                                                                     - monotonicity_tol * observations;
            if(rejected_acceleration)
            {
                stepMax = istep;
                cao.warn(tim.date(), "reset stepMax to 4, normal EM yields better likelihoods than the accelerated EM.",
                         normal_like, " -", accelerated_like, ">", monotonicity_tol * observations);
                loglike = normal_like;
            }
            else
            {
                restore_parameters(faith, accelerated_candidate);
                loglike = accelerated_like;
            }
            preserve_best(loglike);
            recent_rejections.push_back(rejected_acceleration);
            if(recent_rejections.size() > 24) recent_rejections.pop_front();
            const double ordinary_progress = std::max(normal_like - accepted_like,
                                                       monotonicity_tol * observations);
            const double realized_gain = rejected_acceleration
                ? 0.0
                : std::min(10.0, std::max(0.0, accelerated_like - normal_like) / ordinary_progress);
            recent_realized_gains.push_back(realized_gain);
            if(recent_realized_gains.size() > 16) recent_realized_gains.pop_front();
            cao.print(tim.date(), "SqS3 proposal alpha =", alpha, ", accelerated likelihood =",
                      accelerated_like, ", normal likelihood =", normal_like, ", selected =",
                      rejected_acceleration ? "ordinary EM" : "SqS3", ", realized gain =", realized_gain);
        }
    }
    if(ordinary_finish_start >= 0 && SIG_COND && !did_converge)
    {
        reset_convergence_history();
        for(int scan = ordinary_finish_start; SIG_COND && scan <= opts.nimpute; ++scan)
        {
            tim.clock();
            faith.initIteration();
            loglike = evaluate_e_step(false);
            preserve_best(loglike);
            cao.print(tim.date(), "ordinary-EM finishing scan", scan, ", likelihood =", loglike,
                      ", time", tim.reltime(), "sec");
            if(joint_converged(scan, loglike))
            {
                did_converge = true;
                cao.print(tim.date(), "joint model converged during ordinary-EM finishing after",
                          stable_iterations, " consecutive stable iterations");
                break;
            }
            if(scan == opts.nimpute) break;
            faith.updateIteration();
        }
    }
    if(have_best_parameters
       && (!std::isfinite(loglike) || best_like > loglike + monotonicity_tol * observations))
    {
        cao.warn(tim.date(), "restoring best observed-likelihood checkpoint:", best_like, "instead of", loglike);
        restore_parameters(faith, best_parameters);
        loglike = best_like;
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
