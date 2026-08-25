#include "joint_cuda.hpp"
#include "phaseless.hpp"
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
#ifdef USE_FLOAT
using Scalar = float;
#else
using Scalar = double;
#endif

void cuda_check(cudaError_t status, const char * operation)
{
    if(status == cudaSuccess) return;
    std::ostringstream message;
    message << operation << ": " << cudaGetErrorString(status);
    throw std::runtime_error(message.str());
}

template<typename T>
struct DeviceBuffer
{
    T * ptr{nullptr};
    size_t capacity{0};
    DeviceBuffer() = default;
    explicit DeviceBuffer(size_t count)
    {
        ensure(count);
    }
    ~DeviceBuffer()
    {
        if(ptr) cudaFree(ptr);
    }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer & operator=(const DeviceBuffer &) = delete;
    DeviceBuffer(DeviceBuffer && other) noexcept : ptr(other.ptr), capacity(other.capacity)
    {
        other.ptr = nullptr;
        other.capacity = 0;
    }
    DeviceBuffer & operator=(DeviceBuffer && other) noexcept
    {
        if(this == &other) return *this;
        if(ptr) cudaFree(ptr);
        ptr = other.ptr;
        capacity = other.capacity;
        other.ptr = nullptr;
        other.capacity = 0;
        return *this;
    }
    void ensure(size_t count)
    {
        if(count <= capacity) return;
        if(ptr) cuda_check(cudaFree(ptr), "cudaFree while growing persistent buffer");
        ptr = nullptr;
        capacity = 0;
        if(count) cuda_check(cudaMalloc(reinterpret_cast<void **>(&ptr), count * sizeof(T)), "cudaMalloc");
        capacity = count;
    }
    void zero(size_t count)
    {
        ensure(count);
        if(count) cuda_check(cudaMemset(ptr, 0, count * sizeof(T)), "cudaMemset");
    }
    void upload(const T * source, size_t count)
    {
        ensure(count);
        if(count)
            cuda_check(cudaMemcpy(ptr, source, count * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy H2D");
    }
};

__host__ __device__ inline int unordered_state_count(int C)
{
    return C * (C + 1) / 2;
}

__host__ __device__ inline int unordered_state_index(int z1, int z2, int C)
{
    if(z1 > z2)
    {
        const int tmp = z1;
        z1 = z2;
        z2 = tmp;
    }
    return z1 * C - z1 * (z1 - 1) / 2 + (z2 - z1);
}

__device__ inline Scalar calculate_emission(const Scalar * gl,
                                            const Scalar * p,
                                            int local_ind,
                                            int s,
                                            int S,
                                            int M,
                                            int m,
                                            int z1,
                                            int z2)
{
    const Scalar p1 = p[m + M * z1];
    const Scalar p2 = p[m + M * z2];
    const Scalar g0 = gl[local_ind * S * 3 + s];
    const Scalar g1 = gl[local_ind * S * 3 + S + s];
    const Scalar g2 = gl[local_ind * S * 3 + 2 * S + s];
    const Scalar out = g0 * (1 - p1) * (1 - p2) + g1 * ((1 - p1) * p2 + p1 * (1 - p2)) + g2 * p1 * p2;
    return out < Scalar(1e-10) ? Scalar(1e-10) : out;
}

__global__ void cache_inputs_kernel(const Scalar * gl,
                                    const Scalar * p,
                                    const Scalar * q,
                                    const Scalar * f,
                                    int K,
                                    int C,
                                    int M,
                                    int S,
                                    int global_start,
                                    int first_ind,
                                    const int * pair_z1,
                                    const int * pair_z2,
                                    Scalar * cached_h,
                                    Scalar * cached_emit)
{
    const int local_ind = blockIdx.x;
    const int ind = first_ind + local_ind;
    const int tid = threadIdx.x;
    const int U = unordered_state_count(C);
    for(int s = 0; s < S; ++s)
    {
        const int m = global_start + s;
        if(tid < C)
        {
            Scalar value = 0;
            for(int k = 0; k < K; ++k) value += q[k + K * ind] * f[k * C * M + tid + C * m];
            cached_h[(local_ind * S + s) * C + tid] = value;
        }
        if(tid < U)
        {
            const int z1 = pair_z1[tid];
            const int z2 = pair_z2[tid];
            cached_emit[(local_ind * S + s) * U + tid] =
                calculate_emission(gl, p, local_ind, s, S, M, m, z1, z2);
        }
    }
}

__global__ void forward_kernel(const Scalar * cached_h,
                               const Scalar * cached_emit,
                               const Scalar * r,
                               int C,
                               int S,
                               int global_start,
                               int first_ind,
                               const int * pair_z1,
                               const int * pair_z2,
                               Scalar * alpha,
                               Scalar * cs,
                               Scalar * likelihood)
{
    const int local_ind = blockIdx.x;
    const int ind = first_ind + local_ind;
    const int tid = threadIdx.x;
    const int U = unordered_state_count(C);
    extern __shared__ Scalar shared[];
    Scalar * previous = shared;
    Scalar * current = previous + U;
    Scalar * row_sum = current + U;
    Scalar * scale = row_sum + C;
    Scalar local_likelihood = 0;
    for(int s = 0; s < S; ++s)
    {
        const int m = global_start + s;
        const Scalar * h = cached_h + (local_ind * S + s) * C;
        const Scalar * emit = cached_emit + (local_ind * S + s) * U;
        if(tid < C)
        {
            Scalar rows = 0;
            if(s)
                for(int other = 0; other < C; ++other) rows += previous[unordered_state_index(tid, other, C)];
            row_sum[tid] = rows * r[1 + 3 * m];
        }
        __syncthreads();
        if(tid < U)
        {
            const int z1 = pair_z1[tid];
            const int z2 = pair_z2[tid];
            if(s == 0)
                current[tid] = emit[tid] * h[z1] * h[z2];
            else
                current[tid] = emit[tid]
                               * (previous[tid] * r[3 * m] + h[z1] * row_sum[z2] + h[z2] * row_sum[z1]
                                  + h[z1] * h[z2] * r[2 + 3 * m]);
        }
        __syncthreads();
        if(tid == 0)
        {
            Scalar total = 0;
            for(int state = 0; state < U; ++state)
                total += (pair_z1[state] == pair_z2[state] ? Scalar(1) : Scalar(2)) * current[state];
            scale[0] = Scalar(1) / total;
            cs[local_ind * S + s] = scale[0];
            local_likelihood += log(total);
        }
        __syncthreads();
        if(tid < U)
        {
            current[tid] *= scale[0];
            alpha[(local_ind * S + s) * U + tid] = current[tid];
            previous[tid] = current[tid];
        }
        __syncthreads();
    }
    if(tid == 0) likelihood[ind] += local_likelihood;
}

__global__ void backward_posterior_kernel(const Scalar * gl,
                                          const Scalar * p,
                                          const Scalar * q,
                                          const Scalar * f,
                                          const Scalar * r,
                                          const Scalar * cached_h,
                                          const Scalar * cached_emit,
                                          const Scalar * alpha,
                                          const Scalar * cs,
                                          int K,
                                          int C,
                                          int M,
                                          int S,
                                          int global_start,
                                          int first_ind,
                                          bool final_iteration,
                                          const int * pair_z1,
                                          const int * pair_z2,
                                          Scalar * partial_a1,
                                          Scalar * partial_a2,
                                          Scalar * partial_cluster,
                                          Scalar * gp)
{
    const int local_ind = blockIdx.x;
    const int ind = first_ind + local_ind;
    const int tid = threadIdx.x;
    const int U = unordered_state_count(C);
    extern __shared__ Scalar shared[];
    Scalar * beta = shared;
    Scalar * previous_beta = beta + U;
    Scalar * row_sum = previous_beta + U;
    Scalar * alpha_marginal = row_sum + C;
    Scalar * constant = alpha_marginal + C;
    if(tid < U) beta[tid] = 1;
    __syncthreads();
    for(int s = S - 1; s >= 0; --s)
    {
        const int m = global_start + s;
        const Scalar * h = cached_h + (local_ind * S + s) * C;
        const Scalar * emit = cached_emit + (local_ind * S + s) * U;
        if(s > 0 && tid < C)
        {
            Scalar value = 0;
            for(int other = 0; other < C; ++other)
                value += alpha[(local_ind * S + s - 1) * U + unordered_state_index(tid, other, C)];
            alpha_marginal[tid] = value;
        }
        __syncthreads();
        if(tid < C)
        {
            const int z = tid;
            const Scalar pz = p[m + M * z];
            const Scalar g0 = gl[local_ind * S * 3 + s];
            const Scalar g1 = gl[local_ind * S * 3 + S + s];
            const Scalar g2 = gl[local_ind * S * 3 + 2 * S + s];
            Scalar a1 = 0, a2 = 0;
            for(int other = 0; other < C; ++other)
            {
                const int state = unordered_state_index(z, other, C);
                const Scalar gamma_over_emit =
                    alpha[(local_ind * S + s) * U + state] * beta[state] / emit[state];
                const Scalar po = p[m + M * other];
                a1 += gamma_over_emit * (1 - pz) * (g0 * (1 - po) + g1 * po);
                a2 += gamma_over_emit * pz * (g1 * (1 - po) + g2 * po);
            }
            partial_a1[(local_ind * S + s) * C + z] = a1;
            partial_a2[(local_ind * S + s) * C + z] = a2;

            Scalar refresh_weight = 0;
            if(s == 0)
                for(int other = 0; other < C; ++other)
                {
                    const int state = unordered_state_index(z, other, C);
                    refresh_weight += alpha[local_ind * S * U + state] * beta[state];
                }
            Scalar transition_weight = 0;
            if(s > 0)
            {
                for(int z2 = 0; z2 < C; ++z2)
                {
                    const int state = unordered_state_index(z, z2, C);
                    transition_weight += emit[state] * beta[state]
                                         * (r[1 + 3 * m] * alpha_marginal[z2] + r[2 + 3 * m] * h[z2]);
                }
                transition_weight *= cs[local_ind * S + s];
            }
            for(int k = 0; k < K; ++k)
            {
                const Scalar count = s == 0
                                         ? refresh_weight * q[k + K * ind] * f[k * C * M + z + C * m] / h[z]
                                         : transition_weight * q[k + K * ind] * f[k * C * M + z + C * m];
                partial_cluster[((local_ind * S + s) * K + k) * C + z] = count;
            }
        }
        if(final_iteration && tid == 0)
        {
            const Scalar g0 = gl[local_ind * S * 3 + s];
            const Scalar g1 = gl[local_ind * S * 3 + S + s];
            const Scalar g2 = gl[local_ind * S * 3 + 2 * S + s];
            Scalar gp0 = 0, gp1 = 0, gp2 = 0;
            for(int state = 0; state < U; ++state)
            {
                const int z1 = pair_z1[state];
                const int z2 = pair_z2[state];
                const Scalar multiplicity = z1 == z2 ? Scalar(1) : Scalar(2);
                const Scalar gamma_over_emit =
                    multiplicity * alpha[(local_ind * S + s) * U + state] * beta[state] / emit[state];
                const Scalar p1 = p[m + M * z1];
                const Scalar p2 = p[m + M * z2];
                gp0 += gamma_over_emit * g0 * (1 - p1) * (1 - p2);
                gp1 += gamma_over_emit * g1 * ((1 - p1) * p2 + p1 * (1 - p2));
                gp2 += gamma_over_emit * g2 * p1 * p2;
            }
            gp[3 * m + 3 * M * ind] = gp0;
            gp[3 * m + 1 + 3 * M * ind] = gp1;
            gp[3 * m + 2 + 3 * M * ind] = gp2;
        }
        __syncthreads();
        if(s > 0)
        {
            if(tid == 0)
            {
                for(int z1 = 0; z1 < C; ++z1)
                {
                    Scalar sum = 0;
                    for(int z2 = 0; z2 < C; ++z2)
                    {
                        const int state = unordered_state_index(z1, z2, C);
                        sum += emit[state] * beta[state] * h[z2] * r[1 + 3 * m];
                    }
                    row_sum[z1] = sum;
                }
                Scalar total = 0;
                for(int state = 0; state < U; ++state)
                {
                    const int z1 = pair_z1[state];
                    const int z2 = pair_z2[state];
                    const Scalar multiplicity = z1 == z2 ? Scalar(1) : Scalar(2);
                    total += multiplicity * emit[state] * beta[state] * h[z1] * h[z2] * r[2 + 3 * m];
                }
                constant[0] = total;
            }
            __syncthreads();
            if(tid < U)
            {
                const int z1 = pair_z1[tid];
                const int z2 = pair_z2[tid];
                previous_beta[tid] =
                    (emit[tid] * beta[tid] * r[3 * m] + row_sum[z1] + row_sum[z2] + constant[0])
                    * cs[local_ind * S + s];
            }
            __syncthreads();
            if(tid < U) beta[tid] = previous_beta[tid];
            __syncthreads();
        }
    }
}

__global__ void reduce_site_partials_kernel(const Scalar * partial_a1,
                                            const Scalar * partial_a2,
                                            int batch,
                                            int C,
                                            int S,
                                            int global_start,
                                            Scalar * e_a1,
                                            Scalar * e_a2)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index >= C * S) return;
    const int z = index % C;
    const int s = index / C;
    Scalar a1 = 0, a2 = 0;
    for(int local_ind = 0; local_ind < batch; ++local_ind)
    {
        const int partial_index = (local_ind * S + s) * C + z;
        a1 += partial_a1[partial_index];
        a2 += partial_a2[partial_index];
    }
    const int m = global_start + s;
    e_a1[z + C * m] += a1;
    e_a2[z + C * m] += a2;
}

__global__ void reduce_cluster_partials_kernel(const Scalar * partial_cluster,
                                               int batch,
                                               int K,
                                               int C,
                                               int S,
                                               int global_start,
                                               Scalar * e_cluster)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index >= C * K * S) return;
    const int z = index % C;
    const int k = (index / C) % K;
    const int s = index / (C * K);
    Scalar total = 0;
    for(int local_ind = 0; local_ind < batch; ++local_ind)
        total += partial_cluster[((local_ind * S + s) * K + k) * C + z];
    e_cluster[z + C * k + C * K * (global_start + s)] += total;
}

__global__ void reduce_ancestry_partials_kernel(const Scalar * partial_cluster,
                                                int batch,
                                                int K,
                                                int C,
                                                int S,
                                                int first_ind,
                                                Scalar * e_ancestry)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index >= batch * K) return;
    const int k = index % K;
    const int local_ind = index / K;
    Scalar total = 0;
    for(int s = 0; s < S; ++s)
        for(int z = 0; z < C; ++z) total += partial_cluster[((local_ind * S + s) * K + k) * C + z];
    e_ancestry[k + K * (first_ind + local_ind)] += total;
}

struct CudaWorkspace
{
    const Phaseless * model_owner;
    const MyFloat2D * gl_owner;
    int K, C, N, M;
    std::vector<int> chunk_sizes;
    std::vector<int> batch_capacities;
    std::vector<DeviceBuffer<Scalar>> gl_chunks;
    std::vector<Scalar> flat_f;
    DeviceBuffer<Scalar> p, q, f, r;
    DeviceBuffer<int> pair_z1, pair_z2;
    DeviceBuffer<Scalar> e_a1, e_a2, e_ancestry, e_cluster, gp, likelihood;
    DeviceBuffer<Scalar> alpha, cs, cached_h, cached_emit;
    DeviceBuffer<Scalar> partial_a1, partial_a2, partial_cluster;

    CudaWorkspace(const Phaseless & model, const MyFloat2D & gl)
    : model_owner(&model), gl_owner(&gl), K(model.K), C(model.C), N(model.N), M(model.M),
      flat_f(static_cast<size_t>(K) * C * M), p(model.P.size()), q(model.Q.size()), f(flat_f.size()),
      r(model.R.size()), pair_z1(unordered_state_count(C)), pair_z2(unordered_state_count(C)), e_a1(C * M),
      e_a2(C * M), e_ancestry(K * N), e_cluster(C * K * M), likelihood(N)
    {
        std::vector<int> host_z1, host_z2;
        host_z1.reserve(unordered_state_count(C));
        host_z2.reserve(unordered_state_count(C));
        for(int z1 = 0; z1 < C; ++z1)
            for(int z2 = z1; z2 < C; ++z2)
            {
                host_z1.push_back(z1);
                host_z2.push_back(z2);
            }
        pair_z1.upload(host_z1.data(), host_z1.size());
        pair_z2.upload(host_z2.data(), host_z2.size());
        chunk_sizes.reserve(gl.size());
        batch_capacities.assign(gl.size(), 0);
        gl_chunks.reserve(gl.size());
        for(size_t chunk = 0; chunk < gl.size(); ++chunk)
        {
            const int S = model.pos_chunk[chunk + 1] - model.pos_chunk[chunk];
            chunk_sizes.push_back(S);
            gl_chunks.emplace_back(gl[chunk].size());
            gl_chunks.back().upload(gl[chunk].data(), gl[chunk].size());
        }
    }
    bool matches(const Phaseless & model, const MyFloat2D & gl) const
    {
        if(model_owner != &model || gl_owner != &gl || K != model.K || C != model.C || N != model.N
           || M != model.M || gl.size() != chunk_sizes.size())
            return false;
        for(size_t chunk = 0; chunk < gl.size(); ++chunk)
            if(chunk_sizes[chunk] != model.pos_chunk[chunk + 1] - model.pos_chunk[chunk]) return false;
        return true;
    }
    void upload_parameters(const Phaseless & model)
    {
        for(int k = 0; k < K; ++k)
            std::copy(model.F[k].data(), model.F[k].data() + C * M,
                      flat_f.begin() + static_cast<size_t>(k) * C * M);
        p.upload(model.P.data(), model.P.size());
        q.upload(model.Q.data(), model.Q.size());
        f.upload(flat_f.data(), flat_f.size());
        r.upload(model.R.data(), model.R.size());
    }
    void clear_outputs(bool final_iteration)
    {
        e_a1.zero(static_cast<size_t>(C) * M);
        e_a2.zero(static_cast<size_t>(C) * M);
        e_ancestry.zero(static_cast<size_t>(K) * N);
        e_cluster.zero(static_cast<size_t>(C) * K * M);
        likelihood.zero(N);
        if(final_iteration) gp.zero(static_cast<size_t>(M) * 3 * N);
    }
    void ensure_batch(int batch, int S)
    {
        const size_t sites = static_cast<size_t>(batch) * S;
        const int U = unordered_state_count(C);
        alpha.ensure(sites * U);
        cs.ensure(sites);
        cached_h.ensure(sites * C);
        cached_emit.ensure(sites * U);
        partial_a1.ensure(sites * C);
        partial_a2.ensure(sites * C);
        partial_cluster.ensure(sites * K * C);
    }
};

std::unique_ptr<CudaWorkspace> workspace;

CudaWorkspace & get_workspace(const Phaseless & model, const MyFloat2D & gl)
{
    if(!workspace || !workspace->matches(model, gl)) workspace = std::make_unique<CudaWorkspace>(model, gl);
    return *workspace;
}
} // namespace

bool joint_cuda_available(std::string & reason)
{
    int count = 0;
    const cudaError_t status = cudaGetDeviceCount(&count);
    if(status != cudaSuccess)
    {
        reason = cudaGetErrorString(status);
        cudaGetLastError();
        return false;
    }
    if(count == 0)
    {
        reason = "no CUDA-capable device was found";
        return false;
    }
    reason.clear();
    return true;
}

double joint_cuda_e_step(Phaseless & model, const MyFloat2D & gl, bool final_iteration)
{
    std::string reason;
    if(!joint_cuda_available(reason)) throw std::runtime_error("CUDA unavailable: " + reason);
    if(model.C < 1 || model.C > 32)
        throw std::invalid_argument("CUDA joint model requires cluster count between 1 and 32");
    if(model.K < 1 || model.N < 1 || model.M < 1)
        throw std::invalid_argument("CUDA joint model requires positive ancestry, sample, and site counts");
    if(gl.size() + 1 != model.pos_chunk.size())
        throw std::invalid_argument("CUDA joint model received inconsistent chunks");

    CudaWorkspace & w = get_workspace(model, gl);
    w.upload_parameters(model);
    w.clear_outputs(final_iteration);
    size_t free_bytes = 0, total_bytes = 0;
    cuda_check(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
    const int U = unordered_state_count(model.C);
    const size_t forward_shared = static_cast<size_t>(2 * U + model.C + 1) * sizeof(Scalar);
    const size_t backward_shared = static_cast<size_t>(2 * U + 2 * model.C + 1) * sizeof(Scalar);
    constexpr int reduction_threads = 256;

    for(size_t chunk = 0; chunk < gl.size(); ++chunk)
    {
        const int S = model.pos_chunk[chunk + 1] - model.pos_chunk[chunk];
        const size_t elements_per_ind =
            static_cast<size_t>(S) * (2 * U + 3 * model.C + model.K * model.C + 1);
        const size_t budget = std::min(free_bytes / 2, total_bytes / 4);
        if(w.batch_capacities[chunk] == 0)
        {
            const size_t capacity_by_memory = budget / std::max<size_t>(1, elements_per_ind * sizeof(Scalar));
            w.batch_capacities[chunk] =
                std::max(1, static_cast<int>(std::min<size_t>(model.N, capacity_by_memory)));
        }
        const int batch_capacity = w.batch_capacities[chunk];
        for(int first = 0; first < model.N; first += batch_capacity)
        {
            const int batch = std::min(batch_capacity, model.N - first);
            w.ensure_batch(batch, S);
            const Scalar * batch_gl = w.gl_chunks[chunk].ptr + static_cast<size_t>(first) * S * 3;
            cache_inputs_kernel<<<batch, U>>>(batch_gl, w.p.ptr, w.q.ptr, w.f.ptr, model.K, model.C, model.M,
                                              S, model.pos_chunk[chunk], first, w.pair_z1.ptr, w.pair_z2.ptr,
                                              w.cached_h.ptr, w.cached_emit.ptr);
            cuda_check(cudaGetLastError(), "launch CUDA input-cache kernel");
            forward_kernel<<<batch, U, forward_shared>>>(
                w.cached_h.ptr, w.cached_emit.ptr, w.r.ptr, model.C, S, model.pos_chunk[chunk], first,
                w.pair_z1.ptr, w.pair_z2.ptr, w.alpha.ptr, w.cs.ptr, w.likelihood.ptr);
            cuda_check(cudaGetLastError(), "launch CUDA forward kernel");
            backward_posterior_kernel<<<batch, U, backward_shared>>>(
                batch_gl, w.p.ptr, w.q.ptr, w.f.ptr, w.r.ptr, w.cached_h.ptr, w.cached_emit.ptr, w.alpha.ptr,
                w.cs.ptr, model.K, model.C, model.M, S, model.pos_chunk[chunk], first, final_iteration,
                w.pair_z1.ptr, w.pair_z2.ptr, w.partial_a1.ptr, w.partial_a2.ptr, w.partial_cluster.ptr,
                w.gp.ptr);
            cuda_check(cudaGetLastError(), "launch CUDA backward/posterior kernel");

            const int site_values = model.C * S;
            reduce_site_partials_kernel<<<(site_values + reduction_threads - 1) / reduction_threads,
                                          reduction_threads>>>(w.partial_a1.ptr, w.partial_a2.ptr, batch,
                                                               model.C, S, model.pos_chunk[chunk], w.e_a1.ptr,
                                                               w.e_a2.ptr);
            const int cluster_values = model.C * model.K * S;
            reduce_cluster_partials_kernel<<<(cluster_values + reduction_threads - 1) / reduction_threads,
                                             reduction_threads>>>(
                w.partial_cluster.ptr, batch, model.K, model.C, S, model.pos_chunk[chunk], w.e_cluster.ptr);
            const int ancestry_values = batch * model.K;
            reduce_ancestry_partials_kernel<<<(ancestry_values + reduction_threads - 1) / reduction_threads,
                                              reduction_threads>>>(w.partial_cluster.ptr, batch, model.K,
                                                                   model.C, S, first, w.e_ancestry.ptr);
            cuda_check(cudaGetLastError(), "launch CUDA sufficient-statistic reduction kernels");
            cuda_check(cudaDeviceSynchronize(), "execute CUDA joint-model batch");
        }
    }

    cuda_check(cudaMemcpy(model.EclusterA1.data(), w.e_a1.ptr, model.C * model.M * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy EclusterA1");
    cuda_check(cudaMemcpy(model.EclusterA2.data(), w.e_a2.ptr, model.C * model.M * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy EclusterA2");
    cuda_check(cudaMemcpy(model.Eancestry.data(), w.e_ancestry.ptr, model.K * model.N * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy Eancestry");
    cuda_check(cudaMemcpy(model.EclusterK.data(), w.e_cluster.ptr,
                          model.C * model.K * model.M * sizeof(Scalar), cudaMemcpyDeviceToHost),
               "copy EclusterK");
    if(final_iteration)
        cuda_check(cudaMemcpy(model.GP.data(), w.gp.ptr, model.M * 3 * model.N * sizeof(Scalar),
                              cudaMemcpyDeviceToHost),
                   "copy GP");
    std::vector<Scalar> individual_likelihood(model.N);
    cuda_check(cudaMemcpy(individual_likelihood.data(), w.likelihood.ptr, model.N * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy likelihood");
    return std::accumulate(individual_likelihood.begin(), individual_likelihood.end(), 0.0);
}
