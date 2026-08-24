#include "joint_cuda.hpp"
#include "phaseless.hpp"
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
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
    DeviceBuffer() = default;
    explicit DeviceBuffer(size_t count)
    {
        if(count) cuda_check(cudaMalloc(reinterpret_cast<void **>(&ptr), count * sizeof(T)), "cudaMalloc");
    }
    ~DeviceBuffer()
    {
        if(ptr) cudaFree(ptr);
    }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer & operator=(const DeviceBuffer &) = delete;
};

__device__ inline Scalar atomic_add(Scalar * address, Scalar value)
{
    return atomicAdd(address, value);
}

__device__ inline Scalar emission(const Scalar * gl,
                                  const Scalar * p,
                                  int local_ind,
                                  int s,
                                  int S,
                                  int M,
                                  int z1,
                                  int z2)
{
    const Scalar p1 = p[s + M * z1];
    const Scalar p2 = p[s + M * z2];
    const Scalar g0 = gl[local_ind * S * 3 + s];
    const Scalar g1 = gl[local_ind * S * 3 + S + s];
    const Scalar g2 = gl[local_ind * S * 3 + 2 * S + s];
    Scalar out = g0 * (1 - p1) * (1 - p2) + g1 * ((1 - p1) * p2 + p1 * (1 - p2)) + g2 * p1 * p2;
    return out < Scalar(1e-10) ? Scalar(1e-10) : out;
}

__global__ void forward_kernel(const Scalar * gl,
                               const Scalar * p,
                               const Scalar * q,
                               const Scalar * f,
                               const Scalar * r,
                               int K,
                               int C,
                               int M,
                               int S,
                               int global_start,
                               int first_ind,
                               Scalar * alpha,
                               Scalar * cs,
                               Scalar * likelihood)
{
    const int local_ind = blockIdx.x;
    const int ind = first_ind + local_ind;
    const int tid = threadIdx.x;
    const int CC = C * C;
    extern __shared__ Scalar shared[];
    Scalar * previous = shared;
    Scalar * current = previous + CC;
    Scalar * h = current + CC;
    Scalar * row_sum = h + C;
    Scalar * scale = row_sum + C;

    for(int s = 0; s < S; ++s)
    {
        const int m = global_start + s;
        if(tid < C)
        {
            Scalar value = 0;
            for(int k = 0; k < K; ++k) value += q[k + K * ind] * f[k * C * M + tid + C * m];
            h[tid] = value;
            Scalar rows = 0;
            if(s)
                for(int other = 0; other < C; ++other) rows += previous[tid + C * other];
            row_sum[tid] = rows * r[1 + 3 * m];
        }
        __syncthreads();
        if(tid < CC)
        {
            const int z1 = tid / C;
            const int z2 = tid - z1 * C;
            const Scalar e = emission(gl, p + global_start, local_ind, s, S, M, z1, z2);
            if(s == 0)
                current[tid] = e * h[z1] * h[z2];
            else
                current[tid] = e
                               * (previous[tid] * r[3 * m] + h[z1] * row_sum[z2] + h[z2] * row_sum[z1]
                                  + h[z1] * h[z2] * r[2 + 3 * m]);
        }
        __syncthreads();
        if(tid == 0)
        {
            Scalar total = 0;
            for(int state = 0; state < CC; ++state) total += current[state];
            scale[0] = Scalar(1) / total;
            cs[local_ind * S + s] = scale[0];
            atomic_add(likelihood, log(total));
        }
        __syncthreads();
        if(tid < CC)
        {
            current[tid] *= scale[0];
            alpha[(local_ind * S + s) * CC + tid] = current[tid];
            previous[tid] = current[tid];
        }
        __syncthreads();
    }
}

__global__ void backward_posterior_kernel(const Scalar * gl,
                                          const Scalar * p,
                                          const Scalar * q,
                                          const Scalar * f,
                                          const Scalar * r,
                                          const Scalar * alpha,
                                          const Scalar * cs,
                                          int K,
                                          int C,
                                          int M,
                                          int S,
                                          int global_start,
                                          int first_ind,
                                          bool final_iteration,
                                          Scalar * e_a1,
                                          Scalar * e_a2,
                                          Scalar * e_ancestry,
                                          Scalar * e_cluster,
                                          Scalar * gp)
{
    const int local_ind = blockIdx.x;
    const int ind = first_ind + local_ind;
    const int tid = threadIdx.x;
    const int CC = C * C;
    extern __shared__ Scalar shared[];
    Scalar * beta = shared;
    Scalar * previous_beta = beta + CC;
    Scalar * h = previous_beta + CC;
    Scalar * row_sum = h + C;
    Scalar * constant = row_sum + C;

    if(tid < CC) beta[tid] = 1;
    __syncthreads();
    for(int s = S - 1; s >= 0; --s)
    {
        const int m = global_start + s;
        if(tid < C)
        {
            Scalar value = 0;
            for(int k = 0; k < K; ++k) value += q[k + K * ind] * f[k * C * M + tid + C * m];
            h[tid] = value;
        }
        __syncthreads();

        if(tid < CC)
        {
            // tid = z + C*other matches Eigen's strided cluster slice.
            const int z = tid % C;
            const int other = tid / C;
            const Scalar e = emission(gl, p + global_start, local_ind, s, S, M, other, z);
            const Scalar gamma_over_emit = alpha[(local_ind * S + s) * CC + tid] * beta[tid] / e;
            const Scalar pz = p[m + M * z];
            const Scalar po = p[m + M * other];
            const Scalar g0 = gl[local_ind * S * 3 + s];
            const Scalar g1 = gl[local_ind * S * 3 + S + s];
            const Scalar g2 = gl[local_ind * S * 3 + 2 * S + s];
            atomic_add(e_a1 + z + C * m, gamma_over_emit * (1 - pz) * (g0 * (1 - po) + g1 * po));
            atomic_add(e_a2 + z + C * m, gamma_over_emit * pz * (g1 * (1 - po) + g2 * po));
            if(final_iteration)
            {
                atomic_add(gp + 3 * m + 3 * M * ind, gamma_over_emit * g0 * (1 - pz) * (1 - po));
                atomic_add(gp + 3 * m + 1 + 3 * M * ind,
                           gamma_over_emit * g1 * ((1 - pz) * po + pz * (1 - po)));
                atomic_add(gp + 3 * m + 2 + 3 * M * ind, gamma_over_emit * g2 * pz * po);
            }
        }
        __syncthreads();

        if(tid < C)
        {
            Scalar refresh_weight = 0;
            if(s == 0)
            {
                for(int other = 0; other < C; ++other)
                {
                    const int state = tid * C + other;
                    refresh_weight += alpha[(local_ind * S) * CC + state] * beta[state];
                }
                for(int k = 0; k < K; ++k)
                {
                    const Scalar count =
                        refresh_weight * q[k + K * ind] * f[k * C * M + tid + C * m] / h[tid];
                    atomic_add(e_ancestry + k + K * ind, count);
                    atomic_add(e_cluster + k * C + tid + C * K * m, count);
                }
            }
            else
            {
                Scalar transition_weight = 0;
                for(int z2 = 0; z2 < C; ++z2)
                {
                    Scalar alpha_previous_sum = 0;
                    for(int other = 0; other < C; ++other)
                        alpha_previous_sum += alpha[(local_ind * S + s - 1) * CC + z2 + C * other];
                    const int state = tid * C + z2;
                    const Scalar eb =
                        emission(gl, p + global_start, local_ind, s, S, M, tid, z2) * beta[state];
                    transition_weight += eb * (r[1 + 3 * m] * alpha_previous_sum + r[2 + 3 * m] * h[z2]);
                }
                transition_weight *= cs[local_ind * S + s];
                for(int k = 0; k < K; ++k)
                {
                    const Scalar count = transition_weight * q[k + K * ind] * f[k * C * M + tid + C * m];
                    atomic_add(e_ancestry + k + K * ind, count);
                    atomic_add(e_cluster + k * C + tid + C * K * m, count);
                }
            }
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
                        const int state = z1 * C + z2;
                        sum += emission(gl, p + global_start, local_ind, s, S, M, z1, z2) * beta[state]
                               * h[z2] * r[1 + 3 * m];
                    }
                    row_sum[z1] = sum;
                }
                Scalar total = 0;
                for(int z1 = 0; z1 < C; ++z1)
                    for(int z2 = 0; z2 < C; ++z2)
                    {
                        const int state = z1 * C + z2;
                        total += emission(gl, p + global_start, local_ind, s, S, M, z1, z2) * beta[state]
                                 * h[z1] * h[z2] * r[2 + 3 * m];
                    }
                constant[0] = total;
            }
            __syncthreads();
            if(tid < CC)
            {
                const int z1 = tid / C;
                const int z2 = tid - z1 * C;
                const Scalar e = emission(gl, p + global_start, local_ind, s, S, M, z1, z2);
                previous_beta[tid] = (e * beta[tid] * r[3 * m] + row_sum[z1] + row_sum[z2] + constant[0])
                                     * cs[local_ind * S + s];
            }
            __syncthreads();
            if(tid < CC) beta[tid] = previous_beta[tid];
            __syncthreads();
        }
    }
}

template<typename T>
void copy_to_device(DeviceBuffer<T> & destination, const T * source, size_t count)
{
    cuda_check(cudaMemcpy(destination.ptr, source, count * sizeof(T), cudaMemcpyHostToDevice),
               "cudaMemcpy H2D");
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
    if(model.C * model.C > 1024) throw std::invalid_argument("CUDA joint model requires cluster count <= 32");
    if(gl.size() + 1 != model.pos_chunk.size())
        throw std::invalid_argument("CUDA joint model received inconsistent chunks");

    std::vector<Scalar> flat_f(static_cast<size_t>(model.K) * model.C * model.M);
    for(int k = 0; k < model.K; ++k)
        std::copy(model.F[k].data(), model.F[k].data() + model.C * model.M,
                  flat_f.begin() + static_cast<size_t>(k) * model.C * model.M);

    DeviceBuffer<Scalar> d_p(model.P.size()), d_q(model.Q.size()), d_f(flat_f.size()), d_r(model.R.size());
    DeviceBuffer<Scalar> d_a1(model.C * model.M), d_a2(model.C * model.M);
    DeviceBuffer<Scalar> d_ancestry(model.K * model.N), d_cluster(model.C * model.K * model.M);
    DeviceBuffer<Scalar> d_gp(final_iteration ? model.M * 3 * model.N : 0), d_likelihood(1);
    copy_to_device(d_p, model.P.data(), model.P.size());
    copy_to_device(d_q, model.Q.data(), model.Q.size());
    copy_to_device(d_f, flat_f.data(), flat_f.size());
    copy_to_device(d_r, model.R.data(), model.R.size());
    cuda_check(cudaMemset(d_a1.ptr, 0, model.C * model.M * sizeof(Scalar)), "cudaMemset EclusterA1");
    cuda_check(cudaMemset(d_a2.ptr, 0, model.C * model.M * sizeof(Scalar)), "cudaMemset EclusterA2");
    cuda_check(cudaMemset(d_ancestry.ptr, 0, model.K * model.N * sizeof(Scalar)), "cudaMemset Eancestry");
    cuda_check(cudaMemset(d_cluster.ptr, 0, model.C * model.K * model.M * sizeof(Scalar)),
               "cudaMemset EclusterK");
    if(final_iteration)
        cuda_check(cudaMemset(d_gp.ptr, 0, model.M * 3 * model.N * sizeof(Scalar)), "cudaMemset GP");
    cuda_check(cudaMemset(d_likelihood.ptr, 0, sizeof(Scalar)), "cudaMemset likelihood");

    size_t free_bytes = 0, total_bytes = 0;
    cuda_check(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
    const int CC = model.C * model.C;
    const size_t shared_bytes = static_cast<size_t>(2 * CC + 2 * model.C + 1) * sizeof(Scalar);
    for(size_t chunk = 0; chunk < gl.size(); ++chunk)
    {
        const int S = model.pos_chunk[chunk + 1] - model.pos_chunk[chunk];
        const size_t bytes_per_ind = static_cast<size_t>(S) * (CC + 4) * sizeof(Scalar);
        const size_t budget = std::min(free_bytes / 2, total_bytes / 4);
        const int batch_capacity =
            std::max(1, std::min(model.N, static_cast<int>(budget / std::max<size_t>(1, bytes_per_ind))));
        for(int first = 0; first < model.N; first += batch_capacity)
        {
            const int batch = std::min(batch_capacity, model.N - first);
            DeviceBuffer<Scalar> d_gl(static_cast<size_t>(batch) * S * 3);
            DeviceBuffer<Scalar> d_alpha(static_cast<size_t>(batch) * S * CC);
            DeviceBuffer<Scalar> d_cs(static_cast<size_t>(batch) * S);
            copy_to_device(d_gl, gl[chunk].data() + static_cast<size_t>(first) * S * 3,
                           static_cast<size_t>(batch) * S * 3);
            forward_kernel<<<batch, CC, shared_bytes>>>(d_gl.ptr, d_p.ptr, d_q.ptr, d_f.ptr, d_r.ptr, model.K,
                                                        model.C, model.M, S, model.pos_chunk[chunk], first,
                                                        d_alpha.ptr, d_cs.ptr, d_likelihood.ptr);
            cuda_check(cudaGetLastError(), "launch CUDA forward kernel");
            backward_posterior_kernel<<<batch, CC, shared_bytes>>>(
                d_gl.ptr, d_p.ptr, d_q.ptr, d_f.ptr, d_r.ptr, d_alpha.ptr, d_cs.ptr, model.K, model.C,
                model.M, S, model.pos_chunk[chunk], first, final_iteration, d_a1.ptr, d_a2.ptr,
                d_ancestry.ptr, d_cluster.ptr, d_gp.ptr);
            cuda_check(cudaGetLastError(), "launch CUDA backward/posterior kernel");
            cuda_check(cudaDeviceSynchronize(), "execute CUDA joint-model kernels");
        }
    }

    cuda_check(cudaMemcpy(model.EclusterA1.data(), d_a1.ptr, model.C * model.M * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy EclusterA1");
    cuda_check(cudaMemcpy(model.EclusterA2.data(), d_a2.ptr, model.C * model.M * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy EclusterA2");
    cuda_check(cudaMemcpy(model.Eancestry.data(), d_ancestry.ptr, model.K * model.N * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy Eancestry");
    cuda_check(cudaMemcpy(model.EclusterK.data(), d_cluster.ptr, model.C * model.K * model.M * sizeof(Scalar),
                          cudaMemcpyDeviceToHost),
               "copy EclusterK");
    if(final_iteration)
        cuda_check(cudaMemcpy(model.GP.data(), d_gp.ptr, model.M * 3 * model.N * sizeof(Scalar),
                              cudaMemcpyDeviceToHost),
                   "copy GP");
    Scalar likelihood = 0;
    cuda_check(cudaMemcpy(&likelihood, d_likelihood.ptr, sizeof(Scalar), cudaMemcpyDeviceToHost),
               "copy likelihood");
    return likelihood;
}
