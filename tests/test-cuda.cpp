#define _DECLARE_TOOLBOX_HERE

#include "../src/joint_cuda.hpp"
#include "../src/phaseless.hpp"
#include "catch.hh"

namespace
{
void require_cuda_matches_cpu(Phaseless & cpu,
                              Phaseless & gpu,
                              const MyFloat2D & gl,
                              bool final_iteration = false)
{
    cpu.initIteration();
    gpu.initIteration();
    if(final_iteration)
    {
        cpu.GP.setZero(cpu.M * 3, cpu.N);
        gpu.GP.setZero(gpu.M * 3, gpu.N);
    }
    double cpu_likelihood = 0;
    for(int ind = 0; ind < cpu.N; ++ind) cpu_likelihood += cpu.runBigass(ind, gl, final_iteration);
    const double gpu_likelihood = joint_cuda_e_step(gpu, gl, final_iteration);

    REQUIRE(gpu_likelihood == Approx(cpu_likelihood).epsilon(2e-4));
    REQUIRE((gpu.EclusterA1 - cpu.EclusterA1).abs().maxCoeff() < 2e-4);
    REQUIRE((gpu.EclusterA2 - cpu.EclusterA2).abs().maxCoeff() < 2e-4);
    REQUIRE((gpu.Eancestry - cpu.Eancestry).abs().maxCoeff() < 2e-4);
    REQUIRE((gpu.EclusterK - cpu.EclusterK).abs().maxCoeff() < 2e-4);
    if(final_iteration) REQUIRE((gpu.GP - cpu.GP).abs().maxCoeff() < 2e-4);
}
} // namespace

TEST_CASE("persistent CUDA joint E-step matches CPU", "[test-cuda]")
{
    std::string reason;
    if(!joint_cuda_available(reason))
    {
        WARN("CUDA test skipped: " << reason);
        return;
    }

    constexpr int K{2}, C{2}, N{2}, M{3};
    Phaseless cpu(K, C, N, M, 7);
    Phaseless gpu(K, C, N, M, 7);
    cpu.pos_chunk = gpu.pos_chunk = {0, M};
    cpu.Q << 0.7, 0.2, 0.3, 0.8;
    cpu.F[0] << 0.8, 0.6, 0.3, 0.2, 0.4, 0.7;
    cpu.F[1] << 0.25, 0.45, 0.75, 0.75, 0.55, 0.25;
    cpu.P << 0.1, 0.8, 0.4, 0.7, 0.2, 0.6;
    cpu.R << 1.0, 0.81, 0.64, 0.0, 0.18, 0.32, 0.0, 0.01, 0.04;
    gpu.Q = cpu.Q;
    gpu.F = cpu.F;
    gpu.P = cpu.P;
    gpu.R = cpu.R;

    MyFloat2D gl(1);
    gl[0] = {0.8, 0.2, 0.1, 0.15, 0.5, 0.2, 0.05, 0.3, 0.7, 0.1, 0.6, 0.2, 0.7, 0.3, 0.1, 0.2, 0.1, 0.7};
    require_cuda_matches_cpu(cpu, gpu, gl);

    // Change every uploaded parameter and run again with the same model/GL
    // objects. This exercises reuse of persistent allocations and cached GL.
    cpu.Q << 0.55, 0.35, 0.45, 0.65;
    cpu.P << 0.2, 0.75, 0.45, 0.6, 0.3, 0.55;
    cpu.F[0] << 0.7, 0.55, 0.35, 0.3, 0.45, 0.65;
    cpu.F[1] << 0.35, 0.5, 0.7, 0.65, 0.5, 0.3;
    cpu.R << 1.0, 0.7225, 0.5625, 0.0, 0.255, 0.375, 0.0, 0.0225, 0.0625;
    gpu.Q = cpu.Q;
    gpu.P = cpu.P;
    gpu.F = cpu.F;
    gpu.R = cpu.R;
    require_cuda_matches_cpu(cpu, gpu, gl);
    require_cuda_matches_cpu(cpu, gpu, gl, true);
}
