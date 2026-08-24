#include "joint_cuda.hpp"
#include <stdexcept>

bool joint_cuda_available(std::string & reason)
{
    reason = "this binary was built without CUDA support (rebuild with `make CUDA=1`)";
    return false;
}

double joint_cuda_e_step(Phaseless &, const MyFloat2D &, bool)
{
    throw std::runtime_error("CUDA joint-model support is not compiled into this binary");
}
