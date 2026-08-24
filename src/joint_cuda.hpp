#ifndef JOINT_CUDA_H_
#define JOINT_CUDA_H_

#include "common.hpp"

class Phaseless;

// Runs one complete joint-model E step.  The caller must call initIteration()
// first and updateIteration() afterwards, exactly as for the CPU path.
bool joint_cuda_available(std::string & reason);
double joint_cuda_e_step(Phaseless & model, const MyFloat2D & gl, bool final_iteration);

#endif
