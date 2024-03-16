#include "init_q.cuh"
#include <algorithm>

__global__ void init_h
(
	AssembledSolution    d_assem_sol
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	d_assem_sol.h[idx] = C(0.0);
	d_assem_sol.v[idx] = C(0.0);
}