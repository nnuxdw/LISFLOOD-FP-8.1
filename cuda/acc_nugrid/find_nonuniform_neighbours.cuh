#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "AssembledSolution.h"
#include "Directions.h"
#include "MortonCode.h"
#include "get_lvl_idx.cuh"
#include "get_level.cuh"
#include "NonUniformNeighbours.h"
#include "GhostCellTypes.h"
#include "compact.cuh"
#include "Boundaries.h"

__global__ void find_nonuniform_neighbours
(
	AssembledSolution    d_assem_sol,
	AssembledSolution    d_nghbr_assem_sol,
	Pars pars,
	Solver solver,
	NonUniformNeighbours d_non_uniform_nghbrs,
	Boundaries           boundaries
);