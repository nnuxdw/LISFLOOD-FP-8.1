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
#include "NonUniformInterfaces.h"
#include "GhostCellTypes.h"
#include "compact.cuh"
#include "Boundaries.h"
#include "PointSources.h"

__global__ void init_h
(
	AssembledSolution    d_assem_sol
);