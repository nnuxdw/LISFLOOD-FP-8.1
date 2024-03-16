#pragma once
#include "../../lisflood.h"
#include "../cuda_simulate.cuh"

namespace lis
{
namespace cuda
{
namespace fv1
{

class Simulation : public cuda::Simulation
{
public:
	void run
	(
		Fnames& filenames,
		States& states,
		Pars& pars,
		::Solver& solver,
		int verbose
	);
};

}
}
}
