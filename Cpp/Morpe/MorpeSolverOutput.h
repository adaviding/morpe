#pragma once
#include <ctime>
#include <memory>
#include "Poly.h"

namespace Morpe
{
	class MorpeSolverOutput
	{
		public:
			MorpeSolverOutput();
			~MorpeSolverOutput();
			//	The best parameter values.  The parameter values that yielded minimum entropy.
			std::shared_ptr<Tensor<float>> Params;
			//	The best (minimum) entropy found by the solver.
			float Entropy;
			//	The number of times that an objective function was evaluated.
			int Nevals;
			//	Optimization began at this datetime (UTC).
			tm UtcStart;
			//	Optimization finished at this datetime (UTC).
			tm UtcStop;
	};
}