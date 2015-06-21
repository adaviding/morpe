#pragma once
#include <memory>
#include <string>

#include "Tensor.h"

namespace Morpe
{
	class MorpeOptions
	{
		public:
			MorpeOptions();
			MorpeOptions(Tensor<float> &paramScale);
			//	If true, the optimization routine will treat each category as equally probable (despite the priors of the training data).
			bool ForceEqualPriors;
			//	The optimization routine will run this many times.  Each time, it will approach the solution a different way.
			int NumberOfApproaches;
			//	default = 2:  Must be greater than 1.  Ensures that parameters actually shrink during optimization.
			float ParamShrinkFactor;
			//	default = 2:  Sets an upper bound on the extent to which a parameter can be scaled.
			float ParamDiffMax;
			// default = 0.01:  Must be greater than 0.  Optimization cannot change parameters less than the following magnitude for any i-th parameter.
			//		[MAGNITUDE] = this.ParamScale.Get(i) * ParamDiffTol
			float ParamDiffTol;
			//	default = 0.005:  Must be greater than 0.  Optimization halts when the change in entropy was less than this amount.
			float EntropyTol;
			//	Returns true if the options make sense, false otherwise.  If false is returned, the error message can be found in this.ErrorMessage.
			bool Check();
			//	The error message set by method this.Check(), if any.
			std::string ErrorMessage;
	};
}