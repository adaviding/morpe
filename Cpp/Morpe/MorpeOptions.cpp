#include "Stdafx.h"

namespace Morpe
{
	MorpeOptions::MorpeOptions()
	{
		this->ForceEqualPriors = true;
		this->NumberOfApproaches = 1;
		this->ParamShrinkFactor = 2.0f;
		this->ParamDiffMax = 2.0f;
		this->ParamDiffTol = 0.01f;
		this->EntropyTol = 0.005f;
	}
	bool MorpeOptions::Check()
	{
		return false;
	}
}