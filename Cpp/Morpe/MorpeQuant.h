#pragma once
#include <vector>

namespace Morpe
{
	class MorpeQuant
	{
		public:
			MorpeQuant();
			std::vector<float> Z;
			std::vector<float> P;
			std::vector<float> Zsep;
			float Pmin;
			float Pmax;
	};
}