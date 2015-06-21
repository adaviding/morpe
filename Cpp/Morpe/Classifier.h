#pragma once
#include <vector>
#include "Tensor.h"
#include "Poly.h"
#include "MorpeOptions.h"
#include "MorpeQuant.h"
#include "MorpeSolverOutput.h"

namespace Morpe
{
	using namespace System;

	class Classifier
	{
		public:		~Classifier();

		//	Prepares to train a new classifier with the following properties.
		//		ncats:	The number of categories.
		//		ndims:	The dimensionality of the feature space.
		//		rank:	The rank of the polynomial expansion.
		public:		Classifier(int ncats, int ndims, int rank);

		//	The number of categories in the training set.  Must be greater than 1.
		public:		int Ncats;

		//	The number of polynomials.
		//		If Ncats==2, then Npoly==1.
		//		If Ncats>2, then Npoly==Ncats.
		public:		int Npoly;

		//	The tableu of polynomial coefficients.
		public:		Poly Coeffs;

		//	The 2D tensor of Morpe parameters.  [#rows=Npoly, #cols=Coeffs.Size()].
		public:		Tensor<float> Params;

		//	The category weights.  These are specified by equation 15 in the paper by Ing.
		public:		std::vector<float> CatWeights;

		//	The quantization data for each polynomial.
		public:		std::vector<MorpeQuant> Quants;
	};
}