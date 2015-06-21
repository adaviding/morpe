#pragma once
#include <memory>
#include "Tensor.h"

namespace Morpe
{
	class Poly
	{
		//	Constructs empty
		public:		Poly(); 
		//	Constructs a polynomial defined by its rank and the number of spatial dimensions.
		public:		Poly(int nDims, int rank); 
		public:		~Poly();
		//	The number of spatial dimensions over which the polynomial is defined.  (i.e. The number of coordinate axes.)
		public:		int Ndims;
		//	The rank of the polynomial.  (i.e. The maximum power of a term.)
		//
		//	For example...
		//		A linear polynomial would be Rank=1.
		//		A quadratic polynomial would be Rank=2.
		//		A cubic polynomial would be Rank=3.
		public:		int Rank;
		//	The polynomial coefficients.  A 2D array, right-padded with -1 values.  The -1 values signal "empty".  Each row defines a coefficient.
		//
		//	For example...
		//		A row of {0, 0, 2, 5, -1, -1, -1} would correspond to the following polynomial term.
		//			x[0] * x[0] * x[2] * x[5]
		//		A row of {0, 1, 2, -1, -1, -1, -1} would correspond to the following polynomial term.
		//			x[0] * x[1] * x[2]
		public:		std::shared_ptr<Tensor<int>> Coeffs;
		//	Computes the expansion coefficients as a row of Y for each row of X.
		//	-------------------------------------------------------------------------------
		//	OUTPUT
		//	-------------------------------------------------------------------------------
		//	Y	:	2D array
		//		Each row is an expanded feature vector.
		//			The number of rows is pre-allocated to match the number of rows for X.
		//			The number of columns should be equal to this.Coeffs.Size (because the expanded feature vector has this many terms).
		//		The expanded feature vector is computed by applying a polynomial expansion (as defined by the given instance of Poly).
		//	-------------------------------------------------------------------------------
		//	INPUT
		//	-------------------------------------------------------------------------------
		//	X	:	2D array
		//		Each row is a feature vector.
		//			The number of columns should be equal to this.Ndims (because this specifies the spatial dimensionality of a feature space).
		public:		void Expand(_Out_ std::shared_ptr<Tensor<float>> Y, _In_ std::shared_ptr<Tensor<float>> X);
		//	Returns the value for a given row and column of the Pascal matrix.
		//	Output is equal to
		//		Factorial(row+col) / Factorial(row) / Factorial(col)
		public:		static int Pascal(int row, int col);
		//	Returns the number of inhomogeneous polynomial coefficients for a given dimensionality and rank.
		public:		static int Ncoeff(int nDims, int rank);
		//	Returns the number of homogeneous polynomial coefficients for a given dimensionality and rank.
		public:		static int NcoeffAtRank(int nDims, int rank);
	};
}