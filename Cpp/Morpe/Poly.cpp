#include "Stdafx.h"

namespace Morpe
{
	Poly::Poly()
	{
		this->Ndims = 0;
		this->Rank = 0;
		this->Coeffs = std::shared_ptr<Tensor<int>>(new Tensor<int>());
	}
	Poly::Poly(int nDims, int rank)
	{
		int i,j;
		this->Ndims = nDims;
		this->Rank = rank;

		//	How many polynomial coefficients total?
		int nCoeff = Poly::Ncoeff(nDims,rank);

		//	Allocate coefficients.
		this->Coeffs = std::shared_ptr<Tensor<int>>(new Tensor<int>(nCoeff,nDims));

		//	One coefficient (or row) at a time.
		std::vector<int> row = std::vector<int>(nDims,-1);
	
		for(int iCoeff=0; iCoeff < nCoeff; iCoeff++)
		{
			//	Increment the smallest digit and begin to handle "carry-over" arithmetic.
			i=-1;
			while( ++row[++i]==nDims );	//	When a digit is maxed, keep incrementing rightward until we're not maxed out.

			//	Finish the "carry-over" by ensuring that any leftward maxed digits have been reset.
			for( j=i-1; j>=0; j-- )
			{
				if( row[j]==nDims )
					row[j]=row[j+1];
			}

			//	Copy last row to output
			for(j=0; j<nDims; j++)
				this->Coeffs->Set(row[j],iCoeff,j);
		}
	}

	int Poly::Ncoeff(int nDims, int rank)
	{
		return Pascal(nDims,rank);

		//				0	1	2	3	4	5	6	7	8	9
		//							Rank of Tensor
		//	0			1	1	1	1	1	1	1	1	1	1
		//	1			1	2	3	4	5	6	7	8	9	10
		//	2	Ndims	1	3	6	10	15	21	28	36	45	55
		//	3	of		1	4	10	20	35	56	84	120	165	220
		//	4	Space	1	5	15	35	70	126	210	330	495	715
		//	5			1	6	21	56	126	252	462	924	1716
	}

	int Poly::NcoeffAtRank(int nDims, int rank)
	{
		return Pascal(nDims-1,rank);
	}


	int Poly::Pascal(int row, int col)
	{
		if( row<=0 )
			return 0;
		if( col<=0 )
			return 0;
		int output = 1;
		int i;
		for(i=(row+col); i>row; i--)
			output *= i;
		for(i=2; i<=col; i++)
			output /= i;
		return output;
	}
}