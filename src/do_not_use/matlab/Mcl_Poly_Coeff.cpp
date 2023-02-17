//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "PolyCoeff.cpp"
#include "mex.h"

#ifndef null
	#define null 0
#endif

//================================================================================================================================
//function [Ncoeff, N, CoeffDims] = Mcl_Poly_Coeff(Rank, Ndims)
//--------------------------------------------------------------------------------------------------------------------------------
// This function outputs the definition of a polynomial defined by Rank and Ndims.  This definition can be used to reprsent and calculate
//	a polynomial.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT
//--------------------------------------------------------------------------------------------------------------------------------
// Rank (int32 scalar)
//	The polynomial's rank.
//		If Rank=1, the polynomial will only contain first-order  (linear) terms.
//		If Rank=2, the polynomial will also contain second-order (quadratic) terms (in addition to the terms of the lower rank).
//		If Rank=3, ...								third-order  (cubic)
//		etc...
// Ndims (int32 scalar)
//	The dimensionality of the space over which the polynomial is defined.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT
//--------------------------------------------------------------------------------------------------------------------------------
// Ncoeff (int32 scalar)
//	The total number of coefficients of the polynomial.
// N (int32 vector: Rank)
//	Each N[iRank] outputs the number of coefficients for each iRank rank.
// CoeffDims (int32 2D array, Ncoeff * Rank) indexed as:
//			C[iCoeff + Ncoeff*iComp]
//		where
//			iComp == 0 first  component of the polynomial expansion
//			iComp == 1 second component of the polynomial expansion
//			iComp == 2 third  component of the polynomial expansion
//			etc...
//		where
//			iCoeff is bound to the range:  [0, Ncoeff(null,Rank,Ndims)]
//		This array enumerates the components of polynomial expansion for each coefficient.  Each (iCoeff,iComp)-th element of CoeffDims
//		can be an integer {0,...,Ndims-1} which identifies an axis of an Ndims dimensional space; OR can be -1 when the entry is "empty"
//		(i.e. when no axis is indicated).   All -1 ("empty") values are packed to right-most columns so the left-most columns will
//		always contain the significant entries (integers {0,...,Ndims-1}).  The first column (associated with iComp==0) will never
//		contain empty (-1) values.
//================================================================================================================================
void Mcl_Poly_Coeff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Check number of argument(input and output)
	if (nrhs<2)
		mexErrMsgTxt("Not enough input arguments.");
	if (nlhs>3)
		mexErrMsgTxt("Too many output arguments.");
	if (nlhs<1)
		mexErrMsgTxt("Not enough output arguments.");

	//--------------------------------------------------------------------------------------------------------------------------------
	//	INPUT
	//--------------------------------------------------------------------------------
	int Rank = 0;
	if( !mxIsInt32(prhs[0]) || (int)mxGetNumberOfElements(prhs[0])!=1 )
		mexErrMsgTxt("The first argument must be a scalar of type int32.");
	Rank = ((int*)mxGetData(prhs[0]))[0];
	//--------------------------------------------------------------------------------
	int Ndims = 0;
	if( !mxIsInt32(prhs[0]) || (int)mxGetNumberOfElements(prhs[0])!=1 )
		mexErrMsgTxt("The first argument must be a scalar of type int32.");
	Ndims = ((int*)mxGetData(prhs[1]))[0];


	//--------------------------------------------------------------------------------------------------------------------------------
	//	OUTPUT (allocate first two output arguments)
	//--------------------------------------------------------------------------------
	plhs[0] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	int *Nco = (int*)mxGetData(plhs[0]);
	//--------------------------------------------------------------------------------
	int *N = null;
	if (nlhs>=2)
	{
		plhs[1] = mxCreateNumericMatrix(Rank,1,mxINT32_CLASS,mxREAL);
		N = (int*)mxGetData(plhs[1]);
	}
	//--------------------------------------------------------------------------------------------------------------------------------

	//	Determine first two arugmenets of output.
	Nco[0] = Ncoeff(N,Rank,Ndims);

	//--------------------------------------------------------------------------------------------------------------------------------
	//	OUTPUT (allocate third argument)
	//--------------------------------------------------------------------------------
	int* CoeffDims = null;
	if (nlhs>=3)
	{
		plhs[2] = mxCreateNumericMatrix(Nco[0],Rank,mxINT32_CLASS,mxREAL);
		CoeffDims = (int*)mxGetData(plhs[2]);
	}
	//--------------------------------------------------------------------------------------------------------------------------------

	//	Determine third argument
	if( CoeffDims!=null )
		EnumCoeffs(CoeffDims,Nco[0],Rank,Ndims);
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Poly_Coeff(nlhs, plhs, nrhs, prhs);
}