//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "QuickSort.cpp"

#ifndef null
	#define null 0
#endif

//================================================================================================================================
//function Mcl_QuickSort(Ind, Xvec)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function returns integer indexes for ascending Xvec.
// Usage:
//		Xvec = rand(100,1); Ind = zeros(length(Xvec),1,'int32'); Mcl_QuickSort(Ind,Xvec); disp(diff(Xvec(Ind)));
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting subsequent comments)
//----------------------------------------------------
// Nx = The number of values to be sorted, numel(Xvec) and numel(Ind).
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not altered by mex function)
//----------------------------------------------------
// Xvec (double vector: Nx)
//  The vector to be sorted.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (values altered by mex function)
//----------------------------------------------------
// Ind (int32 vector: Nx)
//	The indexes into Xvec that gives the ascending order.
//================================================================================================================================
void Mcl_QuickSort(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<2 )
		mexErrMsgTxt("Not enough input arguments.");

	if (!mxIsInt32(prhs[0]))
		mexErrMsgTxt("The input argument Ind must be int32.");
	if (!mxIsDouble(prhs[1]))
		mexErrMsgTxt("The input argument Xvec must be double.");

	int Nx = (int)mxGetNumberOfElements(prhs[0]);
	if( (int)mxGetNumberOfElements(prhs[1])!=Nx )
		mexErrMsgTxt("The input arguments must be vectors of the same length.");

	int* Ind = (int*)mxGetData(prhs[0]);
	double* Xvec = (double*)mxGetData(prhs[1]);

	//	Initialize Idx
	for( int ix=0; ix<Nx; )
		Ind[ix] = ix++;

	//	Call QuickSort
	QuickSort(Ind, Xvec, 0, Nx-1);

	//	Change to 1-based
	for( int ix=0; ix<Nx; )
		Ind[ix++]++;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_QuickSort(nlhs, plhs, nrhs, prhs);
}