//	Almon David Ing
//	Ctr. Perceptual Systems
//	University of Texas at Austin
//	Compiled July 6, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"

//================================================================================================================================
//function Mcl_RangeLimit(outVec, inVec, outMin, outMax)
//	Warning:  Some data in right-hand-side argument outVec will be altered by this function.
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function copies the input vector inVec to the output vector outVec.  The output is limited to the range [outMin, outMax].
//================================================================================================================================
void Mcl_RangeLimit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<4 )
		mexErrMsgTxt("Not enough input arguments.");
	if (!mxIsDouble (prhs[0]))
		mexErrMsgTxt("The argument outVec must be type double.");
	if (!mxIsDouble (prhs[1]))
		mexErrMsgTxt("The argument inVec must be type double.");
	if (!mxIsDouble (prhs[2]))
		mexErrMsgTxt("The argument outMin must be type double.");
	if (!mxIsDouble (prhs[3]))
		mexErrMsgTxt("The argument outMax must be type double.");
	int nTable = (int)(mxGetNumberOfElements(prhs[0]));
	if( (int)(mxGetNumberOfElements(prhs[1]))<nTable )
		mexErrMsgTxt("The argument outVec must be the same length as the argument inVec.");
	double outMin = mxGetScalar(prhs[2]);
	double outMax = mxGetScalar(prhs[3]);
	if( outMax<outMin )
		mexErrMsgTxt("The argument outMax must be greater than the argument outMin.");

	double* outVec = mxGetPr(prhs[0]);
	double* inVec = mxGetPr(prhs[1]);
	double val;

	for(int i=0; i<nTable; i++)
	{
		val = inVec[i];
		if( val<outMin )
			outVec[i]=outMin;
		else if( val>outMax )
			outVec[i]=outMax;
		else
			outVec[i] = val;
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_RangeLimit(nlhs, plhs, nrhs, prhs);
}