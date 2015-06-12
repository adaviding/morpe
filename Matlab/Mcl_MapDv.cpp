//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "Linterp.cpp"

//================================================================================================================================
//function Mcl_MapDv(P, Dv, Quant)
//	Warning:  Some data in right-hand-side argument P will be altered by this function.
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function assists Mcl_Exemplar and Mcl_Quadratic methods by calculating the probability of membership in each category for
//	each DecisionValue Dv.  The function utilizes linear interpolation through the tabled valuues defined by TableP and TableDv.  The
//	argument TableDv must be a monotonic increasing function.
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting these comments)
//----------------------------------------------------
// nSamples = The total number of training samples.
// nCats = The number of categories.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not changed by this mex function)
//----------------------------------------------------
// Dv (double 2D array: nSamples * nCats)
//	Inputs the values of the Decision Function associated with each category (columns) and each training sample (rows).  Note that the
//	samples (rows) are made of all training samples (from all categories).
// Quant  (matlab structure vector: nCats)
//	Inputs the tabled values of the classifier's probability of category membership as a function of the decision variable.
// Quant(iCat).Dv  (double vector: Quant(iCat).Nquantiles)
//	Inputs the tabled values of the decision values for each iCat decision function.
// Quant(iCat).PcMonoLim  (double vector: Quant(iCat).Nquantiles)
//	Inputs the tabled values of the probability of category membership for each iCat decision function.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (pre-allocated memory will be filled with the output of this function.)
//----------------------------------------------------
// P (double 2D array: nSamples * nCats)
//	Outputs the probability of membership in each category for each sample.  Each sample is a row, each category is a column.
//================================================================================================================================
void Mcl_MapDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<2 )
		mexErrMsgTxt("Not enough input arguments.");

	//--------------------------------------------------------------------------------------------------------------------------------
	//	INPUT
	//--------------------------------------------------------------------------------------------------------------------------------
	if (!mxIsDouble (prhs[1]) || mxGetNumberOfDimensions(prhs[1])<2)
		mexErrMsgTxt("The input argument Dv must be a double 2D array: ntSamp * nCats.");
	int ntSamp = (int)(mxGetDimensions(prhs[1])[0]);
	int nCats = (int)(mxGetDimensions(prhs[1])[1]);
	if( nCats==1 )
		mexErrMsgTxt("There must be more than 1 category.  In other words, the argument Dv must have more than 1 column.");
	double* Dv = mxGetPr(prhs[1]);
	//--------------------------------------------------------------------------------------------------------------------------------
	if (!mxIsStruct (prhs[2]) || (int)mxGetNumberOfElements(prhs[2])<nCats)
		mexErrMsgTxt("The argument Dv must be a matlab structure with nCats elements.  In other words, it must have the same number of elements as Dv has columns.");
	//--------------------------------------------------------------------------------------------------------------------------------
	//	OUTPUT
	//--------------------------------------------------------------------------------------------------------------------------------
	if (!mxIsDouble(prhs[0]) || (int)mxGetNumberOfDimensions(prhs[0])<2 || (int)(mxGetDimensions(prhs[0])[0])!=ntSamp || (int)(mxGetDimensions(prhs[0])[1])!=nCats)
		mexErrMsgTxt("The input argument P must be a double 2D array the same size as Dv: ntSamp * nCats.");
	double* P = mxGetPr(prhs[0]);
	//--------------------------------------------------------------------------------------------------------------------------------
	
	int iCat;
	int* s0Dv = new int[nCats]; // stride category into sample
	for(iCat=0; iCat<nCats; iCat++)
		s0Dv[iCat] = iCat*ntSamp;

	double pSum;
	int nq, ind;
	double* TableP;
	double* TableDv;
	for( int iCat=0; iCat<nCats; iCat++ )
	{
		nq = ((int*)mxGetData(mxGetField(prhs[2],iCat,"Nquantiles")))[0];
		TableP = mxGetPr(mxGetField(prhs[2],iCat,"PcMonoLim"));
		TableDv = mxGetPr(mxGetField(prhs[2],iCat,"Dv"));
		//mexPrintf("Mcl_MapDv [nq, TableP[0], TableDv[0], TableP[nq-1], TableDv[nq-1]]=[%i, %f, %f, %f, %f]\n", nq, TableP[0], TableDv[0], TableP[nq-1], TableDv[nq-1]);
		for(int iSamp=0; iSamp<ntSamp; iSamp++)
		{
			ind = iSamp+s0Dv[iCat];
			P[ind] = Linterp_Increasing( TableDv, TableP, nq, Dv[ind] );
		}
	}

	for(int iSamp=0; iSamp<ntSamp; iSamp++)
	{
		pSum = 0.0;
		for(int iCat=0; iCat<nCats; iCat++)
			pSum += P[iSamp+s0Dv[iCat]];
		for(int iCat=0; iCat<nCats; iCat++)
			P[iSamp+s0Dv[iCat]] /= pSum;
	}
	
	delete s0Dv;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_MapDv(nlhs, plhs, nrhs, prhs);
}