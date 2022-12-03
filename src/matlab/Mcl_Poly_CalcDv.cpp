//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#ifndef null
	#define null 0
#endif

//================================================================================================================================
//function Mcl_Poly_CalcDv(o, Dv, X, wOptimized)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function uses a trained Mcl_Poly solver to classify new (untrained) data.  This function reduces each multivariate
//	stimulus coordinate to Ncats decision values.  Use the Mcl_MapDv function to map these values to probabilities.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not altered by mex function)
//----------------------------------------------------
// o 
//	The polynomial model specified by Mcl_Poly_Ctor.
// X (double 2D array: ntSamp * o.Ncoeff) optional parameter, defaults to o.X
//	The expanded spatial coordinate of each sample.  The polynomial expansion must be already applied.  If the solver internalizes
//	transforms on each dimension (via o.Xtrans), X must have been transformed prior to application of the polynomial expansion.
// wOptimized (double vector: o.Npoly * o.Ncoeff) optional parameter, defaults to o.wOptimized
//	If provided, this vector provides the optimized coefficients of the multivariate polynomial.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (values are copied to pre-allocated arrays provided).
//----------------------------------------------------
// Dv (double 2D array: nTest * Ncats)
//	For each test coordinate, the decision variable associated with each category is provided.
//================================================================================================================================
void Mcl_Poly_CalcDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<1 )
		mexErrMsgTxt("Not enough input arguments.");
	
	// mexPrintf("Mcl_Poly_CalcDv: Entry\n");
	//-----------------------------------------------------------------------------------
	// INPUT from structure and optional parameters
	//------------------------------------------------
	if (!mxIsStruct(prhs[0]))
		mexErrMsgTxt("The input argument o must be a matlab structure.");
	mxArray *fld;
	//------------------------------------------------
	fld = mxGetField(prhs[0],0,"Npoly");
	if( fld==null )
		mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name was not found.");
	int Npoly = ((int*)mxGetData(fld))[0];
	//------------------------------------------------
	fld = mxGetField(prhs[0],0,"Ncats");
	if( fld==null )
		mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name was not found.");
	int Ncats = ((int*)mxGetData(fld))[0];
	//------------------------------------------------
	fld = mxGetField(prhs[0],0,"Ncoeff");
	if( fld==null )
		mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name was not found.");
	int Ncoeff = ((int*)mxGetData(fld))[0];
	//------------------------------------------------
	fld = mxGetField(prhs[0],0,"ForceEqualPriors");
	if( fld==null )
		mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name was not found.");
	bool ForceEqualPriors = ((bool*)mxGetData(fld))[0];
	//------------------------------------------------
	double* X = null;
	int Ntsamp = 0;
	if( nrhs>=3 && prhs[2]!=null && (int)mxGetNumberOfElements(prhs[2])>0 )
	{
		X = (double*)mxGetData(prhs[2]);
		//------------------------------------------------
		Ntsamp = (int)mxGetM(prhs[2]);
		if( (int)mxGetN(prhs[2])<Ncoeff )
			mexErrMsgTxt ("The input argument X must have at least o.Ncoeff columns.  Did you apply the polynomial expansion?");
	}
	else
	{
		fld = mxGetField(prhs[0],0,"X");
		if( fld==null )
			mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name (X) was not found.");
		X = (double*)mxGetData(fld);
		//------------------------------------------------
		fld = mxGetField(prhs[0],0,"Ntsamp");
		if( fld==null )
			mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name (Ntsamp) was not found.");
		Ntsamp = ((int*)mxGetData(fld))[0];
	}
	//------------------------------------------------
	double* wParams = null;
	if( nrhs>=4 && prhs[3]!=null && (int)mxGetNumberOfElements(prhs[3])>0 )
	{
		wParams = (double*)mxGetData(prhs[3]);
		if( (int)mxGetNumberOfElements(prhs[3])<(Npoly*Ncoeff) )
			mexErrMsgTxt("The argument wOptimimzed must have at least o.Npoly*o.Ncoeff elements.");
	}
	else
	{
		fld = mxGetField(prhs[0],0,"wOptimized");
		if( fld==null )
			mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name (wOptimized) was not found.");
		wParams = (double*)mxGetData(fld);
	}
	//-----------------------------------------------------------------------------------
	
	// mexPrintf("Mcl_Poly_CalcDv: Input Completed\n");

	//-----------------------------------------------------------------------------------
	// OUTPUT to Dv
	//------------------------------------------------
	double* Dv = null;
	if( nrhs>=2 && prhs[1]!=null &&(int)mxGetNumberOfElements(prhs[1])>0 )
	{
		if( (int)mxGetNumberOfElements(prhs[1])<(Ntsamp*Ncats) )
			mexErrMsgTxt("The output argument Dv must be pre-allocated to a 2D double array: NtSamp * Ncats.");
		if (!mxIsDouble(prhs[1]))
			mexErrMsgTxt("The output argument Dv must be type double.");
		Dv = (double*)mxGetData(prhs[1]);
	}
	else
	{
		fld = mxGetField(prhs[0],0,"Dv");
		if( fld==null )
			mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Poly structure.  A required field name was not found.");
		Dv = (double*)mxGetData(fld);
	}
	//-----------------------------------------------------------------------------------

	// mexPrintf("Mcl_Poly_CalcDv: Output Completed\n");

	double dv = 0.0;
	int iSamp, iPoly, iCoeff;
	//-----------------------------------------------------------------------------------
	// Stride
	//------------------------------------------------
	int nx = Ncats>Ncoeff ? Ncats : Ncoeff;
	int* s0X = new int[nx];
	for(iCoeff=0; iCoeff<nx; iCoeff++)
		s0X[iCoeff] = Ntsamp * iCoeff;
	int* s0Poly = new int[Npoly];
	for(iPoly=0; iPoly<Npoly; iPoly++)
		s0Poly[iPoly] = Ncoeff*iPoly;
	double* wPoly = null;
	//-----------------------------------------------------------------------------------

	// mexPrintf("Mcl_Poly_CalcDv: Stride Completed\n");

	//-----------------------------------------------------------------------------------
	// Calculate the decision values for each sample and polynomial.
	//------------------------------------------------
	//	For each polynomial
	for(iPoly=0; iPoly<Npoly; iPoly++)
	{
		//	Point to the parameters for this polynomial.
		wPoly = wParams + s0Poly[iPoly];

		//	For each sample
		for(iSamp=0; iSamp<Ntsamp; iSamp++)
		{
			//	Initialize sum of the decision value.
			dv = 0.0;
			//	For each coefficient of the polynomial, accumulate a total for the decision value.
			for(iCoeff=0; iCoeff<Ncoeff; iCoeff++)
				dv += wPoly[iCoeff]*X[iSamp+s0X[iCoeff]];
			//	Store decision value.
			Dv[iSamp+s0X[iPoly]] = dv;
		}
	}
	//	If there is only one polynomial, then a complementary Dv (the negative reflection) is filled into second column.
	if( Npoly==1 )
	{
		for(iSamp=0; iSamp<Ntsamp; iSamp++)
			Dv[iSamp+s0X[1]] = -Dv[iSamp];
	}
	//-----------------------------------------------------------------------------------

	// mexPrintf("Mcl_Poly_CalcDv: Exit\n");

	delete s0X;
	delete s0Poly;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Poly_CalcDv(nlhs, plhs, nrhs, prhs);
}