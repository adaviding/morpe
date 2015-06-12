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
void Mcl_MapDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);