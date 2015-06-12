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
void Mcl_Poly_CalcDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);