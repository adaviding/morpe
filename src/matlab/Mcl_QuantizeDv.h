//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "memory.h"
#include "QuickSort.cpp"

#ifndef null
	#define null 0
#endif
//#ifndef DEBUG_MEXPRINTF
//	#define DEBUG_MEXPRINTF
//#endif

static int* QuantizeDv_Idx=null;
static int nQuantizeDv_Idx=0;
// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void QuantizeDv_DeleteStaticMemory(void)
{
	delete QuantizeDv_Idx;
}

//================================================================================================================================
//function Mcl_QuantizeDv(o)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function assists Mcl_Exemplar and Mcl_Poly methods by quantizing decision values.
//  A "decision value" is a number that can be computed for each training sample using a decision function associated with a particular
//	category.  Each category must be associated with a unique decision function, and each of these decision function is assumed (but not
//	required) to have an isotonic relationship with the posterior probability of memberhsip in the category.
//
//	This procedure makes it possible to determine the posterior probability of category membership given the evaluated decision function
//	using a unique monotonic mapping function for each category.  The output of this procedure provides an interpolation table that
//	defines the mapping of each category's decision value with the probability of membership in that category.
//
//  If you have only one decision function for two categories, then you can always mimic a decision function for the second category by
//	transforming the output of the decision function for the first category.  This transform can usually (but now always) be accomplished
//	by taking the negative, inverse, or 1-complement of the first decision function.  For example, a linear model can be defined by the
//	following linear decision function where "x" is the vector coordinate of a training sample, "B" is a vector defining the slope of a
//	linear decision boundary, and "B0" is a scalar that defines the boundary's intercept.
//		f( category "1" , x ) =   B * x + B0
//	In this example, the linear decision function for the second category is the negative of the same function for the first category.
//		f( category "2" , x ) = -(B * x + B0)
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting these comments)
//----------------------------------------------------
// Ntsamp = The total number of training samples.
// Ncats = The number of categories.
// Nquant = The number of quantiles (bins). 
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not changed by this mex function)
//----------------------------------------------------
// o (mxArray*, Matlab structure)
//	A structure constructed from Mcl_Exemplar_Ctor (refer to its comments) and then run at least once through Mcl_Exemplar_Eval.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (pre-allocated memory will be filled with the output of this function.)
//----------------------------------------------------
// o.Quant(iCat) (mxArray*, Matlab structure)
//	The quantization table for each iCat category is mostly completed by this function.  However, o.PcMono and o.PcMonoLim are not
//		filled by this function.
//--------------------------------------------------------------------------------------------------------------------------------
// After this function terminates, a decision value can be mapped to a posterior probability of assigning the iCat category label by
//	using o.Quant(iCat).Dv and o.Quant(iCat).Pc as an interpolation table.  Since probabilities are computed separately for each
//	category (using separate decision functions), there is no guarantee that the probabilities will sum to 1.0 when a stimulus is
//	classified.  Therefore, subsequent methods should always normalize the probabilities (making them sum to 1.0).
//
// Since posterior probabilities should never reach the boundaries [0.0, 1.0], it is recommended that the user set boundaries for
//	o.Quant(iCat).Pc(iBin) for each iBin quantile to be:
//		[1/Ncats/o.Quant.Weight(iBin), 1.0-(Ncats-1)/Ncats/o.Quant.Weight(iBin)]
//	This algorithm does not do this automatically, therefore o.Quant(iCat).Pc(iBin) can reach the boundaries of [0.0, 1.0].
//
// Since quantization was used to compute the interpolation table, there will be quantization (binning) noise in o.Quant(iCat).Dv
//	and o.Quant(iCat).Pc.
//================================================================================================================================
void Mcl_QuantizeDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);