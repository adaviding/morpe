//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "QuickSort.cpp"
#include "Accuracy.cpp"

#ifndef null
	#define null 0
#endif

static int *Poly_Init_Idx=null;
static int nPoly_Init_Idx=0;
static double *Poly_Init_AccMin=null;
static int nPoly_Init_AccMin=0;

// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void Poly_Init_DeleteStaticMemory(void)
{
	delete Poly_Init_Idx;
	delete Poly_Init_AccMin;
}

//================================================================================================================================
//function Mcl_Poly_Init(o, Xcell)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function initializes a new Mcl_Poly classifier.  It is called from the constructor Mcl_Poly_Ctor and basically
//	serves as a helper for the constructor.
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting subsequent comments)
//----------------------------------------------------
// Ncats = The number of categories (i.e. the length of the input X cell array)
// Nsamp(iCat) = The number of training samples in each category (i.e. the number of rows of the input X{iCat} summed across all iCat).
// Ntsamp = The total number of training samples = sum(Nsamp);
// Ndims = The spatial dimension of the data (i.e. the number of columns for each input X{iCat}, this is the same for each iCat).
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not altered by mex function)
//----------------------------------------------------
// Xcell (cell vector, Ncats)
//  Each iCat element of this cell array, Xcell{iCat}, is a 2D array (Ncats,Ndims) contains samples of training data from category iCat.
//
// [Xtrans] (matlab structure Mcl_Trfm vector, Ndims)
//	Each spatial dimension can be associated with a transform.  The contents of Xcell are transformed accordingly before the classifier
//		is initialized.  You can also leave this argument empty [] and no transforms will be applied.
//
// Rank (int32 scalar)
//	A positive integer identifying the polynomial rank.  Flat=0, Linear=1, Quadratic=2, Cubic=3, ...
//
// ForceEqualPriors (bool scalar) Suggest true.
//	If true, then the categories are assumed to have equal prior probabilities, even if there are different numbers of samples for
//	each category in the training data.
//
// Nquantiles (int32 scalar) Suggest any integer from 10 to 50.
//	The number of quantiles to be used in estimating the probability of category membership based on the decision variable.  In general,
//	a good value for Nquantiles is 10*Ncats or 20*Ncats;
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT is a newly allocated structure filled with the following data
//----------------------------------------------------
// Ncats (int32 scalar)
//	The number of categories.
// Ndims (int32 scalar)
//	The number of spatial dimensions (i.e. the number of columns of every iCat element of input Xcell, Xcell{iCat}).
// Ntsamp  (int32 vector, Ncats)
//	The total number of training samples (i.e. the total  of rows for all iCat elements of input Xcell, Xcell{iCat}).
// Nsamp  (int32 vector, Ncats)
//	The number of samples in each category (i.e. the number of rows for each iCat element of input Xcell, Xcell{iCat}).
// Npoly (int32 scalar)
//	The number of polynomial decision functions.
//		if Ncats==2
//			Npoly=1;
//		else if Ncats>2
//			Npoly=Ncats;
//	When there are two categories, only one decision function (i.e. one polynomial) is needed because the decision function for the second
//		category is merely the negative of the decision function for the first category.
//	When there are more than two categories, each category needs its own distinct decision function.
// Rank (int32 scalar)
//	Same as the input Rank.
// Ncoeff (int32 scalar)
//	The total number of polynomial coefficients, also equal to sum(NcoeffByRank).
// NcoeffByRank (int32 vector)
//	The number of polynomial coefficients at each rank.
// CoeffDims (int32 2D array, Ncoeff * Rank)
//	The components of polynomial expansion (by columns) for each coefficient (row).  Each (iCoeff,iComp)-th element of CoeffDims
//		can be an integer {0,...,Ndims-1} which identifies an axis of an Ndims dimensional space; OR can be -1 when the entry is "empty"
//		(i.e. when no axis is indicated).   All -1 ("empty") values are packed to right-most columns so the left-most columns will
//		always contain the significant entries (integers {0,...,Ndims-1}).  The first column (associated with iComp==0) will never
//		contain empty (-1) values.
// ForceEqualPriors (bool scalar)
//	Same as input ForceEqualPriors.
// CatWeight (double vector, Ncats)
//	If the training sample has equal numbers of samples from each category, then this will be a vector of ones.  Otherwise, if ForceEqualPriors
//		is true, then these weights will scale the influence of training samples so that the total weight is Ntsamp and that each category
//		has equal influence on the solution.
// Cat (int32 vector, Ntsamp)
//	The 0-based category label of each sample after samples from all categories are concatenated into a single list.
// X (double 2D array, Ntsamp * Ncoeff)
//	The unscaled polynomial expansions of each sample.  When the polynomial function (a.k.a. decision function) is evaluated, each column
//		becomes scaled by a polynomial coefficient.  Each row is a sample with its category label denoted by Cat.
// Xmeans (double 2D array, Ncats * Ncoeff)
//	The mean of the columns of X computed separately for each category.
// Xvars (double 2D array, Ncats * Ncoeff)
//	The variance of the columns of X computed separately for each category.
// Xmean (double vector, Ncoeff)
//	The pooled mean (where CatWeight was used to weight each category's influence on the result) of the columns of X.  This is literally the
//		weighted some of Xmeans across categories.  Weights are 1./CatWeight.
// Xvar (double vector, Ncoeff)
//	The pooled variance (where CatWeight was used to weight each category's influence on the result) of the columns of X.  Variance of category
//		means is not included, so this is literally the weighted some of Xvars across categories.  Weights are 1./CatWeight.
// Xscale (double vector, Ncoeff)
//	The scale of data in each column of X (i.e. each dimension of the polynomial expansion).  This is equal to the weighted average of the
//		within-category difference between the 0.25 and 0.75 quantile.  Weights are 1./CatWeight.
// Dv (double 2D array, Ntsamp * Ncats):  This variable is calculated in Mcl_Poly_Train.
//	Lists the decision variable (i.e. the value of the polynomial) for each training sample and each category's polynomial decision function.
//		In general, there are Ncats polynomials except when Ncats==2, in which case there is only one polynomial.
// P (double 2D array, Ntsamp * Ncats):  This variable is calculated in Mcl_Poly_Train.
//	This lists the probability that each category should be assigned to each training sample.  Probabilities are computed by linear interpolation
//		of Dv through the quantization table [Quant.Dv, Quant.PcMonoLim] and then normalizing each row to sum to 1.0.
// H (double vector, Ntsamp):  This variable is calculated in Mcl_Poly_Train.
//	This lists the conditional entropy of each training sample.  H[iSamp] appraches 0.0 when the iSamp sample is categorized perfectly.  Chance perforamnce
//		is associated with H[iSamp] of
//			if ForceEqualPriors, CatWeight[Cat[iSamp]]/Ntsamp; 
//			else, 1.0/Ntsamp
// wInit (double vector, Ndims)
//	The initial value of the free parameters.
// h (double scalar):  This variable is calculated in Mcl_Poly_Train.
//	The conditional entropy of the entire data set in the range [Quant.hMin, 1.0].  This value does not approach 0.0 because the probabilities of
//		category membership are limited by Quant.PcMonoLim, and these probabilities are never allowed to equal 0.0 or 1.0.  Instead they are limited
//		to the range [Quant.pLow, Quant.pHigh].
// hLim (double scalar):  This variable is calculated in Mcl_Poly_Train.
//	The conditional entropy mapped uniformly from the range [Quant.hMin, 1.0] to the range [0.0, 1.0].  This is like "stretching" the range.
// Acc (double scalar):  This variable is calculated in Mcl_Poly_Train.
//	The weighted accuracy.  1.0 is perfect accuracy, 1/Ncats is chance.  Weights are specified in CatWeight.
// Quant (mxArray matlab structure, Ncats):  This variable is calculated in Mcl_Poly_Train.
//	This structure specifies the mapping between the decision function for category iCat and the probability of category membership.  It is
//		calculated based on quantization of each category's decision variable and training sample.
// Quant(iCat).Nquantiles (int32 scalar):  This variable is calculated in Mcl_Poly_Train.
//	The number of quantiles (i.e. entries of the interpolation table).
// Quant(iCat).Dv (double vector, Nquantiles)
//	The mean decision value for each quantile (or row of interpolation table).
// Quant(iCat).DvBinSep (double vector, Nquantiles)
//	The boundary of each quantile.  The first value is -Inf and last value is +Inf.
// Quant(iCat).Weight
//	The total weight of training samples accumulated in each quantile (or row of the interpolation table).  See CatWeight for explanation of weighting.
// Quant(iCat).Pc
//	The observed proportion of correct assignments (for the training data) for each quantile (or row of interpolation table).
// Quant(iCat).PcMono
//	A non-decreasing function fit to Pc.
// Quant(iCat).PcMonoLim
//	The same function as PcMono, but probabilities have been limited not to equal 0.0 or 1.0, instead they are limited to [pLow, pHigh].
// Quant(iCat).pLow (double scalar)
//	The lower limit of probability.  This approaches 0.0 as average Weight approaches infinity and is always close to 0.0.
// Quant(iCat).pHigh (double scalar)
//	The upper limit of probability.  This approaches 1.0 as average Weight size approaches infinity and is always close to 1.0.
// Quant(iCat).hMin (double scalar)
//	The minimum possible entropy based on pHigh.  This is always positive but approaches 0.0 as average Weight approaches infinity.
// SolverOptions (mxArray, Matlab structure):  This variable is empty.
//	This is just a placeholder for storing the structure that configures the optimization algorithm.
// SolverOutput (mxArray, Matlab structure):  This variable is calculated in Mcl_Poly_Train.
//	After Mcl_Poly_Train is executed, this variable contains raw output from the optimization algorithm.
// Mem (mxArray, Matlab structure)
//	This is just memory allocated to make Mcl_Poly_Train run as efficiently as possible.
//================================================================================================================================
void Mcl_Poly_Init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);