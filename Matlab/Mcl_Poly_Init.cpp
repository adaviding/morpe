//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "QuickSort.cpp"
#include "Accuracy.cpp"

//#include "PolyCoeff.cpp"

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
//	The initial value of the free parameters for the blurring or generalization function.  These parameters control the width of the blurring
//		kernel for each spaital dimension.
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
//	The probability that the correct category will be assigned (for the training data) for each quantile (or row of interpolation table).
// Quant(iCat).PcMono
//	A non-decreasing function fit to Pc.
// Quant(iCat).PcMonoLim
//	The same function as PcMono, but probabilities have been limited not to equal 0.0 or 1.0, instead they are limited to [pLow, pHigh].
// Quant(iCat).pLow (double scalar)
//	The lower limit of probability.  This approaches 0.0 as average Weight approaches infinity, otherwise is always close to 0.0.
// Quant(iCat).pHigh (double scalar)
//	The upper limit of probability.  This approaches 1.0 as average Weight size approaches infinity, otherwise is always close to 1.0.
// Quant(iCat).hMin (double scalar)
//	The minimum possible entropy based on pHigh.  This is always positive but approaches 0.0 as average Weight approaches infinity.
// SolverOptions (mxArray, Matlab structure):  This variable is empty.
//	This is just a placeholder for storing the structure that configures the optimization algorithm.
// SolverOutput (mxArray, Matlab structure):  This variable is calculated in Mcl_Poly_Train.
//	After Mcl_Poly_Train is executed, this variable contains raw output from the optimization algorithm.
// Mem (mxArray, Matlab structure)
//	This is just memory allocated to make Mcl_Poly_Train run as efficiently as possible.
//================================================================================================================================
void Mcl_Poly_Init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int iDim, iCat, jCat, i;
	mxArray* arg;
	mxArray* subArg;
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<2 )
		mexErrMsgTxt("Not enough input arguments.");

	//-----------------------------------------------------------------------------------
	// INPUT 0, X
	//-----------------------------------------------------------------------------------
	if (!mxIsStruct(prhs[0]))
		mexErrMsgTxt("The input argument o must be a structure.");
	if (!mxIsCell(prhs[1]))
		mexErrMsgTxt("The input argument Xcell must be a cell array.");
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Ncats"); 
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Ncats");
	int Ncats = ((int*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Ndims");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Ndims");
	int Ndims = ((int*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Ntsamp");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Ntsamp");
	int Ntsamp = ((int*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Nsamp");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Nsamp");
	int* Nsamp = (int*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Npoly");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Npoly");
	int Npoly = ((int*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Rank");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Npoly");
	int Rank = ((int*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Ncoeff");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Ncoeff");
	int Ncoeff = ((int*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"CoeffDims");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field CoeffDims");
	int* CoeffDims = (int*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"ForceEqualPriors");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field ForceEqualPriors");
	bool ForceEqualPriors = ((bool*)mxGetData(arg))[0];
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"CatWeight");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field CatWeight");
	double* CatWeight = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Cat");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Cat");
	int* Cat = (int*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"X");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field X");
	double* X = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Xmeans");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Xmeans");
	double* Xmeans = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Xvars");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Xvars");
	double* Xvars = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Xmean");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Xmean");
	double* Xmean = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Xvar");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Xvar");
	double* Xvar = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Xscale");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Xscale");
	double* Xscale = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"Xacc");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field Xacc");
	double* Xacc = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"wScale");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field wScale");
	double* wScale = (double*)mxGetData(arg);
	//--------------------------------------------------------
	arg = mxGetField(prhs[0],0,"wInit");
	if( arg==null ) mexErrMsgTxt("The input structure must contain a field wInit");
	double* wInit = (double*)mxGetData(arg);
	bool wInitOverride = mxIsNaN(wInit[0]);
	//--------------------------------------------------------
	double** Xcell = new double*[Ncats];
	arg = (mxArray*)mxGetData(prhs[0]);
	for(iCat=0; iCat<Ncats; iCat++)
	{
		arg = mxGetCell(prhs[1], iCat);
		if(!mxIsDouble(arg))
			mexErrMsgTxt("The elements inside the input cell array X must be double arrays.");
		Xcell[iCat] = (double*)mxGetData(arg);
		if( iCat>0 && Ndims!=((int*)mxGetDimensions(arg))[1] )
			mexErrMsgTxt("The elements inside the input cell array X must be double arrays with the same number of columns.");
	}
	//----------------------------------------------------------------------------------------------
	// STATIC
	//--------------------------------------------------------
	//	Register the function that must delete static memory
	if( nPoly_Init_Idx==0 )
		mexAtExit(Poly_Init_DeleteStaticMemory);
	//	Determine whether static memory must be initialized
	if( nPoly_Init_Idx<Ntsamp )
	{
		delete Poly_Init_Idx;
		nPoly_Init_Idx = Ntsamp;
		Poly_Init_Idx = new int[nPoly_Init_Idx];
	}
	if( nPoly_Init_AccMin<Ncats )
	{
		delete Poly_Init_AccMin;
		nPoly_Init_AccMin = Ncats;
		Poly_Init_AccMin = new double[nPoly_Init_AccMin];
	}
	//----------------------------------------------------------------------------------------------
	// OUTPUT
	//-----------------------------------------------------------------------------------
	// OUTUPT CatWeight
	//--------------------------------------------------------
	double cwSum = 0.0;
	if (ForceEqualPriors)
	{
		for( iCat=0; iCat<Ncats; iCat++ )
		{
			CatWeight[iCat]=(double)Ntsamp/(double)Nsamp[iCat]/(double)Ncats;
			cwSum += CatWeight[iCat]*(double)Nsamp[iCat];
		}
	}
	else
	{
		for( iCat=0; iCat<Ncats; iCat++ )
		{
			CatWeight[iCat]=1.0;
			cwSum += (double)Nsamp[iCat];
		}
	}
	//-----------------------------------------------------------------------------------
	// Prep
	//--------------------------------------------------------
	//	Get minimum accuracy for each category.
	for( iCat=0; iCat<Ncats; iCat++ )
	{
		Poly_Init_AccMin[iCat] = (cwSum-CatWeight[iCat]*(double)Nsamp[iCat])/cwSum;
		if( Poly_Init_AccMin[iCat] < 0.5 )
			Poly_Init_AccMin[iCat] = 1.0-Poly_Init_AccMin[iCat];
	}
	//	Declare
	int iSamp,jSamp,iCoeff,iRank,iPoly;
	double x, q025, q975;
	double *pX;
	//	Strides
	int *s0_X = new int[Ncoeff];
	int *s0_Xmeans = new int[Ncoeff];
	for( iCoeff=0; iCoeff<Ncoeff; iCoeff++ )
	{
		s0_X[iCoeff] = Ntsamp*iCoeff;
		s0_Xmeans[iCoeff] = Ncats*iCoeff;
		for(iCat=0; iCat<Ncats; iCat++)
		{
			i = iCat + s0_Xmeans[iCoeff];
			//	Zero out summary statistics here.
			Xmeans[i]=0.0;
			Xvars[i]=0.0;
		}
	}
	int *s0_CoeffDims = new int[Rank];
	for( iRank=0; iRank<Rank; iRank++ )
		s0_CoeffDims[iRank] = Ncoeff*iRank;
	int *s0_Xcell = new int[Ndims]; // computed later (for each iCat category)

	//-----------------------------------------------------------------------------------
	// OUTPUT Cat, X, Xmeans, Xvars, Xmean, Xvar, Xscale
	//--------------------------------------------------------
	double catPrior=1.0/(double)Ncats;
	//	Initialize Xmean, Xvar to 0
	for( iCoeff=0; iCoeff<Ncoeff; iCoeff++ )
		Xmean[iCoeff]=Xvar[iCoeff]=0.0;
	jSamp=0;
	for(iCat=0; iCat<Ncats; iCat++)
	{
		//	Stride for cell array, current category.
		for( iDim=0; iDim<Ndims; iDim++ )
			s0_Xcell[iDim] = iDim*Nsamp[iCat];
		//	Go through each stimulus and dimension for the current category
		for(iSamp=0; iSamp<Nsamp[iCat]; iSamp++,jSamp++)
		{
			Cat[jSamp] = iCat;
			for( iCoeff=0; iCoeff<Ncoeff; iCoeff++ )
			{
				x=1.0;
				for( iRank=0; iRank<Rank; iRank++ )
				{
					iDim = CoeffDims[iCoeff+s0_CoeffDims[iRank]];
					if( iDim<0 )
						break;
					x *= Xcell[iCat][iSamp+s0_Xcell[iDim]];
				}
				X[jSamp+s0_X[iCoeff]] = x;
				i = iCat + s0_Xmeans[iCoeff];
				Xmeans[i] += x;
				Xvars[i]  += x*x;
			}
		}
		//	Finish Xmeans, Xvars.  Do Xmean, Xvar.
		if( !ForceEqualPriors )
			catPrior = (double)Nsamp[iCat]/(double)Ntsamp;
		for( iCoeff=0; iCoeff<Ncoeff; iCoeff++ )
		{
			i = iCat + s0_Xmeans[iCoeff];
			//	Finish Xmeans, Xvars.  Do Xmean, Xvar.
			Xmeans[i] /= (double)Nsamp[iCat];
			Xvars[i]   = (Xvars[i]/(double)Nsamp[iCat] - Xmeans[i]*Xmeans[i]) * ( (double)Nsamp[iCat]/(double)(Nsamp[iCat]-1) );
			Xmean[iCoeff] += Xmeans[i] * catPrior;
			Xvar[iCoeff] += Xvars[i] * catPrior;
			//	Initialize integer index for sorting
			for(iSamp=0; iSamp<Nsamp[iCat]; iSamp++)
				Poly_Init_Idx[iSamp] = iSamp;
			//	Point to the values of X that should be sorted
			pX = X+s0_X[iCoeff]+jSamp-Nsamp[iCat];
			//	Sort to get quantiles
			QuickSort(Poly_Init_Idx, pX, 0, Nsamp[iCat]-1);
			q025 = pX[Poly_Init_Idx[(int)(0.5+0.025*(double)(Nsamp[iCat]-1))]];
			q975 = pX[Poly_Init_Idx[(int)(0.5+0.975*(double)(Nsamp[iCat]-1))]];
			//	Record the scale
			Xscale[i] = q975-q025;
		}
	}
	//-----------------------------------------------------------------------------------
	// OUTPUT Xacc
	//--------------------------------------------------------
	int idxCrit = 0;
	for( iCoeff=0; iCoeff<Ncoeff; iCoeff++ )
	{
		//	Initialize integer index for sorting
		for(iSamp=0; iSamp<Ntsamp; iSamp++)
			Poly_Init_Idx[iSamp] = iSamp;
		//	Point to the values of X that should be sorted
		pX = X+s0_X[iCoeff];
		//	Sort to prepare for accuracy calculation
		QuickSort(Poly_Init_Idx, pX, 0, Ntsamp-1);
		//	For each category
		for(iCat=0; iCat<Ncats; iCat++)
		{
			//	Signed Accuracy
			i = iCat+s0_Xmeans[iCoeff];
			Xacc[i] = Accuracy(&idxCrit, Poly_Init_Idx, pX, Cat, CatWeight, Ntsamp, Ncats, iCat);
			if( idxCrit<0 || idxCrit > Ntsamp )
				Xacc[i] = -Xacc[i];
		}
	}
	//-----------------------------------------------------------------------------------
	// OUTPUT wScale, wInit
	//--------------------------------------------------------
	double tAcc,acc;
	double invNsamp = 1.0/(double)Ntsamp;
	for( iCoeff=0; iCoeff<Ncoeff; iCoeff++ )
	{
		for( iCat=0; iCat<Npoly; iCat++ )
		{
			//	Point to wInit and wScale
			i = iCoeff + iCat*Ncoeff;
			//-----------------------
			//	To start, calculate "x" the weighted Xscale.  Note, (double)Ncats has not been divided out yet.
			//-----------------------
			x = Xscale[s0_Xmeans[iCoeff]];
			for( jCat=1; jCat<Ncats; jCat++ )
				x += Xscale[jCat+s0_Xmeans[iCoeff]];
			//	Now calculate scale of weights (should be the inverse of the scale of the values).
			wScale[i] = (double)Ncats/x;
			//-----------------------
			//	Now set wInit[i] equal to the inverse of this scale and weight by the signed cushioned fisher transform of absolte accuracy
			//-----------------------
			if( wInitOverride )
			{
				//	Get accuracy
				acc = Xacc[iCat+s0_Xmeans[iCoeff]];
				//	Transform it
				tAcc = fabs(acc)-Poly_Init_AccMin[iCat];
				tAcc = tAcc>0.0 ? tAcc : 0.0; // ensure non-negative
				tAcc /= (1.0-Poly_Init_AccMin[iCat]+invNsamp); // one-sided fisher transform with upper limit based on accuracy margin invNsamp
				if( acc>0.0 )
					wInit[i] =  tAcc*wScale[i]; // signed accuracy is positive
				else
					wInit[i] = -tAcc*wScale[i]; // signed accuracy is negative
			}
			//-----------------------
		}
	}
	//-----------------------------------------------------------------------------------
	
	//	Delete allocated memory.
	delete Xcell;
	delete s0_Xcell;
	delete s0_X;
	delete s0_Xmeans;
	delete s0_CoeffDims;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Poly_Init(nlhs, plhs, nrhs, prhs);
}