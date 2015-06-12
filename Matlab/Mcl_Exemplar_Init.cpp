//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "QuickSort.cpp"

#ifndef null
	#define null 0
#endif

//================================================================================================================================
//function Mcl_Exemplar_Init(o, Xcell)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function initializes a new Mcl_Exemplar classifier.  It is called from the constructor Mcl_Exemplar_Ctor and basically
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
// IdMethod (int32 scalar) Suggest 1.
//	Identifies the generalization (blurring) function to be used.  Must be an integer belonging to the set {0,1,2}
//		0	k = 1/dist
//		1	k = exp(-dist)
//		2	k = exp(-dist*dist)
//
// MaxNeighbors (int32 scalar)  Suggest 50.
//	The maximum number of nearest neighbors to consider.  When distance is computed, the optional parameter wInit is used.  If wInit
//	is not provided, then it is estimated based on the summary statistics of the training data.  Thus, each axis is scaled before
//	distance is computed (according to the elements of wInit).
//
// ForceEqualPriors (bool scalar) Suggest true.
//	If true, then the categories are assumed to have equal prior probabilities, even if there are different numbers of samples for
//	each category in the training data.
//
// Nquantiles (int32 scalar) Suggest any integer from 10 to 50.
//	The number of quantiles to be used in estimating the probability of category membership based on the decision variable.  In general,
//	a good value for Nquantiles is 10*Ncats or 20*Ncats;
//
// [wInit] (optional double vector, Ndims) Suggest not providing this argument, let this function compute it automatically.
//	The initial weight vector (initial free parameters) for the exemplar model.  These weights are used to scale each dimension of the 
//	training data before distance is computed.  In other words, this controls the width of the blurring function on each spatial dimension.
//	It is highly reccommended that you allow this input to be automatically computed.  Do not specify this input unless you know what
//	you're doing.  It is possible to screw up drastically if you supply bad values in wInit.  This is because the "nearest neighbor"
//	calculations are based on wInit, and failure to define a sensible shape of the neighborhood (via wInit) leads to failure of
//	neighborhood blurring.
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
// IdMethod (int32 scalar)
//	Same as the input IdMethod.
// MaxNeighbors (int32 scalar)
//	Same as the input MaxNeighbors.
// Nneighbors (int32 scalar)
//	Same as the input MaxNeighbors.  This value can be adjusted before Mcl_Exemplar_Train is executed to reduce the number of nearest
//	neighbors considered by the model.
// ForceEqualPriors (bool scalar)
//	Same as input ForceEqualPriors.
// CatWeight (double vector, Ncats)
//	If the training sample has equal numbers of samples from each category, then this will be a vector of ones; otherwise, if ForceEqualPriors
//	is true, then these weights will scale the influence of training samples so that the total weight is Ntsamp and that each category
//	has equal influence on the solution.
// Cat (int32 vector, Ntsamp)
//	The 0-based category label of each sample after samples from all categories are concatenated into a single list.
// X (double 2D array, Ntsamp * Ndims)
//	The spatial coordinate of each sample after samples from all categories are concatenated into a single list.
// Xmeans (double 2D array, Ncats * Ndims)
//	The mean of the sample from each category.
// Xvars (double 2D array, Ncats * Ndims)
//	The variance of the sample from each category.
// Xmean (double vector, Ndims)
//	The pooled mean (where CatWeight was used to weight each category's influence on the result).
// Xvar (double vector, Ndims)
//	The pooled variance (where CatWeight was used to weight each category's influence on the result).  Variance of category means is not included.
//	This is literally the weighted some of Xvars across categories.
// Cube (float 3D array, Ndims * MaxNeighbors * Ntsamp)
//	This cube lists the squared difference between every sample and its MaxNeighbors nearest neighbors.  Squared differences are listed for each
//	iDim spatial dimension as Cube(1+iDim,:,:).  The entries of Cube(1,:,:) are the category label (1-based) of the nearest neighbor.
// Dv (double 2D array, Ntsamp * Ncats):  This variable is calculated in Mcl_Exemplar_Train.
//	This lists the decision variable calculated for each training sample (i.e. the kernel density of each sample applied to each category, normalized)
//	across categories.  A decision variable is computed separately for each column.  Each column represents the decision values with respect to a
//	particular category (i.e. there are Ncats decision functions).
// P (double 2D array, Ntsamp * Ncats):  This variable is calculated in Mcl_Exemplar_Train.
//	This lists the probability that each category should be assigned to each training sample.  Probabilities are computed by linear interpolation
//	of Dv through the quantization table [Quant.Dv, Quant.PcMonoLim] and then normalizing each row to sum to 1.0.
// H (double vector, Ntsamp):  This variable is calculated in Mcl_Exemplar_Train.
//	This lists the conditional entropy of each training sample.  Entropy appraches 0.0 when the sample is categorized correctly, 1.0 is chance performance.
// wInit (double vector, Ndims)
//	The initial value of the free parameters for the blurring or generalization function.  These parameters control the width of the blurring
//	kernel for each spaital dimension.
// h (double scalar):  This variable is calculated in Mcl_Exemplar_Train.
//	The conditional entropy of the entire data set in the range [Quant.hMin, 1.0].  This value does not approach 0.0 because the probabilities of
//	category membership are limited by Quant.PcMonoLim, and these probabilities are never allowed to equal 0.0 or 1.0.  Instead they are limited
//	to the range [Quant.pLow, Quant.pHigh].
// hLim (double scalar):  This variable is calculated in Mcl_Exemplar_Train.
//	The conditional entropy mapped uniformly from the range [Quant.hMin, 1.0] to the range [0.0, 1.0].  This is like "stretching" the range.
// Acc (double scalar):  This variable is calculated in Mcl_Exemplar_Train.
//	The weighted accuracy.  1.0 is perfect accuracy, 1/Ncats is chance.  Weights are specified in CatWeight.
// Quant (mxArray matlab structure, Ncats):  This variable is calculated in Mcl_Exemplar_Train.
//	This structure specifies the mapping between the decision function for category iCat and the probability of category membership.  It is calculated
//	based on quantization of each category's decision variable and training sample.
// Quant(iCat).Nquantiles (int32 scalar):  This variable is calculated in Mcl_Exemplar_Train.
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
// SolverOutput (mxArray, Matlab structure):  This variable is calculated in Mcl_Exemplar_Train.
//	After Mcl_Exemplar_Train is executed, this variable contains raw output from the optimization algorithm.
// Mem (mxArray, Matlab structure)
//	This is just memory allocated to make Mcl_Exemplar_Train run as efficiently as possible.
//================================================================================================================================
void Mcl_Exemplar_Init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int iDim, iCat;
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

	int Ncats = ((int*)mxGetData(mxGetField(prhs[0],0,"Ncats")))[0];
	int* Nsamp = (int*)mxGetData(mxGetField(prhs[0],0,"Nsamp"));
	int Ntsamp = ((int*)mxGetData(mxGetField(prhs[0],0,"Ntsamp")))[0];
	int Ndims = ((int*)mxGetData(mxGetField(prhs[0],0,"Ndims")))[0];
	int IdMethod = ((int*)mxGetData(mxGetField(prhs[0],0,"IdMethod")))[0];
	int idGnrlz = IdMethod%10; // 0 decays 1/dist, 1 decays exp, 2 decays gauss
	int idPrms = IdMethod/10;  // 0 categories share parameters, 1 unique params for each category
	int MaxNeighbors = ((int*)mxGetData(mxGetField(prhs[0],0,"MaxNeighbors")))[0];
	if (MaxNeighbors < Ncats || MaxNeighbors > Ntsamp-1)
		mexErrMsgTxt("The input argument MaxNeighbors must be an integer less than Ntsamp and cannot be less than Ncats.");
	bool ForceEqualPriors = ((bool*)mxGetData(mxGetField(prhs[0],0,"ForceEqualPriors")))[0];
	int Nquantiles = ((int*)mxGetData(mxGetField(mxGetField(prhs[0],0,"Quant"),0,"Nquantiles")))[0];
	double* wInit = (double*)mxGetData(mxGetField(prhs[0],0,"wInit"));
	bool wInitOverride = mxIsNaN(wInit[0]);
	
	double logBase = 1.0/log((double)Ncats);  //	Multiply by this number to convert logarithm to base Ncats
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
	//-----------------------------------------------------------------------------------

	double* CatWeight = (double*)mxGetData(mxGetField(prhs[0],0,"CatWeight"));
	double* Xmean = (double*)mxGetData(mxGetField(prhs[0],0,"Xmean"));
	double* Xvar = (double*)mxGetData(mxGetField(prhs[0],0,"Xvar"));
	double* Xmeans = (double*)mxGetData(mxGetField(prhs[0],0,"Xmeans"));
	double* Xvars = (double*)mxGetData(mxGetField(prhs[0],0,"Xvars"));
	double* X = (double*)mxGetData(mxGetField(prhs[0],0,"X"));
	float* Cube = (float*)mxGetData(mxGetField(prhs[0],0,"Cube"));
	int* Cat = (int*)mxGetData(mxGetField(prhs[0],0,"Cat"));
	int HalfNeighbors = MaxNeighbors/2;
	int hSkip;

	//	Calculate CatWeight
	if (ForceEqualPriors)
	{
		for( iCat=0; iCat<Ncats; iCat++ )
			CatWeight[iCat]=(double)Ntsamp/(double)Nsamp[iCat]/(double)Ncats;
	}
	else
	{
		for( iCat=0; iCat<Ncats; iCat++ )
			CatWeight[iCat]=1.0;
	}

	//	Pre-compute a stride variable.
	int* dimStride = new int[Ndims+2];
	for( iDim=0; iDim<Ndims+2; iDim++ )
		dimStride[iDim] = iDim*Ntsamp;
	int* dimStrideCat = new int[Ndims];

	//-----------------------------------------------------------------------------------
	//	OUTPUT:  Cat, X, Xmeans, Xvars, Xmean, Xvar.
	//-----------------------------------------------------------------------------------
	double x,nCatsInv=1.0/(double)Ncats;
	int i,j,iSamp,jSamp=0;
	//	Initialize Xmean, Xvar to 0
	for( iDim=0; iDim<Ndims; iDim++ )
		Xmean[iDim]=Xvar[iDim]=0.0;
	for(iCat=0; iCat<Ncats; iCat++)
	{
		//	Pre-compute a stride variable and zero out Xmeans, Xvars
		for( iDim=0; iDim<Ndims; iDim++ )
		{
			i = iCat + iDim*Ncats;
			dimStrideCat[iDim] = iDim*Nsamp[iCat];
			Xmeans[i]=0.0;
			Xvars[i]=0.0;
		}
		//	Go through each stimulus and dimension for the current category
		for(iSamp=0; iSamp<Nsamp[iCat]; iSamp++,jSamp++)
		{
			for( iDim=0; iDim<Ndims; iDim++ )
			{
				Cat[jSamp] = iCat;
				x = Xcell[iCat][iSamp+dimStrideCat[iDim]];
				X[jSamp+dimStride[iDim]] = x;
				i = iCat + iDim*Ncats;
				Xmeans[i] += x;
				Xvars[i]  += x*x;
			}
		}
		//	Finish Xmeans, Xvars.  Do Xmean, Xvar.
		if( !ForceEqualPriors )
			nCatsInv = (double)Nsamp[iCat]/(double)Ntsamp;
		for( iDim=0; iDim<Ndims; iDim++ )
		{
			i = iCat + iDim*Ncats;
			Xmeans[i] /= (double)Nsamp[iCat];
			Xvars[i]   = (Xvars[i]/(double)Nsamp[iCat] - Xmeans[i]*Xmeans[i]) * ( (double)Nsamp[iCat]/(double)(Nsamp[iCat]-1) );
			Xmean[iDim] += Xmeans[i] * nCatsInv;
			Xvar[iDim] += Xvars[i] * nCatsInv;
		}
	}
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//	OUTPUT:  wInit
	//-----------------------------------------------------------------------------------
	if( wInitOverride )
	{
		if (idPrms==0)
		{
			for( iDim=0; iDim<Ndims; iDim++ )
				wInit[iDim] = 1.0/sqrt(Xvar[iDim]);
		}
		else
		{
			for( iDim=0; iDim<Ndims; iDim++ )
			{
				for( iCat=0; iCat<Ncats; iCat++ )
				{
					i = iDim + iCat*Ndims;
					wInit[i] = 1.0/sqrt(Xvars[i]);
				}
			}
		}
	}
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//	OUTPUT:  Cube
	//-----------------------------------------------------------------------------------
	double* xsdInv = new double[Ndims];
	for(iDim=0; iDim<Ndims; iDim++)
		xsdInv[iDim] = 1.0/sqrt(Xvar[iDim]);
	double* square = new double[ (Ndims+2)*Ntsamp ];
	int* idx = new int[Ntsamp];
	double dsq;
	hSkip = (Ntsamp-HalfNeighbors)/MaxNeighbors;
	for( iSamp=0; iSamp<Ntsamp; iSamp++ )
	{
		for( jSamp=0; jSamp<Ntsamp; jSamp++ )
		{
			//	Reset sort index
			idx[jSamp]=jSamp;
			//	Ensure that a data point cannot reference itself
			if(iSamp==jSamp)
				square[jSamp] = mxGetInf(); //	 Ensures that data is automatically omitted from the cube.
			else
			{
				//	Store category label (0-based)
				square[jSamp+dimStride[1]] = (double)Cat[jSamp];
				//	Initialized squared distance
				dsq = 0.0;
				//  Calculate squared distance and store the squared distance on each dimension.
				for( iDim=0; iDim<Ndims; iDim++ )
				{
					x = (X[iSamp+dimStride[iDim]] - X[jSamp+dimStride[iDim]]) * xsdInv[iDim];
					x*=x;
					square[jSamp+dimStride[iDim+2]] = x;
					dsq += x;
				}
				if( ForceEqualPriors )
				{
					//	Squared distance is fine.  Don't waste computer power taking square-root.  Result would be the same.
					square[jSamp] = dsq;
				}
				else
				{
					//	Compute the distance weighted to offset the priors.
					square[jSamp] = sqrt(dsq) / CatWeight[Cat[jSamp]-1];
				}
			}
		}
		//	Sort by generalization
		QuickSort(idx, square, 0, Ntsamp-1);
		//	Put the nearest neighbors into the first half of the cube, remaining into the second half..
		//	Having some distant neighbors fill up a half of the space prevents the kernel from assuming that nothing exists in distant places.
		for(j=0; j<MaxNeighbors; j++)
		{
			if( j<=HalfNeighbors )
				jSamp = idx[j];
			else
				jSamp = idx[HalfNeighbors + (j-HalfNeighbors)*hSkip];
			i = iSamp*(Ndims+1)*MaxNeighbors + j*(Ndims+1);
			Cube[i] = (float)square[jSamp+dimStride[1]];
			for(iDim=0; iDim<Ndims; iDim++)
				Cube[i+iDim+1] = (float)square[jSamp+dimStride[2+iDim]];
		}
	}
	//-----------------------------------------------------------------------------------

	//	Delete allocated memory.
	delete Xcell;
	delete dimStride;
	delete dimStrideCat;
	delete idx;
	delete square;
	delete xsdInv;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Exemplar_Init(nlhs, plhs, nrhs, prhs);
}