//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"

#ifndef null
	#define null 0
#endif

//================================================================================================================================
//function [h, Hcross] = Mcl_ConditionalEntropy(H, P, Cat, ForceEqualPriors)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function calculates the conditional entropy of ntSamp samples based on the probabilities of category membership for each
//	category P, the actual (known) category labels Cat, and the weight of each category CatWeight.
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting subsequent comments)
//----------------------------------------------------
// nCats = The number of categories of the training data (i.e. the number of columns of P)
// ntSamp = The total number of training samples = sum(nSamp);
// p(jCat|iSamp) = The posterior probability assigned (by the classifier) to category jCat.
// p(iCat) = The prior probability assigned (by the classifier) to category jCat.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not altered by mex function)
//----------------------------------------------------
// P (double 2D array: ntSamp * nCats)
//	The probability (assigned by a classifier) that each sample belongs to each category.  The rows must sum to 1.0.
// Cat (int32 vector: ntSamp)
//	Gives the cagtegory label of each sample.  Each category label must be an integer belonging to the set {0,...,nCats-1}.
// ForceEqualPriors (logical scalar)
//	If true, the prior probability of seeing a sample from each category is forced to be equal.  In such a case, the conditional
//		entropy must weight the influence of each sample so that each category has equal influence on the conditional entropy.
//		While the weight of each row is adjusted, the actual probabilities (entries of P) are not adjusted.  It is assumed that P
//		already contains posterior probabilities before h, H, and Hcross are computed.
//	If false, then no weighting is applied, so each row of P has equal influence on the conditional entropy (irrespective of Cat).
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (values are copied to pre-allocated arrays provided).
//----------------------------------------------------
// h
//	A double scalar value that is the conditional entropy of all samples.
//		 0.0, the classifier performs perfectly.
//		 1.0, the classifier performs at chance (as defined by the priors).
//		>1.0, the classifier performs worse than chance (as defined by the priors).
//
// Hcross (double matrix, nCats * nCats)
//	The conditional cross-entropy matrix.  Each entry Hcross(iCat,jCat) provides the expected (mean) value of the following expression
//		for cases when the sample is actually a member of category iCat.
//			-p(jCat,iSamp) * log[p(jCat|iSamp)]
//	If ForceEqualPriors, this expected value is a weighted mean so that each category is weighted to enforce equal presence in the
//		rows of input P.  However, the values of p(jCat|iSamp) are not adjusted to compensate for priors before entropy is computed.
//
// H (double vector: ntSamp)
//	The conditional entropy for each row of P.  The output h is the sum of all elements of H (always).
//	If ForceEqualPriors, each element of H is weighted so that all categories contribute equally to h.
//================================================================================================================================
void Mcl_ConditionalEntropy(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<4 )
		mexErrMsgTxt("Not enough input arguments.");

	//-----------------------------------------------------------------------------------------
	// INPUT
	//-----------------------------------------------------------------------------------------
	if (!mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1])<2)
		mexErrMsgTxt("The input argument P must be a double 2D array: ntSamp * nCats.");
	int ntSamp = (int)(mxGetDimensions(prhs[1])[0]);
	int nCats = (int)(mxGetDimensions(prhs[1])[1]);
	if( nCats==1 )
		mexErrMsgTxt("There must be more than 1 category.  In other words, the argument P must have more than 1 column.");
	double* P = mxGetPr(prhs[1]);
	//-----------------------------------------------------------------------------------------
	if (!mxIsInt32(prhs[2]) || (int)(mxGetNumberOfElements(prhs[2]))<ntSamp)
		mexErrMsgTxt("The input argument Cat must be an int32 vector: ntSamp.");
	int* Cat = (int*)mxGetData(prhs[2]);
	//-----------------------------------------------------------------------------------------
	if (!mxIsLogical(prhs[3]))
		mexErrMsgTxt("The input argument ForceEqualPriors must be a logical scalar:  Suggest true.");
	bool ForceEqualPriors = ((bool*)mxGetData(prhs[3]))[0];

	int iCat, jCat, iSamp;
	//-----------------------------------------------------------------------------------------
	// OUTPUT
	//-----------------------------------------------------------------------------------------
	if (!mxIsDouble(prhs[0]) || (int)(mxGetNumberOfElements(prhs[2]))<ntSamp)
		mexErrMsgTxt("The output argument H must be a double vector: ntSamp.");
	double* H = mxGetPr(prhs[0]);

	//---------------------------------------
	//	Init cross-entropy matrix (if requested).
	//---------------------------------------
	//	The stride of category in Hcross
	int* s0Hcross=null;
	double* Hcross = null;
	if( nlhs>1 )
	{
		s0Hcross = new int[nCats];
		for(iCat=0; iCat<nCats; iCat++)
			s0Hcross[iCat] = nCats*iCat;
		plhs[1] = mxCreateDoubleMatrix(nCats,nCats,mxREAL);
		Hcross = mxGetPr(plhs[1]);
		for(iCat=0; iCat<nCats; iCat++)
		{
			for(jCat=0; jCat<nCats; jCat++)
				Hcross[iCat+s0Hcross[jCat]]=0.0;
		}
	}
	//---------------------------------------

	//	The stride of category in P
	int* s0P = new int[nCats];
	//	The weight assigned to members of each category.
	double* wCat = new double[nCats];
	int* wiCat = new int[nCats];
	for(iCat=0; iCat<nCats; iCat++)
	{
		wiCat[iCat]=0;
		s0P[iCat] = ntSamp*iCat;
	}

	//	Start by counting the number of samples in each category (important for weighting)
	for(iSamp=0; iSamp<ntSamp ;iSamp++)
		wiCat[Cat[iSamp]]++;

	//	Get the weight for each sample (by category)
	for(iCat=0; iCat<nCats; iCat++)
	{
		if( ForceEqualPriors )
			wCat[iCat] = (double)ntSamp/(double)nCats/(double)wiCat[iCat];
		else
			wCat[iCat] = 1.0;
	}

	//	The log-base scalar, normalizer, and entropy sign reversor.
	double hScalar = -1.0/(double)ntSamp/log((double)nCats);
	double hThis;
	double h = 0.0;
	double pcat;
	for(iSamp=0; iSamp<ntSamp ;iSamp++)
	{
		iCat = Cat[iSamp];
		hThis = wCat[iCat] * hScalar * log(P[iSamp + s0P[iCat]]);
		H[iSamp] = hThis;
		h += hThis;
		if(Hcross != null)
		{
			for(jCat=0; jCat<nCats; jCat++)
			{
				pcat = P[iSamp + s0P[jCat]];
				if( pcat>0.0 )
					Hcross[iCat+s0Hcross[jCat]] += wCat[iCat] * pcat * log(pcat);
			}
		}
	}

	plhs[0] = mxCreateDoubleScalar(h);

	//	Normalize the cross-entropy matrix
	if(Hcross!=null)
	{
		//	The log-base scalar, diagonal normalizer, and entropy sign reversor.
		hScalar = -(double)nCats/(double)ntSamp/log((double)nCats);
		for(iCat=0; iCat<nCats; iCat++)
		{
			for(jCat=0; jCat<nCats; jCat++)
				Hcross[iCat+s0Hcross[jCat]] *= hScalar;
		}
		//	Delete allocated memory.
		delete s0Hcross;
	}

	//	Delete allocated memory.
	delete wCat;
	delete wiCat;
	delete s0P;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_ConditionalEntropy(nlhs, plhs, nrhs, prhs);
}