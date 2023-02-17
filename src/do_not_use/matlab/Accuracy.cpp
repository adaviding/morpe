//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Written July 23, 2009
//
//	This code is useful for computing accuracy of a classifier given a univariate 2-category classifier which applies a single
//	decision criterion at one location on the univariate axis.

//#include "mex.h"

#ifndef DBL_EPSILON
	#define DBL_EPSILON 2.220446049250313e-16
#endif

#ifndef null
	#define null 0
#endif

// double Accuracy(int* idxCrit, int* Idx, double* Z, int* Cat, double* CatWeight, int ntSamp, int nCats, int aCat)
//-----------------------------------------------------------------------------------------------------------------------
//	Finds the accuracy and decision criterion for a 1-criterion classifier that operates on univariate data (Z) where each datum
//		is associated with a category label (Cat), where the classifier's goal is to maximize the accuracy of assigning category aCat, 
//		given pre-computed 1-based indices Idx which sort the data (Z) in ascending order.  This method accomodates the enforcement
//		of priors via the optional (nullable) CatWeight parameter.
//-----------------------------------------------------------------------------------------------------------------------
// INPUT
//-----------------------------------------------------------------------------------------------------------------------
// Idx (int32 vector: ntSamp)
//	A vector of zero-based indices into Cat and Z.  This is the sort index that indexes Z into ascending order.
// Z (double vector: ntSamp)
//	The univariate values associated with each sample.
// Cat (int32 vector: ntSamp)
//	The 0-based category ID of each sample.
// CatWeight (double vector: nCats)
//	The weight of a sample from each category.  If priors are not being enforced, this can be a vector of 1.0 in each entry.
//	If priors are being enforced, each iCat entry should be (double)ntSamp/(double)nCats/(double)nSamp[iCat] where nSamp
//		is a vector giving the number of samples in each category.  This can be null.
// ntSamp (int32 scalar)
//	The total number of samples from all categories.
// nCats (int32 scalar)
//	The number of categories.
// aCat (int32 scalar)
//	The category ID (0-based) for which the classifier is designed to detect.  Must be an integer from the set {0,...,nCats-1}.
//-----------------------------------------------------------------------------------------------------------------------
// OUTPUT
//-----------------------------------------------------------------------------------------------------------------------
// Return value (float scalar)
//	Outputs the accuracy of the most accurate univariate decision criterion (idxCrit).  This is the total correctly classified
//		weight divided by the total weight.
// idxCrit (pointer to float scalar)
//	If not null, outputs the univariate criterion which divides the samples into two regions maximizing accuracy.
//		If positive, all samples where Z >= Z[Idx[ *idxCrit-1]] are classified as aCat.
//		If negative, all samples where Z <= Z[Idx[-*idxCrit-1]] are classified as aCat.
//		If zero, all samples are classified as aCat.
//		If *idxCrit < -ntSamp, all indices are classified as NOT aCat.
//-----------------------------------------------------------------------------------------------------------------------
double Accuracy(int* idxCrit, int* Idx, double* Z, int* Cat, double* CatWeight, int ntSamp, int nCats, int aCat)
{
	//mexPrintf("Accuracy 01:  [ntSamp, nCats, aCat ,Idx[0], Idx[ntSamp-1]] = [%i, %i, %i, %i, %i]\n", ntSamp, nCats, aCat ,Idx[0], Idx[ntSamp-1]);

	int i, idCat;
	double w;
	double acc = 0.0;
	//	If CatWeight is null, create one.
	double *catWeight = CatWeight;
	if( catWeight==null )
	{
		catWeight = new double[nCats];
		for(i=0; i<nCats; i++)
			catWeight[i] = 1.0;
	}
	//-----------------------------------------------------------------------------------------------------------------------
	//	Compute the total weight for the category and the not-category samples.
	//--------------------------------------------------------
	double wt0=0.0; // not-category weight
	double wt1=0.0; //     category weight
	for(i=0; i<ntSamp; i++)
	{
		idCat = Cat[Idx[i]];
		if( idCat==aCat )
			wt1 += catWeight[idCat];
		else
			wt0 += catWeight[idCat];
	}
	double wt=wt0+wt1;	// total weight
	//-----------------------------------------------------------------------------------------------------------------------
	//	Simulate the criteria by stepping through Z[Idx[.]].
	//	For each criterion value, two decision rules are actually simulated:  "Pos" and "Neg".
	//		The "Pos" rule assigns aCat to values of Z higher than the criterion.
	//		The "Neg" rule assigns aCat to values of Z lower  than the criterion.  
	//	Maximize accuracy by maximizing the difference between the correctly and incorrectly classified weight.
	//-----------------------------------------------------------------------------------------------------------------------
	//	Initialize variables to i == -1
	double wPos=wt1;		// The correctly classified weight for "Pos" rule: aCat--> Z > Z[Idx[i]], initialize to [0]
	double wNeg=wt0;		// The correctly classified weight for "Neg" rule: aCat--> Z < Z[Idx[i]], initialize to [0];
	
	double wPosMax = wPos;		// The maximum of wPos over all i
	int ilPos = -1;				// The lowest value of i for which wPos==wPosMax.
	int  cPos = 0;				// The number of values of i for which wPos==wPosMax.
	
	double wNegMax = wNeg;		// The maximum of wNeg over all i
	int ilNeg = -1;				// The lowest value of i for which wNeg==wNegMax.
	int  cNeg = 0;				// The number of values of i for which wNeg==wNegMax.

	//mexPrintf("Accuracy 02:  [wPos, wNeg, catWeight[0], catWeight[nCats-1]] = [%f, %f, %f, %f]\n", wPos, wNeg, catWeight[0], catWeight[nCats-1]);
	
	//	Try shifting the criterion at each value of the sample (for samples in ascending order).
	for( i=0; i<ntSamp; i++ )
	{
		idCat = Cat[Idx[i]];
		w = catWeight[idCat];
		//	Increment or decrement the correctly classified weight for each decision rule.
		if( idCat==aCat )
		{
			//	Positive rule declines
			wPos -= w;
			//	Negative rule improves
			wNeg += w;
		}
		else
		{
			//	Positive rule improves
			wPos += w;
			//	Negative rule declines
			wNeg -= w;
		}
		//	Track the best-performing positive and negative rules.
		//	If this Z-value is the same as the last Z-value, then a criterion cannot be placed at that i, so the rule does not get tracked.
		if( i==0 || Z[Idx[i]]!=Z[Idx[i-1]] )
		{
			//	A criterion can be placed here because this Z-value increases.
			if( wPos>wPosMax )
			{
				wPosMax = wPos;
				ilPos = i;
				cPos = 1;
			}
			else if( wPos==wPosMax )
			{
				cPos++;
			}
			if( wNeg>wNegMax )
			{
				wNegMax = wNeg;
				ilNeg = i;
				cNeg = 1;
			}
			else if( wNeg==wNegMax )
			{
				cNeg++;
			}
		}
	}

	//mexPrintf("Accuracy 03:  [wPosMax, ilPos, cPos, wNegMax, ilNeg, cNeg] = [%f, %i, %i, %f, %i, %i]\n", wPosMax, ilPos, cPos, wNegMax, ilNeg, cNeg);

	//-----------------------------------------------------------------------------------------------------------------------
	//	Determine the OUTPUT
	//-----------------------------------------------------------------------------------------------------------------------
	if( wPosMax >= wNegMax )
	{
		//	Compute the accuracy
		acc = wPosMax/wt;
		//	Return if index is not required.
		if( idxCrit==null )
			return acc;
		else if( fabs(wPos-wPosMax)<=DBL_EPSILON )
			*idxCrit = ntSamp;
		else if( cPos<=2 || ilPos==-1 )
			*idxCrit = ilPos+1; // Index is known easily, add 1 to shift from Z > Z[Idx[i]] to Z >= Z[Idx[i]]
		else
		{
			//------------------------------------------------------------
			//	Search for the cPos/2 maximum somewhere after ilPos
			//------------------------------------------------------------
			cPos /= 2;
			wPos=wPosMax;
			wPosMax-=DBL_EPSILON; // Allow for truncation error 
			i=ilPos;
			while( cPos>0 )
			{
				idCat = Cat[Idx[i++]];
				w = catWeight[idCat];
				if( idCat==aCat )
					wPos -= w;
				else
					wPos += w;
				//	Assert
				if( i==ntSamp )
				{
					//	Assertion failed.  This should never happen.
					//	Set accuracy that makes no sense.
					acc = -9.9999;
					i=ilPos;
					//	Exit loop
					cPos = 0;
				}
				//	If this Z-value is the same as the last Z-value, then a criterion cannot be placed at that i.
				if( i==0 || Z[Idx[i]]!=Z[Idx[i-1]] )
				{
					//	A criterion can be placed here because this Z-value increases.
					if( wPos>=wPosMax )
						cPos--;
				}
			}
			*idxCrit = i+1;  //	Add 1 to shift from Z > Z[Idx[i]] to Z >= Z[Idx[i]]
			//------------------------------------------------------------
		}
	}
	else
	{
		//	Compute the accuracy.
		acc = wNegMax/wt;
		//	Return if index is not required.
		if( idxCrit==null )
			return acc;
		else if( fabs(wNeg-wNegMax)<=DBL_EPSILON )
			*idxCrit = -ntSamp-1;
		else if( ilNeg==-1 )
			*idxCrit = ilNeg;
		else if( cNeg<=2 )
			*idxCrit = -ilNeg-1; // Index is known easily, flip sign and subtract 1 to shift from Z < Z[Idx[i]] to Z >= Z[Idx[i]].
		else
		{
			//------------------------------------------------------------
			//	Search for the cNeg/2 minimum somewhere after ilNeg
			//------------------------------------------------------------
			cNeg /= 2;
			wNeg = wNegMax;
			wNegMax-=DBL_EPSILON; // Allow for truncation error 
			i=ilNeg;
			while( cNeg>0 )
			{
				idCat = Cat[Idx[i++]];
				w = catWeight[idCat];
				if( idCat==aCat )
					wNeg += w;
				else
					wNeg -= w;
				//	Assert
				if( i==ntSamp )
				{
					//	Assertion failed.  This should never happen.
					//	Set accuracy that makes no sense.
					acc = -9.9999;
					i=ilNeg;
					//	Exit loop
					cNeg = 0;
				}
				//	If this Z-value is the same as the last Z-value, then a criterion cannot be placed at that i.
				if( i==0 || Z[Idx[i]]!=Z[Idx[i-1]] )
				{
					//	A criterion can be placed here because this Z-value increases.
					if( wNeg>=wNegMax )
						cNeg--;
				}
			}
			*idxCrit = -i-1;  //	Flip sign and subtract 1 to shift from Z < Z[Idx[i]] to Z >= Z[Idx[i]].
			//------------------------------------------------------------
		}
		//	Indicate that all samples are to be categorized as NOT aCat.
		if( *idxCrit==0 )
			*idxCrit = -ntSamp-1;
	}
	//mexPrintf("Accuracy 04:  [acc, *idxCrit] = [%f, %i]\n", acc, *idxCrit);
	//-----------------------------------------------------------------------------------------------------------------------
	if( catWeight != CatWeight )
		delete catWeight;
	return acc;
}