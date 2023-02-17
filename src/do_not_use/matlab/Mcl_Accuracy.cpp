//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "math.h"
#include "mex.h"
#include "Accuracy.cpp"
#include "QuickSort.cpp"

#ifndef null
	#define null 0
#endif

static int *Mcl_Accuracy_Idx=null;
static int nMcl_Accuracy_Idx=0;
// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void Mcl_Accuracy_DeleteStaticMemory(void)
{
	delete Mcl_Accuracy_Idx;
}

//================================================================================================================================
//function [Accuracy, Criterion] = Mcl_Accuracy(X, Cat, Ncats, CatWeight)
//--------------------------------------------------------------------------------------------------------------------------------
// Finds the accuracy and decision criterion for a 1-criterion univariate classifier.  The user supplies data as X and Cat and 
//	the category for which to compute the univariate criterion as idCat.  If there are many columns of X, a separate classifier
//	is computed for each column.  Each row of X is a data sample associated with a category label supplied by Cat.  The samples
//	of each categories are permitted to have differential influence on the solution via the CatWeight parameter which is useful
//	for enforcing priors.
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE
//----------------------------------------------------
// nSamp (int32 vector: nCats)
//	The number of samples in each category.  For each iCat category, nSamp(1+iCat) is equal to sum(Cat==int32(iCat)).
// ntSamp (int32 scalar)
//	The total number of samples (in all categories).  Equal to sum(nSamp) or numel(Cat).
// nCols
//	The number of columns of X.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not altered by mex function)
//----------------------------------------------------
// X (double 2D array: ntSamp * nCols)
//	Each row of X gives a the coordinate of a sample of data.  Each row is a separate sample, each column is a spatial dimension.
//	The univariate classifiers are applied separately to each column of X.
// Cat (int32 vector: ntSamp)
//	The 0-based category ID.  The elements of this vector must be integers from the set {0,...,Ncats-1}.
// Ncats (int32 scalar)
//	The number of categories.
// CatWeight (optional double vector: nCats)
//	For category iCat, the weighted accuracy is computed using weights for each category's samples, CatWeight[iCat];  Default = ones(nCats,1);
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (values altered by mex function)
//----------------------------------------------------
// Accuracy (double 2D array: Ncats * nCols)
//	The absolute value of Accuracy(iCat,iCol) gives the weighted accuracy for category iCat in X(:,iCol).  The sign of accuracy
//		disambiguates the decision rule used.
//	Accuracy(iCat,iCol) is negative when the decision rule is: Assign iCat if X(:,iCol)<Criterion(iCat,iCol).
//	Accuracy(iCat,iCol) is positive when the decision rule is: Assign iCat if X(:,iCol)>Criterion(iCat,iCol).
// Criterion (double 2D array: Ncats * nCols)
//================================================================================================================================
void Mcl_Accuracy(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<3 )
		mexErrMsgTxt("Not enough input arguments.");
	if( nlhs==0 )
		mexErrMsgTxt("Not enough output argument.");
	if( nlhs>2 )
		mexErrMsgTxt("Too many output arguments.");
	//----------------------------------------------------------------------------------------------
	// INPUT
	//----------------------------------------------------------------------------------------------
	double *X = null;
	if( !mxIsDouble(prhs[0]) )
		mexErrMsgTxt("X must be type double.");
	X = (double*)mxGetData(prhs[0]);
	int ntSamp = (int)mxGetNumberOfElements(prhs[0]);
	int nCols = 1;
	if( (int)mxGetNumberOfDimensions(prhs[0])>1 )
	{
		ntSamp = (int)mxGetM(prhs[0]);
		nCols = (int)mxGetN(prhs[0]);
	}
	//----------------------------------------------------------------------------------------------
	int *Cat = null;
	if( !mxIsInt32(prhs[1]) )
		mexErrMsgTxt("Cat must be type int32.");
	if( (int)mxGetNumberOfElements(prhs[1]) < ntSamp )
		mexErrMsgTxt("Cat must have ntSamp rows.  Usually, ntSamp = size(X,1).");
	Cat = (int*)mxGetData(prhs[1]);
	//----------------------------------------------------------------------------------------------
	int Ncats = 1;
	if( !mxIsInt32(prhs[2]) )
		mexErrMsgTxt("Ncats must be type int32.");
	if( mxGetNumberOfElements(prhs[2])!=1 )
		mexErrMsgTxt("Ncats must be a scalar.");
	Ncats = ((int*)mxGetData(prhs[2]))[0];
	//----------------------------------------------------------------------------------------------
	double *CatWeight = null;
	if( nrhs>=4 )
	{
		if( !mxIsDouble(prhs[3]) || (int)mxGetNumberOfElements(prhs[3])<Ncats )
			mexErrMsgTxt("CatWeight must be a double vector of length Ncats.");
		CatWeight = (double*)mxGetData(prhs[3]);
	}
	//----------------------------------------------------------------------------------------------
	// OUTPUT
	//----------------------------------------------------------------------------------------------
	double* Acc = null;
	plhs[0] = mxCreateDoubleMatrix(Ncats,nCols,mxREAL);
	Acc = (double*)mxGetData(plhs[0]);
	//----------------------------------------------------------------------------------------------
	double* Crit = null;
	if( nlhs>=2 )
	{
		plhs[1] = mxCreateDoubleMatrix(Ncats,nCols,mxREAL);
		Crit = (double*)mxGetData(plhs[1]);
	}
	//----------------------------------------------------------------------------------------------
	// STATIC
	//----------------------------------------------------------------------------------------------
	//	Register the function that must delete static memory
	if( nMcl_Accuracy_Idx==0 )
		mexAtExit(Mcl_Accuracy_DeleteStaticMemory);
	//	Determine whether static memory must be initialized
	if( nMcl_Accuracy_Idx<ntSamp )
	{
		delete Mcl_Accuracy_Idx;
		nMcl_Accuracy_Idx = 2*ntSamp;
		Mcl_Accuracy_Idx = new int[nMcl_Accuracy_Idx];
	}
	//----------------------------------------------------------------------------------------------

	//mexPrintf("Mcl_Accuracy:  [ntSamp, nCols, Ncats] = [%i, %i, %i]\n", ntSamp, nCols, Ncats);

	int idxCrit = 0;
	int iCol, iSamp, iCat, iOut;
	double *Xcol;
	//	For each column of data
	for(iCol=0; iCol<nCols; iCol++)
	{
		//	Initialize the index
		for(iSamp=0; iSamp<ntSamp; iSamp++)
			Mcl_Accuracy_Idx[iSamp] = iSamp;

		//mexPrintf("Mcl_Accuracy 02:  [iCol] = [%i]\n", iCol);

		//	Sort
		Xcol = X+iCol*ntSamp;
		QuickSort(Mcl_Accuracy_Idx, Xcol, 0, ntSamp-1);

		//mexPrintf("Mcl_Accuracy 03:  [iCol] = [%i]\n", iCol);

		//	For each category
		for(iCat=0; iCat<Ncats; iCat++)
		{
			
			//	Accuracy
			iOut = iCat+Ncats*iCol;

			//mexPrintf("-------------------------------------------------------------------------\n");
			//mexPrintf("Mcl_Accuracy 04:  [iCol, iCat, iOut] = [%i, %i, %i]\n", iCol, iCat, iOut);
			//mexPrintf("-----------------\n");

			Acc[iOut] = Accuracy(&idxCrit, Mcl_Accuracy_Idx, Xcol, Cat, CatWeight, ntSamp, Ncats, iCat);
			if( idxCrit<0 || idxCrit > ntSamp )
				Acc[iOut] = -Acc[iOut];

			//mexPrintf("Mcl_Accuracy 05:  [iCol, iCat, iOut, Acc[iOut], idxCrit] = [%i, %i, %i, %f, %i]\n", iCol, iCat, iOut, Acc[iOut], idxCrit);

			//	Get the crit
			if(Crit!=null)
			{
				if(idxCrit<-ntSamp || idxCrit>ntSamp)
					Crit[iOut] = Xcol[Mcl_Accuracy_Idx[ntSamp-1]]+1.0;
				else if(idxCrit==0 )
					Crit[iOut] = Xcol[Mcl_Accuracy_Idx[0]]-1.0;
				else if( idxCrit>0 )
					Crit[iOut] = 0.5 * (Xcol[Mcl_Accuracy_Idx[ idxCrit]] + Xcol[Mcl_Accuracy_Idx[ idxCrit-1]]);
				else
					Crit[iOut] = 0.5 * (Xcol[Mcl_Accuracy_Idx[-idxCrit]] + Xcol[Mcl_Accuracy_Idx[-idxCrit-1]]);
			}

			//mexPrintf("Mcl_Accuracy 06:  [iCol, iCat, iOut, Acc[iOut], Crit[iOut]] = [%i, %i, %i, %f, %f]\n", iCol, iCat, iOut, Acc[iOut], Crit[iOut]);
		}
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Accuracy(nlhs, plhs, nrhs, prhs);
}