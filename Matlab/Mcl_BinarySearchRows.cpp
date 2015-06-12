//	Almon David Ing
//	Ctr. Perceptual Systems
//	University of Texas at Austin
//	Compiled July 6, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"

//================================================================
//function OUT = Mcl_BinarySearchRows(M, ICOL, TARGET, [ROWMIN], [ROWMAX])
//----------------------------------------------------------------
// This mex function searches the rows of M for a TARGET value.
//----------------------------------------------------------------
// M  A matrix (double) that has already been sorted (ascending) by M(:,ICOL).
//
// ICOL  An integer (as double) specifying the column to search in.  This is a 1-based column index.
//
// TARGET  The target value to search for.
//
// [ROWMIN]  An optional integer (as double) specifying the lowest row to search in.  This is a 1-based row index.
//	Default ROWMIN = 1.
//
// [ROWMAX]  An optional integer (as double) specifying the lowest row to search in.  This is a 1-based row index.
//	Default ROWMAX = size(M,1).
//----------------------------------------------------------------
// OUT  (double) The floating-point index of the TARGET.
//	If all M(:,ICOL) are below TARGET, then OUT is -Inf;
//	If all M(:,ICOL) are above TARGET, then OUT is +Inf;
//	If M(:,ICOL) contains many target values, then OUT is a floating-point index (1-based) centered in that range.
//	If M(:,ICOL) contains one target value, then OUT is the 1-based index of the target in M(:,ICOL).
//	If M(:,ICOL) contains no target values, then OUT is a floating-point value between the indexes of lower-adjacent and upper-adjacent values in M(:,ICOL).
//================================================================
void Mcl_BinarySearchRows(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<3 )
		mexErrMsgTxt ("There must be three inputs.");
	if( mxGetNumberOfDimensions(prhs[0])<2 )
		mexErrMsgTxt ("The argument M must be a 2-Dimensional array.");
	if (!mxIsDouble (prhs[0]))
		mexErrMsgTxt ("The argument M must be type double.");
	if (!mxIsDouble (prhs[1]))
		mexErrMsgTxt ("The argument ICOL must be type double.");
	if (!mxIsDouble (prhs[2]))
		mexErrMsgTxt ("The argument TARGET must be type double.");
	if (!mxIsDouble (prhs[3]))
		mexErrMsgTxt ("The argument ROWMIN must be type double.");
	if (!mxIsDouble (prhs[4]))
		mexErrMsgTxt ("The argument ROWMAX must be type double.");

	//	Determine dimensionality of M
	int nRows = (int)((mxGetDimensions(prhs[0]))[0]);
	int nCols = (int)((mxGetDimensions(prhs[0]))[1]);
	double target = mxGetScalar(prhs[2]);
	double* M = mxGetPr(prhs[0]);

	//	Get argument ICOL and then turn it into a stride.
	int sCol = 0;
	if ((int)(mxGetNumberOfElements(prhs[1]))>=1)
		sCol = (int)(0.5+mxGetScalar(prhs[1]))-1;
	if( sCol<0 || sCol>=nCols )  //	Ensure that ICOL is in bounds
		mexErrMsgTxt("The column you specified for ICOL does not exist in the matrix M.");

	//	Make sCol the stride.
	sCol *= nRows;

	//	Get row minimum and maximum
	int rowMin = 0;
	int rowMax = nRows-1;

	//	Get argument ROWMIN
	if( nrhs>=4 && ((int)(mxGetNumberOfElements(prhs[3]))>=1) )
		rowMin = (int)(0.5+mxGetScalar(prhs[3]))-1;
	if( rowMin<0 || rowMin>=nRows )
		mexErrMsgTxt("The row index you specified for ROWMIN does not exist in the matrix M.");
	//	Get argument ROWMAX
	if( nrhs>=5 && ((int)(mxGetNumberOfElements(prhs[4]))>=1) )
		rowMax = (int)(0.5+mxGetScalar(prhs[4]))-1;
	if( rowMax<0 || rowMax>=nRows )
		mexErrMsgTxt("The row index you specified for ROWMAX does not exist in the matrix M.");

	//	Start off the range delimiters [iLow, iHigh]
	int iLow=rowMin;
	int iHigh=rowMax;
	int iMid;
	if( M[sCol+iLow]<target )
	{
		plhs[0] = mxCreateDoubleScalar( -mxGetInf() );
		return;
	}
	if( M[sCol+iHigh]>target )
	{
		plhs[0] = mxCreateDoubleScalar(  mxGetInf() );
		return;
	}

	//	Get the middle coordinate
	while (true)
	{
		iMid=(iHigh+iLow)/2;

		//-----------------------------------------------------------------------------
		//	Finish
		//-----------------------------------------------------------------------------
		if( M[sCol+iMid]==target )
		{
			//	Put both indexes in range and rely on execution of next "if" statement to finish.
			iLow=iHigh=iMid;
		}
		if( iHigh-iLow<5 )
		{
			//	Get iLow to be the lowest index of M(:,ICOL) equal to the target.
			//	Otherwise it should be just below the target.
			while( iLow>rowMin && M[sCol+iLow-1]>=target ) iLow--;
			while( iLow<rowMax && M[sCol+iLow+1]<target ) iLow++;
			if( iLow<rowMax && M[sCol+iLow+1]==target) iLow++;
			//	Get iHigh to be the highest index of M(:,ICOL) equal to the target.
			//	Otherwise it should be just above the target.
			while( iHigh<rowMax && M[sCol+iHigh+1]<=target ) iHigh++;
			while( iHigh>rowMin && M[sCol+iHigh-1]<target ) iHigh--;
			if( iHigh>rowMin && M[sCol+iHigh-1]==target) iHigh--;
			
			if( iHigh==iLow )
			{
				// One unique target value.  Return its index.
				plhs[0] = mxCreateDoubleScalar( (double)(iLow+1) );
				return;
			}
			if( M[sCol+iLow]==target )
			{
				// Many unique target values exist.  Return the center of the range.
				plhs[0] = mxCreateDoubleScalar( 0.5*(double)(iLow+iHigh+2) );
				return;
			}
			if( iHigh-iLow != 1 )
				mexErrMsgTxt("Error in user input:  The input M(ROWMIN:ROWMAX,ICOL) must be sorted in ascending order.  In other words, it must be non-decreasing.");
			plhs[0] = mxCreateDoubleScalar( (double)(iLow+1) + (target-M[sCol+iLow])/(M[sCol+iHigh]-M[sCol+iLow]) );
			return;
		}
		//-----------------------------------------------------------------------------

		//  Reduce possible range by 1/2
		if( M[sCol+iMid] > target )
			iLow=iMid;
		else // Must be < target
			iHigh=iMid;
	}		
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_BinarySearchRows(nlhs, plhs, nrhs, prhs);
}