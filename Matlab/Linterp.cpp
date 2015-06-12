//	Almon David Ing.
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Written July 6, 2009

//	If there is a tabled increasing function x(i) tabulated at each integer value of i:{0,...,iMax}, this function
//	returns the value of the output iTarget associated with given input xTarget which is guaranteed to be in the
//	range [iMin, iMax].  The values iMin and iMax are restricted to 0 <= iMin < iMax <= N-1 where N is the length of the
//	input table xTable.  The values in xTable must be non-decreasing.  The mapping from xTarget (input) to iTarget (output)
//	is determined by linear interpolation of the table xTable which lists the values of x(i) for all i:{0,...,iMax}.
//  This algorithm's running time is log_2( iMax-iMin ).
double BinarySearch_Increasing(int iMin, int iMax, double* xTable, double xTarget)
{
	int iHigh=iMax;
	int iLow=iMin;
	int iMid;
	while (true)
	{
		iMid=(iHigh+iLow)/2;

		//-----------------------------------------------------------------------------
		//	Finish
		//-----------------------------------------------------------------------------
		if( xTable[iMid]==xTarget )
		{
			//	Put both indexes in range and rely on execution of next "if" statement to finish.
			iLow=iHigh=iMid;
		}
		if( iHigh-iLow<5 )
		{
			//	Get iLow to be the lowest index of xTable possibly equal to the xTarget.
			//	Otherwise it should be just below the xTarget.
			while( iLow>iMin && xTable[iLow-1]>=xTarget ) iLow--;
			while( iLow<iMax && xTable[iLow+1]<xTarget ) iLow++;
			if( iLow<iMax && xTable[iLow+1]==xTarget) iLow++;
			//	Get iHigh to be the highest index of xTable possibly equal to the xTarget.
			//	Otherwise it should be just above the xTarget.
			while( iHigh<iMax && xTable[iHigh+1]<=xTarget ) iHigh++;
			while( iHigh>iMin && xTable[iHigh-1]>xTarget ) iHigh--;
			if( iHigh>iMin && xTable[iHigh-1]==xTarget) iHigh--;
			
			if( iHigh==iLow ) // One unique xTarget value.  Return its index.
				return (double)iLow;
			if( xTable[iLow]==xTarget ) // Many unique xTarget values exist.  Return the center of the range.
				return 0.5*(double)(iLow+iHigh);
			if( iHigh-iLow != 1 )
				//	Return linear interpolant.
				return (double)iLow + (xTarget-xTable[iLow])/(xTable[iHigh]-xTable[iLow]);
		}
		//-----------------------------------------------------------------------------

		//  Reduce possible range by 1/2
		if( xTable[iMid] > xTarget )
			iLow=iMid;
		else // Must be < xTarget
			iHigh=iMid;
	}
}

//	If there is a table of [xTable, yTable] with two rows and nTable columns, this uses 1D linear interpolation to estimate
//	yTarget given xTarget.  It is required that the input xTarget be non-decreasing.
//  This algorithm's running time is log_2( iMax-iMin ) because it employs binary search.
double Linterp_Increasing( double* xTable, double* yTable, int nTable, double xTarget )
{
	double i = BinarySearch_Increasing(0, nTable-1, xTable, xTarget);
	double iMod = i-(double)(int)i;
	if( iMod<1e-10 || iMod>0.9999999999 )
		return yTable[(int)(i+0.5)];
	else
		return (1.0-iMod)*yTable[(int)i] + iMod*yTable[1+(int)i];
}