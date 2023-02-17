//	Almon David Ing
//	Ctr. Perceptual Systems
//	University of Texas at Austin
//	First compiled July 1, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"

#ifndef null
	#define null 0
#endif

static double *Mcl_ForceMonotonic_Dy=null;
static double *Mcl_ForceMonotonic_Ynd=null;
static int nMcl_ForceMonotonic_Dy=0;

// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void DeleteStaticMemory(void)
{
	delete Mcl_ForceMonotonic_Dy;
	delete Mcl_ForceMonotonic_Ynd;
}

//===========================================================================================================================
//function nTrips = Mcl_ForceMonotonic(Ymono, Y, [IdMethod])
// This function performs monotonic regression on Y.
// Warning:  The values of Ymono and MEM_DY will be changed inside this function.  Any data in these vectors will be destroyed.
//	This function outputs Ymono (a monotonic increasing function approximating Y).
//---------------------------------------------------------------------------------------------------------------------------
// OUTPUT
//----------------------------------------------------------------
// Ymono (double vector: nY)
//	Outputs a non-decreasing 1-d vector the same size as Y and is a monotonic function fit to Y.
//	Guaranteed:
//		mean(Ymono) == mean(Y)
//		min(Ymono) >= min(Y)
//		max(Ymono) <= max(Y)
//
// nTrips (int32 scalar)
//	Outputs the number of trips through a refining loop.
//---------------------------------------------------------------------------------------------------------------------------
// INPUT
//----------------------------------------------------------------
// Y:	(double) A 1-d vector.  This is the signal that will be forced monotonic (on output).  The input Y is not changed.
//
// IdMethod:  An integer specifying the method to be used.  This argument is optional.  If no argument is supplied, default is 1.
//	(0):  A non-decreasing function is returned.  The function may be flat in regions.
//  (1):  The non-decreasing function is interpolated through it's flat spots.
//			Each flat spot receives a proportion of positive derivative energy from immediately adjacent spots which are not flat.
//			The amount of energy received by the flat spots is proportional to the length (number of samples) of each flat spot.
//			If the non-decreasing function is flat everywhere (i.e. totally non-increasing), then a flat function will be returned.
//  (2):  A blending of the first two methods.
//===========================================================================================================================
void Mcl_ForceMonotonic (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i,j;

	//----------------------------------------------------------------------------------------------
	// INPUT
	//---------------------------------------------
	if( nrhs<2 || nrhs>3 )
		mexErrMsgTxt("There must be 2 or 3 input arguments.");
	//---------------------------------------------
	double* y = null;
	if (!mxIsDouble (prhs[1]))
        mexErrMsgTxt("The input vector Y must be type double.");
	int n = (int)(mxGetNumberOfElements(prhs[1]));
	if (n<5)
        mexErrMsgTxt("The input vector Y must have a length greater than 4.");
	y = mxGetPr(prhs[1]);
	//---------------------------------------------
	int IdMethod = 2;
	if (nrhs>=3)
	{
		if (!mxIsInt32 (prhs[2]))
			mexErrMsgTxt ("The input scalar IdMethod must be a 32 bit integer.");
		IdMethod = ((int*)mxGetData(prhs[2]))[0];
		if (IdMethod<0 || IdMethod>2)
			mexErrMsgTxt ("The input scalar IdMethod must belong to the set {0,1,2}.");
	}
	//----------------------------------------------------------------------------------------------
	// STATIC
	//---------------------------------------------
	//	Register the function that must delete static memory
	if( nMcl_ForceMonotonic_Dy==0 )
		mexAtExit(DeleteStaticMemory);
	//	Determine whether static memory must be initialized
	if( nMcl_ForceMonotonic_Dy<n)
	{
		delete Mcl_ForceMonotonic_Dy;
		delete Mcl_ForceMonotonic_Ynd;
		nMcl_ForceMonotonic_Dy = 2*n;
		Mcl_ForceMonotonic_Dy = new double[ nMcl_ForceMonotonic_Dy];
		Mcl_ForceMonotonic_Ynd = new double[nMcl_ForceMonotonic_Dy];
	}
	//----------------------------------------------------------------------------------------------
	// OUTPUT
	//---------------------------------------------
	double* out = null;
	if (!mxIsDouble (prhs[0]))
		mexErrMsgTxt("The argument Ymono must be type double.");
	if( (int)(mxGetNumberOfElements(prhs[0]))<n )
		mexErrMsgTxt("The argument Ymono must have at least as many elements as Y.");
	out = mxGetPr(prhs[0]);
	//---------------------------------------------
	int nTrips = 0;
	//----------------------------------------------------------------------------------------------
	
	//	Limit the number of trips
	int tripLimit = 100 + (int)(30.0*log( (double)n ));
	int nm1 = n-1;

	//	Compute min, max, mean, etc.
	double yMin=y[0];
	double yMax=y[0];
	double yMean=y[0];
	double dyThis;
	double dyMin=y[1]-y[0];
	out[0]=0.0;
	for(i=1; i<n; i++)
	{
		if (y[i] < yMin) yMin=y[i];
		if (y[i] > yMax) yMax=y[i];
		yMean += y[i];
		dyThis = y[i]-y[i-1];
		Mcl_ForceMonotonic_Dy[i-1] = dyThis;
		if (dyThis<dyMin)
			dyMin = dyThis;
		out[i]=0.0;
	}
	yMean/=(double)n;
	double rngY = yMax-yMin;
	double tolY = -rngY*0.001;
	double halfThis;
	double aCarry;
	bool ndInput = dyMin>=0.0;

	//----------------------------------------------------------------------------------
	//	Repeat non-decreasing to become within tolerance.
	//----------------------------------------------------------------------------------
	while (dyMin<tolY)
	{
		i=0;  // Point to beginning of dy
		j=n-2;  // Point to end of dy.
		//	Go through all of dy
		while (j>0)
		{
			halfThis = Mcl_ForceMonotonic_Dy[i]/2.0;
			if (halfThis<0.0)
			{
				//	Element y(i+1) must be incremented and y(i) must be decremented.
				out[i+1] -= halfThis; // increment y(i+1)
				out[i] += halfThis; // decrement y(i)
				//	Now, 'Mcl_ForceMonotonic_Dy' must be adjusted to account for this new change.
				if (i>0)
					Mcl_ForceMonotonic_Dy[i-1] += halfThis; // The difference between y(i-1) and y(i) goes down.
				if (i<n-2)
					Mcl_ForceMonotonic_Dy[i+1] += halfThis; // The difference between y(i+1) and y(i+2) goes down.
				Mcl_ForceMonotonic_Dy[i] = 0.0; // The difference between y(i) and y(i+1) goes to zero.
			}
			halfThis = Mcl_ForceMonotonic_Dy[j]/2.0;
			if (halfThis<0.0)
			{
				//	Element y(j+1) must be incremented and y(j) must be decremented.
				out[j+1] -= halfThis; // increment y(j+1)
				out[j] += halfThis; // decrement y(j)
				//	Now, 'Mcl_ForceMonotonic_Dy' must be adjusted to account for this new change.
				if (j>0)
					Mcl_ForceMonotonic_Dy[j-1] += halfThis; // The difference between y(j-1) and y(j) goes down.
				if (j<n-2)
					Mcl_ForceMonotonic_Dy[j+1] += halfThis; // The difference between y(j+1) and y(j+2) goes down.
				Mcl_ForceMonotonic_Dy[j] = 0.0; // The difference between y(j) and y(j+1) goes to zero.
			}
			j--;
			i++;
		}

		//	Update minimum of Mcl_ForceMonotonic_Dy
		dyMin = Mcl_ForceMonotonic_Dy[0];
		for (i=1; i<nm1; i++)
			if (Mcl_ForceMonotonic_Dy[i]<dyMin) dyMin=Mcl_ForceMonotonic_Dy[i];
		//	Increment number of trips.
		nTrips++;
		if (nTrips>=tripLimit)
			break;
	}
	//----------------------------------------------------------------------------------

	double sTarget=0.0, sBig=0.0, sSmall=0.0;

	if (dyMin>=0.0) // out is non-decreasing.
	{
		for(i=0; i<n; i++)
			out[i] += y[i];
	}
	else // out is not yet non-decreasing.  Take extra care.
	{
		//----------------------------------------------------------------------------------
		//	Fix small Mcl_ForceMonotonic_Dy values to be non-negative (so that "out" is non-decreasing).
		//	This operation will bias the output in an undesirable way:  The function will increase too rapidly, so there is a "clean-up" procedure after.
		//----------------------------------------------------------------------------------
		out[0] += y[0];
		aCarry=0.0;
		for(i=0; i<nm1; )
		{
			if( Mcl_ForceMonotonic_Dy[i++]<0.0 )
				aCarry -= Mcl_ForceMonotonic_Dy[i-1];
			// Set the output
			out[i] += y[i]+aCarry;
			//	Handle built up truncation error
			if( out[i]<out[i-1] )
			{
				aCarry += out[i-1]-out[i];
				out[i]=out[i-1];
			}
		}
		//---------------------------
		//	The prior operation  will have (1) expanded the range of our function and (2) increased the mean of the function.
		//---------------------------
		rngY = out[nm1]-out[0];
		aCarry = (rngY-aCarry)/rngY;  // To counter range expansion
		for(i=0; i<n; i++)
		{
			out[i] = (out[i]-out[0])*aCarry + out[0];  // Counter the expansion of our function's range
			//	Track what happened to the mean of the function
			if( out[i]<yMean )
				sSmall += yMean-out[i];
			else
				sBig += out[i]-yMean;
		}
		//---------------------------
		
		//---------------------------
		//	Counter the increase to the function's mean.
		//---------------------------
		if (sBig!=sSmall) // Signals the mean is off-center
		{
			if( sSmall==0.0 || sBig==0.0 )
			{
				//	This part of the loop should never be encountered, but it's possible (perhaps if the function is totally flat).
				if( sSmall>0.0 )
				{
					sTarget = sSmall/(double)n; // The extent to which the actual mean is below the desired mean.
					for(i=0; i<n; i++)
						out[i] += sTarget;
				}
				else
				{
					sTarget = sBig/(double)n; // The extent to which the actual mean is above the desired mean.
					for(i=0; i<n; i++)
						out[i] -= sTarget;
				}
			}
			else
			{
				//	The distance of all values from the mean will be scaled so that sSmall==sBig (i.e. so that the mean is what we intended).
				sTarget = (sBig+sSmall)/2.0;
				sSmall = sTarget/sSmall;
				sBig = sTarget/sBig;
				//	Ensure that the range is not violated
				aCarry=1.0;
				if(sSmall>1.0)
					aCarry = (yMean-yMin)/(yMean-out[0])/sSmall;
				else if( sBig>1.0 )
					aCarry = (yMax-yMean)/(out[nm1]-yMean)/sBig;
				if(aCarry<1.0)
				{
					sSmall*=aCarry;
					sBig*=aCarry;
				}
				for(i=0; i<n; i++)
					if (out[i]<yMean)
						out[i] = (out[i]-yMean)*sSmall + yMean;
					else if(out[i]>yMean)
						out[i] = (out[i]-yMean)*sBig + yMean;
			}
		}
		//----------------------------------------------------------------------------------
	}

	//----------------------------------------------------------------------------------
	//	Range violation:  At this point, we definitely have a non-decreasing function that has the same mean as the input function,
	//	but the output function could span a greater range (which is undesirable).  We must limit the range of the ouptut function
	//	to the range of the input function.  This can be done without changing the mean.
	//----------------------------------------------------------------------------------
	if( ndInput )
	{
		for(i=0; i<n; i++)
			out[i] = y[i];
	}
	else
	{
		//	Limit small side
		aCarry=0.0;
		for(i=0; i<n; i++)
		{
			if(out[i]<yMin)
			{
				aCarry += yMin-out[i];
				out[i]=yMin;
			}
			else if(aCarry>0.0)
			{
				dyThis = out[i]-yMin;
				if(aCarry<dyThis)
				{
					out[i]-=aCarry;
					break;
				}
				else
				{
					out[i]=yMin;
					aCarry-=dyThis;
				}
			}
			else
				break;
		}
		//	Limit big side
		aCarry=0.0;
		for(i=nm1; i>=0; i--)
		{
			if(out[i]>yMax)
			{
				aCarry += out[i]-yMax;
				out[i]=yMax;
			}
			else if(aCarry>0.0)
			{
				dyThis = yMax-out[i];
				if(aCarry<dyThis)
				{
					out[i]+=aCarry;
					break;
				}
				else
				{
					out[i]=yMax;
					aCarry-=dyThis;
				}
			}
			else
				break;
		}
	}
	//----------------------------------------------------------------------------------

	if( nlhs>0 )
		plhs[0] = mxCreateDoubleScalar( (double)nTrips );

	//	If the blending method is selected, store the output in the static variable Ynd
	int maxLen = 0;
	int len = 0;
	double blendFactor=1.0;
	if (IdMethod>=2)
	{
		Mcl_ForceMonotonic_Ynd[0] = out[0];
		for(i=1; i<n; i++)
		{
			Mcl_ForceMonotonic_Ynd[i] = out[i];
			if( out[i]!=out[i-1] )
				len = 0;
			else
				len++;
			if( len>maxLen )
				maxLen = len;
		}
		blendFactor = 1.0/(1.0+   pow( (double)(4*maxLen)/(double)(n-1), 3.0 )   );
	}

	//mexPrintf("Mcl_ForceMonotonic:  [maxLen, blendFactor] = [%i, %f]\n", maxLen, blendFactor);

	//	Change a non-decreasing function into a monotonic function.
	if (IdMethod>=1)
	{
		//	The function is now guaranteed flat and must be converted to a monotonic function.
		//	Compute min, max, and recompute derivative.
		dyMin=out[1]-out[0];
		for(i=1; i<n; i++)
		{
			dyThis = out[i]-out[i-1];
			Mcl_ForceMonotonic_Dy[i-1] = dyThis;
			if (dyThis<dyMin)
				dyMin = dyThis;
		}
		sSmall=(out[nm1]-out[0])*0.00001/(double)n;
		//	Only continue if the derivative is zero and if the function is not totally flat.
		if (out[nm1]>out[0] && dyMin<=sSmall)
		{
			//----------------------------------------------------------------------------------
			//	Alter the derivative to be positive everywhere, then rebuild the function from this altered derivative.
			//----------------------------------------------------------------------------------

			//---------------------------
			//	Init indexes at a location where [iLow, iMid] straddles the first flat region.
			//---------------------------
			int iLow, iMid=0, iHigh;
			//	iMid:  The index of a non-flat spot
			//	iLow:  The preceding non-flat index before iMid
			//	iHigh: The subsequent non-flat index following iMid
			iLow=-1;
			while( iMid<nm1 && Mcl_ForceMonotonic_Dy[iMid]<=sSmall ) iMid++;
			//---------------------------

			//	Advance from iMid's initial value up through iMid = nm-1.
			while( iMid<nm1 )
			{
				//	Advance index of iHigh to the subsequent non-flat index following iMid
				iHigh=iMid+1;
				while( iHigh<nm1 && Mcl_ForceMonotonic_Dy[iHigh]<=sSmall ) iHigh++;

				if( iHigh-iMid>1 )
				{
					if( iMid-iLow>1 )
					{
						//---------------------------
						//	Distribute mass to both sides
						//---------------------------
						//	The mass to each index (center index counts twice and receives double mass).
						dyThis = Mcl_ForceMonotonic_Dy[iMid]/(double)(iHigh-iLow);
						//	Mass to iMid (double mass)
						Mcl_ForceMonotonic_Dy[iMid] = dyThis*2.0;
						//	Mass below iMid
						for(i=iLow+1; i<iMid; i++)
							Mcl_ForceMonotonic_Dy[i] += dyThis;
						//	Mass above iMid
						for(i=iMid+1; i<iHigh; i++)
							Mcl_ForceMonotonic_Dy[i] += dyThis;
						//---------------------------
					}
					else
					{
						//---------------------------
						//	Distribute mass to high side
						//---------------------------
						//	The mass to each index (center index will get single mass).
						dyThis = Mcl_ForceMonotonic_Dy[iMid]/(double)(iHigh-iMid);
						Mcl_ForceMonotonic_Dy[iMid] = dyThis;
						//	Mass from iMid to iHigh.
						for(i=iMid+1; i<iHigh; i++)
							Mcl_ForceMonotonic_Dy[i] += dyThis;
						//---------------------------
					}
					
				}
				else if( iMid-iLow>1)
				{
					//---------------------------
					//	Distribute mass to low side
					//---------------------------
					//	The mass to each index (center index will get single mass).
					dyThis = Mcl_ForceMonotonic_Dy[iMid]/(iMid-iLow);
					//	Mass to iMid
					Mcl_ForceMonotonic_Dy[iMid] = dyThis;
					//	Mass from iLow to iMid (must be added to mass that may already exist there).
					for(i=iLow+1; i<iMid; i++)
						Mcl_ForceMonotonic_Dy[i] += dyThis;
					//---------------------------
				}
				//	Advance indexes for next trip
				iLow=iMid;		// The "subseqent flat region" on the next trip is set to the "preceding flat region" from the last trip.
				iMid=iHigh;		// We center on last trip's high index.
			}
			//----------------------------------------------------------------------------------

			//----------------------------------------------------------------------------------
			//	Integrate derivative. Ensure the mean of the original input.
			//----------------------------------------------------------------------------------
			sSmall=yMean-out[0];
			sBig=0.0;
			for(i=0; i<nm1; )
			{
				out[i+1] = out[i] + Mcl_ForceMonotonic_Dy[i];
				if(out[++i]>yMean)
					sBig+=out[i]-yMean;
				else
					sSmall+=yMean-out[i];
			}
			//	All values will be shifted by a factor of 0.5*(sBig-sSmall)/(sBig+sSmall)
			sTarget = (sBig+sSmall)/2.0;
			sSmall = sTarget/sSmall;
			sBig = sTarget/sBig;
			//	Ensure that the range is not violated
			aCarry=1.0;
			if(sSmall>1.0)
				aCarry = (yMean-yMin)/(yMean-out[0])/sSmall;
			else if( sBig>1.0 )
				aCarry = (yMax-yMean)/(out[nm1]-yMean)/sBig;
			if(aCarry<1.0)
			{
				sSmall*=aCarry;
				sBig*=aCarry;
			}
			//	Adjust elements on big and small sides so that mean is preserved.
			for(i=0; i<n; i++)
				if (out[i]<yMean)
					out[i] = (out[i]-yMean)*sSmall + yMean;
				else if(out[i]>yMean)
					out[i] = (out[i]-yMean)*sBig + yMean;
			//----------------------------------------------------------------------------------
		}
	}

	//	If the output is NaN (which is theoretically possible, but generally should not happen), change it all to be flat (yMean)
	if( mxIsNaN(out[0]) )
	{
		for(i=0; i<n; i++)
			out[i] = yMean;
	}
	else if( IdMethod>=2 )
	{
		//	Blend the non-decreasing and monotonic functions according to blendFactor (a 1.0-->0.0 monotonic sigmoid transform of the length of the longest flat portion).
		aCarry = 1.0-blendFactor;
		for(i=0; i<n; i++)
		{
			out[i] = blendFactor*out[i] + aCarry*Mcl_ForceMonotonic_Ynd[i];
		}
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_ForceMonotonic(nlhs, plhs, nrhs, prhs);
}