using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Morpe.F
{
	public class MonotonicRegressor
	{
		/// <summary>
		/// The type of monotonic regression performed.  Suggest MonotonicRegressionType.Blended.
		/// </summary>
		public MonotonicRegressionType Type = MonotonicRegressionType.Blended;
		/// <summary>
		/// Performs a monotonic regression of a tabulated function.  This method calculates a monotonic function that is approximately equal
		/// to the tabulated function provided.
		/// </summary>
		/// <param name="output">The monotonic function.  This should be supplied as vector of length identical to input.
		/// When the function has finished executing, the following are guaranteed to be true:
		///		mean(output) == mean(input)
		///		min(output) >= min(input)
		///		max(output) <= max(input)
		/// </param>
		/// <param name="input">The tabulated function.</param>
		/// <returns>The number of repetitions through a loop.</returns>
		public int Run(float[] output, float[] input)
		{
			if (input == null || output == null || output.Length < input.Length)
				throw new ArgumentException("The arguments cannot be null, and the length of the output must be at least the length of the input.");
			if (dy == null || dy.Length < input.Length)
			{
				dy = new float[input.Length];
				ynd = new float[input.Length];
			}
			//	Count the number of trips through the loop.
			int nTrips = 0;
			int tripLimit = 100 + (int)(30.0 * Math.Log((float)input.Length));
			int n = input.Length;
			int nm1 = input.Length - 1;
			int i,j;
			//	Compute min, max, mean, etc.
			float yMin=input[0];
			float yMax=input[0];
			float yMean=input[0];
			float dyThis;
			float dyMin=input[1]-input[0];
			output[0]=0.0f;
			for(i=1; i<n; i++)
			{
				if (input[i] < yMin) yMin=input[i];
				if (input[i] > yMax) yMax=input[i];
				yMean += input[i];
				dyThis = input[i]-input[i-1];
				dy[i-1] = dyThis;
				if (dyThis<dyMin)
					dyMin = dyThis;
				output[i]=0.0f;
			}
			yMean/=(float)n;
			float rngY = yMax-yMin;
			float tolY = -rngY*0.001f;
			float halfThis;
			float aCarry;
			bool ndInput = dyMin>=0.0f;

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
					halfThis = dy[i]/2.0f;
					if (halfThis<0.0)
					{
						//	Element y(i+1) must be incremented and y(i) must be decremented.
						output[i+1] -= halfThis; // increment y(i+1)
						output[i] += halfThis; // decrement y(i)
						//	Now, 'dy' must be adjusted to account for this new change.
						if (i>0)
							dy[i-1] += halfThis; // The difference between y(i-1) and y(i) goes down.
						if (i<n-2)
							dy[i+1] += halfThis; // The difference between y(i+1) and y(i+2) goes down.
						dy[i] = 0.0f; // The difference between y(i) and y(i+1) goes to zero.
					}
					halfThis = dy[j]/2.0f;
					if (halfThis<0.0)
					{
						//	Element y(j+1) must be incremented and y(j) must be decremented.
						output[j+1] -= halfThis; // increment y(j+1)
						output[j] += halfThis; // decrement y(j)
						//	Now, 'dy' must be adjusted to account for this new change.
						if (j>0)
							dy[j-1] += halfThis; // The difference between y(j-1) and y(j) goes down.
						if (j<n-2)
							dy[j+1] += halfThis; // The difference between y(j+1) and y(j+2) goes down.
						dy[j] = 0.0f; // The difference between y(j) and y(j+1) goes to zero.
					}
					j--;
					i++;
				}

				//	Update minimum of dy
				dyMin = dy[0];
				for (i=1; i<nm1; i++)
					if (dy[i]<dyMin) dyMin=dy[i];
				//	Increment number of trips.
				nTrips++;
				if (nTrips>=tripLimit)
					break;
			}
			//----------------------------------------------------------------------------------

			float sTarget=0.0f, sBig=0.0f, sSmall=0.0f;

			if (dyMin>=0.0) // out is non-decreasing.
			{
				for(i=0; i<n; i++)
					output[i] += input[i];
			}
			else // out is not yet non-decreasing.  Take extra care.
			{
				//----------------------------------------------------------------------------------
				//	Fix small dy values to be non-negative (so that "out" is non-decreasing).
				//	This operation will bias the output in an undesirable way:  The function will increase too rapidly, so there is a "clean-up" procedure after.
				//----------------------------------------------------------------------------------
				output[0] += input[0];
				aCarry=0.0f;
				for(i=0; i<nm1; )
				{
					if( dy[i++]<0.0 )
						aCarry -= dy[i-1];
					// Set the output
					output[i] += input[i]+aCarry;
					//	Handle built up truncation error
					if( output[i]<output[i-1] )
					{
						aCarry += output[i-1]-output[i];
						output[i]=output[i-1];
					}
				}
				//---------------------------
				//	The prior operation  will have (1) expanded the range of our function and (2) increased the mean of the function.
				//---------------------------
				rngY = output[nm1]-output[0];
				aCarry = (rngY-aCarry)/rngY;  // To counter range expansion
				for(i=0; i<n; i++)
				{
					output[i] = (output[i]-output[0])*aCarry + output[0];  // Counter the expansion of our function's range
					//	Track what happened to the mean of the function
					if( output[i]<yMean )
						sSmall += yMean-output[i];
					else
						sBig += output[i]-yMean;
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
							sTarget = sSmall/(float)n; // The extent to which the actual mean is below the desired mean.
							for(i=0; i<n; i++)
								output[i] += sTarget;
						}
						else
						{
							sTarget = sBig/(float)n; // The extent to which the actual mean is above the desired mean.
							for(i=0; i<n; i++)
								output[i] -= sTarget;
						}
					}
					else
					{
						//	The distance of all values from the mean will be scaled so that sSmall==sBig (i.e. so that the mean is what we intended).
						sTarget = (sBig+sSmall)/2.0f;
						sSmall = sTarget/sSmall;
						sBig = sTarget/sBig;
						//	Ensure that the range is not violated
						aCarry=1.0f;
						if(sSmall>1.0)
							aCarry = (yMean-yMin)/(yMean-output[0])/sSmall;
						else if( sBig>1.0 )
							aCarry = (yMax-yMean)/(output[nm1]-yMean)/sBig;
						if(aCarry<1.0)
						{
							sSmall*=aCarry;
							sBig*=aCarry;
						}
						for(i=0; i<n; i++)
							if (output[i]<yMean)
								output[i] = (output[i]-yMean)*sSmall + yMean;
							else if(output[i]>yMean)
								output[i] = (output[i]-yMean)*sBig + yMean;
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
					output[i] = input[i];
			}
			else
			{
				//	Limit small side
				aCarry=0.0f;
				for(i=0; i<n; i++)
				{
					if(output[i]<yMin)
					{
						aCarry += yMin-output[i];
						output[i]=yMin;
					}
					else if(aCarry>0.0)
					{
						dyThis = output[i]-yMin;
						if(aCarry<dyThis)
						{
							output[i]-=aCarry;
							break;
						}
						else
						{
							output[i]=yMin;
							aCarry-=dyThis;
						}
					}
					else
						break;
				}
				//	Limit big side
				aCarry=0.0f;
				for(i=nm1; i>=0; i--)
				{
					if(output[i]>yMax)
					{
						aCarry += output[i]-yMax;
						output[i]=yMax;
					}
					else if(aCarry>0.0)
					{
						dyThis = yMax-output[i];
						if(aCarry<dyThis)
						{
							output[i]+=aCarry;
							break;
						}
						else
						{
							output[i]=yMax;
							aCarry-=dyThis;
						}
					}
					else
						break;
				}
			}
			//----------------------------------------------------------------------------------

			//	If the blending method is selected, store the output in the static variable Ynd
			int maxLen = 0;
			int len = 0;
			float blendFactor=1.0f;
			if (this.Type==MonotonicRegressionType.Blended)
			{
				ynd[0] = output[0];
				for(i=1; i<n; i++)
				{
					ynd[i] = output[i];
					if( output[i]!=output[i-1] )
						len = 0;
					else
						len++;
					if( len>maxLen )
						maxLen = len;
				}
				blendFactor = (float)(1.0/(1.0+   Math.Pow( (double)(4*maxLen)/(double)(n-1), 3.0 )   ));
			}

			//	Change a non-decreasing function into a monotonic function.
			if (this.Type == MonotonicRegressionType.Increasing || this.Type == MonotonicRegressionType.Blended)
			{
				//	The function is now guaranteed flat and must be converted to a monotonic function.
				//	Compute min, max, and recompute derivative.
				dyMin=output[1]-output[0];
				for(i=1; i<n; i++)
				{
					dyThis = output[i]-output[i-1];
					dy[i-1] = dyThis;
					if (dyThis<dyMin)
						dyMin = dyThis;
				}
				sSmall=(output[nm1]-output[0])*0.00001f/(float)n;
				//	Only continue if the derivative is zero and if the function is not totally flat.
				if (output[nm1]>output[0] && dyMin<=sSmall)
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
					while( iMid<nm1 && dy[iMid]<=sSmall ) iMid++;
					//---------------------------

					//	Advance from iMid's initial value up through iMid = nm-1.
					while( iMid<nm1 )
					{
						//	Advance index of iHigh to the subsequent non-flat index following iMid
						iHigh=iMid+1;
						while( iHigh<nm1 && dy[iHigh]<=sSmall ) iHigh++;

						if( iHigh-iMid>1 )
						{
							if( iMid-iLow>1 )
							{
								//---------------------------
								//	Distribute mass to both sides
								//---------------------------
								//	The mass to each index (center index counts twice and receives double mass).
								dyThis = dy[iMid]/(float)(iHigh-iLow);
								//	Mass to iMid (double mass)
								dy[iMid] = dyThis*2.0f;
								//	Mass below iMid
								for(i=iLow+1; i<iMid; i++)
									dy[i] += dyThis;
								//	Mass above iMid
								for(i=iMid+1; i<iHigh; i++)
									dy[i] += dyThis;
								//---------------------------
							}
							else
							{
								//---------------------------
								//	Distribute mass to high side
								//---------------------------
								//	The mass to each index (center index will get single mass).
								dyThis = dy[iMid]/(float)(iHigh-iMid);
								dy[iMid] = dyThis;
								//	Mass from iMid to iHigh.
								for(i=iMid+1; i<iHigh; i++)
									dy[i] += dyThis;
								//---------------------------
							}
					
						}
						else if( iMid-iLow>1)
						{
							//---------------------------
							//	Distribute mass to low side
							//---------------------------
							//	The mass to each index (center index will get single mass).
							dyThis = dy[iMid]/(iMid-iLow);
							//	Mass to iMid
							dy[iMid] = dyThis;
							//	Mass from iLow to iMid (must be added to mass that may already exist there).
							for(i=iLow+1; i<iMid; i++)
								dy[i] += dyThis;
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
					sSmall=yMean-output[0];
					sBig=0.0f;
					for(i=0; i<nm1; )
					{
						output[i+1] = output[i] + dy[i];
						if(output[++i]>yMean)
							sBig+=output[i]-yMean;
						else
							sSmall+=yMean-output[i];
					}
					//	All values will be shifted by a factor of 0.5*(sBig-sSmall)/(sBig+sSmall)
					sTarget = (sBig+sSmall)/2.0f;
					sSmall = sTarget/sSmall;
					sBig = sTarget/sBig;
					//	Ensure that the range is not violated
					aCarry=1.0f;
					if(sSmall>1.0)
						aCarry = (yMean-yMin)/(yMean-output[0])/sSmall;
					else if( sBig>1.0 )
						aCarry = (yMax-yMean)/(output[nm1]-yMean)/sBig;
					if(aCarry<1.0)
					{
						sSmall*=aCarry;
						sBig*=aCarry;
					}
					//	Adjust elements on big and small sides so that mean is preserved.
					for(i=0; i<n; i++)
						if (output[i]<yMean)
							output[i] = (output[i]-yMean)*sSmall + yMean;
						else if(output[i]>yMean)
							output[i] = (output[i]-yMean)*sBig + yMean;
					//----------------------------------------------------------------------------------
				}
			}

			//	If the output is NaN (which is theoretically possible, but generally should not happen), change it all to be flat (yMean)
			if( Double.IsNaN(output[0]) )
			{
				for(i=0; i<n; i++)
					output[i] = yMean;
			}
			else if( this.Type == MonotonicRegressionType.Blended )
			{
				//	Blend the non-decreasing and monotonic functions according to blendFactor (a 1.0-->0.0 monotonic sigmoid transform of the length of the longest flat portion).
				aCarry = 1.0f-blendFactor;
				for(i=0; i<n; i++)
				{
					output[i] = blendFactor*output[i] + aCarry*ynd[i];
				}
			}
			return nTrips;
		}
		protected float[] dy;
		protected float[] ynd;
	}
}
