//	Almon David Ing.
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Written July 22, 2009

//#include <limits>
#include <math.h>
#include "mex.h"

#ifndef NULL
	#define NULL 0
#endif

//	Creates the table necessary for evaluating FastExp_Eval.
double* FastExp_CreateTable()
{
	double* table = new double[6826];
	int i;
	for(i=0; i<=1000; i++)
		table[i] = exp( (double)i*0.001 );
	for(i=0; i<=900; i++)
	{
		table[i+1001] = exp( double(i+100)*0.01 );
		table[i+1902] = exp( double(i+100)*0.1 );
	}
	for(i=0; i<=609; i++)
		table[i+2803] = exp( (double)(i+100) );
	for(i=0; i<3413; i++)
		table[i+3413] = 1.0/table[i];
	return table;
}
//	Fast evaluation of Exp using linear interpolation through a pre-tabled exp values.  
double FastExp_Eval(double* Table, double z)
{
	double wHigh;
	int iLow;

	if( Table==NULL )
		return exp(z);

	if(z<0.0)
	{
		if(z>-10.0)
		{
			if(z>-1.0)
			{
				// z is [0, 1.0)
				z*=-1000.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[3413+iLow]*(1.0-wHigh) + Table[3414+iLow]*wHigh;
			}
			else
			{
				// z is [1.0, 10.0)
				z = (-z-1.0)*100.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[4414+iLow]*(1.0-wHigh) + Table[4415+iLow]*wHigh;
			}
		}
		else
		{
			if(z<=-709.0)
				return 0.0;
			if(z>-100.0)
			{
				// z is [10.0, 100.0)
				z = (-z-10.0)*10.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[5315+iLow]*(1.0-wHigh) + Table[5316+iLow]*wHigh;
			}
			else
			{
				//	z is [100.0, 709.0)
				z = -z-100.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[6216+iLow]*(1.0-wHigh) + Table[6217+iLow]*wHigh;
			}
		}
	}
	else if(z>0.0)
	{
		if(z<10.0)
		{
			if(z<1.0)
			{
				// z is [0, 1.0)
				z*=1000.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[iLow]*(1.0-wHigh) + Table[1+iLow]*wHigh;
			}
			else
			{
				// z is [1.0, 10.0)
				z = (z-1.0)*100.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[1001+iLow]*(1.0-wHigh) + Table[1002+iLow]*wHigh;
			}
		}
		else
		{
			if(z>=709.0)
				return mxGetInf();
			if(z<100.0)
			{
				// z is [10.0, 100.0)
				z = (z-10.0)*10.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[1902+iLow]*(1.0-wHigh) + Table[1903+iLow]*wHigh;
			}
			else
			{
				//	z is [100.0, 709.0)
				z -= 100.0;
				iLow = (int)z;
				wHigh = fmod(z,1.0);
				return Table[2803+iLow]*(1.0-wHigh) + Table[2804+iLow]*wHigh;
			}
		}
	}
	return 1.0;
}