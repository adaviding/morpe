//	Almon Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef null
	#define null 0
#endif
#ifndef M_2PI
	#define M_2PI 6.28318530717959;
#endif

static int* Mcl_RandRotate_Stride = null;
static int Mcl_RandRotate_nStride = 0;
// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void DeleteStaticMemory(void)
{
	delete Mcl_RandRotate_Stride;
}

//================================================================================================================================
//function M = Mcl_RandRotate(n)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function returns a random n * n rotation matrix.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT
//----------------------------------------------------
// n (int32 scalar)
//	The rank of the desired rotation matrix.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT
//----------------------------------------------------
// M (double 2D matrix: n * n)
//	User can supply any rotation matrix (e.g. eye(n)).  On output, the matrix will be altered to 
//================================================================================================================================
void Mcl_RandRotate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs!=1 )
		mexErrMsgTxt("There must be 1 input argument.");
	if( nlhs!=1 )
		mexErrMsgTxt("There must be 1 output argument.");

	//----------------------------------------------------------------------------------------------
	// INPUT
	//---------------------------------------------
	if (!mxIsInt32(prhs[0]) && (int)mxGetNumberOfElements(prhs[0])!=1)
		mexErrMsgTxt("The input argument n must be an int32 scalar.");
	int n = ((int*)mxGetData(prhs[0]))[0];
	if( n<1 )
		mexErrMsgTxt("The input argument n must be a positive integer.");
	//----------------------------------------------------------------------------------------------

	//----------------------------------------------------------------------------------------------
	// OUTPUT
	//---------------------------------------------
	if( plhs[0]==null || !mxIsDouble(plhs[0]) || mxGetM(plhs[0])!=n || mxGetN(plhs[0])!=n )
		plhs[0] = mxCreateNumericMatrix(n,n,mxDOUBLE_CLASS,mxREAL);
	//	Point to output.
	double* M = (double*)mxGetData(plhs[0]);
	//----------------------------------------------------------------------------------------------

	//----------------------------------------------------------------------------------------------
	// STATIC
	//---------------------------------------------
	//	Register the function that must delete static memory
	if( Mcl_RandRotate_nStride==0 )
	{
		//	Initializa random number generator to new seed, but only when the mex function is first opened.
		//	Time is returned to the nearest second, so avoid calling this more than once.
		srand(time(0));
		//	Register the mex function's shut down.
		mexAtExit(DeleteStaticMemory);
	}
	//	Determine whether static memory must be initialized
	if( Mcl_RandRotate_nStride<n)
	{
		delete Mcl_RandRotate_Stride;
		Mcl_RandRotate_nStride = 2*n;
		Mcl_RandRotate_Stride = new int[Mcl_RandRotate_nStride];
	}		
	//----------------------------------------------------------------------------------------------

	double theta;
	double rScale = 1.0/(double)RAND_MAX;
	//	Handle special case
	if( n==1 )
	{
		theta = (double)rand()*rScale;
		if( theta>0.5 )
			M[0] = 1.0;
		else
			M[0] = -1.0;
		return;
	}
	int i,j,k,ii,jj;
	rScale *= M_2PI;
	double c=1.0, s=0.0;
	double z;

	//	Initialize output to identity matrix and fill Mcl_RandRotate_Stride
	for(i=0; i<n; i++)
	{
		//	Fill Mcl_RandRotate_Stride
		Mcl_RandRotate_Stride[i] = n*i;
		//	Initialize output as identity matrix
		for(j=0; j<n; j++)
			M[j+Mcl_RandRotate_Stride[i]]= 0.0;
		M[i+Mcl_RandRotate_Stride[i]]= 1.0;
	}
	
	//	For each random pair (i,j) of rows in the rotation matrix, for all i!=j
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			if(i!=j)
			{
				//	Random rotation angle.
				theta = (double)rand()*rScale; // planar rotation counterclockwise by theta
				c = cos(theta);	// R(i,i) and  R(j,j)
				s = sin(theta);	// R(j,i) and -R(i,j)
				
				//	For each column R(:,i) and R(:,j), multiply each row M(i,:) and M(j,:).  In-place storage (in M) is possible.
				for(k=0; k<n; k++)
				{
					ii = i+Mcl_RandRotate_Stride[k];
					jj = j+Mcl_RandRotate_Stride[k];
					//	X(i,k)	=	R(i,i)*M(i,k)	+	R(i,j)*M(j,k)
					z			=	c*M[ii]			-	s*M[jj];
					//	X(j,k)	=	R(j,i)*M(i,k)	+	R(j,j)*M(j,k)
					M[jj]		=	s*M[ii]			+	c*M[jj];
					//	In place, M=X
					M[ii] = z;
				}
			}
		}
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_RandRotate(nlhs, plhs, nrhs, prhs);
}