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
//	(0):  A non-decreasing function is returned.  The function may may be flat in regions.
//  (1):  The non-decreasing function is interpolated through it's flat spots.
//			Each flat spot receives a proportion of positive derivative energy from immediately adjacent spots which are not flat.
//			The amount of energy received by the flat spots is proportional to the length (number of samples) of each flat spot.
//			If the non-decreasing function is flat everywhere (i.e. totally non-increasing), then a flat function will be returned.
//===========================================================================================================================
void Mcl_ForceMonotonic (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);