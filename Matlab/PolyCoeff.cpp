//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Written July 20, 2009
//
//	This code is useful for constructing or defining the parameters of a multivariate polynomial.
//	A polynomial is defined by a set of parameters.  If the polynomial is defined over an Ndims-dimensional space, then the
//		number of polynomial parameters is always Ndims (for a linear polynomial) or greater (for more complex polynomials like
//		quadratic, cubic, etc...).  In addition, a polynomial always has one extra parameter (an additive constant), so the
//		parameter space (i.e. the space of polynomial coefficients) always has a higher dimensionality than Ndims.  A polynomial
//		is a function defined by two things: (1) It's rank and (2) the vector of polynomial coefficients.
//
//	The number of polynomial coefficients is a function of its Rank.  For example, if the polynomial coefficients are indexes
//		of "a" and the spatial dimensions are indexes of "x", then the polynomial is a function like:
//	Rank 0 = Flat		-->  E.g. When Ndims=?, f(x) = a0
//	Rank 1 = Linear		-->  E.g. When Ndims=2, f(x) = a0 + a1*x1 + a2*x2
//	Rank 2 = Quadratic  -->  E.g. When Ndims=2, f(x) = a0 + a1*x1 + a2*x2 + a11*x1*x1 + a22*x2*x2 + a12*x1*x2
//	Rank 3 = Cubic		-->  E.g. When Ndims=2, f(x) = a0 + a1*x1 + a2*x2 + a11*x1*x1 + a22*x2*x2 + a12*x1*x2 + a111*x1*x1*x1 + a112*x1*x1*x2 + a122*x1*x2*x2 + a222*x2*x2*x2
//		etc...												
//	Note that Flat (Rank=0) and Linear (Rank=1) are typically not called polynomials, but they are special cases of polynomials
//		and their functionality is encapsulated by Mcl polynomial methods.
//
//	The number polynomial coefficients for a rank Rank assuuming an Ndims dimensional space is entry (Ndims,Rank) of the symmetric Pascal Matrix.
//		This is also equal to "n choose k" or n!/k!/(n-k)! where n = Ndims+Rank-2 and k = Ndims-1

//#include "mex.h"

#ifndef null
	#define null 0
#endif

//	Returns the entry (a,b) of the Pascal matrix (which is symmetric).  Output is equal to (a+b)! / a / b
int Pascal(int a, int b)
{
	int aa = a;
	int bb = b;
	if(a<b)
	{
		aa = b;  // aa = max(a,b)
		bb = a;  // bb = min(a,b)
	}
	int c = a+b;
	int output = 1;
	int i;
	for (i = aa + 1; i <= c; i++)
		output *= i;
	for (i = 2; i <= bb; i++)
		output /= i;
	return output;
}

//	Returns the total number of polynomial terms (coefficients) not including the additive constant for a polynomial of rank
//	Rank defined over an Ndims dimensional space.  If N is not null, it must point to a pre-allocated vector of length Rank.
//	On output, each N[i] will contain the number of (i+1)-order terms for a polynomial defined over an Ndims-dimensional space.
//		e.g. N[0]  Returns the number of first order ("linear") terms when Rank>=1 (this is always equal to Ndims).
//		e.g. N[1]  Returns the number of second order ("quadratic") terms when Rank>=2
//		e.g. N[2]  Returns the number of third order ("cubic") terms when Rank>=3
//		etc...
//	This function's return value is equal to the sum of elements in N.
int Ncoeff(int* N, int Rank, int Ndims)
{	
	//	Indexors
	int i;
	//	If 0-dimensionality, set everything to empty.
	if( Ndims<=0 )
	{
		if( N!=null )
		{
			for(i=0; i<Rank; i++)
				N[i] = 0;
		}
		return 0;
	}
	//	If there is no rank, return 0.
	if (Rank<=0)
		return 0;
	//	Return linear terms.
	if( Rank==1)
	{
		//	The number of linear terms is always Ndims
		if (N != null)
			N[0] = Ndims;
		return Ndims;
	}
	//	Compute the number of inhomogeneous polynomial terms for the given rank and dimensionality.
	int output = Pascal(Rank,Ndims);
	
	//	Subtract 1 to remove the 0-th order term.
	output--;
	
	//mexPrintf("Pascal(%i,$i)=%i\n",Rank,Ndims,output);
	
	//	Fill output vector N and sum to compute output.
	if (N != null)
	{
		N[0] = Ndims;  //	There are Ndims linear terms
		//mexPrintf("N[0]=%i\n", Ndims);
		for( i=1; i<Rank; i++)
		{
			N[i] = Pascal(i+1,Ndims-1);
			//mexPrintf("N[%i]=%i\n", i, N[i]);
		}
	}

	//	Return value.
	return output;
}

// Enumerates the components of the polynomial expansion associated with each coefficient of the polynomial.  User must provide an
//	allocated 2D (Ncoeff * Rank) array CoeffDims to contain the output.  Basically, this function outputs the "definition" of a
//	polynomial defined by Rank and Ndims.
// INPUT
// Ncoeff (int32 scalar)
//	The number of coefficients of the polynomial.  User can call the function Ncoeff(null,Rank,Ndims) to obtain this value.
// Rank (int32 scalar)
//	The polynomial's rank.
//		If Rank=1, the polynomial will only contain first-order  (linear) terms.
//		If Rank=2, the polynomial will also contain second-order (quadratic) terms.
//		If Rank=3, ...								third-order  (cubic)
//	...
// Ndims (int32 scalar)
//	The dimensionality of the space over which the polynomial is defined.
// OUTPUT
//	CoeffDims (int32 2D array, Ncoeff * Rank) indexed as:
//		C[iCoeff + Ncoeff*iComp]
//	where
//		iComp == 0 first  component of the polynomial expansion
//		iComp == 1 second component of the polynomial expansion
//		iComp == 2 third  component of the polynomial expansion
//		etc...
//	where
//		iCoeff is bound to the range:  [0, Ncoeff(null,Rank,Ndims)]
//	CoeffDims enumerates the components of polynomial expansion for each coefficient.  Each (iCoeff,iComp)-th element of CoeffDims
//	can be an integer {0,...,Ndims-1} which identifies an axis of an Ndims dimensional space; OR can be -1 when the entry is "empty"
//	(i.e. when no axis is indicated).   All -1 ("empty") values are packed to right-most columns so the left-most columns will
//	always contain the significant entries (integers {0,...,Ndims-1}).  The first column (associated with iComp==0) will never
//	contain empty (-1) values.
void EnumCoeffs(int* CoeffDims, int Ncoeff, int Rank, int Ndims)
{
	int i,j;
	//	Init digits to first term of linear tensor.
	int* ind = new int[Rank];
	for(i=0; i<Rank; ) ind[i++] = -1;
	ind[0]=0;
	//	For each coefficient
	int iCoeff=0;
	//mexPrintf("EnumCoeffs\n");
	while(true)
	{

		//	Copy digits to output
		for(i=0; i<Rank; i++)
			CoeffDims[iCoeff + Ncoeff*i] = ind[i];

		//mexPrintf("EnumCoeffs iCoeff=%i, Ncoeff=%i, CoeffDims(iCoeff+1,:)=[%i", iCoeff,Ncoeff,CoeffDims[iCoeff]);
		//for(i=1; i<Rank; i++)
		//	mexPrintf(", %i", CoeffDims[iCoeff + Ncoeff*i]);
		//mexPrintf("]\n");

		if (++iCoeff >= Ncoeff)
			break;

		//	Increment the smallest digit and begin to handle "carry-over" arithmetic.
		i=-1;
		while( ++ind[++i]==Ndims );

		//	Finish the "carry-over" by ensuring that leftward digits have been cleared from Ndims.
		for( j=i-1; j>=0; j-- )
		{
			if( ind[j]==Ndims )
				ind[j]=ind[j+1];
		}
	}

	delete ind;
}