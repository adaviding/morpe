//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "Exemplar_Include.cpp"

#ifndef null
	#define null 0
#endif

//================================================================================================================================
//function Mcl_Exemplar_CalcDv(o, Dv, Xtest, [AllowSame], [wOptimized])
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function uses a trained Mcl_Exemplar solver to classify new (untrained) data.  This function reduces each multivariate
//	stimulus coordinate to Ncats decision values.  Use the Mcl_MapDv function to map these values to probabilities.
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting subsequent comments)
//----------------------------------------------------
// Ncats = The number of categories of the training data (i.e. the length of the input Xtest cell array)
// Ntsamp = The total number of training samples = sum(nSamp);
// Ndims = The spatial dimension of the data (i.e. the number of columns of Xtest).
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not altered by mex function)
//----------------------------------------------------
// o 
//	The exemplar model specified by Mcl_Exemplar_Ctor.
// Xtest (double 2D array: nTest * Ndims)
//	The spatial coordinate of each testing sample.  Note:  If the solver internalizes transforms on each dimension, Xtest must
//	be already transformed.
// [AllowSame] (logical bool scalar) defaults to false
//	If false, then a training coordinate which is exactly equal to a testing coordinate (i.e. the same exact spatial coordinate)
//	will not be used in the blurring or generalization operation.  This is automatically false when mod(o.IdMethod,10)==0.
//	If true, an exactly equal training coordinate from the training set (if it exists) can be used.
// [wOptimized] (double vector: Ndims) defaults to o.wOptimized
//	If provided, this vector provides the optimized parameters of the exemplar model.  Each parameter scales a dimension.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (values are copied to pre-allocated arrays provided).
//----------------------------------------------------
// Dv (double 2D array: nTest * Ncats)
//	For each test coordinate, the decision variable associated with each category is provided.
//================================================================================================================================
void Mcl_Exemplar_CalcDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<3 )
		mexErrMsgTxt("Not enough input arguments.");
	
	if (!mxIsDouble(prhs[1]))
		mexErrMsgTxt("The output argument Dv must be type double.");

	//-----------------------------------------------------------------------------------
	// INPUT from structure
	//-----------------------------------------------------------------------------------
	if (!mxIsStruct(prhs[0]))
		mexErrMsgTxt("The input argument o must be a matlab structure.");
	mxArray *fld = mxGetField(prhs[0],0,"IdMethod");
	if( fld==null )
		mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Exemplar structure.  A required field name was not found.");
	int IdMethod = ((int*)mxGetData(fld))[0];
	//------------------------------------------------
	int idGnrlz = IdMethod%10; // 0 decays 1/dist, 1 decays exp, 2 decays gauss
	int idPrms = IdMethod/10;  // 0 categories share parameters, 1 unique params for each category
	int Ncats = ((int*)mxGetData(mxGetField(prhs[0],0,"Ncats")))[0];
	int Ntsamp = ((int*)mxGetData(mxGetField(prhs[0],0,"Ntsamp")))[0];
	int Ndims = ((int*)mxGetData(mxGetField(prhs[0],0,"Ndims")))[0];
	bool ForceEqualPriors = ((bool*)mxGetData(mxGetField(prhs[0],0,"ForceEqualPriors")))[0];
	double* CatWeight = (double*)mxGetData(mxGetField(prhs[0],0,"CatWeight"));
	double* Xtrain = (double*)mxGetData(mxGetField(prhs[0],0,"X"));
	int* Cat = (int*)mxGetData(mxGetField(prhs[0],0,"Cat"));
	double* wOpt = (double*)mxGetData(mxGetField(prhs[0],0,"wOptimized"));
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//  INPUT from Xtest
	//-----------------------------------------------------------------------------------
	if ( !mxIsDouble(prhs[2]) || Ndims>1 && (mxGetNumberOfDimensions(prhs[2])!=2 || Ndims!=mxGetDimensions(prhs[2])[1]) )
		mexErrMsgTxt("The input argument Xtest must be a 2D double array: nTest * Ndims.");
	int nTest = ((int*)mxGetDimensions(prhs[2]))[0];
	double* Xtest = (double*)mxGetData(prhs[2]);
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//  INPUT from AllowSame
	//-----------------------------------------------------------------------------------
	bool AllowSame = false;
	if( nrhs>3 && mxGetNumberOfElements(prhs[3])>0 )
	{
		if (!mxIsLogical(prhs[3]))
			mexErrMsgTxt("The input argument AllowSame must be a logical scalar.");
		AllowSame = ((bool*)mxGetData(prhs[3]))[0];
	}
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//  INPUT from wOptimized
	//-----------------------------------------------------------------------------------
	if( nrhs>4 && mxGetNumberOfElements(prhs[4])>0 )
	{
		if (!mxIsDouble(prhs[4]))
			mexErrMsgTxt("The input argument wOptimized must be type double.");
		if( idPrms==0 && (int)mxGetNumberOfElements(prhs[4])<Ndims )
			mexErrMsgTxt("The input argument wOptimized must be a double vector of length Ndims.");
		if( idPrms==1 && (int)mxGetNumberOfElements(prhs[4])<(Ndims * Ncats) )
			mexErrMsgTxt("The input argument wOptimized must be a double vector (or 2D array) of length Ndims * Ncats.");
		wOpt = mxGetPr(prhs[4]);
	}
	//-----------------------------------------------------------------------------------

	//-----------------------------------------------------------------------------------
	//  OUTPUT to Dv
	//-----------------------------------------------------------------------------------
	if( (int)mxGetNumberOfElements(prhs[1])<nTest*Ncats )
		mexErrMsgTxt("The output argument Dv must be pre-allocated to a 2D double array: nTest * Ncats.");
	if( Ncats!=((int*)mxGetDimensions(prhs[1]))[1] )
		mexErrMsgTxt("The output argument Dv must be pre-allocated to a 2D double array: nTest * Ncats.");
	if( nTest!=((int*)mxGetDimensions(prhs[1]))[0] )
		mexErrMsgTxt("The output argument Dv must be pre-allocated to a 2D double array: nTest * Ncats.");
	double* Dv = (double*)mxGetData(prhs[1]);
	//-----------------------------------------------------------------------------------

	int iSamp,iTest,iDim,iCat;

	//-----------------------------------------------------------------------------------
	//	Pre-calculate strides for faster indexing.
	//-----------------------------------------------------------------------------------
	int* sDimTest = new int[Ndims];
	int* sDimSamp = new int[Ndims];
	for(iDim=0; iDim<Ndims; iDim++)
	{
		sDimSamp[iDim] = iDim*Ntsamp;
		sDimTest[iDim] = iDim*nTest;
	}
	int* sCatDv = new int[Ncats];
	int* sDimCat = new int[Ncats];
	for(iCat=0; iCat<Ncats; iCat++)
	{
		sCatDv[iCat] = iCat*nTest;
		if( idPrms==0 )
			sDimCat[iCat] = 0;  // Do not stride through weight vector
		else
			sDimCat[iCat] = iCat*Ndims; // Stride through weight vector.
	}
	//-----------------------------------------------------------------------------------

	//	Ensure that static memory has been initialized.
	if( idGnrlz>0 )
		Exemplar_EnsureStaticMemory();

	double dx;
	double sqdist;

	//	Generalize (blur) each iTest test coordinate to every training sample.
	for(iTest=0; iTest<nTest; iTest++)
	{
		//	Init Dv for this row to zero (Dv sums generaliztion or blurring).
		for(iCat=0; iCat<Ncats; iCat++)
			Dv[iTest+sCatDv[iCat]]=0.0;
		//-----------------------------------------------------------------------------------
		//	Generalize this test coordinate to all training samples.
		//-----------------------------------------------------------------------------------
		for(iSamp=0; iSamp<Ntsamp; iSamp++)
		{
			//	Category id
			iCat = Cat[iSamp];
			//	Init zero distance.
			sqdist=0.0;
			for(iDim=0; iDim<Ndims; iDim++)
			{
				dx = wOpt[iDim+sDimCat[iCat]]*( Xtest[iTest+sDimTest[iDim]] - Xtrain[iSamp+sDimSamp[iDim]] );
				sqdist += dx*dx;
			}
			//	Ensure this test coordinate is not identical to the training sample.
			if(AllowSame || sqdist>0.0)
			{
				if( ForceEqualPriors )
				{
					if( idGnrlz==0 )
						Dv[iTest+sCatDv[iCat]] += CatWeight[iCat]/sqrt(sqdist);
					else if( idGnrlz==1 )
						Dv[iTest+sCatDv[iCat]] += CatWeight[iCat]*FastExp_Eval(FastExp_Table, -sqrt(sqdist));
					else if( idGnrlz==2 )
						Dv[iTest+sCatDv[iCat]] += CatWeight[iCat]*FastExp_Eval(FastExp_Table, -sqdist);
				}
				else
				{
					if( idGnrlz==0 )
						Dv[iTest+sCatDv[iCat]] += 1.0/sqrt(sqdist);
					else if( idGnrlz==1 )
						Dv[iTest+sCatDv[iCat]] += FastExp_Eval(FastExp_Table, -sqrt(sqdist));
					else if( idGnrlz==2 )
						Dv[iTest+sCatDv[iCat]] += FastExp_Eval(FastExp_Table, -sqdist);
				}
			}
		}
		//-----------------------------------------------------------------------------------

		//-----------------------------------------------------------------------------------
		//	Normalize across dims, use sqdist because it is not being used anymore.
		//-----------------------------------------------------------------------------------
		sqdist=0.0;
		for(iCat=0; iCat<Ncats; iCat++)
			sqdist += Dv[iTest+sCatDv[iCat]];
		
		if( sqdist>0.0 )
		{
			for(iCat=0; iCat<Ncats; iCat++)
				Dv[iTest+sCatDv[iCat]] /= sqdist;
		}
		else
		{
			//	It is unlikely (but possible) that this block will be executed.
			sqdist = 1.0/(double)Ncats;
			if( ForceEqualPriors )
				for( iCat=0; iCat<Ncats; iCat++ )
					Dv[iTest+sCatDv[iCat]] = sqdist;
			else
				for( iCat=0; iCat<Ncats; iCat++ )
					Dv[iTest+sCatDv[iCat]] = sqdist/CatWeight[iCat];
		}
		//-----------------------------------------------------------------------------------
	}

	//	Delete allocated memory.
	delete sDimTest;
	delete sDimSamp;
	delete sCatDv;
	delete sDimCat;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Exemplar_CalcDv(nlhs, plhs, nrhs, prhs);
}