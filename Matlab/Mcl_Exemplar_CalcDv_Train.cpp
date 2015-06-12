//	Almon David Ing
//	Ctr. Perceptual Systems
//	University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "Exemplar_Include.cpp"

#ifndef null
	#define null 0
#endif
//================================================================
//function Mcl_Exemplar_Train_CalcDv(o, wSq)
//----------------------------------------------------------------
// This processes training data in the structure o with respect to the classifier parameters (squared) represented in wSq.
//	Each training sample is reduced to nCats multivariate coordinates.  Only generalization to o.Nneighbors is considered.
//	The pre-computed cube facilitates fast generalization.  For slower generalization (without neighborhood limitation), 
//	use the following code...
//		Mcl_Exemplar_CalcDv(o, o.Dv, o.X, false, sqrt(wSq));
//----------------------------------------------------------------
// o
//	An initialized instance of the exemplar model.
// wSq
//	A 1-d vector of squared weights for each dimension.  The weights scale the distance of each dimension before the
//	kernel is applied.  Therefore, these weights control the size and shape of the blurring (generalization) kernel.
//================================================================
void Mcl_Exemplar_Train_CalcDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<2 )
		mexErrMsgTxt ("There must be two inputs.");

	//------------------------------------------------------------------------------------------------------------------
	//	Read from intput o
	//------------------------------------------------------------------------------------------------------------------
	if( !mxIsStruct(prhs[0]) )
		mexErrMsgTxt ("The argument o must be a structure.");
	mxArray *fld = mxGetField(prhs[0],0,"IdMethod");
	if( fld==null )
		mexErrMsgTxt ("The first argument is not a properly formatted Mcl_Exemplar structure.  A required field name was not found.");
	int IdMethod = ((int*)mxGetData(fld))[0];
	int idGnrlz = IdMethod%10; // 0 decays 1/dist, 1 decays exp, 2 decays gauss
	int idPrms = IdMethod/10;  // 0 categories share parameters, 1 unique params for each category
	bool ForceEqualPriors = ((bool*)mxGetData(mxGetField(prhs[0],0,"ForceEqualPriors")))[0];
	int Ndims = ((int*)mxGetData(mxGetField(prhs[0],0,"Ndims")))[0];
	int NdimsP1 = Ndims+1;
	int Ncats = ((int*)mxGetData(mxGetField(prhs[0],0,"Ncats")))[0];
	int MaxNeigh = ((int*)mxGetData(mxGetField(prhs[0],0,"MaxNeighbors")))[0];
	int Nneigh = ((int*)mxGetData(mxGetField(prhs[0],0,"Nneighbors")))[0];
	int Ntsamp = ((int*)mxGetData(mxGetField(prhs[0],0,"Ntsamp")))[0];
	double* CatWeight = (double*)mxGetData(mxGetField(prhs[0],0,"CatWeight"));
	float* Cube = (float*)mxGetData(mxGetField(prhs[0],0,"Cube"));
	//------------------------------------------------------------------------------------------------------------------

	//	Read input wSq
	if (!mxIsDouble (prhs[1]))
		mexErrMsgTxt ("The argument wSq must be type double.");
	double* wSq = mxGetPr(prhs[1]);

	//------------------------------------------------------------------------------------------------------------------
	//	Read from output o
	//------------------------------------------------------------------------------------------------------------------
	double* Dv = (double*)mxGetData(mxGetField(prhs[0],0,"Dv"));
	//------------------------------------------------------------------------------------------------------------------

	//	Ensure that static memory has been initialized.
	if( idGnrlz>0 )
		Exemplar_EnsureStaticMemory();

	int iCat;
	int iNeigh;
	int iSamp;
	int iDim;
	int iCube;
	double sqdist;

	int* sDimCat = new int[Ncats];
	for(iCat=0; iCat<Ncats; iCat++)
	{
		if( idPrms==0 )
			sDimCat[iCat] = 0;  // Do not stride through weight vector
		else
			sDimCat[iCat] = iCat*Ndims; // Stride through weight vector.
	}

	//	For each datum
	for( iSamp=0; iSamp<Ntsamp; iSamp++ )
	{
		//	Clear out Dv
		for( iCat=0; iCat<Ncats; iCat++ )
			Dv[Ntsamp*iCat+iSamp]=0.0;

		//	For each nearest neighbor
		for( iNeigh=0; iNeigh<Nneigh; iNeigh++ )
		{
			//	Init zero
			sqdist=0.0;
			//	Set strideD
			iCube = (  iNeigh + iSamp*MaxNeigh  ) * NdimsP1;
			//	For each spatial dimension
			iDim=0;
			iCat = (int)(Cube[iCube]+0.49999);
			for(iDim=0;iDim<Ndims;iDim++)
				sqdist += (double)Cube[++iCube]*wSq[iDim+sDimCat[iCat]];
			//	Evaluate without self-reference.
			if( sqdist>0.0 )
			{
				//	Multiply weights assuming equal priors are being forced
				if( ForceEqualPriors)
				{
					//	Accumulate the generalization function a.k.a. blurred kerenel density, augment height to force equal priors.
					if( idGnrlz==0 )
						Dv[Ntsamp*iCat+iSamp] += CatWeight[iCat] / sqrt(sqdist);
					else if( idGnrlz==1)
						Dv[Ntsamp*iCat+iSamp] += CatWeight[iCat] * FastExp_Eval(FastExp_Table, -sqrt(sqdist));
					else if( idGnrlz==2)
						Dv[Ntsamp*iCat+iSamp] += CatWeight[iCat] * FastExp_Eval(FastExp_Table, -sqdist);
				}
				else
				{
						//	Accumulate the generalization function a.k.a. blurred kerenel density.
					if( idGnrlz==0 )
						Dv[Ntsamp*iCat+iSamp] += 1.0 / sqrt(sqdist);
					else if( idGnrlz==1)
						Dv[Ntsamp*iCat+iSamp] += FastExp_Eval(FastExp_Table, -sqrt(sqdist));
					else if( idGnrlz==2)
						Dv[Ntsamp*iCat+iSamp] += FastExp_Eval(FastExp_Table, -sqdist);
				}
			}
		}
		//	Normalize the generalization for that sample (rather than declaring a new variable, just use sqdist)
		sqdist=0.0;
		for( iCat=0; iCat<Ncats; iCat++ )
			sqdist += Dv[Ntsamp*iCat+iSamp];
		if( sqdist>0.0 )
		{
			for( iCat=0; iCat<Ncats; iCat++ )
				Dv[Ntsamp*iCat+iSamp] /= sqdist;
		}
		else
		{
			//	It is unlikely (but possible) that this block will be executed.
			sqdist = 1.0/(double)Ncats;
			if( ForceEqualPriors )
				for( iCat=0; iCat<Ncats; iCat++ )
					Dv[Ntsamp*iCat+iSamp] = sqdist;
			else
				for( iCat=0; iCat<Ncats; iCat++ )
					Dv[Ntsamp*iCat+iSamp] = sqdist/CatWeight[iCat];
		}
	}
	//	Delete allocated memory
	delete sDimCat;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_Exemplar_Train_CalcDv(nlhs, plhs, nrhs, prhs);
}