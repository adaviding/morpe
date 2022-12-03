//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "memory.h"
#include "QuickSort.cpp"

#ifndef null
	#define null 0
#endif
//#ifndef DEBUG_MEXPRINTF
//	#define DEBUG_MEXPRINTF
//#endif

static int* QuantizeDv_Idx=null;
static int nQuantizeDv_Idx=0;
// This function is registered in mexAtExit() when static memory is used.  It is called when Matlab exits or when a user types
//	"clear" or "clear mex" at the command prompt.  This is necessary to clean up static memory resources.
static void QuantizeDv_DeleteStaticMemory(void)
{
	delete QuantizeDv_Idx;
}

//================================================================================================================================
//function Mcl_QuantizeDv(o)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function assists Mcl_Exemplar and Mcl_Poly methods by quantizing decision values.
//  A "decision value" is a number that can be computed for each training sample using a decision function associated with a particular
//	category.  Each category must be associated with a unique decision function, and each of these decision function is assumed (but not
//	required) to have an isotonic relationship with the posterior probability of memberhsip in the category.
//
//	This procedure makes it possible to determine the posterior probability of category membership given the evaluated decision function
//	using a unique monotonic mapping function for each category.  The output of this procedure provides an interpolation table that
//	defines the mapping of each category's decision value with the probability of membership in that category.
//
//  If you have only one decision function for two categories, then you can always mimic a decision function for the second category by
//	transforming the output of the decision function for the first category.  This transform can usually (but now always) be accomplished
//	by taking the negative, inverse, or 1-complement of the first decision function.  For example, a linear model can be defined by the
//	following linear decision function where "x" is the vector coordinate of a training sample, "B" is a vector defining the slope of a
//	linear decision boundary, and "B0" is a scalar that defines the boundary's intercept.
//		f( category "1" , x ) =   B * x + B0
//	In this example, the linear decision function for the second category is the negative of the same function for the first category.
//		f( category "2" , x ) = -(B * x + B0)
//--------------------------------------------------------------------------------------------------------------------------------
// NOMENCLATURE (for interpreting these comments)
//----------------------------------------------------
// Ntsamp = The total number of training samples.
// Ncats = The number of categories.
// Nquant = The number of quantiles (bins). 
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not changed by this mex function)
//----------------------------------------------------
// o (mxArray*, Matlab structure)
//	A structure constructed from Mcl_Exemplar_Ctor (refer to its comments) and then run at least once through Mcl_Exemplar_Eval.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (pre-allocated memory will be filled with the output of this function.)
//----------------------------------------------------
// o.Quant(iCat) (mxArray*, Matlab structure)
//	The quantization table for each iCat category is mostly completed by this function.  However, o.PcMono and o.PcMonoLim are not
//		filled by this function.
//--------------------------------------------------------------------------------------------------------------------------------
// After this function terminates, a decision value can be mapped to a posterior probability of assigning the iCat category label by
//	using o.Quant(iCat).Dv and o.Quant(iCat).Pc as an interpolation table.  Since probabilities are computed separately for each
//	category (using separate decision functions), there is no guarantee that the probabilities will sum to 1.0 when a stimulus is
//	classified.  Therefore, subsequent methods should always normalize the probabilities (making them sum to 1.0).
//
// Since posterior probabilities should never reach the boundaries [0.0, 1.0], it is recommended that the user set boundaries for
//	o.Quant(iCat).Pc(iBin) for each iBin quantile to be:
//		[1/Ncats/o.Quant.Weight(iBin), 1.0-(Ncats-1)/Ncats/o.Quant.Weight(iBin)]
//	This algorithm does not do this automatically, therefore o.Quant(iCat).Pc(iBin) can reach the boundaries of [0.0, 1.0].
//
// Since quantization was used to compute the interpolation table, there will be quantization (binning) noise in o.Quant(iCat).Dv
//	and o.Quant(iCat).Pc.
//================================================================================================================================
void Mcl_QuantizeDv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<1 )
		mexErrMsgTxt("Not enough input arguments.");
	if( !mxIsStruct(prhs[0]) )
		mexErrMsgTxt("The first argument of this function (o) must be a matlab structure.");
	//--------------------------------------------------------
	mxArray* Q = mxGetField(prhs[0], 0, "Quant");
	if( Q==null )
		mexErrMsgTxt("o.Quant must be created before calling this function.");
	if( !mxIsStruct(Q) )
		mexErrMsgTxt("o.Quant must be a matlab structure.");
	//--------------------------------------------------------
	int Ncats = ((int*)mxGetData(mxGetField(prhs[0], 0, "Ncats")))[0];
	int Ntsamp = ((int*)mxGetData(mxGetField(prhs[0], 0, "Ntsamp")))[0];
	//--------------------------------------------------------
	mxArray* dum = mxGetField(prhs[0], 0, "Dv");
	if( dum==null )
		mexErrMsgTxt("o.Dv must be created before calling this function.");
	if( !mxIsDouble(dum) || mxGetNumberOfDimensions(dum)<2)
		mexErrMsgTxt("o.Dv must be allocated as a two-dimensional array of type double.");
	double* Dv = mxGetPr(dum);
	//--------------------------------------------------------
	dum = mxGetField(prhs[0], 0, "Cat");
	if( dum==null )
		mexErrMsgTxt("o.Cat must be created before calling this function.");
	if( !mxIsInt32(dum) || (int)mxGetNumberOfElements(dum)<Ntsamp)
		mexErrMsgTxt("o.Cat must be allocated as a vector of type int32 and of length Ntsamp where Ntsamp is the number of rows of o.Dv.");
	int* Cat = (int*)mxGetData(dum);
	//--------------------------------------------------------
	dum = mxGetField(prhs[0], 0, "CatWeight");
	if( dum==null )
		mexErrMsgTxt("o.CatWeight must be created before calling this function.");
	if( !mxIsDouble(dum) || (int)mxGetNumberOfElements(dum)<Ncats)
		mexErrMsgTxt("o.CatWeight must be allocated as a vector of type double and of length Ncats where Ncats is the number of columns of o.Dv.");
	double* CatWeight = mxGetPr(dum);
	//----------------------------------------------------------------------------------------------
	// STATIC
	//--------------------------------------------------------
	//	Register the function that must delete static memory
	if( nQuantizeDv_Idx==0 )
		mexAtExit(QuantizeDv_DeleteStaticMemory);
	//	Determine whether static memory must be initialized
	if( nQuantizeDv_Idx<Ntsamp )
	{
		delete QuantizeDv_Idx;
		nQuantizeDv_Idx = Ntsamp;
		QuantizeDv_Idx = new int[nQuantizeDv_Idx];
	}
	//--------------------------------------------------------

	double *QdvSep, *Qdv, *Qpc, *Qweight;
	int iCat, jCat, iSamp, iEntry, jEntry, Nquant, iBin, iSorted;

	//--------------------------------------------------------
	//	Initialize the weights and check to see that the category labels supplied by the user are valid.
	//--------------------------------------------------------
	double wTotalBins = 0.0;
	for(iSamp=0; iSamp<Ntsamp; iSamp++)
	{
		if( Cat[iSamp]<0 || Cat[iSamp]>=Ncats )
			mexErrMsgTxt ("All elements of o.Cat must be integers in the set {0,...,Ncats-1}.");
		wTotalBins+=CatWeight[Cat[iSamp]];
	}
	//--------------------------------------------------------

	//	Keep track of the cumulative weight 
	double wPerBin,wNextBin,w,wLastEntry,wLastBin,dwThis,pcBin,dvBin;

	//	Quantize decision values for each jCat column of Dv
	for( jCat=0; jCat<Ncats; jCat++ )
	{
		//--------------------------------------------------------
		// OUTPUT Each jCat decision function has separate
		//--------------------------------------------------------
		dum = mxGetField(Q, jCat, "Nquantiles");
		if( !mxIsInt32(dum) )
			mexErrMsgTxt("Nquantiles must be type Int32.");
		Nquant = ((int*)mxGetData(dum))[0];
		if( Nquant < Ncats )
			mexErrMsgTxt("There must be at least Ncats quantiles.");
		//--------------------------------------------------------
		dum = mxGetField(Q, jCat, "Dv");
		if( dum==null )
			mexErrMsgTxt("o.Quant(iCat).Dv must be created before calling this function.");
		if( !mxIsDouble(dum) )
			mexErrMsgTxt("o.Quant(iCat).Dv must be allocated as a vector of type double.");
		Qdv = mxGetPr(dum);
		//--------------------------------------------------------
		dum = mxGetField(Q, jCat, "Pc");
		if( dum==null )
			mexErrMsgTxt("o.Quant(iCat).Pc must be created before calling this function.");
		if( !mxIsDouble(dum) )
			mexErrMsgTxt("o.Quant(iCat).Pc must be allocated as a vector of type double.");
		if( (int)mxGetNumberOfElements(dum) < Nquant )
			mexErrMsgTxt("The length of o.Quant(iCat).Pc must be at least Nquant.  In other words, the length of o.Quant(iCat).Pc must be at least the length of o.Quant(iCat).Dv.");
		Qpc = mxGetPr(dum);
		//--------------------------------------------------------
		dum = mxGetField(Q, jCat, "Weight");
		if( dum==null )
			mexErrMsgTxt("o.Quant(iCat).Weight must be created before calling this function.");
		if( !mxIsDouble(dum) )
			mexErrMsgTxt("o.Quant(iCat).Weight must be allocated as a vector of type double.");
		if( (int)mxGetNumberOfElements(dum) < Nquant )
			mexErrMsgTxt("The length of o.Quant(iCat).Weight must be at least Nquant.  In other words, the length of o.Quant(iCat).Weight must be at least the length of o.Quant(iCat).Dv.");
		Qweight = mxGetPr(dum);
		//--------------------------------------------------------
		dum = mxGetField(Q, jCat, "DvBinSep");
		if( dum==null )
			mexErrMsgTxt("o.Quant(iCat).DvBinSep must be created before calling this function.");
		if( !mxIsDouble(dum) )
			mexErrMsgTxt("o.Quant(iCat).DvBinSep must be allocated as a vector of type double.");
		if( (int)mxGetNumberOfElements(dum) < Nquant+1 )
			mexErrMsgTxt("The length of o.Quant(iCat).DvBinSep must be at least Nquant+1.  In other words, the length of o.Quant(iCat).DvBinSep must be at least one greater than the length of o.Quant(iCat).Dv.");
		QdvSep = mxGetPr(dum);
		QdvSep[0] = -mxGetInf();
		QdvSep[Nquant] = mxGetInf();
		//--------------------------------------------------------

		//	Compute the amount of weight per bin.
		wPerBin = wTotalBins/(double)(Nquant+0.01);

		//	Initialize the sort indexes.  Normally, it would be straightforward to simply create ordered indexes as [0..Ntsamp-1].
		//	However, the quicksort algorithm will proceed faster if we make an informed guess about how these elements would likely
		//	be sorted.  So assume correct entries will have a higher decision value, incorrect will have a lower decision value when
		//	assigning indexes.  This has no effect on the outcome, it just makes QuickSort a bit faster.
		iEntry=0;
		jEntry=Ntsamp-1;
		for(iSamp=0; iSamp<Ntsamp; iSamp++)
		{
			//	Correct
			if( Cat[iSamp]==jCat )
			{
				QuantizeDv_Idx[jEntry]=jEntry;
				jEntry--;
			}
			//	Incorrect
			else
			{
				QuantizeDv_Idx[iEntry]=iEntry;
				iEntry++;
			}
		}
		//	Sort the decision values in Dv.  The pointer to the jCat column of Dv is Dv+Ntsamp*jCat.
		QuickSort(QuantizeDv_Idx,Dv+Ntsamp*jCat,0,Ntsamp-1);

		//	Initialize quantization of sorted list by zeroing the cumulative weight through the list.
		wNextBin = wPerBin;
		w=0.0,wLastEntry=0.0,wLastBin=0.0,dwThis=0.0;
		pcBin=0.0;
		dvBin=0.0;
		//	The current quantile (bin) number.
		iBin=0;
		//	Look ahead to subsequent entries (after iEntry) to see when the current quantile should end.
		//	Then continue advancing through subsequent quantiles until all of the entries have been quantized.
		for(iSamp=1; iSamp<Ntsamp; iSamp++)
		{
			//	Get the index for ascending sort.
			iSorted=QuantizeDv_Idx[iSamp];
			//	Offset iSorted to the correct column of Dv
			iEntry = iSorted + Ntsamp*jCat;
			//	Figure out the sample's category.
			iCat = Cat[iSorted]; // The category of this sample (shift from 1-based to 0-based).
			//	Current weight of this sample
			dwThis = CatWeight[iCat];
			//	Remember accumulated weight of previous entry
			wLastEntry=w;
			//	Accumulate weight for current entry.
			w += dwThis;
			//	If we have accumulated enough to advance to the next bin
			if( w >= wNextBin )
			{
				//	But if the previous entry (iSamp-1) was closer to the accumulated weight we were looking for
				if( wNextBin-wLastEntry < w-wNextBin )
				{
					//	Then rewind to the previous entry before defining the bin.
					iSorted=QuantizeDv_Idx[--iSamp];
					iEntry = iSorted + Ntsamp*jCat;
					w = wLastEntry;
				}
				else
				{
					//	Accumulate totals -- don't rewind
					dvBin += dwThis*Dv[iEntry];
					if( iCat==jCat ) pcBin += dwThis;
				}
				//	Finish the current bin
				Qweight[iBin] = w-wLastBin;
				Qpc[iBin] = pcBin/Qweight[iBin];
				Qdv[iBin] = dvBin/Qweight[iBin];
				//	The bin should never be finished on the last entry, unless the largest decision value has a weight bigger than a bin,
				//	which should never happen, but could cause a runtime crash... so we have the if statement.
				if( iSamp<Ntsamp-1 )
					// The separator between bins is the average of the entries adjacent to the bin separator.
					QdvSep[iBin+1] = 0.5 * (Dv[iEntry]+Dv[QuantizeDv_Idx[iSamp+1]+Ntsamp*jCat]);
				else
					// This should never be executed unless the weight of the last entry is too big (i.e. the user supplied bad input).
					QdvSep[iBin+1] = Dv[iEntry] + 1.0e-4 * (Dv[iEntry]-Dv[QuantizeDv_Idx[iSamp-1]+Ntsamp*jCat]);
				//	Advance bins
				iBin++;
				pcBin=0.0;
				dvBin=0.0;
				wLastBin=w;
				wNextBin+=wPerBin;
				//	If we made it to the last bin, go through and accumulate the rest.
				if( iBin==Nquant-1 )
				{
					//	For each subsequent entry (until finished)
					while(++iSamp<Ntsamp)
					{
						//	Get the index for ascending sort.
						iSorted=QuantizeDv_Idx[iSamp];
						iEntry = iSorted + Ntsamp*jCat;
						//	Figure out which sample this is (and the sample's category).
						iCat = Cat[iSorted]; // The category of this sample (shift from 1-based to 0-based).
						//	Current weight
						dwThis = CatWeight[iCat];
						//	Accumulate weight for current entry.
						w += dwThis;
						//	Accumulate totals
						dvBin += dwThis*Dv[iEntry];
						if( iCat==jCat ) pcBin += dwThis;
					}
					//	Finalize the last bin.
					Qweight[iBin] = w-wLastBin;
					Qpc[iBin] = pcBin/Qweight[iBin];
					Qdv[iBin] = dvBin/Qweight[iBin];
				}
			}
			else
			{
				//	Accumulate totals
				dvBin += dwThis*Dv[iEntry];
				if( iCat==jCat ) pcBin += dwThis;
			}
		}
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_QuantizeDv(nlhs, plhs, nrhs, prhs);
}