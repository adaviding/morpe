//	Almon David Ing
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Compiled July 11, 2009 using Microsoft Visual C++ 2008 on 32-bit Windows Vista

#include "mex.h"
#include "math.h"
#include "Linterp.cpp"

#ifndef null
	#define null (void*)0
#endif
//#ifndef DEBUG_MEXPRINTF
//	#define DEBUG_MEXPRINTF
//#endif

//================================================================================================================================
//function Mcl_QuantizedDv_Complement(Qout, Qin)
//--------------------------------------------------------------------------------------------------------------------------------
// This mex function creates a quantization table Qout that is the complement (negative reflection) of another table Qin.  This is
//	useful for the Mcl_Poly classifier.  When there are two categories, the Mcl_Poly classifier only trains a single polynomial.
//	A second polynomial is implied as the negative decision value of the first classifier.
//--------------------------------------------------------------------------------------------------------------------------------
// INPUT (values are not changed by this mex function)
//----------------------------------------------------
// Qin (mxArray*, Matlab structure, Mcl_QuantizedDv scalar)
//	An Mcl_QuantizedDv structure that has already been set.
//--------------------------------------------------------------------------------------------------------------------------------
// OUTPUT (pre-allocated memory will be filled with the output of this function.)
//----------------------------------------------------
// Qin (mxArray*, Matlab structure, Mcl_QuantizedDv scalar)
//	An Mcl_QuantizedDv structure that has already been initialized (all fields exist and all memory is allocated in those fields).
//================================================================================================================================
void Mcl_QuantizedDv_Complement(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	Basic error checking of arguments doesn't take very long and helps insure against major disasters.
	if( nrhs<2 )
		mexErrMsgTxt("Not enough input arguments.");

	mxArray *fld;
	//-------------------------------------------------------------------------------------------------------
	// INPUT
	//----------------------------------------------------
	if( !mxIsStruct(prhs[1]) )
		mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	const mxArray* Qin = prhs[1];
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"Nquantiles");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsInt32(fld) ) mexErrMsgTxt("Qin.Nquantiles must be type int32.");
	int Nquant = ((int*)mxGetData(fld))[0];
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"Dv");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qin.Dv must be type double.");
	double* QiDv = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"DvBinSep");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qin.DvBinSep must be type double.");
	double* QiDvBinSep = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"Weight");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qin.Weight must be type double.");
	double* QiWeight = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"Pc");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qin.Pc must be type double.");
	double* QiPc = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"PcMono");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qin.PcMono must be type double.");
	double* QiPcMono = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qin,0,"PcMonoLim");
	if( fld==null ) mexErrMsgTxt("The Qin argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qin.PcMonoLim must be type double.");
	double* QiPcMonoLim = (double*)mxGetData(fld);
	//-------------------------------------------------------------------------------------------------------
	// OUTPUT
	//----------------------------------------------------
	if( !mxIsStruct(prhs[0]) )
		mexErrMsgTxt("The first argument of this function (o) must be a matlab structure of type Mcl_QuantizedDv.");
	const mxArray* Qout = prhs[0];
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"Nquantiles");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsInt32(fld) ) mexErrMsgTxt("Qout.Nquantiles must be type int32.");
	if( Nquant != ((int*)mxGetData(fld))[0] ) mexErrMsgTxt("Qout.Nquantiles must be equal to Qin.Nquantiles.");
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"Dv");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qout.Dv must be type double.");
	double* QoDv = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"DvBinSep");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qout.DvBinSep must be type double.");
	double* QoDvBinSep = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"Weight");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qout.Weight must be type double.");
	double* QoWeight = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"Pc");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qout.Pc must be type double.");
	double* QoPc = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"PcMono");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qout.PcMono must be type double.");
	double* QoPcMono = (double*)mxGetData(fld);
	//----------------------------------------------------
	fld = mxGetField(Qout,0,"PcMonoLim");
	if( fld==null ) mexErrMsgTxt("The Qout argument must be a matlab structure of type Mcl_QuantizedDv.");
	if( !mxIsDouble(fld) ) mexErrMsgTxt("Qout.PcMonoLim must be type double.");
	double* QoPcMonoLim = (double*)mxGetData(fld);
	//----------------------------------------------------
	mxGetPr(mxGetField(Qout,0, "pLow"))[0] = mxGetScalar(mxGetField(Qin,0, "pLow"));
	mxGetPr(mxGetField(Qout,0,"pHigh"))[0] = mxGetScalar(mxGetField(Qin,0,"pHigh"));
	mxGetPr(mxGetField(Qout,0, "hMin"))[0] = mxGetScalar(mxGetField(Qin,0, "hMin"));
	//-------------------------------------------------------------------------------------------------------
	int ii, io;
	//	Bin limits are always -Inf, +Inf.
	QoDvBinSep[0] = -mxGetInf();
	QoDvBinSep[Nquant-1] =  mxGetInf();
	//	For each quantile
	for(ii=0; ii<Nquant; ii++)
	{
		io = Nquant-ii-1;
		//	 Decision values are an off-center reflection
		QoDv[io] = -QiDv[ii];
		if( ii>0 )
			QoDvBinSep[io+1] = -QiDvBinSep[ii];
		//	Weight is the same
		QoWeight[io] = QiWeight[ii];
		//	Pc are 1.0 complements.
		QoPc[io]		= 1.0-QiPc[ii];
		QoPcMono[io]	= 1.0-QiPcMono[ii];
		QoPcMonoLim[io]	= 1.0-QiPcMonoLim[ii];
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Mcl_QuantizedDv_Complement(nlhs, plhs, nrhs, prhs);
}