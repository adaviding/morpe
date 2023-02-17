function varargout = Mcl_Poly_Hfunc(o, wParams)

% function h = Mcl_Poly_Hfunc(o, wParams)
%	Classifies each training sample and returns the entropy for the entire training set.  Useful during optimization.
% ------------------------------------------------------------------
% INPUT
% o
%	A structure constructed from Mcl_Poly_Ctor (refer to its comments).
% wParams
%	A vector of free parameters for the exemplar model, each parameter is in the range [-Inf, +Inf].  2.^wParams represents the
%	scale of each spatial dimension and therefore determines the size and shape of the blurring (generalization) kernel.
% ------------------------------------------------------------------
% OUTPUT
% h (double scalar)
%	The entropy of the training data in (0, 1].  Entropy approaching 0 represents excellent performance, approaching 1
%	represents poor performance.
% Hcross (double matrix, o.Ncats * o.Ncats)
%	If nargout>=2, this function also outputs the conditional cross-entropy matrix as a second return value.
%		See Mcl_ConditionalEntropy for more info about Hcross.
% ------------------------------------------------------------------

%	Calculate the decision values
Mcl_Poly_CalcDv(o, o.Dv, o.X, wParams);

%	Quick check to see if there were any problems
if any(isnan(o.Dv(:))) || any(isinf(o.Dv(:)))
	dbstop if error;
	error(['Decision values evaluated to NaN or Inf when wParams = [' num2str(wParams') ']']);
end

%	Quantize the decision values
Mcl_QuantizeDv(o);

%	Monotonic correction of quantized decision values.
for iCat = 1:o.Ncats
	
	%	The probability of category membership must be fit with a monotonic function.
	Mcl_ForceMonotonic(o.Quant(iCat).PcMono, o.Quant(iCat).Pc);
	
	%	Limit the range of each probability.
	Mcl_RangeLimit(o.Quant(iCat).PcMonoLim, o.Quant(iCat).PcMono, o.Quant(iCat).pLow, o.Quant(iCat).pHigh);
	
end

%	Map the decision values to probabilities.
Mcl_MapDv(o.P, o.Dv, o.Quant);

%	Output conditional entropy
if nargout<2
	h = Mcl_ConditionalEntropy(o.H, o.P, o.Cat, o.ForceEqualPriors);
	varargout = {h};
else
	[h, Hcross] = Mcl_ConditionalEntropy(o.H, o.P, o.Cat, o.ForceEqualPriors);
	varargout = {h, Hcross};
end