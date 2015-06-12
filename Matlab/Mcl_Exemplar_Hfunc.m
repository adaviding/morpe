function varargout = Mcl_Exemplar_Hfunc(o, goFast, wParams)

% function h = Mcl_Exemplar_Hfunc(o, isTraining, wParams)
%	Classifies each training sample and returns the entropy for the entire training set.  Useful during optimization.
% ------------------------------------------------------------------
% INPUT
% o
%	A structure constructed from Mcl_Exemplar_Ctor (refer to its comments).
% goFast
%	If true, the neighborhood cube is used for generalization (blurring) because this provides a fast and nearly accurate method
%	for evaluating entropy. If false, the slower method of utilizing the full set of training samples is used.
% wParams
%	A vector of free parameters for the exemplar model, each parameter is in the range [-Inf, +Inf].  2.^wParams represents the
%	scale of each spatial dimension and therefore determines the size and shape of the blurring (generalization) kernel.
% ------------------------------------------------------------------
% OUTPUT
% h (double scalar)
%	The entropy of the training data in (0, 1].  Entropy approaching 0 represents excellent performance, approaching 1
%	represents poor performance.
% Hcross (double matrix, o.Ncats * o.Ncats)
%	If nargout>2, this function also outputs the conditional cross-entropy matrix as a second return value.
%		See Mcl_ConditionalEntropy for more info about Hcross.
% ------------------------------------------------------------------


%	Evaluate the decision function of each category for each training sample.
if goFast
	%	Compute the squared weights from the free parameters
	wSq = 2.^(2*wParams);
	Mcl_Exemplar_CalcDv_Train(o, wSq);
else
	%	Compute the weights from the free parameters
	w = 2.^wParams;
	Mcl_Exemplar_CalcDv(o, o.Dv, o.X, false, w);
end

%disp('Hfunc 2');

%	Quantize the decision values to make it possible to compute the probability of category membership
Mcl_QuantizeDv(o);

%disp('Hfunc 3');

for iCat=1:o.Ncats
	%	If the function is a hard step (threshold), then only solve for a non-decreasing function.
	doMono = int32(0);
	if 0.2>(  mean(o.Quant(iCat).PcMono<(3*o.Quant(iCat).pLow))+mean(o.Quant(iCat).PcMono>(1-3*(1-o.Quant(iCat).pHigh)))  )
		%	At least 20% of the p-values are significantly away from 0 and 1, so it's not such a hard step (threshold).
		%	Monotonic smoothing is okay.
		doMono = int32(1);
	end
	%	The probability of category membership must be fit with a monotonic function.
	Mcl_ForceMonotonic(o.Quant(iCat).PcMono, o.Quant(iCat).Pc, doMono);
	
	%disp('Hfunc 4');

	%	Limit the range of each probability.
	Mcl_RangeLimit(o.Quant(iCat).PcMonoLim, o.Quant(iCat).PcMono, o.Quant(iCat).pLow, o.Quant(iCat).pHigh);
	
	%disp('Hfunc 5');
end

%	Map the decision variables Dv to probabilities of category membership P
Mcl_MapDv(o.P, o.Dv, o.Quant);

%disp('Hfunc 6');

%	Calculate the conditional entropy
if nargout<2
	%disp('Hfunc 7');
	h = Mcl_ConditionalEntropy(o.H, o.P, o.Cat, o.ForceEqualPriors);
	varargout = {h};
else
	%disp('Hfunc 8');
	[h, Hcross] = Mcl_ConditionalEntropy(o.H, o.P, o.Cat, o.ForceEqualPriors);
	varargout = {h, Hcross};
end