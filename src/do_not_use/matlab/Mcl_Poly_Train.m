function o = Mcl_Poly_Train(o)

% function o = Mcl_Exemplar_Train(o)
% Trains the polynomial model based on the training specifications contained in the input.  Returns optimization results and updates the input
%	"o" then returns it as output.
% Suggested usage:
%	o = Mcl_Poly_Train(o);
% ------------------------------------------------------------------
% INPUT and OUTPUT
% ------------------------------------------------------------------
% o
%	The structure created by Mcl_Exemplar_Ctor.  It's contents are altered by many of the mex functions called by this function.
%	When this function exits, it will return this structure with updated values (some of which are described below).  Caution:  Mex-functions
%	will change some of the values in the input structure, so the input will be changed.
% ------------------------------------------------------------------
% OUTPUT
% ------------------------------------------------------------------
% o.SolverOutput
%	The output of the minimization solver.
% o.h
%	The conditional entropy
% o.Hcross
%	The conditional cross-entropy matrix.
% o.hLim
%	The limiting entropy
% o.wOptimized
%	The optimized weights
% o.Acc
%	The classifier accuracy
% o.ConfusionMatrix
%	The proportion of times when the sample was from category iCat and was classified as jCat is given as o.ConfusionMatrix(iCat,jCat).
% ------------------------------------------------------------------

%	Get the default solver options and pass them through.
if( isempty(o.SolverOptions) )
	solvOptions = Mcl_MinimizeEntropy('OPTIONS');
else
	solvOptions = o.SolverOptions;
end

%	Set the scale.
solvOptions.wScale = o.wScale;

%	Initialize approach countdown.
if isfield(solvOptions, 'Napproaches') && ~isempty(solvOptions.Napproaches)
	iApproach = solvOptions.Napproaches;
else
	iApproach = 1;
end

%	Set the entropy function (for usage during optimization)
hFunc = @(w) Mcl_Poly_Hfunc(o,w);

%	Make several approaches and take the one that maximizes accuracy.
o.h = +Inf;
bestApproach=iApproach;

%	Continue until iApproach drops below 1.
%	Each approach is from the same starting point, but along a different path.
while iApproach>0
	
	%	Attempt an optimized solution of the free parameters
	slnTry = Mcl_MinimizeEntropy(o.wInit,hFunc,solvOptions);
	
	if (slnTry.h < o.h)
		% ------------------------------------------------------------------
		% Improvement detected.  Calculate performance and save.
		% ------------------------------------------------------------------
		%	Evaluate the model using the optimal parameters.
		[hTry, HcrossTry] = Mcl_Poly_Hfunc(o,slnTry.wOptimized);
		%	Compute the accuracy and confusion matrix.
		accTry = 0;
		cfnTry = zeros(o.Ncats,o.Ncats);
		%	Decision rule is based on this probability cutoff.
		pCut = 1/double(o.Ncats);
		[maxResp, indResp] = max(o.P,[],2);
		if o.ForceEqualPriors
			for iCat=1:o.Ncats
				%	0-based
				indCat = o.Cat==int32(iCat-1);
				for jCat=1:o.Ncats
					%	1-based
					cfnTry(iCat,jCat) = mean(indResp(indCat)==jCat)*pCut;
				end
				accTry = accTry + cfnTry(iCat,iCat);
			end
		else
			for iCat=1:o.Ncats
				%	0-based
				indCat = o.Cat==int32(iCat-1);
				for jCat=1:o.Ncats
					%	1-based
					cfnTry(iCat,jCat) = sum(indResp(indCat)==jCat)*length(indCat)/o.Ntsamp;
				end
				accTry = accTry + cfnTry(iCat,iCat);
			end
		end
		o.SolverOutput = slnTry;
		o.Acc = accTry;
		o.ConfusionMatrix = cfnTry;
		o.h = hTry;
		o.Hcross = HcrossTry;
		bestApproach = iApproach;  % save best approach
		% ------------------------------------------------------------------
	end
	
	%	Decrement the remaining number of approaches
	iApproach = iApproach - 1;
	
	%	Display result
	if o.DisplayModulus>0
		disp(['Mcl_Poly_Train: ' num2str(iApproach) ' approaches left.  [h, Acc] = [' num2str(slnTry.h) ',' num2str(accTry) ']']);
		disp('--------------------------------------------------------------------------------');
	end
end

%	If necessary, restore the object (o) to its optimal state (by evaluating Hfunc with optimal parameters).
if bestApproach > 1
	Mcl_Poly_Hfunc(o,o.SolverOutput.wOptimized);
end
%	Store optimized weights
o.wOptimized = o.SolverOutput.wOptimized;