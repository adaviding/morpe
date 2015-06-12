function o = Mcl_Exemplar_Train(o, clipMem)

% function o = Mcl_Exemplar_Train(o, clipMem)
% Trains the exemplar model based on the training specifications contained in the input.  Returns optimization results and updates fields
%	of the input structure "o" (using mex functions).
% Suggested usage:
%	o = Mcl_Exemplar_Train(o, true);
% Memory clipping:
%	When you are done optimizing "o", you can reduce its size (in memory or on disk) by clearing some fields as follows.
%		o.Cube = []; % The Cube occupies a huge amount of memory
%		o.Mem = []; % Totally meaningless memory
% ------------------------------------------------------------------
% INPUT and OUTPUT
% o
%	The structure created by Mcl_Exemplar_Ctor.  It's contents are altered by many of the mex functions called by this function.
%	When this function exits, it will return many values that should also be incorporated into the structure o. This can be accomplished
%	by using this function as suggested above.
% [clipMem] (optional logical scalar, default = true)
%	If true, this algorithm will "clip" some memory that is only required for optimization but not afterwards.  This results in a reduction
%		of the solver's size in memory or on disk by deleting o.Cube and o.Mem.
% ------------------------------------------------------------------
% OUTPUT
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

if isfield(solvOptions, 'Napproaches') && ~isempty(solvOptions.Napproaches) && isnumeric(solvOptions.Napproaches)
	iApproach = solvOptions.Napproaches;
else
	iApproach = 1;
end

%	Use the fast evaluation method which utilizes a fixed set of nearby neighbors.
hFunc = @(w) Mcl_Exemplar_Hfunc(o,true,w);

%	Make several approaches and take the one that maximizes accuracy.
o.h = +Inf;
bestApproach=iApproach;

%	Continue until iApproach drops below 1, but if these approaches yield worse than chance performance (o.h>1), go at least for two more special approaches.
while iApproach>0 || (o.h>1 && iApproach>-2)
	
	if iApproach<1 && o.h>1
		% -------------------------------------------------
		% One of two special approaches is needed
		% -------------------------------------------------
		if iApproach>=0
			% First special approach:
			%	Solve the exemplar parameters starting from a very blurry, almost totally flat kernel.
			slnTry = Mcl_MinimizeEntropy( log2(o.wInit)-16,hFunc,solvOptions );
		else
			% Second special approach:
			%	Set the exemplar parameters to a very blurry, almost totally flat kernel.
			slnTry.wOptimized = log2(o.wInit)-16;
		end
		% -------------------------------------------------
	else
		% -------------------------------------------------
		% Proceed normally
		% -------------------------------------------------
		%	Solve the exemplar parameters normally
		slnTry = Mcl_MinimizeEntropy(    log2(o.wInit),hFunc,solvOptions );
		% -------------------------------------------------
	end
	
	%	Evaluate the exemplar model using the optimal parameters on the full set of training data (not using the fast shortcut).
	[hTry, HcrossTry] = Mcl_Exemplar_Hfunc(o,false,slnTry.wOptimized);
	
	if (hTry < o.h)
		% ------------------------------------------------------------------
		% Improvement detected while optimization was used.  Calculate performance and save.
		% ------------------------------------------------------------------
		%	Compute the accuracy and confusion matrix.
		accTry = 0;
		cfnTry = zeros(o.Ncats,o.Ncats);
		%	Decision rule is based on this probability cutoff.
		pCut = 1/double(o.Ncats);
		[maxResp, indResp] = max(o.P,[],2);
		if o.ForceEqualPriors
			for iCat=1:o.Ncats
				%	0-based
				indCat = o.Cat==(iCat-1);
				for jCat=1:o.Ncats
					%	1-based
					cfnTry(iCat,jCat) = mean(indResp(indCat)==jCat)*pCut;
				end
				accTry = accTry + cfnTry(iCat,iCat);
			end
		else
			for iCat=1:o.Ncats
				%	0-based
				indCat = o.Cat==(iCat-1);
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
	elseif iApproach==-1
		% ------------------------------------------------------------------
		% Flat model must be enforced (all optimization failed to obtain at least chance performance).
		% ------------------------------------------------------------------
		if o.ForceEqualPriors
			o.Acc = 1/double(o.Ncats);
			o.ConfusionMatrix(:) = 1/double(o.Ncats*o.Ncats);
			o.h = 1;
			o.Hcross(:) = 1/double(o.Ncats);
			for iCat=1:o.Ncats
				o.Quant(iCat).Pc(:) = 1/double(o.Ncats);
			end
		else
			[cwMin, iMin] = min(o.CatWeight);
			o.Acc = 1/cwMin/double(o.Ncats);
			o.ConfusionMatrix(:) = 0;
			o.ConfusionMatrix(iMin,:)=1./o.CatWeight/double(o.Ncats);
			o.Hcross(:)=1;
			for iCat=1:o.Ncats
				o.Quant(iCat).Pc(:) = 1/o.CatWeight(iCat)/double(o.Ncats);
				o.Hcross(iCat,:) = o.Hcross(iCat,:)/o.CatWeight(iCat)/double(o.Ncats);
				o.Hcross(iCat,:) = o.Hcross(:,iCat)* log(o.CatWeight(iCat)*double(o.Ncats)) / log(double(o.Ncats)); % log base o.Ncats
			end
		end
		for iCat=1:o.Ncats
			o.Quant(iCat).PcMono(:) = o.Quant(iCat).Pc(:);
			o.Quant(iCat).Pc(:) = o.Quant(iCat).PcMonoLim(:);
		end
		bestApproach = 1; % save final approach (this one).
		iApproach = -2;   % exit
		% ------------------------------------------------------------------
	end
	iApproach = iApproach - 1; % cycle until exit
end

%	If necessary, return with the object in its correct state by evaluating the winning attempt on the
%	full set of training data (not using the fast shortcut).
if bestApproach > 1
	Mcl_Exemplar_Hfunc(o,false,o.SolverOutput.wOptimized);
end
%	Store optimized weights
o.wOptimized = 2.^o.SolverOutput.wOptimized;
if o.IdMethod>=10
	o.wOptimized = reshape(o.wOptimized, [o.Ndims,o.Ncats]);
end
%	Clip memory
if nargin<2 || clipMem
	o.Cube = [];
end