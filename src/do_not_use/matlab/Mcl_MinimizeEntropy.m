function out = Mcl_MinimizeEntropy( wInit, Func, Options )

% function out = Mcl_MinimizeEntropy( wInit, Func, Options ) returns a structure of optimization results.
%
% function out = Mcl_MinimizeEntropy('OPTIONS') returns the defalt Options input parameter (a structure).
%-----------------------------------------------------------------------------------------------------------------
%	This function optimizes the free parameters w of the user-specified anonymous function Func(w) by taking large
%	steps in random directions through the parameter space.  Gradient descent is also utilized by this algorithm.
%	The Jacobian and Hessian matrices are not utilized in any way.  The algorithm begins by employing large
%	step-sizes and then gradually reduces the step size over time.  Finite differencing of the gradient is based on
%	that step-size.
%
%	This optimization function can be used to find parameters associated with minimum entropy when the hyper-surface
%	of entropy in the parameter space is "rippled".  In such cases, small finite parameter differences cannot be
%	utilized to perform gradient descent towards a global entropy minimum.  This method solves the problem by using
%	large finite differences to perform simulated annealing while measuring large-scale gradients (for gradient descent).
%	The step size is gradually reduced until a solution is found.  The solution is not guaranteed to be optimal, but
%	this algorithm is likely to avoid most local minima assuming the amount of rippling of the entropy hyper-surface
%	is not too strong.
%-----------------------------------------------------------------------------------------------------------------
% wInit:    Initial parameter values
% 
% Func:     User specified anonymous function that acceps a double vector (like wInit) and returns the entropy.
%
% Options:  The options structure.  For default solver options, evaluate Mcl_MinimizeEntropy('OPTIONS').
%-----------------------------------------------------------------------------------------------------------------
% out.Nevals
%	Scalar, The number of times the function Func() was evaluated.
%
% out.History
%	2D array, The entropy (column 1) and input parameters (subsequent columns) for each time Func() was evaluated (rows).
%
% out.h
%	Scalar, the optimized entropy.
%
% out.wOptimized
%	The optimized parameters.
%-----------------------------------------------------------------------------------------------------------------
%	Optimization is halted after two criteria have been met:
%	(1)	The step size falls below Options.wDiffTol.
%	(2) The entropy changes by less than Options.hTol.
%
%	This function can only be used for unconstrained parameter spaces (i.e. parameters can range in [-inf, +inf]).
%	This function can only be used when a range of maximum and minimum step sizes (for adjusting parameters)
%	is known beforehand.
%
%	Assumes each free parameter is unconstrained (it can exist in the range [-inf, inf]).
%	Assumes parameter iW can adjusted by a maximum amount of
%		Options.wScale(iW) * Options.wDiffMax    OR if Options.wScale is a scalar:  Options.wScale * Option.wDiffMax
%	and a minimum amount of 
%		Options.wScale(iW) * Options.wDiffTol    OR if Options.wScale is a scalar:  Options.wScale * Option.wDiffTol
%-----------------------------------------------------------------------------------------------------------------

%	Return solver options
if nargin==1
	out = struct(...
		'wScale', 1, ...  The scale for all parameters.  Vectorize to set scale of each parameter.
		'wShrinkFactor', 2, ...  The factor by which the parameter adjustment size 
		'wDiffMax', 2, ...	The maximum scale.
		'wDiffTol', 0.01, ... The minimum scale and parameter adjustment tolerance.
		'hTol', 0.0005, ... The precision of entropy.
		'DisplayModulus', 1);  % Display the status every Options.DisplayModulus evaluations of Func.  Set to <1 for no display.
	return;
end

%	Options were not supplied.  Get default options.
if nargin < 3
	Options = Mcl_MinimizeEntropy('OPTIONS');
end

%	Record input shape of wInit
wInitSize = size(wInit);

%	Reshape the initial parameters to be a tall vector.
wInit = wInit(:);

%	Check input options.
if numel(Options.wScale)==1
	Options.wScale = Options.wScale * ones(size(wInit));
elseif numel(Options.wScale) ~= numel(wInit)
	error('Options.wScale must be a scalar or have the same number of elements as wInit.');
end
if Options.wDiffTol<0 || Options.wDiffMax<=Options.wDiffTol
	error('Options.wDiffTol must be greater than 0.  Options.wDiffMax must be greater than Options.wDiffTol.');
end
if Options.wShrinkFactor<=1
	error('Options.wShrinkFactor must be greater than 1.  Suggest 2.');
end
wGrowFactor = Options.wShrinkFactor^(1/3);  %	Grow slower than shrink.

%	Initialize output.
nDims = length(wInit);
out = struct(...
	'Nevals', 1, ...
	'h', Func(wInit), ...
	'wOptimized', wInit, ...
	'ClockBegin', clock, ...
	'ClockEnd', []);

bestStr = ' BEST';

% ---------------------------------------------------------
%	Store history, and display status
% ---------------------------------------------------------
%out.History(out.Nevals,1          ) = out.h;
%out.History(out.Nevals,2:(nDims+1)) = out.wOptimized;
if Options.DisplayModulus>0
	disp(['Init [Neval, EstFinished, h | w] = [' num2str(out.Nevals) ', 0%, ' num2str(out.h) '  |  ' num2str((wInit./Options.wScale(:))') ']' bestStr]);
end
% ---------------------------------------------------------


%	Set the initial step size
diffSize = Options.wDiffMax;

%	There are two modes:  Gradient and Ortho.  When a full orthogonal section has been explored, we can switch to gradient
%	for a while, but not until then.
modeGradient = false;

%	The number 128 is hard-coded as the maximum number of orthogonal bases.  If there are more than 128, than 
%	simple random orthogonal vectors (assumed to be orthonormal) are used in place.
nOrtho = min(128, nDims);
dhOrtho = zeros(nOrtho,1);
iOrtho = 0;
dhOrthoNorm = 0;
bOrtho = f_NewOrthonormalBasis( diffSize * Options.wScale(:) );
dw_dh = Options.wScale(:);
dh_TotalForThisDiffSize = 0;  %	The change in entropy accumulated for the current difference size.

optimizing = true;

while optimizing
	
	%	Try another step.
	if modeGradient
		%	Keep track of total dh for this gradient.
		dh_TotalForThisGradient = 0;
		%	Check to see if the gradient can be computed (or would be worth computing)
		gradLim = Options.hTol*sqrt(diffSize/Options.wDiffTol)/3;
		if dhOrthoNorm>gradLim
			%	Gradient line search
			gradSize = 1;
			%	Try for a decrease in h by setting the direction negative.  This begins a line search along dw_dh.
			%	Later, we will reverse direction on the line search by switching the sign of this variable.
			iDirection = -1;
			while gradSize>gradLim
				%	Try a gradient that will produce a change in h by a factor of gradSize*out.h
				wTry = out.wOptimized + iDirection*max(0.1,out.h)*gradSize*dw_dh;
				% ---------------------------------------------------------
				%	Evaluate function, store history, and display status
				% ---------------------------------------------------------
				hTry = Func(wTry);
				out.Nevals = out.Nevals + 1;
				%if size(out.History,1)<out.Nevals
				%	out.History = [out.History; zeros(size(out.History))];
				%end
				%out.History(out.Nevals,1) = hTry;
				%out.History(out.Nevals,2:(nDims+1)) = wTry;
				if Options.DisplayModulus>0 && mod(out.Nevals,Options.DisplayModulus)==0
					minFin = log(Options.wDiffTol);
					estFin = max(0,min(100, 100-floor(100*(log(diffSize)-minFin)/(log(Options.wDiffMax)-minFin))));
					if hTry < out.h
						bestStr = ' BEST';
					else
						bestStr = '';
					end
					disp(['Grad [Neval, EstFinished, h | w] = [' num2str(out.Nevals) ', ' num2str(estFin) '%, ' num2str(hTry) '  |  ' num2str((wTry./Options.wScale(:))') ']'  bestStr]);
				end
				% ---------------------------------------------------------
				if isinf(hTry) || isnan(hTry)
					disp(['WARNING:  Entropy was NaN or Inf at w:  ' num2str(wTry')]);
					gradSize = gradSize / 2;
					if dh_TotalForThisGradient<0
						%	Only switch direction of line-search if an improvement was previously made.
						iDirection = -iDirection;
					end
				else
					dhGrad = hTry-out.h; %	The change in entropy
					if dhGrad < 0
						%	Accumulate the improvement in entropy for this orthogonal basis.
						dh_TotalForThisGradient = dh_TotalForThisGradient + dhGrad;
						out.h = hTry;
						out.wOptimized = wTry;
					else
						if dh_TotalForThisGradient<0
							%	Only switch direction of line-search if an improvement was previously made.
							iDirection = -iDirection;
						end
						gradSize = gradSize / 2;
					end
				end
			end
		end
		dh_TotalForThisDiffSize = dh_TotalForThisDiffSize + dh_TotalForThisGradient;
		%	Terminate algorithm or continue with a new orthogonal basis?
		if dh_TotalForThisDiffSize > -Options.hTol
			%	Shrink the size of the finite difference.
			diffSize = diffSize / Options.wShrinkFactor;
			%	Terminate algorithm if finite differences are small enough.
			if diffSize <= Options.wDiffTol
				%	Termination criteria has been reached
				optimizing = false;
			end
		else
			%	This orthogonal basis yielded an improvement in entropy.  Therefore, grow the finite differencer.
			diffSize = min(Options.wDiffMax, diffSize*wGrowFactor);
		end
		%	Zero the entropy improvement for the next accumulation (next orthogonal basis).
		dh_TotalForThisDiffSize = 0;
		%	Gradient mode has stopped working for the current step size and orthogonal basis.
		%	Calculate a new orthogonal basis and switch away from gradient mode.
		bOrtho = f_NewOrthonormalBasis( diffSize * Options.wScale(:) );
		dhOrtho(:)=0;
		modeGradient = false;
	else
		%	Increment the count of steps taken inside the orthogonal basis.
		iOrtho = iOrtho + 1;
		%	If we have stepped all the way through...
		if iOrtho>nOrtho
			%	An orthogonal basis was just completed.
			%	Calculate the weighted average of dw (weighted by dh) then normalize to convert to dw / dh, then enter gradient mode.
			dhOrthoNorm = sqrt(sum(dhOrtho.^2));
			dw_dh = bOrtho * dhOrtho / dhOrthoNorm;
			modeGradient = true;
			%	Reset iOrtho
			iOrtho=0;
		else
			%	Try both directions
			for iDirection = [-1, 1]
				wTry = out.wOptimized + iDirection * bOrtho(:,iOrtho);
				% ---------------------------------------------------------
				%	Evaluate function and store history
				% ---------------------------------------------------------
				hTry = Func(wTry);
				out.Nevals = out.Nevals + 1;
				%if size(out.History,1)<out.Nevals
				%	out.History = [out.History; zeros(size(out.History))];
				%end
				%out.History(out.Nevals,1) = hTry;
				%out.History(out.Nevals,2:(nDims+1)) = wTry;
				if Options.DisplayModulus>0 && mod(out.Nevals,Options.DisplayModulus)==0
					minFin = log(Options.wDiffTol);
					estFin = max(0,min(100, 100-floor(100*(log(diffSize)-minFin)/(log(Options.wDiffMax)-minFin))));
					if hTry < out.h
						bestStr = ' BEST';
					else
						bestStr = '';
					end
					disp(['Orth [Neval, EstFinished, h | w] = [' num2str(out.Nevals) ', ' num2str(estFin) '%, ' num2str(hTry) '  |  ' num2str((wTry./Options.wScale(:))') ']'  bestStr]);
				end
				% ---------------------------------------------------------
				%	The change in entropy for each member of the orthogonal basis
				if ~isnan(hTry) && ~isinf(hTry)
					%	Track the change in h for this orthonormal basis.
					dhOrtho(iOrtho) = dhOrtho(iOrtho) + iDirection * (hTry-out.h)/2;
					%	Check for improvement
					if hTry < out.h
						%	Accumulate the improvement in entropy for this orthogonal basis.
						dh_TotalForThisDiffSize = dh_TotalForThisDiffSize + hTry-out.h;
						%	Store any improvement
						out.wOptimized = wTry;
						out.h = hTry;
						%	Do not try the opposite direction
						if iDirection==-1
							dhOrtho(iOrtho) = 2*dhOrtho(iOrtho);  % Only a single difference (weighted 1/2) contributed, so double it.
							break;
						end
					end
				else
					disp(['WARNING:  Entropy was NaN or Inf at w:  ' num2str(wTry')]);
					%	dhOrtho(iOrtho) is already zero (unless it was set by this basis in the other direction).
				end
			end
		end
	end
end

%	Clip extra memory from solver's history
%out.History = out.History(1:out.Nevals,:);

%	Reshape the free parameters to match the original shape of the input.
out.wOptimized = reshape(out.wOptimized, wInitSize);

%	Note the stop time.
out.ClockEnd = clock;

return;

%	Get a new orthogonal basis vectors D(:,iDiff) for searching the parameter space in each iDiff direction.
%	The argument del is a tall vector representing the scale of desired adjustments for each parameter.
function D = f_NewOrthonormalBasis( del )
nOrtho = min(128, numel(del));
if nOrtho<128
	%	Get a real orthonormal basis.
	D = Mcl_RandRotate( int32(nOrtho) );
	for iDiff=1:nOrtho
		D(:,iDiff) = D(:,iDiff).*del;
	end
else
	%	Get an approximate orthonormal basis by creating 128 random unit vectors.
	D = rand(numel(del),nOrtho)-0.5;
	for iDiff=1:nOrtho
		D(:,iDiff) = D(:,iDiff).*del/sqrt(sum(D(:,iDiff).^2));
	end
end