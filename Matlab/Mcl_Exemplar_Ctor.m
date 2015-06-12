function out = Mcl_Exemplar_Ctor(Xcell, Xtrans, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles, SolverOptions)

% function Output = Mcl_Exemplar_Ctor(Xcell, [Xtrans], IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles, SolverOptions)
%--------------------------------------------------------------------------------------------------------------------------------
% This constructs a new Mcl_Exemplar classifier.  This function initializes the classifier and utilizes the mex function
%	Mcl_Exemplar_Init.  After construction, the classifier can be trained by calling the Mcl_Train function.  After training,
%	the classifier can be used to classify new data using Mcl_Classify.
%--------------------------------------------------------------------------------------------------------------------------------
% NOMENCLATURE (for interpreting subsequent comments)
%----------------------------------------------------
% nCats = The number of categories (i.e. the length of the input X cell array)
% nSamp(iCat) = The number of training samples in each category (i.e. the number of rows of the input X{iCat} summed across all iCat).
% ntSamp = The total number of training samples = sum(nSamp);
% nDims = The spatial dimension of the data (i.e. the number of columns for each input X{iCat}, this is the same for each iCat).
%--------------------------------------------------------------------------------------------------------------------------------
% INPUT (values are not altered by mex function)
%----------------------------------------------------
% Xcell (cell vector: nCats)
%  Each iCal element of this cell array, Xcell{iCat}, is a 2D array (nCats,nDims) contains samples of training data from category iCat.
%
% [Xtrans] Optional (matlab structure Mcl_Trfm array: nDims)
%  Each spatial dimension can be associated with a transform.  The contents of Xcell are transformed accordingly before the classifier
%	is initialized.  You can also leave this argument empty [] and no transforms will be applied.
%
% IdMethod (int32 scalar) Suggest 1.
%	Identifies the generalization (blurring) function to be used.  Must be an integer belonging to the set {0,1,2,10,11,12}
%		0	k = 1/dist
%		1	k = exp(-dist)
%		2	k = exp(-dist*dist)
%		------------------------
%		+0	Use nDims weight parameters (each category uses the same weight parameters)
%		+10 Use nDims*nCats weight parameters (each category uses unique weight parameters).
%
% MaxNeighbors (int32 scalar)  Suggest 100.
%	The maximum number of nearest neighbors to consider while the solver is being optimized.  The shape of the neighborhood is defined
%	by the optional parameter wInit.  If wInit is not provided, then it is estimated based on the summary statistics (means and variances)
%	of the training data.  wInit defines the shape of the neighborhood because it defines how each axis is scaled.  As this number is
%	reduced, the solver becomes faster, but the accuracy of the solution will also decrease.  In general, a neighborhood size of 200
%	is sufficient to get accurate results when the number of spatial dimensions (D) is less than 20.
%
% ForceEqualPriors (bool scalar) Suggest true.
%	If true, then the categories are assumed to have equal prior probabilities, even if there are different numbers of samples for
%	each category in the training data.
%
% Nquantiles (int32 scalar) Suggest any integers from 10 to 50.
%	The number of quantiles to be used in estimating the probability of category membership based on each cateogory's decision variable.
%	In general, a good value for Nquantiles is int32(ntSamp/30) because it yields ~30 samples per quantile.
%
% [SolverOptions] (optional matlab structure Mcl_MinimizeEntropy('OPTIONS') scalar)
%	Optimization options.  If provided, the solver Mcl_Exemplar_Train is called.
%--------------------------------------------------------------------------------------------------------------------------------
% OUTPUT is a newly allocated structure filled with the following data
%----------------------------------------------------
% Ncats (int32 scalar)
%	The number of categories.
% Ndims (int32 scalar)
%	The number of spatial dimensions (i.e. the number of columns of every iCat element of input Xcell, Xcell{iCat}).
% Ntsamp  (int32 vector, nCats)
%	The total number of training samples (i.e. the total  of rows for all iCat elements of input Xcell, Xcell{iCat}).
% Nsamp  (int32 vector, nCats)
%	The number of samples in each category (i.e. the number of rows for each iCat element of input Xcell, Xcell{iCat}).
% IdMethod (int32 scalar)
%	Same as the input IdMethod.
% MaxNeighbors (int32 scalar)
%	Same as the input MaxNeighbors.
% Nneighbors (int32 scalar)
%	Same as the input MaxNeighbors.  This value can be adjusted before Mcl_Exemplar_Train is executed to reduce the number of nearest
%	neighbors considered by the model.
% ForceEqualPriors (bool scalar)
%	Same as input ForceEqualPriors.
% CatWeight (double vector, Ncats)
%	If the training sample has equal numbers of samples from each category, then this will be a vector of ones; otherwise, if ForceEqualPriors
%	is true, then these weights will scale the influence of training samples so that the total weight is ntSamp and that each category
%	has equal influence on the solution.
% Cat (int32 vector, ntSamp)
%	The 0-based category label of each sample after samples from all categories are concatenated into a single list.
% Xtrans (structure, nDims)
%	The transforms applied to the column of X.  All subsequent output varaibles assume X is in the transformed space already.
% X (double 2D array, ntSamp * nDims)
%	The spatial coordinate of each sample after samples from all categories are concatenated into a single list.  Also, dimensions have
%	been transformed according to Xtrans.
% Xmeans (double 2D array, nCats * nDims)
%	The mean of the sample from each category
% Xvars (double 2D array, nCats * nDims)
%	The variance of the sample from each category.
% Xmean (double vector, nDims)
%	The pooled mean (where CatWeight was used to weight each category's influence on the result).
% Xvar (double vector, nDims)
%	The pooled variance (where CatWeight was used to weight each category's influence on the result).  Variance of category means is not included.
%	This is literally the weighted some of Xvars across categories.
% Cube (double 3D array, (1+nDims) * MaxNeighbors * ntSamp)
%	This cube lists the squared difference between every sample and its MaxNeighbors nearest neighbors.  Squared differences are listed for each
%	iDim spatial dimension as Cube(1+iDim,:,:).  The entries of Cube(1,:,:) are the category label (1-based) of the nearest neighbor.
% Dv (double 2D array, ntSamp * nCats):  This variable is calculated in Mcl_Exemplar_Train.
%	This lists the decision variable (i.e. the total kernel density) calculated for each training sample relative to all other training samples.
%	A decision variable is computed using each category's decision function (i.e. the blurring or generalization function) applied to each sample.
% P (double 2D array, ntSamp * nCats):  This variable is calculated in Mcl_Exemplar_Train.
%	This lists the probability that each category should be assigned to each training sample.  Probabilities are computed by linear interpolation
%	of Dv through the quantization table [Quant.Dv, Quant.PcMonoLim] and then normalizing each row to sum to 1.0.
% H (double vector, ntSamp):  This variable is calculated in Mcl_Exemplar_Train.
%	This lists the conditional entropy of each training sample.  Entropy appraches 0.0 when the sample is categorized correctly, 1.0 is chance performance.
% wInit (double vector, nDims)
%	The initial value of the free parameters for the blurring or generalization function.  These parameters control the width of the blurring
%	kernel for each spaital dimension.
% h (double scalar):  This variable is calculated in Mcl_Exemplar_Train.
%	The conditional entropy of the entire data set in the range [Quant.hMin, 1.0].  This value does not approach 0.0 because the probabilities of
%	category membership are limited by Quant.PcMonoLim, and these probabilities are never allowed to equal 0.0 or 1.0.  Instead they are limited
%	to the range [Quant.pLow, Quant.pHigh].
% hLim (double scalar):  This variable is calculated in Mcl_Exemplar_Train.
%	The conditional entropy mapped uniformly from the range [Quant.hMin, 1.0] to the range [0.0, 1.0].  This is like "stretching" the range.
% Acc (double scalar):  This variable is calculated in Mcl_Exemplar_Train.
%	The weighted accuracy.  1.0 is perfect accuracy, 1/nCats is chance.  Weights are specified in CatWeight.
% Quant (mxArray matlab structure, nCats):  This variable is calculated in Mcl_Exemplar_Train.
%	This structure specifies the mapping between the decision function for category iCat and the probability of category membership.  It is calculated
%	based on quantization of each category's decision variable and training sample.
% Quant(iCat).Nquantiles (int32 scalar):  This variable is calculated in Mcl_Exemplar_Train.
%	The number of quantiles (i.e. entries of the interpolation table).
% Quant(iCat).Dv (double vector, Nquantiles)
%	The mean decision value for each quantile (or row of interpolation table).
% Quant(iCat).DvBinSep (double vector, Nquantiles)
%	The boundary of each quantile.  The first value is -Inf and last value is +Inf.
% Quant(iCat).Weight
%	The total weight of training samples accumulated in each quantile (or row of the interpolation table).  See CatWeight for explanation of weighting.
% Quant(iCat).Pc
%	The probability that the correct category will be assigned (for the training data) for each quantile (or row of interpolation table).
% Quant(iCat).PcMono
%	A non-decreasing function fit to Pc.
% Quant(iCat).PcMonoLim
%	The same function as PcMono, but probabilities have been limited not to equal 0.0 or 1.0, instead they are limited to [pLow, pHigh].
% Quant(iCat).pLow (double scalar)
%	The lower limit of probability.  This approaches 0.0 as average Weight approaches infinity, otherwise is always close to 0.0.
% Quant(iCat).pHigh (double scalar)
%	The upper limit of probability.  This approaches 1.0 as average Weight size approaches infinity, otherwise is always close to 1.0.
% Quant(iCat).hMin (double scalar)
%	The minimum possible entropy based on pHigh.  This is always positive but approaches 0.0 as average Weight approaches infinity.
% SolverOptions (mxArray, Matlab structure):  This variable is empty.
%	This is just a placeholder for storing the structure that configures the optimization algorithm.
% Napproaches (int32 scalar):  This can also be an integer as double.
%	The number of times in which the solver makes a unique approach towards the solution from the initial parameter values.  As this number increases,
%	so does the probability that you found the real global optimum.  Suggest Napproaches = 3.
% SolverOutput (mxArray, Matlab structure):  This variable is calculated in Mcl_Exemplar_Train.
%	After Mcl_Exemplar_Train is executed, this variable contains raw output from the optimization algorithm.
% Mem (mxArray, Matlab structure)
%	This is just memory allocated to make Mcl_Exemplar_Train run as efficiently as possible.

nCats = numel(Xcell);
nDims = size(Xcell{1},2);
nSamp = zeros(1,nCats);
for iCat=1:nCats
	nSamp(iCat) = size(Xcell{iCat},1);
end
ntSamp = sum(nSamp);
if ~exist('wInit', 'var') || isempty(wInit)
	if IdMethod<10
		wInit = NaN(nDims,1);
	elseif IdMethod<20
		wInit = NaN(nDims,nCats);
	end
end
if ~exist('DisplayMessages', 'var') || isempty(DisplayMessages)
	DisplayMessages = true;
end
if ~exist('Nquantiles', 'var') || isempty(Nquantiles)
	Nquantiles = min(20,ceil(ntSamp/20));
end
if ~exist('ForceEqualPriors', 'var') || isempty(ForceEqualPriors)
	ForceEqualPriors = true;
end
if ~exist('SolverOptions', 'var')
	SolverOptions = [];
end

%	All memory must be allocated in M-file, not Mex-File.  This protects against funny business with Matlab's memory manager.
out = struct(...
	'Class', 'Mcl_Exemplar', ...
	'Ncats', int32(nCats), ...
	'Ndims', int32(nDims), ...
	'Ntsamp', int32(ntSamp), ...
	'Nsamp', int32(nSamp), ...
	'IdDims', zeros(1,nDims,'int32'), ... % Globally unique identifiers for each stimulus dimension.
	'IdMethod', int32(IdMethod), ...
	'MaxNeighbors', int32(MaxNeighbors), ...
	'Nneighbors', int32(MaxNeighbors), ...
	'ForceEqualPriors', logical(ForceEqualPriors), ...
	'CatWeight', zeros(1,nCats), ...
	'Cat', zeros(ntSamp,1,'int32'), ...
	'Xtrans', Xtrans, ...
	'X', zeros(ntSamp,nDims), ...
	'Xmeans', zeros(nCats,nDims), ...
	'Xvars', zeros(nCats,nDims), ...
	'Xmean', zeros(1,nDims), ...
	'Xvar', zeros(1,nDims), ...
	'Cube', zeros(1+nDims,MaxNeighbors,ntSamp,'single'), ...
	'Dv', zeros(ntSamp,nCats), ...
	'P', zeros(ntSamp,nCats), ...
	'H', zeros(ntSamp,1), ...
	'wInit', wInit, ...
	'wOptimized', zeros(size(wInit)), ...
	'h', 0, ...
	'Hcross', zeros(nCats,nCats), ...
	'Acc', 0, ...
	'ConfusionMatrix', zeros(nCats,nCats), ...
	'Quant', [], ...
	'SolverOptions', SolverOptions, ...
	'SolverOutput', [] );
pLow = 0.25/nCats/(ntSamp/Nquantiles);
pHigh = 1-(nCats-1)*pLow;
for iCat=1:nCats
	out.Quant = [out.Quant, ...
		struct(...
			'Class', 'Mcl_QuantizedDv', ...
			'Nquantiles', int32(Nquantiles), ...
			'Dv', zeros(Nquantiles,1), ...
			'DvBinSep', zeros(Nquantiles+1,1), ...
			'Weight', zeros(Nquantiles,1), ...
			'Pc', zeros(Nquantiles,1), ...
			'PcMono', zeros(Nquantiles,1), ...
			'PcMonoLim', zeros(Nquantiles,1), ...
			'pLow', pLow, ...
			'pHigh', pHigh, ...
			'hMin', -log(1-(nCats-1)*pLow) / log(nCats)...
		) ];
end
%	Apply transforms of the spatial dimensions
if ~isempty(Xtrans)
	for iCat=1:nCats
		Xcell{iCat} = Mcl_Trfm_Cols(Xtrans,true,Xcell{iCat});
	end
end
%	Calculate initial representation in the structure.
Mcl_Exemplar_Init(out, Xcell);
if ~isempty(SolverOptions)
	out = Mcl_Exemplar_Train(out, true);
end