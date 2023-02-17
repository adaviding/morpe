function out = Mcl_AccuracyN(x, qSize, forceEqualPriors, pThreshold)

% function OUT = Mcl_AccuracyN(X, [QSIZE], [FORCEEQUALPRIORS], [PTHRESHOLD])
%	Gives the accuracy and lamda for length(X) categories of samples.  This is based on quantiling
%	along the unidimensional measure.  Note, each X has only 1 column (i.e. it represents a unidimensional coordinate).
%	The quantiling uses quantiles no smaller than QSIZE (but as small as possible given this restriction).
%	Default QSIZE = 50.
%	Default FORCEEQUALPRIORS = true;
%	Default PTHRESHOLD = 0.5;
%--------------------------------------------------------------
% OUT.Acc	The accuracy
% OUT.Rsq	The value of R^2
% OUT.Lamda	The value of lamda
% OUT.qTween  The value that separates quantiles.  (Note, if there are nq quantiles, this vector's length is nq+1)
% OUT.qMid  The central value for this quantile.  (Note, if there are nq quantiles, this vector's length is nq)
% OUT.p_q  The probability of each category given the quantile.  (Note, if there are nq quantiles, this matrix's size is [nq, length(X)])
%	Note, probabilities are never 0 or 1 because an event occurring 1/2 time per quantile was added with equal probability for each category.
%---------------------------------------------------------------

out = struct(...
	'Acc', [], ...
	'Rsq', [], ...
	'Lamda', [], ...
	'qTween', [], ...
	'qMid', [], ...
	'p_q', [] );

nCats = length(x);
nStim = zeros(nCats,1);
for i=1:nCats
	[ns, nd] = size(x{i});
	nStim(i) = ns;
end

if nargin<2
	qSize = 50;
end
if nargin < 3
	forceEqualPriors = true;
end
if nargin < 4
	pThreshold = 0.5;
end
if pThreshold < 0.5 || pThreshold > 1
	error('pThreshold must be in the range [0.5, 1]');
end
	

%	Create weights to equalize the priors
if forceEqualPriors
	catWeight = mean(nStim)./nStim;
else
	catWeight = ones(size(nStim));
end

%	A sortable table
%	1 row per sample
%	column 1 = value of unidimensional variable 'x'
%	column 2 = category label
%	column 3 = catWeight
%	column 4 = cumulative catWeight (after sort)
xx = [];
for icat = 1:nCats
	xx = [xx; x{icat}, repmat([icat catWeight(icat) 0], [nStim(icat), 1])];
end

%	Sort by unidimensional value 'x'
xx = sortrows(xx, 1);

%	Compute cumulative catWeight
xx(:,4) = cumsum(xx(:,3));

%	'il' is the last row index
[il, dum] = size(xx);
wsum = xx(il,4);

%	Now, column 4 will be switched from a cumulative catWeight to something where ceil(xx(:,4)) will give the quantile number
nq = ceil(wsum/qSize);%	the number of quantiles
xx(:,4) = ceil(xx(:,4)/(wsum+0.001)*nq);

%	The total catWeight for samples from each category in each quantile
qLims = zeros(nq,2);
wq = zeros(nq, nCats);
ib = 0;
for iq = 1:nq
	iiq = find(xx(:,4)==iq);
	qLims(iq,1) = xx(iiq(1),1);
	qLims(iq,2) = xx(iiq(numel(iiq)),1);
	%	Count the total catWeight for each category
	for icat = 1:nCats
		wq(iq,icat) = sum( (xx(iiq,2)==icat).*(xx(iiq,3)) );
	end
end

%	The catWeight matrix is augmented to count events that occur 1/2 catWeighted times per quantile with probability of each category of 1/nCats.
%	The event is made to have equal probability for each category (so that no probabilities can be 0 or 1).
wq = wq + 1/2/nCats;
W = sum(wq, 2);  %	The total catWeight in each quatile

%	Now a probability of each category given the quantile is computed.
out.p_q = wq;
for iq = 1:nq
	%	Calculate the probabilities for this quantile
	pq = out.p_q(iq,:)/W(iq);
	out.p_q(iq,:) = pq;
	%	Get the second-highest probability, store it in wq(iq,1).
	pq = sort(pq);
	wq(iq,1) = pq(nCats-1);
end

%----------------------------------------------------------------------------------------------------
%	Now check for pThreshold (remove binning noise from performance estimates).
%----------------------------------------------------------------------------------------------------
if pThreshold>0.5
	pqPassedNoise = false(nq,1);
	if forceEqualPriors
		priors = ones(nCats,1)/nCats;
	else
		priors = nStim/(sum(nStim));
	end
	%	Use wq(:,1) as memory, the total weight for the second-most-popular category
	wq(:,1) = W.*wq(:,1);
	for iCat=1:nCats
		%	Use wq(:,2) as memory, the total weight for the current category.
		wq(:,2) = W.*out.p_q(:,iCat);
		biProb = betainc(0.5, wq(:,1), wq(:,2));
		pqPassedNoise = pqPassedNoise | (biProb>pThreshold);
	end
	for iCat=1:nCats
		out.p_q(~pqPassedNoise,iCat) = priors(iCat);
	end
end
%----------------------------------------------------------------------------------------------------

% The proportion correct in each quantile
P = max(out.p_q, [], 2);

out.Acc = sum(W.*P)/sum(W);
out.Rsq = sum(sum(wq.*out.p_q)) / sum(W);
out.Lamda = exp(  sum(sum(wq.*log(out.p_q))) / sum(W)  );
out.qTween = [qLims(1,1); (qLims(1:(nq-1),2)+qLims(2:nq,1))/2; qLims(nq,2)];
out.qMid = (qLims(:,1)+qLims(:,2))/2;