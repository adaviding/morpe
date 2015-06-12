function varargout = Mcl_Lamda(ahx, bhx, minPerBin, forceEqualPriors)

% function [lamda, rsq, hxBins , paBins, nBins, hxBinCenters] = Mcl_Lamda(ahx, bhx, [minPerBin], [forceEqualPriors])
%	Returns a histogram analysis of h(x) with respect to the relative densities of categories 'a' and 'b'.
%	The calculation is based on bin sizes no smaller than 'minPerBin' and assumes the categorization rule that
%	category 'a' is assigned when h(x)>0, otherwise category 'b'.
% ahx:	Samples of h(x) from category 'a' (a 1-d vector)
% bhx:  Samples of h(x) from category 'b' (a 1-d vector)
% [minPerBin]:  The smallest number of samples allowed in each bin.  Must be at least 10.  default minPerBin=50.
% -------------
% lamda:  The value of lamda (a single double-precision number) representing the 
%	geometric mean likelihood of correctly determining the probability that a given
%	value of h(x) belongs to category "a" or "b".
% hxBins:  A 1-d vector giving the h(x) value in-between bins.
%	The first bin is unbounded, therefore hxBins(1) always equals -inf.
%	The last bin is also unbounded, therefore hxBins(length(hxBins)) always equals +inf.
%	The number of bins is equal to length(hxBins)-1.
%	The value h(x)=0 is guaranteed to be one of the bin separators, therefore, one member of hxBins will always equal zero (this will be near the center).
% paBins:  The probability of category 'a' for each bin.
% nBins:  The number of h(x) samples (from both categories) in each bin.
% hxBinCenters:  Gives the value of h(x) at the center of each bin.  This is ideal for executing -->  plot(hxBinCenters, paBins);
% rsq:  The value of R^2 (the proportion of variance accounted for)

sz = size(ahx);
if sz(2)>sz(1)
	ahx = ahx';
end
sz = size(bhx);
if sz(2)>sz(1)
	bhx = bhx';
end

sz = size(ahx);
if sz(2)>1
	error('ahx must be a 1-d vector');
end
sz = size(bhx);
if sz(2)>1
	error('bhx must be a 1-d vector');
end

if nargin<3
	minPerBin=50;
elseif isempty(minPerBin) | isnan(minPerBin)
	minPerBin=50;
else
	if length(minPerBin)>1
		error('minPerBin must be a positive integer greater than 10.  e.g. minPerBin=100.');
	end
	if minPerBin<10
		error('minPerBin must be a positive integer greater than 10.  e.g. minPerBin=100.');
	end
	minPerBin = round(minPerBin);
end
if length(ahx)<minPerBin | length(bhx)<minPerBin
	error('ahx and/or bhx are not big enough.');
end
if nargin<4
	forceEqualPriors = true;
elseif isempty(forceEqualPriors)
	forceEqualPriors = true;
end

%	Count how many samples exist from each category
na = length(ahx);
nb = length(bhx);

%	Weights equalize the priors when na~=nb
if forceEqualPriors
	aWeight = 0.5*(nb+na)/na;
	bWeight = 0.5*(nb+na)/nb;
else
	aWeight = 1;
	bWeight = 1;
end

%	Mix and sort the samples by h(x)
data = sortrows( [...
	ahx, ones(size(ahx)); ...
	bhx, zeros(size(bhx))], 1 );
[len, dum] = size(data);

%	Find the first sample greater than zero
imin = 1;
imax = len;
if data(imin,1)>0
	iZero = imin;
elseif data(imax,1)<=0
	iZero = imax;
else
	while (imax-imin)>0
		itry = round((imax+imin)/2);
		if data(itry,1)>0
			imax = itry;
		else
			imin = itry;
		end
		if (imax-imin)==1
			if data(imin,1)>0
				itry = imin;
			else
				itry = imax;
			end
			imin = itry;
			imax = itry;
		end
	end
	iZero = itry;
end

nSmall = iZero-1;
nBig = len-iZero;

hxSmall = data(1:nSmall,:); %	will be categorized as 'b'
hxBig  = data(iZero:len,:);  %	will be categorized as 'a'

%	Specify the number of each kind of bin (done separately for bins on the small and large side of the criterion)
nSmallBins = floor(nSmall/minPerBin);
nBigBins = floor(nBig/minPerBin);

if nSmallBins>0
	%	Specify the size of each bin
	nPerSmall = repmat( floor(nSmall/nSmallBins), [1 nSmallBins]);
	nPerSmall(1) = nSmall-sum(nPerSmall(2:nSmallBins)); % The end bin picks up any extra samples (if necessary)
	iSmallBins = cumsum(nPerSmall);  %	the index of small bin endings.
	if nSmallBins == 1
		iSmallBinCenters = round(nSmall/2);
	else
		iSmallBinCenters = [iSmallBins(1:(nSmallBins-1))-round(diff(iSmallBins)/2), round(nSmall- nPerSmall(length(nPerSmall))/2 )];  % the index of small bin centers
	end
	iSmallBins = iSmallBins(1:(nSmallBins-1)); % the last bin ending is removed
	xSmallBins = [(hxSmall(iSmallBins,1)+hxSmall(iSmallBins+1,1))/2]; % the hx coordinate after each iSmallBins (halfway to the next bin)
	xSmallBinCenters = hxSmall(iSmallBinCenters,1);
	if length(iSmallBinCenters)==1 & iSmallBinCenters==0
		xSmallBinCenters = [];
	else
		xSmallBinCenters = hxSmall(iSmallBinCenters,1);
	end

	if nSmallBins==1
		ctaSmallBins = sum( (hxSmall(:,2)==1) * aWeight);
		ctbSmallBins = sum( (hxSmall(:,2)==0) * bWeight);
	else
		ctaSmallBins = zeros(size(nPerSmall));
		ctbSmallBins = zeros(size(nPerSmall));

		ctaSmallBins(1) = sum( (hxSmall(1:iSmallBins(1),2)==1) * aWeight );
		ctbSmallBins(1) = sum( (hxSmall(1:iSmallBins(1),2)==0) * bWeight );

		for i=2:(nSmallBins-1)
			ctaSmallBins(i) = sum( (hxSmall((iSmallBins(i-1)+1):iSmallBins(i),2)==1) * aWeight );
			ctbSmallBins(i) = sum( (hxSmall((iSmallBins(i-1)+1):iSmallBins(i),2)==0) * bWeight );
		end

		ctaSmallBins(nSmallBins) = sum( (hxSmall((iSmallBins(nSmallBins-1)+1):nSmall,2)==1) * aWeight );
		ctbSmallBins(nSmallBins) = sum( (hxSmall((iSmallBins(nSmallBins-1)+1):nSmall,2)==0) * bWeight );
	end

	paSmallBins = (0.25+ctaSmallBins)./(0.5+ctaSmallBins+ctbSmallBins);
end

if nBigBins>0
	nPerBig = repmat( floor(nBig/nBigBins), [1 nBigBins] );
	nPerBig(nBigBins) = nBig-sum(nPerBig(1:(nBigBins-1))); % The end bin picks up any extra samples (if necessary)

	iBigBins = cumsum(nPerBig);  %	the index of big bin endings
	if nBigBins==1
		iBigBinCenters = round(nBig/2);
	else
		iBigBinCenters = [iBigBins(1:(nBigBins-1))-round(diff(iBigBins)/2), round(nBig- nPerBig(length(nPerBig))/2 )];  % the index of big bin centers
	end
	iBigBins = [iBigBins(1:(nBigBins-1))];
	xBigBins = [(hxBig(iBigBins,1)+hxBig(iBigBins+1,1))/2];
	if length(iBigBinCenters)==1 & iBigBinCenters==0
		xBigBinCenters = [];
	else
		xBigBinCenters = hxBig(iBigBinCenters,1);
	end

	if nBigBins==1
		ctaBigBins = sum( (hxBig(:,2)==1) * aWeight);
		ctbBigBins = sum( (hxBig(:,2)==0) * bWeight);
	else
		ctaBigBins = zeros(size(nPerBig));
		ctbBigBins = zeros(size(nPerBig));

		ctaBigBins(1) = sum( (hxBig(1:iBigBins(1),2)==1) * aWeight );
		ctbBigBins(1) = sum( (hxBig(1:iBigBins(1),2)==0) * bWeight );

		for i=2:(nBigBins-1)
			ctaBigBins(i) = sum( (hxBig((iBigBins(i-1)+1):iBigBins(i),2)==1) * aWeight );
			ctbBigBins(i) = sum( (hxBig((iBigBins(i-1)+1):iBigBins(i),2)==0) * bWeight );
		end

		ctaBigBins(nBigBins) = sum( (hxBig((iBigBins(nBigBins-1)+1):nBig,2)==1) * aWeight );
		ctbBigBins(nBigBins) = sum( (hxBig((iBigBins(nBigBins-1)+1):nBig,2)==0) * bWeight );
	end

	paBigBins = (0.25+ctaBigBins)./(0.5+ctaBigBins+ctbBigBins);
end
%	Make sure that no two bins subtend the same range.
if nSmallBins > 0
	iSmall = 1;
	while iSmall < length(xSmallBinCenters)
		if xSmallBinCenters(iSmall)==(xSmallBinCenters(iSmall+1))
			%	Add information from bin iSmall to bin iSmall+1.
			paSmallBins(iSmall+1) = (paSmallBins(iSmall)*nPerSmall(iSmall)+paSmallBins(iSmall+1)*nPerSmall(iSmall+1))/(nPerSmall(iSmall)+nPerSmall(iSmall+1));
			nPerSmall(iSmall+1) = (nPerSmall(iSmall)+nPerSmall(iSmall+1));
			%	Combine the bins by removing iSmall
			paSmallBins = [paSmallBins(1:(iSmall-1)), paSmallBins( (iSmall+1):length(paSmallBins) )];
			nPerSmall = [nPerSmall(1:(iSmall-1)), nPerSmall( (iSmall+1):length(nPerSmall) )];
			xSmallBinCenters = [xSmallBinCenters(1:(iSmall-1)); xSmallBinCenters( (iSmall+1):length(xSmallBinCenters) )];
			%	Remove the separation between iSmall and iSmall + 1
			xSmallBins = [xSmallBins(1:(iSmall-1)); xSmallBins( (iSmall+1):length(xSmallBins) )];
		else
			iSmall = iSmall + 1;
		end
	end
end
%	Make sure that no two bins subtend the same range.
if nBigBins > 0
	iBig = 1;
	while iBig < length(xBigBinCenters)
		if xBigBinCenters(iBig)==(xBigBinCenters(iBig+1))
			%	Add information from bin iBig to bin iBig+1.
			paBigBins(iBig+1) = (paBigBins(iBig)*nPerBig(iBig)+paBigBins(iBig+1)*nPerBig(iBig+1))/(nPerBig(iBig)+nPerBig(iBig+1));
			nPerBig(iBig+1) = (nPerBig(iBig)+nPerBig(iBig+1));
			%	Combine the bins by removing iBig
			paBigBins = [paBigBins(1:(iBig-1)), paBigBins( (iBig+1):length(paBigBins) )];
			nPerBig = [nPerBig(1:(iBig-1)), nPerBig( (iBig+1):length(nPerBig) )];
			xBigBinCenters = [xBigBinCenters(1:(iBig-1)); xBigBinCenters( (iBig+1):length(xBigBinCenters) )];
			%	Remove the separation between iBig and iBig + 1
			xBigBins = [xBigBins(1:(iBig-1)); xBigBins( (iBig+1):length(xBigBins) )];
		else
			iBig = iBig + 1;
		end
	end
end

if nSmallBins==0
	varargout = {...
		[], ... lamda
		[], ... rsq
		[-Inf; xBigBins ; +Inf], ... hxBins
		Mcl_ForceNondecreasing(paBigBins'), ...paBins
		nPerBig', ... nBins
		xBigBinCenters }; % hxBinCenters
elseif nBigBins==0
	varargout = {...
		[], ... lamda
		[], ... rsq
		[-Inf; xSmallBins; +Inf], ... hxBins
		Mcl_ForceNondecreasing(paSmallBins'), ...paBins
		nPerSmall', ... nBins
		xSmallBinCenters }; % hxBinCenters
else
	varargout = {...
		[], ... lamda
		[], ... rsq
		[-Inf; xSmallBins; 0; xBigBins ; +Inf], ... hxBins
		Mcl_ForceNondecreasing([paSmallBins, paBigBins]'), ...paBins
		[nPerSmall, nPerBig]', ... nBins
		[xSmallBinCenters; xBigBinCenters] }; % hxBinCenters
end
%	Compute Rsq and Lamda
minBin = varargout{6}(1);
maxBin = varargout{6}(length(varargout{6}));
ahx(ahx<minBin) = minBin;
ahx(ahx>maxBin) = maxBin;
bhx(bhx<minBin) = minBin;
bhx(bhx>maxBin) = maxBin;
apa = interp1(varargout{6}, varargout{4}, ahx);
bpa = interp1(varargout{6}, varargout{4}, bhx);
varargout{2} = sum([apa*aWeight; (1-bpa)*bWeight])/(na+nb); %Rsq
varargout{1} = exp(sum([log(apa)*aWeight; log(1-bpa)*bWeight])/(na+nb)); %Lamda