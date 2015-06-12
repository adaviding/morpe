function varargout = Mcl_TruncHist1(x,ax,sz)

% function [binCenters, binHeights] = Mcl_TruncHist2(x,y,ax,sz)
% Computes a 1-D truncated histogram
% x:  The 1-D data point
% [ax]:  The axes of the histogram [xmin, xmax], defaults to the result of function Mcl_Truncation.
% [sz]:  The size of the histogram or number of bins.
%-------------------
% out: The 1D histogram with integer values counting the samples in each bin.
% ax:  The axes of the histogram [xmin, xmax]

if nargin < 2
	ax = [];
	sz = 50;
elseif nargin < 3
	sz = 50;
end

if isempty(ax)
	mu = mean(x);
	sd = std(x);
	skew = sum( ((x-mu(1))/sd(1)).^3 )/(length(x)-3);
	kurt = sum( ((x-mu(1))/sd(1)).^4 )/(length(x)-3);
	[xmin, xmax] = Mcl_Truncation( skew, kurt, Mcl_Quantiles(x,[.01, .05, .1, .9, .95, .99]) );
	ax = [xmin, xmax];
end

xrng = ax(2)-ax(1);
xmin = ax(1);
binCenters = ((0.5:sz)*xrng/sz+xmin)';

binHeights = zeros([sz, 1], 'int32');
for i=1:length(x)
	ix = round(min(sz,max(1, 1+(x(i)-xmin)/xrng*(sz-1) )));
	binHeights(ix) = binHeights(ix) + 1;
end

varargout = {binCenters, binHeights};