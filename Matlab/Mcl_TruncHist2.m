function varargout = Mcl_TruncHist2(x,y,ax,sz,blurLevel)

% function [out, ax] = Mcl_TruncHist2(x,y,ax,sz)
% Computes a 2-D truncated histogram
% x:  The x-coordinate of each data point
% y:  The y-coordinate of each data point
% [ax]:  The axes of the histogram [xmin, xmax, ymin, ymax], defaults to the result of function Mcl_Truncation.
% [sz]:  The size of the histogram [height,width] or equivalently [row,col], defaults to [50, 50].
% [blurLevel]:  (0 to 1)  The approximate proportion of high-pass energy removed by the low-pass filter.  A Gaussian blur
%	is applied by a 2D Gaussian function with two standard deviations (one for x, the other for y).  The standard
%	deviations are equal to:
%		blurLevel / 6 * [ax(2)-ax(1), ax(4)-ax(3)]
%-------------------
% out: The 2D histogram counting the number of samples in each bin.
% ax:  The axes of the histogram [xmin, xmax, ymin, ymax]

if nargin < 3
	ax = [];
	sz = [50,50];
	blurLevel = 0;
elseif nargin < 4
	sz = [50,50];
	blurLevel = 0;
elseif nargin < 5
	blurLevel = 0;
end

if ischar(blurLevel)
	blurLevel = 0;
end
if blurLevel < 0
	blurLevel = 0;
elseif blurLevel >= 1
	blurLevel = 0.15;
end

if isempty(ax)
	mu = mean([x,y]);
	sd = std([x,y]);
	skew = sum([(x-mu(1))/sd(1), (y-mu(2))/sd(2)].^3)/(length(x)-3);
	kurt = sum([(x-mu(1))/sd(1), (y-mu(2))/sd(2)].^4)/(length(x)-3);
	[xmin, xmax] = Mcl_Truncation( skew(1), kurt(1), Mcl_Quantiles(x,[.01, .05, .1, .9, .95, .99]) );
	[ymin, ymax] = Mcl_Truncation( skew(2), kurt(2), Mcl_Quantiles(y,[.01, .05, .1, .9, .95, .99]) );
	ax = [xmin, xmax, ymin, ymax];
end

xrng = ax(2)-ax(1);
yrng = ax(4)-ax(3);
xmin = ax(1);
ymin = ax(3);

out = zeros(sz);
for i=1:length(x)
	ix = round(min(sz(2),max(1, 1+(x(i)-xmin)/xrng*(sz(2)-1) )));
	iy = round(min(sz(1),max(1, 1+(y(i)-ymin)/yrng*(sz(1)-1) )));
    iy = sz(2)-iy+1;
	out(iy,ix) = out(iy,ix) + 1;
end

%	Blur the histogram with a low-pass raised cosine window.
if blurLevel > 0
	xLen = sz(2);
	yLen = sz(1);
	
	%	Create a blur kernel with an even number of rows and columns
	bxLen = 2*floor(0.5*xLen*blurLevel);
	byLen = 2*floor(0.5*xLen*blurLevel);
	bxc = (bxLen+1)/2;
	byc = (byLen+1)/2;
	bxs = bxLen/6;
	bys = bxLen/6;
	[bxMesh, byMesh] = meshgrid(...
		(((1:bxLen)-bxc)/bxs).^2, ...
		(((1:byLen)-byc)/bys).^2 );
	bxMesh = exp( -0.5* (bxMesh + byMesh) );
	
	%	The output is normalized to have the mass that was originally output.
	outSum = sum(out(:));
	out = conv2(out,bxMesh,'same');
	out = out * (outSum / sum(out(:)));
end

varargout = {out, ax};