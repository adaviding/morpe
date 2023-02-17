function varargout = Mcl_Truncation(xskew, xkurt, xquantiles)

%function [xmin, xmax] = Mcl_Truncation(xskew, xkurt, xquantiles)
%	This function comes in handy when you need to make a truncated histogram.  This function picks a good range for the histogram's axis.
%---------------
% xskew: The skewness of a variable "x"
% xkurt:  The kurtosis of the variable "x".  This number is positive and sampling a normal distribution results in a mean of 3.
% xquantiles:  A vector length 6, the values of a variable "x" at six quantiles:  1%, 5%, 10%, 90%, 95%, 99%.
%---------------
% xmin:  The minimum of the truncated axis.
% xmax:  The maximum of the truncated axis.

w = exp(-abs(xskew)^(1/3)*xkurt^(1/4));
if xskew>0 % positively skewed, truncate the positive side
	xmin = xquantiles(1);
	if w<1
		xmax = xquantiles(6)*w + xquantiles(5)*(1-w);
	else
		xmax = xquantiles(5)*(1/w) + xquantiles(4)*(1-1/w);
	end
else % negatively skewed, truncate the negative side
	xmax = xquantiles(6);
	if w<1
		xmin = xquantiles(1)*w + xquantiles(2)*(1-w);
	else
		xmin = xquantiles(2)*(1/w) + xquantiles(3)*(1-1/w);
	end
end

varargout = {xmin, xmax};