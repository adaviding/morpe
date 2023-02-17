function out = Mcl_FracDist_Sample(o, n)

% function out = Mcl_FracDist_Sample(o, n)
%	Creates n random samples from the fractured Gaussian distribution specified by o.
%-----------------------
% INPUT
%-----------------------
% o (scalar struct Mcl_FracDist)
%	An Mcl_FracDist structure specifying the parameters of the fractured Gaussian distribution.
% n (scalar int32)
%	The number of random samples desired.
%-----------------------
% OUTPUT
%-----------------------
% out (double 2D array, n * o.Ndims)
%	A list of n random samples from the fractured Gaussian distribution.
%-----------------------

rFrac = max(1,ceil(rand(n,1)*o.Nfrac));
out = randn(n,o.Ndims) + o.MuFrac(rFrac,:);