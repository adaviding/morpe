function x = Mcl_Trfm_Vec(o,goForward,x)

% function x = Mcl_Trfm_Vec(o,goForward,x)
% This function transforms the vector of input using the transform specified by the input Mcl_Trfm structure.
%
% The Mcl_Trfm structure is used by the Mcl classifiers to label and transform spatial dimensions that require
%	transformation between two domains.  One domain is usually based on the data's raw unit, but such data is sometimes
%	skewed or distributed in a way that is troublesome for Mcl classifiers.  These classifiers work best when the data
%	is distributed as Gaussian as possible.  Thus, the Mcl_Trfm object specifies a univariate and usually isotonic
%	transform that can be used to map data from one domain into another.  Thesecond domain is defined by the transform
%	(specified by an Mcl_Trfm structure) and a reverse transform (also specified by the Mcl_Trfm structure).
%-------------------------------------------------------------------------------------------------------
% INPUT
%-------------------------------------------------------------------------------------------------------
% o (matlab structure Mcl_Trfm, scalar)
%	A structure that defines the type of transform to apply.  This structure can be constructed using Mcl_Trfm_Ctor.m.
%	Alternatively, this can be an Mcl_Axis structure constructed by Mcl_Axis_Ctor.
% goForward
%	If true, the transform is applied; if false, the inverse transform is applied.
% x
%	The vector to be transformed.
%-------------------------------------------------------------------------------------------------------
% OUTPUT
%-------------------------------------------------------------------------------------------------------
% x
%	The transformed vector.
%-------------------------------------------------------------------------------------------------------

if ~isstruct(o) || ~isfield(o,  'Class')
	error('The input o must be a structure of type Mcl_Trfm or Mcl_Axis.');
end
if isequal(o.Class, 'Mcl_Axis')
	x = Mcl_Trfm_Vec(o.Transform,goForward,x);
elseif isequal(o.Class, 'Mcl_Trfm')
    if isequal(o.Name, 'Affine')
		if goForward
			%	Subtract a value, then divide a value.
            x = (x - o.Vals(1)) / o.Vals(2);
		else
			%	Multiply a value, then add a value.
			x = x * o.Vals(2) + o.Vals(1);
		end
    elseif isequal(o.Name, 'FisherZ')
		if goForward
			%	Go from [-1, 1] to Z
			x = min(1,max(-1,x/o.Vals(2)));
			x = 0.5*log( (o.Vals(1)+x)./(o.Vals(1)-x) );
		else
			%	Go from Z to [-1, 1]
			x = 2*x*o.Vals(2);
			x = (exp(x)-o.Vals(1))./(exp(x)+o.Vals(1));
		end
	elseif isequal(o.Name, 'LogAbs')
		if goForward
			%	Go from [0, inf] to Z
			x = log10(abs(x+o.Vals(2))+o.Vals(1)); % Loss of sign
		else
			%	Go from Z to [0, inf]
			x = 10.^x-o.Vals(1)-o.Vals(2); % Loss of sign
		end
	elseif isequal(o.Name, 'SgnPow')
		if goForward
			%	Go from original domain to power domain
			x = sign(x).*abs(x).^o.Vals(1);
		else
			%	Go from power domain to original domain
			x = sign(x).*abs(x).^(1/o.Vals(1));
		end
	end
else
	error('The input o must be a structure of type Mcl_Trfm or Mcl_Axis.');
end