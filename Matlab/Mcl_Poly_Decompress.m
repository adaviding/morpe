function o = Mcl_Poly_Decompress(o)

% function o = Mcl_Poly_Decompress(o)
% This reverses the Mcl_Poly_Compress function.
% -------------------------------------------------------------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------
% o
%	An Mcl_Poly structure.
%--------------------------------------------------------------------------------------------------------------------------------
% OUTPUT
%----------------------------------------------------
% o.X (double 2D array, o.Ncoeff * o.Ndims)
%	Polynomial expansion columns are clipped.  Cast to double (64-bit) precision from single (32-bit) precision.
% o.Dv
%	Cast to double (64-bit) precision from single (32-bit) precision.
% o.P
%	Cast to double (64-bit) precision from single (32-bit) precision.
% o.H
%	Cast to double (64-bit) precision from single (32-bit) precision.
%--------------------------------------------------------------------------------------------------------------------------------

if ~isstruct(o) || ~isfield(o, 'Class')
	error('You must provide an Mcl_Poly or Mcl_Poly_Tngl structure as input.');
end
if isequal(o.Class, 'Mcl_Poly')
	o.X		= Mcl_Poly_Expand(double(o.X), o.CoeffDims);
	o.Dv	= double(o.Dv);
	o.P		= double(o.P);
	o.H		= double(o.H);
elseif isequal(o.Class, 'Mcl_Poly_Tngl')
	for i = 1:numel(o.UniClassifier)
		o.UniClassifier(i) = Mcl_Poly_Decompress(o.UniClassifier(i));
	end
	for i = 1:numel(o.PairClassifier)
		o.PairClassifier(i) = Mcl_Poly_Decompress(o.PairClassifier(i));
	end
else
	error('You must provide an Mcl_Poly or Mcl_Poly_Tngl structure as input.');
end