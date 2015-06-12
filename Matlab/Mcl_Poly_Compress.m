function o = Mcl_Poly_Compress(o)

% function o = Mcl_Poly_Compress(o)
% This compresses some information inside an Mcl_Poly object (to consume less storage area).  This reduces some 64-bit floating point numbers
%	down to 32-bit and clips o.X to only include columns for the linear terms.
% -------------------------------------------------------------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------
% o
%	An Mcl_Poly structure.
%--------------------------------------------------------------------------------------------------------------------------------
% OUTPUT
%----------------------------------------------------
% o.X (single 2D array, o.Ntsamp * o.Ndims)
%	Polynomial expansion columns are clipped.  Cast from double (64-bit) precision to single (32-bit) precision.
% o.Dv
%	Cast from double (64-bit) precision to single (32-bit) precision.
% o.P
%	Cast from double (64-bit) precision to single (32-bit) precision.
% o.H
%	Cast from double (64-bit) precision to single (32-bit) precision.
%--------------------------------------------------------------------------------------------------------------------------------

if ~isstruct(o) || ~isfield(o, 'Class')
	error('You must provide an Mcl_Poly or Mcl_Poly_Tngl structure as input.');
end
if isequal(o.Class, 'Mcl_Poly')
	o.X		= single(o.X(:,1:o.Ndims));
	o.Dv	= single(o.Dv);
	o.P		= single(o.P);
	o.H		= single(o.H);
elseif isequal(o.Class, 'Mcl_Poly_Tngl')
	for i = 1:numel(o.UniClassifier)
		o.UniClassifier(i) = Mcl_Poly_Compress(o.UniClassifier(i));
	end
	for i = 1:numel(o.PairClassifier)
		o.PairClassifier(i) = Mcl_Poly_Compress(o.PairClassifier(i));
	end
else
	error('You must provide an Mcl_Poly or Mcl_Poly_Tngl structure as input.');
end