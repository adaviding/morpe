function X = Mcl_Trfm_Cols(o,goForward,X)

% function x = Mcl_Trfm_Cols(o,goForward,X)
% This function transforms the columns of input X using the transform specified by the input Mcl_Trfm structures.
%
% The Mcl_Trfm structure is used by the Mcl classifiers to label and transform spatial dimensions that require
%	transformation between two domains.  One domain is usually based on the data's raw unit, but such data is sometimes
%	skewed or distributed in a way that is troublesome for Mcl classifiers.  These classifiers work best when the data
%	is distributed as Gaussian as possible.  Thus, the Mcl_Trfm object specifies a univariate and usually isotonic
%	transform that can be used to map data from one domain into another.  Thesecond domain is defined by the transform
%	(specified by an Mcl_Trfm structure) and a reverse transform (also specified by the Mcl_Trfm structure).
%-------------------------------------------------------------------------------------------------------
% NOMENCLATURE
%-------------------------------------------------------------------------------------------------------
% nDims
%	The number of spatial dimensions in the data.  Equal to size(X,2) and numel(o).
% nSamp
%	The number of datum to be transformed.  Equal to size(X,1).
%-------------------------------------------------------------------------------------------------------
% INPUT
%-------------------------------------------------------------------------------------------------------
% o (matlab structure Mcl_Trfm or Mcl_Axis vector: nDims)
%	A structure that defines the type of transform to apply to each column of the input X.  This structure can be constructed
%	using Mcl_Trfm_Ctor.m.
% goForward (logical scalar)
%	If true, the transform is applied; if false, the inverse transform is applied.
% X (double 2D array: nSamp, nDims)
%	The dataset to be transformed.  Each iDim column is transformed according to o(iDim).
%-------------------------------------------------------------------------------------------------------
% OUTPUT
%-------------------------------------------------------------------------------------------------------
% X
%	The transformed dataset.
%-------------------------------------------------------------------------------------------------------

if iscell(X)
	for iCell=1:numel(X)
		X{iCell} = Mcl_Trfm_Cols(o,goForward,X{iCell});
	end
else
	if numel(o) < size(X,2)
		error('The number of transforms, numel(o), must equal the number of columns in X.');
	end
	for iCol=1:numel(o)
		if numel(o)>=iCol && ~isempty(o(iCol))
			X(:,iCol) = Mcl_Trfm_Vec(o(iCol),goForward,X(:,iCol));
		end
	end
end