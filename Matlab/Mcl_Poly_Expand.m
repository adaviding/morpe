function Xhigh = Mcl_Poly_Expand(Xlow, CoeffDims)

% function Xhigh = Mcl_Poly_Expand(Xlow, CoeffDims)
% Outputs the polynomial expansion Xhigh of data samples Xlow for the polynomial defined by CoeffDims.
% ------------------------------------------------------------------------------------------------------------------------------------
% INPUT
% ------------------------------------------------------------------
% Xlow (double 2D array: ntSamp * nDims)
%	Data samples.  Fits specifications of Mcl_Poly.X(:,1:Mcl_Poly.Ndims).  See the comments of Mcl_Poly constructor for more information.
% CoeffDims (int32 2D array: Ncoeff * Rank)
%	Coefficient dimensions.  Fits specifications of Mcl_Poly.CoeffDims.  See the comments of Mcl_Poly constructor for more information.
% ------------------------------------------------------------------------------------------------------------------------------------
% OUTPUT
% ------------------------------------------------------------------
% Xhigh (double 2D array: ntSamp * Ncoeff)
%	Polynomial expansion of data samples.  Data samples.  Fits specifications of Mcl_Poly.X.  See the comments of Mcl_Poly constructor
%		for more information.
% ------------------------------------------------------------------------------------------------------------------------------------

Ncoeff = size(CoeffDims,1);
Rank = size(CoeffDims,2);
Xhigh = ones(size(Xlow,1),Ncoeff);
for iCoeff=1:Ncoeff
	for iRank=1:Rank
		iDim = CoeffDims(iCoeff,iRank)+1;
		if iDim>0
			Xhigh(:,iCoeff) = Xhigh(:,iCoeff) .* Xlow(:,iDim);
		else
			break;
		end
	end
end