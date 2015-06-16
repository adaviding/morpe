function out = Mcl_Poly_Tngl_Ctor(Xcell, Xtrans, Rank, ForceEqualPriors, Nquantiles, DisplayModulus, OptNapproaches)

% function out = Mcl_Poly_Tngl_Ctor(Xcell, Xtrans, Rank, ForceEqualPriors, Nquantiles, DisplayModulus, OptNapproaches)
%	This function constructs a "Triangle" of Ndims*(Ndims-1)/2 classifiers.  Each classifier considers only two spatial dimensions.
%	There are Ndims*(Ndims-1)/2 unique pairs of dimensions, hence that determines the number of classifiers in the triangle.
%	Note that if SolverOptions are not supplied, then this method will merely construct the classifiers without training them.
%	This triangle is useful for generating triangular scatterplots including classifier performance.
%--------------------------------------------------------------------------------------------------------------------------
% NOMENCLATURE
% Ndims
%	The number of spatial dimensions of the data.  Each iCat element of Xcell, Xcell{iCat}, has Ndims columns.
% Npairs
%	The number of dimension pairs, equal to Ndims*(Ndims-1)/2.
%--------------------------------------------------------------------------------------------------------------------------
% INPUT:  See Mcl_Poly_Ctor
%--------------------------------------------------------------------------------------------------------------------------
% OUTPUT
% out.Npairs(int32 scalar)
%	The number of dimension pairs.
% out.DimPairs (int32 2D array: Npairs * 2)
%	Each iPair row out.DimPairs(iPair,:) lists the pair of dimensions considered.
% out.Classifiers (Mcl_Poly matlab structure vector: Npairs)
%	A classifier constructed for each pair of dimensions by calling the constructor Mcl_Poly_Ctor.
%--------------------------------------------------------------------------------------------------------------------------

Ncats = numel(Xcell);
Ndims = size(Xcell{1},2);

%	Implement default behavior
if ~exist('DisplayModulus', 'var') || isempty(DisplayModulus)
	DisplayModulus = 250;
end
if ~exist('Nquantiles', 'var') || isempty(Nquantiles)
	Nquantiles = min(20,ceil(ntSamp/20));
end
if ~exist('ForceEqualPriors', 'var') || isempty(ForceEqualPriors)
	ForceEqualPriors = true;
end
if ~exist('Rank', 'var') || isempty(Rank)
	Rank = 1;
end
%	Handle solver options
if ~exist('OptNapproaches', 'var') || isempty(OptNapproaches)
	OptNapproaches = 0;
end
if OptNapproaches<=0
	trainStr = '.  Classifier is not yet trained.';
else
	trainStr = '.  Training...';
end

out = struct(...
	'Class', 'Mcl_Poly_Tngl', ...
	'Ndims', int32(0.001+Ndims), ...
	'IdDims', [], ...
	'UniClassifier', [], ...
	'Npairs', int32(0.001+Ndims*(Ndims-1)/2), ...
	'DimPairs', [], ...
	'PairClassifier', [], ...
	'DisplayModulus', DisplayModulus );

XsubTrans = [];
Xsub = cell(size(Xcell));
out.DimPairs = zeros(out.Npairs,2,'int32');
iPair = int32(0);
for iDim=int32(1:Ndims)
	%-------------------------------------------------------------------------------------------------------------------------------
	% Univariate
	%-------------------------------------------------------------------------------------------------------------------------------
	for iCat=1:Ncats
		Xsub{iCat} = Xcell{iCat}(:,iDim);
	end
	if ~isempty(Xtrans)
		XsubTrans = Xtrans(iDim);
	end
	if out.DisplayModulus>0
		disp('--------------------------------------------------------------------------------');
		disp(['Mcl_Poly_Tngl_Ctor: Uni ' num2str(iDim) ' classifier being constructed' trainStr]);
	end
	if isempty(out.UniClassifier)
		out.UniClassifier = Mcl_Poly_Ctor(Xsub, XsubTrans, Rank, ForceEqualPriors, Nquantiles, DisplayModulus, OptNapproaches);
	else
		out.UniClassifier(iDim) = Mcl_Poly_Ctor(Xsub, XsubTrans, Rank, ForceEqualPriors, Nquantiles, DisplayModulus, OptNapproaches);
	end
	if OptNapproaches>0
		if out.DisplayModulus>0
			disp(['Mcl_Poly_Tngl_Ctor: Uni ' num2str(iDim) ' classifier trained.  [h, Acc] = [' num2str(out.UniClassifier(iDim).h) ',' num2str(out.UniClassifier(iDim).Acc) ']']);
			disp('--------------------------------------------------------------------------------');
		end
	end
	%-------------------------------------------------------------------------------------------------------------------------------
end
for iDim=int32(1:(Ndims-1))
	for jDim= (iDim+1):int32(Ndims)
		%-------------------------------------------------------------------------------------------------------------------------------
		% Bivariate
		%-------------------------------------------------------------------------------------------------------------------------------
		iPair = iPair+1;
		out.DimPairs(iPair,:) = [iDim,jDim];
		for iCat=1:Ncats
			Xsub{iCat} = Xcell{iCat}(:,out.DimPairs(iPair,:));
		end
		if ~isempty(Xtrans)
			XsubTrans = Xtrans([iDim, jDim]);
		end
		if out.DisplayModulus>0
			disp('--------------------------------------------------------------------------------');
			disp(['Mcl_Poly_Tngl_Ctor: Pair ' num2str(iPair)  ' (' num2str(iDim) ',' num2str(jDim) ') classifier being constructed' trainStr]);
		end
		if isempty(out.PairClassifier)
			out.PairClassifier = Mcl_Poly_Ctor(Xsub, XsubTrans, Rank, ForceEqualPriors, Nquantiles, DisplayModulus, OptNapproaches);
		else
			out.PairClassifier(iPair) = Mcl_Poly_Ctor(Xsub, XsubTrans, Rank, ForceEqualPriors, Nquantiles, DisplayModulus, OptNapproaches);
		end
		if OptNapproaches>0
			%	Display messages to command window
			if out.DisplayModulus>0
				disp(['Mcl_Poly_Tngl_Ctor: Pair ' num2str(iPair)  ' (' num2str(out.DimPairs(iPair,1)) ',' num2str(out.DimPairs(iPair,2)) ') classifier trained.  [h, Acc] = [' num2str(out.PairClassifier(iPair).h) ',' num2str(out.PairClassifier(iPair).Acc) ']']);
				disp('--------------------------------------------------------------------------------');
			end
			%--------------------------------------------
			% Check to see if a unidimensional classifier has outperformed a bivariate classifier.  This is theoretically impossible, but it
			%	can happen if (1) the Cube parameter led to poor optimization or (2) the optimization algorithm got stuck in a local minimum.
			%--------------------------------------------
			iUni = iDim; % iUni is the better performing of the univariate classifiers.
			jUni = jDim;
			ipwc = 0;
			if out.UniClassifier(jDim).h < out.UniClassifier(iDim).h
				iUni = jDim;
				jUni = iDim;
				ipwc = 1;
			end
			if out.UniClassifier(iUni).h < out.PairClassifier(iPair).h
				%--------------------------------------------------------------------------------------------------
				% A univariate model has out-performed the bivariate.
				%--------------------------------------------------------------------------------------------------
				unc = out.UniClassifier(iUni).Ncoeff;
				%	Express the dimensions of each coefficient in terms of the pairwise model
				udims = out.UniClassifier(iUni).CoeffDims;
				udims(udims(:)==0)=ipwc;
				%	A map from the bivariate coefficients to the coefficients of the full model.
				cMap = zeros(unc,1,'int32');
				for iMap=1:unc
					for iCoeff=1:out.PairClassifier(iPair).Ncoeff
						if all(udims(iMap,:)==out.PairClassifier(iPair).CoeffDims(iCoeff,:))
							cMap(iMap)=iCoeff;
							break;
						end
					end
				end
				if any(cMap==0)
					error('Assertion failed.  This code needs some work.');
				end
				%	Map the coefficients from univariate to bivariate.
				out.PairClassifier(iPair).wInit(:) = 0;
				for iPoly=1:out.PairClassifier(iPair).Npoly
					out.PairClassifier(iPair).wInit(cMap,iPoly) = out.UniClassifier(iUni).wOptimized(:,iPoly);
				end
				%	Retrain bivaraite classifier.
				if out.DisplayModulus>0
					disp(['Bivariate classifier ('  num2str(out.DimPairs(iPair,1)) ',' num2str(out.DimPairs(iPair,2)) ') is less powerful than a univariate child (' num2str(iUni) ').  Re-initialized as child.  Retraining...']);
				end
				out.PairClassifier(iPair) = Mcl_Poly_Train(out.PairClassifier(iPair));
				if out.DisplayModulus>0
					disp(['Mcl_Poly_Tngl_Ctor: Pair ' num2str(iPair)  ' (' num2str(out.DimPairs(iPair,1)) ',' num2str(out.DimPairs(iPair,2)) ') classifier re-trained.  [h, Acc] = [' num2str(out.PairClassifier(iPair).h) ',' num2str(out.PairClassifier(iPair).Acc) ']']);
					disp('--------------------------------------------------------------------------------');
				end
			end
			%--------------------------------------------
		end
		%-------------------------------------------------------------------------------------------------------------------------------
	end
end