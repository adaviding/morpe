function out = Mcl_Exemplar_Tngl_Ctor(Xcell, Xtrans, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles, DisplayMessages, SolverOptions)

% function out = Mcl_Exemplar_Tngl_Ctor(Xcell, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles, [DisplayMessages])
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
% INPUT:  See Mcl_Exemplar_Ctor
%--------------------------------------------------------------------------------------------------------------------------
% OUTPUT
% out.Npairs(int32 scalar)
%	The number of dimension pairs.
% out.DimPairs (int32 2D array: Npairs * 2)
%	Each iPair row out.DimPairs(iPair,:) lists the pair of dimensions considered.
% out.Classifiers (Mcl_Exemplar matlab structure vector: Npairs)
%	A classifier constructed for each pair of dimensions by calling the constructor Mcl_Exemplar_Ctor.
%--------------------------------------------------------------------------------------------------------------------------

Ncats = size(Xcell);
Ndims = size(Xcell{1},2);

if ~exist('DisplayMessages', 'var')
	DisplayMessages = true;
end
if ~exist('SolverOptions', 'var')
	SolverOptions = [];
	trainStr = '.  Classifier is not yet trained.';
else
	trainStr = '.  Training...';
end

out = struct(...
	'Class', 'Mcl_Exemplar_Tngl', ...
	'Ndims', int32(0.001+Ndims), ...
	'UniClassifier', [], ...
	'Npairs', int32(0.001+Ndims*(Ndims-1)/2), ...
	'DimPairs', [], ...
	'PairClassifier', [], ...
	'DisplayMessages', DisplayMessages );

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
	if isempty(out.UniClassifier)
		out.UniClassifier = Mcl_Exemplar_Ctor(Xsub, XsubTrans, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles);
	else
		out.UniClassifier(iDim) = Mcl_Exemplar_Ctor(Xsub, XsubTrans, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles);
	end
	if out.DisplayMessages
		disp(['Mcl_Exemplar_Tngl_Ctor: Uni ' num2str(iDim) ' classifier constructed' trainStr]);
	end
	if ~isempty(SolverOptions)
		out.UniClassifier(iDim).SolverOptions = SolverOptions;
		out.UniClassifier(iDim).SolverOptions.wShrinkFactor = max(2,out.UniClassifier(iDim).SolverOptions.wShrinkFactor);
		out.UniClassifier(iDim).SolverOptions.wDiffMax = max(2,out.UniClassifier(iDim).SolverOptions.wDiffMax);
		out.UniClassifier(iDim).SolverOptions.Napproaches = 1;
		out.UniClassifier(iDim) = Mcl_Exemplar_Train(out.UniClassifier(iDim),true);
		if out.DisplayMessages
			disp(['Mcl_Exemplar_Tngl_Ctor: Uni ' num2str(iDim) ' classifier trained.  [h, Acc] = [' num2str(out.UniClassifier(iDim).h) ',' num2str(out.UniClassifier(iDim).Acc) ']']);
		end
	end
	%-------------------------------------------------------------------------------------------------------------------------------
end
for iDim=int32(1:(Ndims-1))
	%-------------------------------------------------------------------------------------------------------------------------------
	% Bivariate
	%-------------------------------------------------------------------------------------------------------------------------------
	for jDim= (iDim+1):int32(Ndims)
		iPair = iPair+1;
		out.DimPairs(iPair,:) = [iDim,jDim];
		for iCat=1:Ncats
			Xsub{iCat} = Xcell{iCat}(:,out.DimPairs(iPair,:));
		end
		if ~isempty(Xtrans)
			XsubTrans = Xtrans([iDim, jDim]);
		end
		if isempty(out.PairClassifier)
			out.PairClassifier = Mcl_Exemplar_Ctor(Xsub, XsubTrans, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles);
		else
			out.PairClassifier(iPair) = Mcl_Exemplar_Ctor(Xsub, XsubTrans, IdMethod, MaxNeighbors, ForceEqualPriors, Nquantiles);
		end
		if out.DisplayMessages
			disp(['Mcl_Exemplar_Tngl_Ctor: Pair ' num2str(iPair)  ' (' num2str(iDim) ',' num2str(jDim) ') classifier constructed' trainStr]);
		end
		if ~isempty(SolverOptions)
			%	Assign solver options to classifier
			out.PairClassifier(iPair).SolverOptions = SolverOptions;
			%	Train classifier
			out.PairClassifier(iPair) = Mcl_Exemplar_Train(out.PairClassifier(iPair),false);
			%	Display messages to command window
			if out.DisplayMessages
				disp(['Mcl_Exemplar_Tngl_Ctor: Pair ' num2str(iPair)  ' (' num2str(out.DimPairs(iPair,1)) ',' num2str(out.DimPairs(iPair,2)) ') classifier trained.  [h, Acc] = [' num2str(out.PairClassifier(iPair).h) ',' num2str(out.PairClassifier(iPair).Acc) ']']);
			end
			%--------------------------------------------
			% Check to see if a unidimensional classifier has outperformed a bivariate classifier.  This is theoretically impossible, but it
			%	can happen if (1) the Cube parameter led to poor optimization or (2) the optimization algorithm got stuck in a local minimum.
			%--------------------------------------------
			iUni = iDim; % iUni is the better performing of the univariate classifiers.
			if out.UniClassifier(jDim).h < out.UniClassifier(iDim).h
				iUni = jDim;
			end
			if out.UniClassifier(iUni).h < out.PairClassifier(iPair).h
				% Yes, the bivariate has underperformed the univariate.  Now we must re-optimize the bivariate starting from the univariate case.
				if iUni==iDim
					%	1 = Relevant dimension, 2 = Irrelevant
					out.PairClassifier(iPair).wInit(1,:) = out.UniClassifier(iUni).wOptimized; % relevant dimension
					out.PairClassifier(iPair).wInit(2,:) = out.PairClassifier(iPair).wInit(2,:) / 2^16;  % irrelevant dimension
				else
					%	2 = Relevant dimension, 1 = Irrelevant
					out.PairClassifier(iPair).wInit(1,:) = out.PairClassifier(iPair).wInit(2,:) / 2^16;  % irrelevant dimension
					out.PairClassifier(iPair).wInit(2,:) = out.UniClassifier(iUni).wOptimized; % relevant dimension
				end
				%	Retrain bivariate classifier.
				out.PairClassifier(iPair).SolverOptions.wDiffMax = 4;
				out.PairClassifier(iPair) = Mcl_Exemplar_Train(out.PairClassifier(iPair),false);
				if out.UniClassifier(iUni).h < out.PairClassifier(iPair).h
					%--------------------------------------------
					% The univariate case is still better.  Now we know it is the fault of the Cube, not the optimization procedure.  Most
					%	likely, this means that the other dimension is totally irrelevant and just adds noise.  Therefore, simply evaluate
					%	the univariate version of the bivariate classifier once and skip optimization (which relies on the Cube).
					%--------------------------------------------
					out.PairClassifier(iPair).wOptimized = out.PairClassifier(iPair).wInit;
					out.PairClassifier(iPair).Dv = out.UniClassifier(iUni).Dv;
					out.PairClassifier(iPair).P = out.UniClassifier(iUni).P;
					out.PairClassifier(iPair).H = out.UniClassifier(iUni).H;
					out.PairClassifier(iPair).h = out.UniClassifier(iUni).h;
					out.PairClassifier(iPair).Hcross = out.UniClassifier(iUni).Hcross;
					out.PairClassifier(iPair).Acc = out.UniClassifier(iUni).Acc;
					out.PairClassifier(iPair).ConfusionMatrix = out.UniClassifier(iUni).ConfusionMatrix;
					out.PairClassifier(iPair).Quant = out.UniClassifier(iUni).Quant;
					%--------------------------------------------
				end
				%	Display messages to the command window
				if out.DisplayMessages
					disp(['Mcl_Exemplar_Tngl_Ctor: Pair ' num2str(iPair)  ' (' num2str(out.DimPairs(iPair,1)) ',' num2str(out.DimPairs(iPair,2)) ') classifier trained.  [h, Acc] = [' num2str(out.PairClassifier(iPair).h) ',' num2str(out.PairClassifier(iPair).Acc) '] and it fell into the univariate case']);
				end
			end
			%--------------------------------------------
			%	Clip unneeded memory
			out.PairClassifier(iPair).Cube = [];
		end
	end
	%-------------------------------------------------------------------------------------------------------------------------------
end