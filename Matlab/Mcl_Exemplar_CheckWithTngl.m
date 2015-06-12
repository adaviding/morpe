function oxFull = Mcl_Exemplar_CheckWithTngl(oxFull, oxTngl)

% function oxFull = Mcl_Exemplar_CheckWithTngl(oxFull, oxTngl)
% Consider two exemplar models:  oxFull utilizes multiple (>2) sample dimensions while oxPair only utilizes
%	a subset (2) of sample dimensions.  It is theoretically impossible for oxFull to perform worse than oxPair
%	assuming the exemplar models are equal in all other ways.  However, in practice (using Mcl), oxFull can
%	underperform for two reaosns:
%		(1) The optimization procedure failed to find the global entropy minumum for oxFull
%		(2) The optimization Cube misdirected the optimization procedure towards suboptimal parameters.
%	This function ensures that oxFull NEVER underperforms any of the inferior Univariate or Bivariate models,
%	all of which are contained within oxTngl.  If oxFull is found to have a higher entropy than any of these
%	models, then it is re-optimized with initial parameters that mimic the most powerful sub-dimensional
%	model in oxTngl.  If the model is still suboptimal, then the Cube is bad, and it is probably true that
%	the model oxFull should be approximately equivalent to the most powerful sub-dimensional model, so it
%	is just set to the sub-dimensional model.
%--------------------------------------------------------------------------------------------------------------------------
% INPUT and OUTPUT
% oxFull (matlab structure Mcl_Exemplar scalar)
%	An exemplar classifier constructed on more than 2 sample dimensions.
%--------------------------------------------------------------------------------------------------------------------------
% INPUT
% oxTngl (matlab structure Mcl_Exemplar_Tngl scalar)
%	This must be constructed from the same input dimensions as oxFull with all options equal to the options used for oxFull.
%--------------------------------------------------------------------------------------------------------------------------

%	Find the best performing sub-dimensional model.
%	(Note that bivariate models are already forced to outperform univarite ones, so don't search univariate models.)
iPair = 1;
for jPair=2:oxTngl.Npairs
	if oxTngl.PairClassifier(jPair).h < oxTngl.PairClassifier(iPair).h
		iPair = jPair;
	end
end

%	Compare
if oxFull.h > oxTngl.PairClassifier(iPair).h
	%--------------------------------------------------------------------------------------------------
	% A bivariate model has out-performed the full.
	%--------------------------------------------------------------------------------------------------
	% Use the starting parameters that effectively make the full classifier into the most powerful bivariate classifier.
	for iDim=1:oxFull.Ndims
		if iDim==oxTngl.DimPairs(iPair,1)
			oxFull.wInit(iDim,:) = oxTngl.PairClassifier(iPair).wInit(1,:);
		elseif iDim==oxTngl.DimPairs(iPair,2)
			oxFull.wInit(iDim,:) = oxTngl.PairClassifier(iPair).wInit(2,:);
		else
			oxFull.wInit(iDim,:) = oxFull.wInit(iDim,:) / 2^16;
		end
	end
	if ~isempty(oxFull.Cube)
		%-----------------------------------------------------------------
		% Memory has not been clipped.  We can re-initialize and re-optimize.
		%-----------------------------------------------------------------
		oxFull.SolverOptions.wDiffMax = 4;
		Mcl_Exemplar_Train(oxFull,false);
		%-----------------------------------------------------------------
	end
	%	Compare again
	if oxFull.h > oxTngl.PairClassifier(iPair).h
		%-----------------------------------------------------------------
		% The bivariate case is still better.  Now we know it is the fault of the Cube, not the optimization procedure.  Most
		%	likely, this means that the dimensions not included in the bivariate model are totally irrelevant and just add noise.
		%	Therefore, simply evaluate the univariate version of the bivariate classifier once and skip optimization (which relies on the Cube).
		%-----------------------------------------------------------------
		oxFull.wOptimized = oxFull.wInit;
		oxFull.Dv = oxTngl.PairClassifier(iPair).Dv;
		oxFull.P = oxTngl.PairClassifier(iPair).P;
		oxFull.H = oxTngl.PairClassifier(iPair).H;
		oxFull.h = oxTngl.PairClassifier(iPair).h;
		oxFull.Hcross = oxTngl.PairClassifier(iPair).Hcross;
		oxFull.Acc = oxTngl.PairClassifier(iPair).Acc;
		oxFull.ConfusionMatrix = oxTngl.PairClassifier(iPair).ConfusionMatrix;
		oxFull.Quant = oxTngl.PairClassifier(iPair).Quant;
		%-----------------------------------------------------------------
	end
	%--------------------------------------------------------------------------------------------------
end
	