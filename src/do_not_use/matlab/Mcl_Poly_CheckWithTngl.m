function oxFull = Mcl_Poly_CheckWithTngl(oxFull, oxTngl)

% function oxFull = Mcl_Poly_CheckWithTngl(oxFull, oxTngl)
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
% oxFull (matlab structure Mcl_Poly scalar)
%	An exemplar classifier constructed on more than 2 sample dimensions.
%--------------------------------------------------------------------------------------------------------------------------
% INPUT
% oxTngl (matlab structure Mcl_Poly_Tngl scalar)
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
	iDim=oxTngl.DimPairs(iPair,1)-1;
	jDim=oxTngl.DimPairs(iPair,2)-1;
	pnc = oxTngl.PairClassifier(iPair).Ncoeff;
	%	Express the dimensions of each coefficient in terms of the full model
	pdims = oxTngl.PairClassifier(iPair).CoeffDims;
	pdims(pdims(:)==1)=jDim;
	pdims(pdims(:)==0)=iDim;
	%	A map from the bivariate coefficients to the coefficients of the full model.
	cMap = zeros(pnc,1,'int32');
	for iMap=1:pnc
		for iCoeff=1:oxFull.Ncoeff
			if all(pdims(iMap,:)==oxFull.CoeffDims(iCoeff,:))
				cMap(iMap)=iCoeff;
				break;
			end
		end
	end
	if any(cMap==0)
		error('Assertion failed.  This code needs some work.');
	end
	%	Initialize the coefficients of the full model to match the optimized coefficients of the pairwise model.
	oxFull.wInit(:) = 0;
	for iPoly=1:oxFull.Npoly
		oxFull.wInit(cMap,iPoly) = oxTngl.PairClassifier(iPair).wOptimized(:,iPoly);
	end
	%	Retrain full classifier.
	if oxFull.DisplayModulus>0
		disp(['Full classifier is less powerful than bivariate child (' num2str(iDim+1) ',' num2str(jDim+1) ').  Re-initialized as child.  Retraining...']);
	end
	oxFull = Mcl_Poly_Train(oxFull);
	%	Display result
	if oxFull.DisplayModulus>0
		disp(['Mcl_Poly_Tngl_Ctor: Full classifier re-trained.  [h, Acc] = [' num2str(oxFull.h) ',' num2str(oxFull.Acc) ']']);
		disp('--------------------------------------------------------------------------------');
	end
	%--------------------------------------------------------------------------------------------------
end
	