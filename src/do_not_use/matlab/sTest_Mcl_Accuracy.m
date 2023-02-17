%	Test the triangular exemplar models, the full exemplar model, and the plotting function Mcl_TriPlot

%-------------------------------------------------------------------------------------------------
%	Set simulation parameters and simulate training samples (i.e. make up N-category multivariate data).
%-------------------------------------------------------------------------------------------------
% The number of spatial dimensions
NDIMS = 3;
% The number of samples from each category

NSAMP = round(3000*[1 1]);  % Simple 2-cat, equal priors
%NSAMP = round(2000*[1, .25]); % Simple 2-cat, unequal priors, but equal priors will be enforced.

%NSAMP = round(1000*[1 1 1]); % 3-cat, equal priors.
%NSAMP = round(1000*[1, 0.5, 0.25]); % 3-cat, unequal priors, but equal priors will be enforced.

%NSAMP = round(1000*[1 1 1 1]); % 4-cat, equal priors.
%NSAMP = round(1000*[.1, 1, .25, .5]); % 4-cat, unequal priors, but equal priors will be enforced.

% The number of standard deviations separating the Gaussian means.
DPRIME = 5;
%	Force equal priors
FORCE_EQUAL_PRIORS = true;

% The number of categories
nCats = numel(NSAMP);
% The total number of samples
ntSamp = sum(NSAMP);

%	Generate the means of each population
if( NDIMS==1 )
	if( nCats==2 )
		mu = [-1; 1] * DPRIME / 2;
	else
		mu = ((0.5:nCats)'-(nCats/2)) * DPRIME / 2;
	end
else
	mu = [];
	while size(mu,1)<nCats
		mu = [mu; Mcl_RandomRotationMatrix(NDIMS)];
	end
	mu = mu(1:nCats,:) * DPRIME / sqrt(2);
end
%	Generate random samples from each population
Xcell = cell(nCats,1);
X = zeros(ntSamp,NDIMS);
Cat = zeros(ntSamp,1,'int32');
iSamp = 0;
for iCat=1:nCats
	Xcell{iCat} = randn(NSAMP(iCat),NDIMS) + repmat(mu(iCat,:), [NSAMP(iCat),1]);
	ind = iSamp+(1:NSAMP(iCat));
	X(ind,:) = Xcell{iCat};
	Cat(ind) = int32(iCat-1);
	iSamp = iSamp + NSAMP(iCat);
end
%-------------------------------------------------------------------------------------------------

%	Check d'
for iCat=1:(nCats-1)
	for jCat=(iCat+1):nCats
		dist = sqrt(sum( (mean(Xcell{iCat})-mean(Xcell{jCat})).^2 ));
		disp(['Check:  dPrime [' num2str(iCat) ',' num2str(jCat) '] = ' num2str(dist)]);
	end
end

disp('=========================================');
disp([num2str(normcdf(DPRIME/2)) ' is the theoretical accuracy for a dPrime of ' num2str(DPRIME)]);
%	Construct theoretical accuracy for each dimension
SgnAccTheory = zeros(nCats,NDIMS);
for iDim=1:NDIMS
	for iCat=1:nCats
		mui = mu(iCat,iDim);
		muj = mean(mu( (1:nCats)~=iCat, iDim ));
		SgnAccTheory(iCat,iDim) = normcdf(abs(mui-muj)/2);
		if mui<muj
			SgnAccTheory(iCat,iDim) = -SgnAccTheory(iCat,iDim);
		end
	end
end
disp(' ');
disp('The signed accuracy gives the accuracy and the direction of the univariate decision rule.');
if nCats>2
	disp(' ');
	disp('For more than two categories, each univariate criterion is obtained by simulating a 2-category classifier with equal priors (so chance is always 50%).');
end
disp(' ');
disp('In theory, we would expect the signed accuracy to be approximately:');
disp(SgnAccTheory);
%	Construct CatWeightsEach
if FORCE_EQUAL_PRIORS
	CatWeight = mean(NSAMP)./NSAMP;
else
	CatWeight = ones(size(NSAMP));
end
[SgnAcc, Crit] = Mcl_Accuracy(X, Cat, int32(nCats), CatWeight);
disp(' ');
disp('Signed Accuracy (nCats * NDIMS) = ');
disp(num2str(SgnAcc));
disp(' ');
disp('Criterion (nCats * NDIMS) = ');
disp(num2str(Crit));

figure; plot([-1 1], [-1 1], 'k');
hold on;
xlabel('Theory');
ylabel('Measured');
plot(SgnAccTheory(:),SgnAcc(:), 'bo');

ops = Mcl_TnglPlot('OPTIONS');
ops.FigureSizeXy = [800 800];
ops.HistSize_2d = [200 200];
ops.AxisMin = -5;
ops.AxisMax = +5;
ops.AxisMid = 0;
ops.Gamma = 1+log10(nCats);
figure;
tngl = Mcl_TnglPlot(Xcell, ops);
for iDim=1:NDIMS
	subplot(tngl.SubplotHandles(iDim,iDim));
	dy = 18;
	yMax = 1.15;
	ylPos = yMax;
	yrPos = yMax;
	pixHeight = ops.FigureSizeXy(2) * tngl.SubplotRects(iDim,iDim,4);
	for iCat=1:nCats
		plot(Crit(iCat,iDim)+[0 0], [0 1], 'Color', tngl.Colors{iCat});
		accStr = [num2str(round(abs(SgnAcc(iCat,iDim)*100))) '%'];
		if SgnAcc(iCat,iDim)<0
			text(ops.AxisMin+0.25, ylPos, ['-' accStr], 'Color', tngl.Colors{iCat}, 'FontWeight', 'bold');
			ylPos = ylPos - dy/pixHeight;
		else
			text(ops.AxisMax-0.05, yrPos, ['+' accStr], 'Color', tngl.Colors{iCat}, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
			yrPos = yrPos - dy/pixHeight;
		end
	end
	for iCat=nCats:-1:1
		plot(Crit(iCat,iDim)+[0 0], [0.5 1], 'Color', tngl.Colors{iCat});
	end
end
disp('=========================================');