%	Test the polynomial solver (Mcl_Poly) and the triangle plotting function (Mcl_TnglPlot).

%-------------------------------------------------------------------------------------------------
%	Set parameters of simulated training samples (i.e. make up N-category multivariate data).
%-----------------------------------------
% The number of spatial dimensions
NDIMS = 6;
% The number of samples from each category
%NSAMP = round(500*[1 1]);
%NSAMP = round(2000*[1 1]);
%NSAMP = round(1000*[1 0.25]);
%NSAMP = round(1000*[1 .5 .25]);
NSAMP = round(2000*[1 1 1]);
% The number of standard deviations separating the Gaussian means.
DPRIME = 3;
%-------------------------------------------------------------------------------------------------
%	Set the options for the polynomial solver.
%-----------------------------------------
% The polynomial rank
RANK = 3;
% The number of quantiles used to quantize the decision values.
NQUANT = min(30,floor(min(NSAMP)/20));
%	Force equal priors (this should be true for most practical situations that I can imagine.).
FORCE_EQUAL_PRIORS = true;
%	Display modulus (controls frequency of status messages to command window).
TNGL_DISPLAY_MOD = 100;
%	Number of optimization approaches
TNGL_NOPT = 1;
%	Display modulus (controls frequency of status messages to command window).
FULL_DISPLAY_MOD = 500;
%	Number of optimization approaches
FULL_NOPT = 1;
%-------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------
%	Create stimuli
%-----------------------------------------
% The number of categories
nCats = numel(NSAMP);
%	Generate the means of each population
mu = [];
while size(mu,2)<nCats
	mu = [mu, Mcl_RandRotate(int32(NDIMS))];
end
mu = mu(:,1:nCats) * DPRIME / sqrt(2);
%	Generate random samples from each population
Xcell = cell(nCats,1);
for iCat=1:nCats
	Xcell{iCat} = randn(NSAMP(iCat),NDIMS) + repmat(mu(:,iCat)', [NSAMP(iCat),1]);
end
%-------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------
%	Create stimulus dimension transforms and then transform stimulus dimensions
%-----------------------------------------
if true
	Xtrans = [Mcl_Trfm_Ctor('FisherZ',1.0025); Mcl_Trfm_Ctor('LogAbs',0.0001)];
	for iDim=3:NDIMS
		Xtrans(iDim) = Mcl_Trfm_Ctor('None');
	end
	Xcell = Mcl_Trfm_Cols(Xtrans,false,Xcell);
else
	Xtrans = [];
end
%-------------------------------------------------------------------------------------------------

for iCat=1:(nCats-1)
	for jCat=(iCat+1):nCats
		dist = sqrt(sum( (mean(Xcell{iCat})-mean(Xcell{jCat})).^2 ));
		disp(['Check:  dPrime [' num2str(iCat) ',' num2str(jCat) '] = ' num2str(dist)]);
	end
end

disp('=========================================');


disp('--------------------------------------------');
tTngl = tic;
disp(['Constructing and optimizing the traingle of ' num2str(NDIMS) ' univariate and ' num2str(NDIMS*(NDIMS-1)/2) ' bivariate classifiers...']);
%	Construct triangle of MCL solvers (one solver per pair of spatial dimensions).
sOps.Napproaches = 2;
oxTngl = Mcl_Poly_Tngl_Ctor(Xcell,Xtrans,RANK,FORCE_EQUAL_PRIORS,NQUANT,TNGL_DISPLAY_MOD,TNGL_NOPT);
disp(['Total time spent in Tngl (triangle) constructor and optimization:  ' num2str(toc(tTngl)) ' seconds']);
disp(' ');

disp('--------------------------------------------');
tFull = tic;
nParams = Mcl_Poly_Coeff(int32(RANK),int32(NDIMS));
if nCats>2
	nParams = nParams * nCats;
end
disp(['Constructing the full ' num2str(NDIMS) ' dimensional classifier and optimizing its ' num2str(nParams) ' parameters...']);
tConstruct = tic;
%	Construct and optimize an Mcl_Poly classifier.  The structure output is conceptually like an instance of a class.
ox = Mcl_Poly_Ctor(Xcell,Xtrans,RANK,FORCE_EQUAL_PRIORS,NQUANT,FULL_DISPLAY_MOD,FULL_NOPT);
disp(['  ' num2str(toc(tConstruct)) ' seconds were spent constructing and optimizing the full model... Now checking Tngl...']);
%	Check the full optimization with the Tngl to ensure sub-dimensional spaces achieved higher entropy than the full space.
ox = Mcl_Poly_CheckWithTngl(ox, oxTngl);
%	Display status.
disp('    Tngl check finished.');
disp(' ');


disp('--------------------------------------------');
tPlot = tic;
disp('Generating trainglular scatterplot using TnglPlot...');
%	Get the plot options
plotOps = Mcl_TnglPlot('OPTIONS');
%	Prepare to pass the classifier into the plotting function.
plotOps.FullClassifier = ox;
plotOps.TnglClassifiers = oxTngl;
%	Set the axis limits for the plotting function.
rng = ceil(DPRIME/sqrt(NDIMS));
if isempty(Xtrans)
	rng = rng+zeros(NDIMS,1);
	plotOps.Axes = Mcl_Axis_Ctor([-2-rng, +2+rng], [-2-rng, zeros(size(rng)), 2+rng]);
else
	plotOps.Axes = [];
	for iDim=1:NDIMS
		if isequal('None', Xtrans(iDim).Name)
			plotOps.Axes = [plotOps.Axes, Mcl_Axis_Ctor([-2-rng, +2+rng], [-2-rng, zeros(size(rng)), 2+rng], Xtrans(iDim))];
		elseif isequal('LogAbs', Xtrans(iDim).Name)
			plotOps.Axes = [plotOps.Axes, Mcl_Axis_Ctor([0.0001, 10000], [0.001,0.01,0.1,1,10,100,1000], Xtrans(iDim))];
		elseif isequal('FisherZ', Xtrans(iDim).Name)
			plotOps.Axes = [plotOps.Axes, Mcl_Axis_Ctor([-1,1], -1:0.5:1, Xtrans(iDim))];
		end
		plotOps.Axes(iDim).TextLabel = ['${\it{x}\rm{_{' num2str(iDim) '}}}$'];
	end
end
plotOps.FullLabel_TextInterpreter = 'latex';
plotOps.HistSize_2d = [200 200];
plotOps.FigureSizeXy = [800 800];
%plotOps.Gamma = 1;
plotOps.BlurLevel = 0.1;
%	Open a figure
figure; hold on;
%	Plot the data
plotOutput = Mcl_TnglPlot(Xcell, plotOps);
print(gcf, '-dtiff', ['-r' num2str(plotOutput.Dpi)], '-painters', ['c:/users/ing/desktop/sMcl_Poly_Test[Ncats=' num2str(nCats) ',Rank=' num2str(RANK) ',Ndims=' num2str(NDIMS) '].tiff']);
%imwrite(frame2im(getframe(gcf)), ['c:/users/ing/desktop/sMcl_Poly_Test[Ncats=' num2str(nCats) ',Rank=' num2str(RANK) ',Ndims=' num2str(NDIMS) '].jpg'], 'jpg', 'Quality', 100);
disp('=========================================');