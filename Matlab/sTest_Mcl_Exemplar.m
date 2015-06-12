%	Test the triangular exemplar models, the full exemplar model, and the plotting function Mcl_TriPlot

%-------------------------------------------------------------------------------------------------
%	Set simulation parameters and simulate training samples (i.e. make up N-category multivariate data).
%-------------------------------------------------------------------------------------------------
% The number of spatial dimensions
NDIMS = 4;
% The number of samples from each category
%NSAMP = round(1000*[1 .5 .25]);
NSAMP = round(500*[1 1 1]);
% The number of neighbors considered in the training function.  Half of these are "nearest neighbors",
%	the other half are random samples from the remaining area.
NNEIGH = max(60,30*NDIMS);
% The number of standard deviations separating the Gaussian means.
DPRIME = 2.5;
% The number of quantiles used to quantize the decision values.
NQUANT = min(25,floor(min(NSAMP)/20));
%	Which methods should we run?  These need to be taken from the set {0,1,2,10,11,12}
IDMETHOD = 1;

% The number of categories
nCats = numel(NSAMP);
%	Generate the means of each population
mu = [];
while size(mu,2)<nCats
	mu = [mu, Mcl_RandomRotationMatrix(NDIMS)];
end
mu = mu(:,1:nCats) * DPRIME / sqrt(2);
%	Generate random samples from each population
Xcell = cell(nCats,1);
for iCat=1:nCats
	Xcell{iCat} = randn(NSAMP(iCat),NDIMS) + repmat(mu(:,iCat)', [NSAMP(iCat),1]);
end
%-------------------------------------------------------------------------------------------------
for iCat=1:(nCats-1)
	for jCat=(iCat+1):nCats
		dist = sqrt(sum( (mean(Xcell{iCat})-mean(Xcell{jCat})).^2 ));
		disp(['Check:  dPrime [' num2str(iCat) ',' num2str(jCat) '] = ' num2str(dist)]);
	end
end

disp('=========================================');
disp(['Method ' num2str(IDMETHOD)]);

% Get the default optimization options.  If options are not gotten, default options will be used later.
sOps = Mcl_MinimizeEntropy('OPTIONS');
% To suppress output from the solver into the command window, set the DisplayModulus to any number less than 1
sOps.DisplayModulus = 1;

disp('--------------------------------------------');
tTngl = tic;
disp('Constructing traingle of NDIMS univariate and NDIMS*(NDIMS-1)/2 bivariate classifiers...');
%	Construct triangle of MCL solvers (one solver per pair of spatial dimensions).
sOps.Napproaches = 2;
oxTngl = Mcl_Exemplar_Tngl_Ctor(Xcell,[],IDMETHOD,NNEIGH,true,NQUANT,true,sOps);
disp(['Total time spent in Tngl (triangle) constructor and optimization:  ' num2str(toc(tTngl)) ' seconds']);
disp(' ');

disp('--------------------------------------------');
tFull = tic;
disp('Constructing full NDIMS dimensional classifier...');
tConstruct = tic;
%	Construct an MCL solver.  The structure output is conceptually like an instance of a class.
ox = Mcl_Exemplar_Ctor(Xcell,[],IDMETHOD,NNEIGH,true,NQUANT);
disp(['  ' num2str(toc(tConstruct)) ' seconds were spent in constructor.']);
% The number of approaches for the optimization procedure.
sOps.Napproaches = 3;
ox.SolverOptions = sOps;
disp('Optimizing full NDIMS dimensional classifier...');
tSolve = tic;
%	Run optimization.
ox = Mcl_Exemplar_Train(ox,false);
ox = Mcl_Exemplar_CheckWithTngl(ox, oxTngl);
%	Clear memory that is needed only for optimization.
ox.Cube = [];
ox.Mem = [];
%	Display status.
disp(['  ' num2str(toc(tSolve)) ' seconds were spent in solver.']);
disp(['Total time spent in constructor and solver:  ' num2str(toc(tFull)) ' seconds']);
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
plotOps.AxisMin = -3-ceil(DPRIME/sqrt(2*NDIMS));
plotOps.AxisMid = 0;
plotOps.AxisMax = +3+ceil(DPRIME/sqrt(2*NDIMS));
plotOps.Labels = cell(NDIMS,1);
plotOps.TextInterpreter = 'latex';
for iDim=1:NDIMS
	plotOps.Labels{iDim} = ['${\it{x}\rm{_{' num2str(iDim) '}}}$'];
end
%plotOps.HistSize_2d = [200 200];
%	Open a figure
figure; hold on;
%	Plot the data
plotOutput = Mcl_TnglPlot(Xcell, plotOps);
disp(['  ' num2str(toc(tSolve)) ' seconds were spent generating the plot.']);
imwrite(frame2im(getframe(gcf)), 'c:/users/ing/desktop/sMcl_Exemplar_Test.jpg', 'jpg', 'Quality', 100);
disp('=========================================');