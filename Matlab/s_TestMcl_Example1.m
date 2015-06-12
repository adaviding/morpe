%	This is a script for testing / demonstrating
%		Mcl_FitExemplar (and associated functions)
%		Mcl_FitQuad (and associated functions)
%		Mcl_TriangularScatterplot
%	In this script, multivariate data is made up (spherical Gaussian 2-category data),
%		fit with quadratic and exemplar models, and then displayed.

DO_EXEMPLAR = false;
DO_QUADRATIC = DO_EXEMPLAR;
WRITE_PERFORMANCE = DO_EXEMPLAR;  % If true, Accuracy (top-right), Rsq (top-left), and Lamda (bottom-right) are drawn in corners of each plot

PLOT_SCALE = 2;

rand('twister',sum(100*clock));

%-------------------------------------------------------------------------------------------------
%	Generate test samples (i.e. 2-category multivariate data is made up).
%-------------------------------------------------------------------------------------------------
if true
	dPrimeApprox = 2.5;		% Any positive number (the greater this number, the more discriminable the categories are)
	nDims = 5;		% The dimensionality of the space.
	nSamples = 200;	% The number of samples from each category

	%	Compute the mean for categories A and B (two category problem)
	oVec = (floor(4.99999999*rand(1,nDims))-2);
	if all(oVec==0)
		oVec(1) = 1;
	end
	muB = dPrimeApprox./sqrt(sum(abs(oVec)))*oVec;
	muA = -muB/2;
	muB = muB/2;

	%	Generate random samples
	x = {...
		randn(nSamples,nDims)+repmat(muA, [nSamples 1]),...
		randn(nSamples,nDims)+repmat(muB, [nSamples 1]) };
end
%-------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------
%	Fit a quadratic categorization model and pairwise sub-models
%-------------------------------------------------------------------------------------------------
if DO_QUADRATIC
	qops = Mcl_FitQuad('OPTIONS');	%	Get default options

	qops.Criterion = 'l';		% 'l' = lamda maximization (i.e. maximum likelihood, minimum entropy), 'a' = accuracy maximization, 'v' = maximize R^2
	qops.ApproachStyle = 2;		%	ApproachStyle 2 is more reliable (but takes longer) than ApproachStyle 1
	qops.Napproaches = 4;		%	The number of unique approaches (each taken with a unique order of parameters).

	disp('Fitting main quadratic model...  (This takes a few minutes.)');

	qfit = Mcl_FitQuad(x,qops)		%	The full quadratic categorization models

	disp('Fitting pairwise quadratic models...  (This takes a few minutes.)');

	qtfits = Mcl_TriangularFitQuad(x,qops);		%	The pairwise sub-models
end
%-------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------
%	Fit an exemplar model in 3 different ways, compare the fits
%-------------------------------------------------------------------------------------------------
if DO_EXEMPLAR
	xops = Mcl_FitExemplar('OPTIONS');

	xops.Criterion = 'l';  % 'l' = lamda maximization (i.e. maximum likelihood, minimum entropy), 'a' = accuracy maximization, 'v' = maximize R^2
	xops.ApproachStyle = 2;
	xops.Generalization = 'e'; % 'h', 'e', or 'g'

	disp('Fitting main exemplar models with exponential generalization...  (This takes a few minutes.)');

	xefit = Mcl_FitExemplar_Fast(x,xops)
	
	disp('Fitting pairwise exemplar model with exponential generalization...  (This takes a few minutes.)');

	xetfits = Mcl_TriangularFitExemplar_Fast(x,xops);
end
%-------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------
%	Create a figure of triangular scatterplots
%--------------------------------------------------------------------------------------------------
pops = Mcl_TriangularScatterplot('OPTIONS');
pops.AxisMin = -4;
pops.AxisMax = +6;
pops.AxisMid = (pops.AxisMax+pops.AxisMin)/2;
pops.BlurLevel = 0.25;  % For display purposes, can be a number in the range [0, 1), suggest 0.2 and no greater than 0.5.
pops.Colors = {[0 0 1], [1 0 0]}; % The color of each category
pops.DrawAccuracy = WRITE_PERFORMANCE;
pops.DrawRsq = false;
pops.DrawLamda = WRITE_PERFORMANCE;
pops.FillHistograms = false;
% pops.HistSize_1d = 25;
pops.HistSize_2d = [100 100];
pops.Labels = cell(nDims,1); for i=1:nDims, pops.Labels{i} = ['$\it{s}\rm{_' num2str(i) '}$']; end;
pops.LineWidth = 1*PLOT_SCALE;
pops.TextInterpreter = 'latex';
pops.FontSize_Label = 18*PLOT_SCALE;
pops.FontSize_Axis = 10*PLOT_SCALE;
pops.Size_Plot = [80 80]*(1+nDims)*PLOT_SCALE;
if DO_EXEMPLAR
	pops.mFullExem = xefit;
	pops.mPairsExem = xetfits;
end
if DO_QUADRATIC
	pops.mFullQuad = qfit;
	pops.mPairsQuad = qtfits;
end
disp('Generating triangular scatterplot... (This takes less than a minute.)');

figure;
xeResult = Mcl_TriangularScatterplot(x,pops);
%-------------------------------------------------------------------------------------------------

%	You could save the image of the plotting area like this:
if DO_QUADRATIC
	imwrite(frame2im(getframe(gcf)), 'C:\Users\ing\Desktop\Mcl_Fig2.jpg', 'jpg', 'Quality', 100);
else
	imwrite(frame2im(getframe(gcf)), 'C:\Users\ing\Desktop\Mcl_Fig1.jpg', 'jpg', 'Quality', 100);
end