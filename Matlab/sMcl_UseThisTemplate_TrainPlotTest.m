%=================================================================================
% WHAT DOES THIS SCRIPT FILE DO?
%	Run this script file to train a MoRPE classifier, plot your data (with classifier-relevant information), and test the classifier's performance.
%---------------------------------------------------------------------------------

%	If the "Mcl" directory is on the same path as your working directory, you can add the Mcl directory to the working path by
%	uncommenting the following line of code.  This allows Mcl functionality from your working directory.
%addpath([cd, '\..\Mcl\']);

%	A figure will be saved in the following directory (one above your working directory).
PATH_OUT = [cd, '\..\'];

%=================================================================================
% IN THIS SECTION, READ IN YOUR DATA.
%	Your job is to create two cell arrays that hold your (1) training and (2) testing data.
%	You also need to create an "Axes" variable that determines how each feature dimension is plotted.
%		1.  Xtrain
%		2.  Xtest
%	Next you want to create an "Mcl_Axis" structure for each feature dimension.
%		3.  Axes
%	
% IN PLACE OF THE CODE THAT YOU WILL INSERT, the following code generates random data for training and testing.
%	The following simulated data is 3 categories which is determined by the length of the cell arrays Xtrain and Xtest.
%	If your data has a different number of categories, your variables Xtrain and Xtest will reflect that, and the subsequent
%	code will execute fine.  (Just pay attention to the comments near COLORS, a variable below.)
%---------------------------------------------------------------------------------

%	Parameters of your dataset.
N_FEATURE_DIMENSIONS = 5;
N_TRAIN = 1000;
N_TEST = 1000;

%	Distribution means.
means = {
	rand(1,N_FEATURE_DIMENSIONS), ...
	rand(1,N_FEATURE_DIMENSIONS), ...
	rand(1,N_FEATURE_DIMENSIONS) };
for i=1:3
	means{i} = (means{i}-mean(means{i}))/std(means{i});
end

%	Simulate Gaussian data with specified means.
Xtrain = {...
	randn(N_TRAIN, N_FEATURE_DIMENSIONS) + repmat(means{1},[N_TRAIN,1]), ...
	randn(N_TRAIN, N_FEATURE_DIMENSIONS) + repmat(means{2},[N_TRAIN,1]), ...
	randn(N_TRAIN, N_FEATURE_DIMENSIONS) + repmat(means{3},[N_TRAIN,1]) };
Xtest = {...
	randn(N_TEST, N_FEATURE_DIMENSIONS) + repmat(means{1},[N_TEST,1]), ...
	randn(N_TEST, N_FEATURE_DIMENSIONS) + repmat(means{2},[N_TEST,1]), ...
	randn(N_TEST, N_FEATURE_DIMENSIONS) + repmat(means{3},[N_TEST,1]) };

% Specify the feature dimensions to be included in the plot and the classifier.
%	(Sometimes your data has more dimensions than you want to consider here, so this variable indexes only those dimensions you want.)
DIMS = 1:N_FEATURE_DIMENSIONS;

%	Clean up the samples (remove NaN and Inf)
%Xtrain = Mcl_CleanSamples(Xtrain, true);
%Xtest = Mcl_CleanSamples(Xtest, true);

%=================================================================================
% SET A FEW PLOT OPTIONS
%---------------------------------------------------------------------------------

%	When you plot figures, each category will be associated with a color.  For a 3-category problem, these are the colors.
%	If the number of categories is more than 3, then add more colors here.
COLORS = {[1 0 0], [0,0,1], [0 1 0]};

% Compute the axes (for use by TnglPlot and the classifier itself).
%	For display purposes (TnglPlot), this will set the range of each features and a nonlinear transform to make the data scatter in a visually pleasing way.
%	As it turns out, the same "visually pleasing" style will yield good classifier performance.
%	Your goal should be to make the data scatter in a way that appears Gaussian.
Axes = Mcl_Axis_Ctor([-3,3], -3:3, Mcl_Trfm_Ctor('None'), [-2,0,2]);
Axes = Axes(ones(numel(DIMS),1));

%=================================================================================
% SET CLASSIFIER OPTIONS
%---------------------------------------------------------------------------------

RANK = 2;  %	1=linear classifier, 2=quadratic,  3=cubic
NQUANT=25;
FORCE_EQUAL_PRIORS = true;
TNGL_DISPMOD = 500;
FULL_DISPMOD = 2500;
TNGL_NTRY =	1; %1+RANK;
FULL_NTRY = 2; %2+RANK;
nPoly = 3;
nCats = numel(Xtrain);

%=================================================================================
% THE FOLLOWING CODE IS WHERE THE CLASSIFIER IS ACTUALLY FIT
%	Actually, the "full classifier" is the one you're interested.  I also fit bivariate classifiers for all possible possible 2-way combinations of
%	features.
%	
%	Note, I start by creating the options structure for the "TnglPlot" which is part of the data visualization (the visual figure is created later).
%	I use the fields of this "TnglPlot" options structure to store classifier output (from training).  After training, these same parameters are
%	required to implement the classifier on new data.
%---------------------------------------------------------------------------------

%	Get a new options structure for the TNGL plot.  This structure will also hold the result of training the classifier.
TnglPlot = Mcl_TnglPlot('OPTIONS');

disp(['Constructing and optimizing traingle of ' num2str(numel(DIMS)) ' univariate and  ' num2str(numel(DIMS)*(numel(DIMS)-1)/2) ' bivariate polynomial(' num2str(RANK) ') classifiers...']);
%	Construct triangle of MCL solvers (one solver per pair of spatial dimensions).
tTngl = tic;	%	 Measure how long this will take.
TnglPlot.TnglClassifiers = Mcl_Poly_Tngl_Ctor(Xtrain,Axes,RANK,FORCE_EQUAL_PRIORS,NQUANT,TNGL_DISPMOD,TNGL_NTRY);
TnglPlot.TnglClassifiers.IdDims = DIMS;  % For plotting
disp(['Total time spent in pairwise constructor and optimiztion:  ' num2str(toc(tTngl)) ' seconds']);

disp(' ');
disp('--------------------------------------------');
tFull = tic;  tConstruct = tic;   % Measure how long this will take.
disp(['Constructing full ' num2str(numel(DIMS)) ' dimensional polynomial(' num2str(RANK) ') classifier and optimizing ' num2str(nPoly*Mcl_Poly_Coeff(int32(RANK),int32(numel(DIMS)))) ' parameters...']);
%	Construct an MCL solver.  The structure output is conceptually like an instance of a class.
TnglPlot.FullClassifier = Mcl_Poly_Ctor(Xtrain,Axes,RANK,FORCE_EQUAL_PRIORS,NQUANT,FULL_DISPMOD,FULL_NTRY);
TnglPlot.FullClassifier = Mcl_Poly_CheckWithTngl(TnglPlot.FullClassifier, TnglPlot.TnglClassifiers);
TnglPlot.FullClassifier.IdDims = DIMS;
disp(['  ' num2str(toc(tConstruct)) ' seconds were spent in constructor and solver.']);
disp(' ');
disp('--------------------------------------------------------------------------------------');
%	Reduce fields of the TnglPlot structure to a size that can be stored on disk without unnecessary space being taken up.
TnglPlot.TnglClassifiers = Mcl_Poly_Compress(TnglPlot.TnglClassifiers);		%	The data is de-expanded so that you don't store full polynomial expansions (which can be easily re-created).
TnglPlot.FullClassifier  = Mcl_Poly_Compress(TnglPlot.FullClassifier);		%	The data is de-expanded so that you don't store full polynomial expansions (which can be easily re-created).

%=================================================================================
% DRAW FIGURE (using training data).  All reports of "accuracy" or "entropy" are training values.
%---------------------------------------------------------------------------------

%	Configure some figure properties.
TnglPlot.FigureSizeXy = [900 900];
TnglPlot.Dpi = 300;
TnglPlot.FontName_Acc = 'Cambria';
TnglPlot.FontName_Label = 'Cambria';
TnglPlot.Axes = Axes;
TnglPlot.BlurLevel = 0.12;  % For display purposes, can be a number in the range [0, 1), suggest 0.12.
TnglPlot.Colors = {[1 0 0], [0 0.8 0], [0 0 1]}; % The color of each category
TnglPlot.DrawH = false;
TnglPlot.FontName_Acc = 'Cambria';
TnglPlot.FontName_Label = 'Cambria';
TnglPlot.FontSize_Axis = ceil(150/double(numel(DIMS)+1)*min(TnglPlot.FigureSizeXy)/1000);
TnglPlot.FontSize_Acc = ceil(165/double(numel(DIMS)+1)*min(TnglPlot.FigureSizeXy)/1000);
TnglPlot.FontSize_Label = ceil(175/double(numel(DIMS)+1)*min(TnglPlot.FigureSizeXy)/1000); % Latex draws larger than normal font size.
TnglPlot.FullLabel = '$\rho(c|\bf{\rm{x}})$';
TnglPlot.HistSize_1d = 20;
TnglPlot.HistSize_2d = [200 200];
TnglPlot.LineWidth = 2;
TnglPlot.FullLabel_TextInterpreter = 'latex';

disp('Generating TnglPlot...');
figure;
TnglPlot = Mcl_TnglPlot(Xtrain, TnglPlot);

%	Save TnglPlot and related information.
drawnow;
if ~exist(PATH_OUT, 'dir')
	mkdir(PATH_OUT);
end

fname = ['TnglPlot_' num2str(numel(DIMS)) '_' num2str(RANK)];

disp('Saving to cache...');
save([fname '.mat'], 'TnglPlot');

disp('Printing...');
print(gcf, '-dtiff', ['-r' num2str(TnglPlot.Dpi)], '-painters', [fname '.tiff']);

close gcf;
pause(1);

%=================================================================================
% CLASSIFY NEW DATA (i.e. Xtest) AND MEASURE TESTING PERFORMANCE (i.e. testing accuracy, entropy).
%	Only the full classifier is tested here.  We don't test the pairwise classifiers because they are just intended
%	as a check on the convergence of the full classifier, plus they enhance the visual figure.
%---------------------------------------------------------------------------------

%	The number of testing samples in each category
nTest = zeros(nCats,1);
for iCat=1:nCats
	nTest(iCat) = size(Xtest{iCat},1);
end

%	Initialize memory (copy the full classifier)
TnglPlot.FullClassifier_Test = TnglPlot.FullClassifier;
TnglPlot.FullClassifier_Test.Ntest = nTest;
TnglPlot.FullClassifier_Test.hTest = nan(nCats,1);
TnglPlot.FullClassifier_Test.AccTest = nan(nCats,1);
TnglPlot.FullClassifier_Test.CfnMatrixTest = nan(nCats,nCats);
TnglPlot.FullClassifier_Test.HcrossTest = nan(nCats,nCats);

disp(' ');
disp('--------------------------------------------');
tFull = tic;  tConstruct = tic;   % Measure how long this will take.
disp(['Testing full ' num2str(numel(DIMS)) ' dimensional polynomial(' num2str(RANK) ') classifier and optimizing ' num2str(nPoly*Mcl_Poly_Coeff(int32(RANK),int32(numel(DIMS)))) ' parameters...']);

%	Test for each category.
for iCat=1:nCats
	if nTest(iCat)>0
		%	Allocate memory
		Dv = zeros(nTest(iCat),nCats);
		P  = zeros(nTest(iCat),nCats);
		H  = zeros(nTest(iCat),1);
		%	Compute model response to test set.
		Mcl_Poly_CalcDv(TnglPlot.FullClassifier_Test, Dv, Mcl_Poly_Expand(Mcl_Trfm_Cols(Axes,true,Xtest{iCat}), TnglPlot.FullClassifier_Test.CoeffDims));
		Mcl_MapDv(P, Dv, TnglPlot.FullClassifier_Test.Quant);
		%	Compute accuracy and confusion matrix.
		[dum, resp] = max(P, [], 2);
		TnglPlot.FullClassifier_Test.AccTest(iCat) = mean(resp==iCat);
		%	Compute confusion matrix
		for jCat=1:nCats
			if FORCE_EQUAL_PRIORS
				TnglPlot.FullClassifier_Test.CfnMatrixTest(iCat,jCat) = mean(resp==jCat) / sum(nTest>0);
			else
				TnglPlot.FullClassifier_Test.CfnMatrixTest(iCat,jCat) = sum( resp==jCat) / sum(nTest);
			end
		end
		%	Measure the conditional entropy
		[h, hc] = Mcl_ConditionalEntropy(H, P(1:nTest(iCat),:), repmat(int32(iCat-1), [nTest(iCat),1]), false);
		TnglPlot.FullClassifier_Test.hTest(iCat) = h;
		if FORCE_EQUAL_PRIORS
			TnglPlot.FullClassifier_Test.HcrossTest(iCat,:) = hc(iCat,:) / sum(nTest>0);
		else
			TnglPlot.FullClassifier_Test.HcrossTest(iCat,:) = hc(iCat,:) * NTEST/sum(nTest);
		end
	end
end
disp(['Testing  [h | Acc] = [' num2str(TnglPlot.FullClassifier_Test.hTest') '   |   ' num2str(TnglPlot.FullClassifier_Test.AccTest') ']']);

%	The data is de-expanded so that you don't store full polynomial expansions (which can be easily re-created).
TnglPlot.FullClassifier_Test  = Mcl_Poly_Compress(TnglPlot.FullClassifier_Test);
% --------------------------------------------------------------------------------

%=================================================================================
% SAVE Everything (so you don't have to train and test again.  Quickly create a new plot by loading this structure and
%	calling the funciton "TnglPlot".
% --------------------------------------------------------------------------------
save TnglPlot;