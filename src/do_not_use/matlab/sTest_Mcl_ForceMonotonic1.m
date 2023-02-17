%	A script for testing Mcl_ForceMonotonic.
%	Run it repeatedly (N_SIMS=1000) to see how the function works.

PLOT_SCALE = 1; %6;
N_SIMS = 1;
IDMETHOD = int32(2);

% Scenarios
SEPARATION = 1;		%	This is equivalent to d-prime when both functions are Gaussian
RATES = [1 1];		%	This is equivalent to "standard deviation" when both functions are Gaussian
EXPONENTS = [2 2];	%	[1 1] --> Exponentials,  [2 2] --> Gaussian, [4 4] --> Very Platykurtic
if EXPONENTS(1)~=EXPONENTS(2)
	warning('WARNING:  When the distribution exponents are not equal, this script is not an accurate simulation of a monotonically corrected quadratic model.  This scripts assumes that the optimal quadratic model is linear (true only when exponents are equal).  This script does not optimize non-linear parameters and therefore the gray, blue, and red lines are not accurate, nor are any simulation results.  Only the black line is accurate because all of the quadratic parameters were optimized to generate it.');
end

%	Ensure that RATES(2) is always higher than RATES(1).  If rates are the same, then ensure that the first exponent is higher than the second.
if RATES(1)>RATES(2)
	%	Swap
	dum = RATES(2);
	RATES(2) = RATES(1);
	RATES(1) = dum;
	
	dum = EXPONENTS(2);
	EXPONENTS(2) = EXPONENTS(1);
	EXPONENTS(1) = dum;
elseif RATES(1)==RATES(2) && EXPONENTS(2)>EXPONENTS(1)
	%	Swap
	dum = RATES(2);
	RATES(2) = RATES(1);
	RATES(1) = dum;
	
	dum = EXPONENTS(2);
	EXPONENTS(2) = EXPONENTS(1);
	EXPONENTS(1) = dum;
end

nBins = 100;
nStim = [1000, 1000];

%	Compute the quantile size
quantSize = ceil(sum(nStim)/nBins);
if quantSize*nBins ~= sum(nStim)
	error('The number of bins must be an integer factor of the total number of samples.');
end
%	Compute quantile boundaries.
iqLow = 1:quantSize:sum(nStim);
iqHigh = quantSize:quantSize:sum(nStim);

%	Compute the theoretical (actual) p-value for each quantile.
dx=0.01;
x = (-5^sqrt(1/max(EXPONENTS))-SEPARATION/2):dx:(5^sqrt(1/max(EXPONENTS))+SEPARATION/2);
uTerp = (0.5:nBins)/nBins;

%	Density (for each category).
dCat = {...
		exp(-(abs(x-SEPARATION/2)/RATES(1)).^EXPONENTS(1)/EXPONENTS(1))	, ...
		exp(-(abs(x+SEPARATION/2)/RATES(2)).^EXPONENTS(2)/EXPONENTS(2))	};
%	Cumulative probability (for each category).
pCat = {[], []};
for iCat=1:2
	pCat{iCat} = cumtrapz(dCat{iCat});
	dCat{iCat} = dCat{iCat} / pCat{iCat}(length(pCat{iCat})) * dx;
	pCat{iCat} = pCat{iCat} / pCat{iCat}(length(pCat{iCat}));
end
pTotal = 0.5*(pCat{1}+pCat{2});
p2 = dCat{1}./(dCat{1}+dCat{2});
pActual = interp1(pTotal,p2,uTerp);

%	The horizontal axis of the plots
binNumber = 1:length(pActual);

%	Init memory for monotonic regression
pMonotonic = zeros(size(iqLow));
pNond = zeros(size(iqLow));
pMono = zeros(size(iqLow));

%------------------------------------------------------------------------------------------
% Run simulations
%------------------------------------------------------------------------------------------
N_SIMS = max(1,N_SIMS);
MseReduction = zeros(N_SIMS,1);
for iSim=1:N_SIMS
	%	Simulate stimuli from each category
	stim = {...
		[zeros(nStim(1),1), -interp1(pCat{1},x,rand(nStim(1),1))], ...
		[ ones(nStim(2),1), -interp1(pCat{2},x,rand(nStim(2),1))] };

	%	Sort the stimuli by the decision variable.
	stimSorted = sortrows([stim{1}; stim{2}], 2);

	%	Compute the uncorrected p-value for each quantile (based on this simulation).
	pRaw = zeros(size(iqLow));
	for iq = 1:length(pRaw)
		pRaw(iq) = mean(stimSorted(iqLow(iq):iqHigh(iq),1));
	end
	%	Compute the monotonically corrected p-value for each quantile (based on this simulation).
	%pMonotonic = Mcl_ForceNondecreasing(pRaw);
	Mcl_ForceMonotonic(pNond, pRaw, int32(0));
	Mcl_ForceMonotonic(pMono, pRaw, int32(1));
	Mcl_ForceMonotonic(pMonotonic, pRaw, IDMETHOD);
	
	%	Store the MseReduction
	MseReduction(iSim) = mean( (pRaw-pActual).^2 ) / mean( (pMonotonic-pActual).^2 );
end
%------------------------------------------------------------------------------------------

qParams = fminsearch(@(a) f_ConditionalEntropy_UnidimensionalQuadratic(a, stim), [0 1 0]);
pFit = interp1(pTotal,1-1./(1+exp(x.^2*qParams(1)+x*qParams(2)+qParams(3))),uTerp);

%------------------------------------------------------------------------------------------
%	Make monotonic regression plot
%------------------------------------------------------------------------------------------
SZ_FIGURE = [300 200];

blueLineColor = [0.3 0.3 1];

figure;
set(gcf, 'Position', [10, 40, SZ_FIGURE*PLOT_SCALE]);
set(gcf, 'RendererMode' ,'manual');
set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0, 0, SZ_FIGURE*PLOT_SCALE]); % [left, bottom, width, height]
szPaper = SZ_FIGURE * 396 / 72;
axFontName = 'Times';
axFontSize = 16*PLOT_SCALE;

plot(binNumber,pRaw,		'Color', blueLineColor,		'LineWidth', PLOT_SCALE*1.0);	% Blue
axis([1, nBins, 0 1]);
hold on;
if any(EXPONENTS~=2)
	plot(binNumber,pFit,		'Color', [1 1 1]*0.6,		'LineWidth', PLOT_SCALE*4.0);	% Grey
end
%plot(binNumber,pActualGaussian,		'Color', [1 1 1]*0.5,		'LineWidth', PLOT_SCALE*4.0);	% Grey
plot(binNumber,pActual,		'Color', [1 1 1]*0.0,		'LineWidth', PLOT_SCALE*4.0);	% Black
plot(binNumber,pMonotonic,	'Color', [1 0 0],			'LineWidth', PLOT_SCALE*2.25);	% Red = chosen method (IDMETHOD)
%plot(binNumber,pMono,		'Color', [0.8 0.5 0],		'LineWidth', PLOT_SCALE*1.25);	% Orange = Monotonic
%plot(binNumber,pNond,		'Color', [0 1 0],			'LineWidth', PLOT_SCALE*0.75);	% Green = Non-decreasing
plot(binNumber,pRaw,		'Color', blueLineColor,		'LineWidth', PLOT_SCALE*1.0);	% Blue = Raw
set(gca, 'FontName', axFontName);
set(gca, 'FontSize', axFontSize);
grid on;
xlabel('Quantile of Kernel Value');
ylabel('Conditional Probability');
print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_ForceMonotonic' num2str(SEPARATION) '_Fig1_' num2str(max(EXPONENTS)) '_.tiff']);
close gcf;
%------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------
%	Make density plot
%------------------------------------------------------------------------------------------
SZ_FIGURE = [300 200];

figure;
set(gcf, 'Position', [10, 40, SZ_FIGURE*PLOT_SCALE]);
set(gcf, 'RendererMode' ,'manual');
set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0, 0, SZ_FIGURE*PLOT_SCALE]); % [left, bottom, width, height]
szPaper = SZ_FIGURE * 396 / 72;
axFontName = 'Times';
axFontSize = 24*PLOT_SCALE;
axis([min(x), max(x), 0 1]);
hold on;
plot(x,dCat{1}/max(dCat{1})*0.9,	'Color', [1 0 0],	'LineWidth', PLOT_SCALE*2.25);
plot(x,dCat{2}/max(dCat{2})*0.9,	'Color', [0 0 1],	'LineWidth', PLOT_SCALE*2.25);
set(gca, 'FontName', axFontName);
set(gca, 'FontSize', axFontSize);
xlabel('\it{s}\rm{_1}');
ylabel('Density');
set(gca, 'XTickMode', 'manual');
set(gca, 'XTick', []);
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', []);
print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_ForceMonotonic' num2str(SEPARATION) '_Fig2_' num2str(max(EXPONENTS)) '_.tiff']);
close gcf;
%------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------
%	Plot simulation results
%------------------------------------------------------------------------------------------
if N_SIMS > 1
	figure; hist(MseReduction,20);
	meanMse = mean(MseReduction);
	stdMse = std(MseReduction);
	title(['mean MSE ratio = ' num2str(single(meanMse)) ' +/- ' num2str(single(stdMse/sqrt(N_SIMS-1)))]);
end
%------------------------------------------------------------------------------------------