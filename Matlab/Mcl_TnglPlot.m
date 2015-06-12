function out = Mcl_TnglPlot(Xcell, ops)

% function out =  Mcl_TnglPlot(Xcell, ops)
% Plots a triangular scatterplot in the currently open figure window (using the entire window).  A triangular scatterplot
%	is an arrangement of plots arranged a nDims*nDims grid.  The purpose of the plots is to illustrate the multivariate
%	distribution of data composed of one or multiple categories.  The upper-right triangle of the grid is composed of
%	nDims*(nDims-1)/2 bivariate scatterplots.  The diagonal of the grid is composed of univariate histograms.  You can
%	supply an optimized classifier (e.g. Mcl_Exemplar) and the classifier's output will be overlaid.  To see an example
%	of how this function is used and how to incorporate a trained classifier, see the script file sMcl_Exemplar_Test2.m.
% The default options for this function can be obtained by calling:
%	Mcl_TngPlot('OPTIONS')
%-------------------------------------------------------
% NOMENCLATURE
% nCats
%	The number of categories, also equal to the length of input Xcell.
% nDims
%	The number of spatial dimensions or the dimensionality of the multivariate space, i.e. size(Xcell{iCat},2).
% nSamp(iCat)
%	The number of samples provided (in Xcell) for each iCat, i.e. size(Xcell{iCat},1).
%-------------------------------------------------------
% INPUT
% Xcell (cell array: nCats)
%	A cell array with a length equal to the number of categories, nCats.  Each iCat element of Xcell, Xcell{iCat}, contains
%	samples from a particular category {1,...,nCats}.
% Xcell{iCat} (double 2D array: nSamp(iCat), nDims)
%	Each data sample is entered as a row in this array.  Each row is a multivariate coordinate.
% ops (matlab structure)
%	A structure of options that controls the appearance of the figure.  You can obtain default options by calling:
%		Mcl_TngPlot('OPTIONS')
%-------------------------------------------------------
% OUTPUT
% The output is similar to the options structure.  In addition, out.SubplotHandles contains a handle to each of the subplots
%	in the grid.
%-------------------------------------------------------

RAISE_LOWER_LEFT = false;

%	Compute screen size and determine the screen rectangle for capturing the image.
screenSize = get(0,'ScreenSize');

%	Return default options.
if nargin==1 && ischar(Xcell)
	out = struct(...
		'Class', 'Mcl_TnglPlot', ...
		'Axes', 'Drawing instructions for the axes.  This can be empty (for automatic behavior), an Mcl_Axis scalar, or an Mcl_Axis vector of length size(Xcell{1},2)).', ...
		'BlurLevel', 'See usage in Mcl_TruncHist2.  Default is 0.12.', ...  See usage in Mcl_TruncHist2.m
		'Colors', 'A cell array the same length as X giving the color for each category.  Each color is a vector of length 3 with rgb values in the range [0,1].', ...
		'Dpi', 300, ...
		'DrawAcc', true, ...
		'DrawH', true, ...
		'DrawBounds', true, ...
		'DrawUniTape', true, ...
		'FigureSizeXy', [0 0]+(min(screenSize(3:4))-200), ... [width, height]
		'FontName_Acc', 'Cambria', ...
		'FontName_Axis', 'Times New Roman', ...
		'FontName_Label', 'Cambria', ...
		'FontSize_Acc', 12, ...
		'FontSize_Axis', 11, ...
		'FontSize_Label', 14, ...
		'ForceEqualPriors',  true, ...
		'FullClassifier', 'An Mcl classifier trained to the full dimensional space.  For example, see Mcl_Poly or Mcl_Exemplar.', ...
		'FullLabel', '$\hat{p}(c|\bf{x})$', ...
		'FullLabel_TextInterpreter', 'latex', ...
		'Gamma', 'This is the gamma correction exponent for the scatterplots in the upper-right triangle.  If length(x)==1, default is 1; otherwise, default is 2.', ...
		'HistSize_1d', 20, ...
		'HistSize_2d', [100, 100], ...
		'LineWidth', 2, ...
		'SubplotHandles', 'On output, returns the handles to each subplot in the 2-D grid arrangement.', ...
		'SubplotRects', 'Defaults to the value of Mcl_SubplotRects(nd,nd) where "nd" is the number of columns in the data, (a.k.a. the number of spatial dimensions).', ...
		'TnglClassifiers', 'A set of Mcl classifieres, each trained to a unique pair of spatial dimensions.  For example, see Mcl_Poly or Mcl_Exemplar_Tngl.' );
	return;
end
%	Retrieve default options if necessary.
if nargin < 2
	ops = Mcl_TnglPlot('OPTIONS');
else
	if ~isstruct(ops)
		error('The second input argument (if provided) must be a Matlab structure. Evaluate Mcl_TnglPlot(''OPTIONS'') for the default structure.');
	end
	if isfield(ops, 'Size_Plot')
		error('The Size_Plot field of the options structure is deprecated. Use the FigureSizeXy field instead.');
	end
end

if ~iscell(Xcell)
	error('The input Xcell must be a cell array.');
end

%	Determine how many categories, how many stimuli in each category, and the number of dimensions.
nCats = length(Xcell);
nStim = zeros(nCats,1);
nDims = size(Xcell{1},2);
for iCat=1:nCats
	nStim(iCat) = size(Xcell{iCat},1);
	if nDims~=size(Xcell{iCat},2)
		error('Samples from each category must have the same number of dimensions (i.e. all elements of Xcell must have the same number of columns).');
	end
end

if nDims<=3
	ADDEND_FONTSIZE_LABEL_LOWERLEFT = 0;
elseif nDims==4
	ADDEND_FONTSIZE_LABEL_LOWERLEFT = 2;
else
	ADDEND_FONTSIZE_LABEL_LOWERLEFT = 4;
end

%	Create weights to equalize the priors for the accuracy calculations
if ops.ForceEqualPriors
	catWeight = mean(nStim)./nStim;
else
	catWeight = ones(size(nStim));
end

Xcell = Mcl_CleanSamples(Xcell, true);
for iCat=1:nCats
	nStim(iCat) = size(Xcell{iCat},1);
	if nStim(iCat)<nDims
		error(['Not enough rows samples for category ' num2str(iCat)]);
	end
end

%	Default subplot rectangles
if isempty(ops.SubplotRects) | ischar(ops.SubplotRects)
    % last argument is [Left, Top, Right, Bottom]
	if nDims<=3
		ops.SubplotRects = Mcl_SubplotRects(nDims,nDims,0.07,0.07,[8/ops.FigureSizeXy(1), 8/ops.FigureSizeXy(2), 10/ops.FigureSizeXy(1), 24/ops.FigureSizeXy(2)]);
	else
		ops.SubplotRects = Mcl_SubplotRects(nDims,nDims,0.035,0.035,[4/ops.FigureSizeXy(1), 4/ops.FigureSizeXy(2), 10/ops.FigureSizeXy(1), 24/ops.FigureSizeXy(2)]);
	end
end
%	Default classifiers
if ischar(ops.TnglClassifiers)
	ops.TnglClassifiers = [];
end
if ischar(ops.FullClassifier)
	ops.FullClassifier = [];
end
if ischar(ops.Axes)
	ops.Axes = [];
end

%	Handles denote the output of each subplot
ops.SubplotHandles = zeros(nDims,nDims);

%	Gamma correction
if ischar(ops.Gamma)
	if nCats>=2
		ops.Gamma = 2;
	else
		ops.Gamma = 1;
	end
end
%	Category colors.
if ischar(ops.Colors)
	ops.Colors = Mcl_ColorSet(nCats);
end
darks = zeros(size(ops.Colors));
for iCat = 1:nCats
	darks(iCat) = 1-max(ops.Colors{iCat});
end
%	Deprecation
if isfield(ops, 'Labels')
	error('ops.Labels is deprecated.  Use ops.Axes.TextLabel field instead.');
end
if isfield(ops, 'TextInterpreter')
	error('ops.TextInterpreter is deprecated.  Use ops.Axes.TextInterpreter or ops.FullLabel_TextInterpreter instead.');
end
%	Blur
if ischar(ops.BlurLevel)
	ops.BlurLevel = 0.12;
end
%	For stroking of text drawn on the bivariate scatterplots.
STROKE_SET = [[-4:-1,1:4]*72/ops.Dpi, 0];

%--------------------------------------------------------------------
%	Axes (automatic default)
%--------------------------------------------------------------------
if numel(ops.Axes)==1
	ops.Axes = ops.Axes(ones(nDims,1));
elseif numel(ops.Axes)>1 && numel(ops.Axes)~=nDims
	error('The length of ops.Axes must match the number of dimensions, or the length can be 1 or empty.');
elseif isempty(ops.Axes)
	%	Measure the limits of each dimension
	xSmall	= min(Xcell{1})';
	xBig	= max(Xcell{1})';
	for iCat=2:nCats
		xSmall	= min(xSmall, min(Xcell{iCat})');
		xBig	= max(xBig,   max(Xcell{iCat})');
	end
	%	Constructs all the axes at once
	ops.Axes = Mcl_Axis_Ctor([xSmall,xBig], repmat(xSmall, [1 3]) + (xBig-xSmall)*[0.1, 0.5, 0.9]);
else
	if ~isstruct(ops.Axes) || ~isfield(ops.Axes, 'Class') || ~isequal('Mcl_Axis', ops.Axes(1).Class)
		error('The members of ops.Axes must be Mcl_Axis structures.');
	end
end
%	Transform Xcell
Xcell = Mcl_Trfm_Cols(ops.Axes,true,Xcell);
%--------------------------------------------------------------------

set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [10, 40, ops.FigureSizeXy]);
set(gcf, 'RendererMode' ,'manual');
set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0, 0, ops.FigureSizeXy]); % [left, bottom, width, height]
szPaper = ops.FigureSizeXy * ops.Dpi / 72;
%set(gcf, 'PaperSize', ops.FigureSizeXy);

%--------------------------------------------------------------------
%	Diagonal histograms
%--------------------------------------------------------------------
clsfrDv = [];
clsfrP = [];
for iRow=1:nDims
	%	Compute the bin centers
	bins = (0.5:ops.HistSize_1d)'/ops.HistSize_1d*diff(ops.Axes(iRow).T_Lims)+ops.Axes(iRow).T_Lims(1);
	%	Compute the histogram for Xcell and make it plottable.
	hx = cell(nCats,1);
	hmax = 0;
	for iCat = 1:nCats
		hx{iCat} = hist(Xcell{iCat}(:,iRow),bins)'*catWeight(iCat);
		hmax = max(hmax,max(hx{iCat}));
	end
	%	Normalize the histogram to a max peak of 1
	for iCat = 1:nCats
		hx{iCat} = hx{iCat}/hmax;
		hx{iCat} = [hx{iCat}(1); hx{iCat}; hx{iCat}(length(hx{iCat}))];
	end
	%	Bins
	bins = [ops.Axes(iRow).T_Lims(1); bins; ops.Axes(iRow).T_Lims(2)];
	%	Begin the subplot for the diagonal
	lbwh = reshape(ops.SubplotRects(iRow,iRow,:),1,4);
	ops.SubplotHandles(iRow,iRow) = subplot('position', lbwh);
	hold on;
	
	accStr = '';
	hStr = '';
	
	%	Make finer sampling.
	nTape = max(100,2*ops.HistSize_1d);
	xTape = (0.5:nTape)'/nTape*diff(ops.Axes(iRow).T_Lims)+ops.Axes(iRow).T_Lims(1);
		
	if nCats>1 && ops.DrawUniTape && ~isempty(ops.TnglClassifiers)
		%	Ensure that memory is allocated for classification
		if isempty(clsfrDv) || size(clsfrDv,1)~=numel(xTape) || size(clsfrDv,2)~=nCats
			clsfrDv = zeros(numel(xTape),nCats);
			clsfrP  = zeros(numel(xTape),nCats);
		end
		%	Calculate decision values from tape.
		if isequal('Mcl_Exemplar', ops.TnglClassifiers.UniClassifier(iRow).Class)
			Mcl_Exemplar_CalcDv(ops.TnglClassifiers.UniClassifier(iRow), clsfrDv, xTape, false);
		elseif isequal('Mcl_Poly', ops.TnglClassifiers.UniClassifier(iRow).Class)
			Mcl_Poly_CalcDv(ops.TnglClassifiers.UniClassifier(iRow), clsfrDv, ...
				Mcl_Poly_Expand(xTape,ops.TnglClassifiers.UniClassifier(iRow).CoeffDims) );
		else
			error(['Unrecognized class of univariate classifier:  ' ops.TnglClassifiers.UniClassifier(iRow).Class]);
		end
		%	Map the decision values.
		Mcl_MapDv(clsfrP,clsfrDv,ops.TnglClassifiers.UniClassifier(iRow).Quant);
		%	Determine the category of each bin (iMax)
		[pMax, iMax] = max(clsfrP,[],2);
		%	Draw the color of the category associated with each quantile
		for iTape=2:numel(xTape)
			if pMax(iTape)-0.001 > 1/double(ops.TnglClassifiers.UniClassifier(iRow).Ncats)
				tapeColor = ops.Colors{iMax(iTape)};
			else
				tapeColor = [1 1 1];
			end
			patch( xTape([iTape-1, iTape, iTape, iTape-1]), [1.1, 1.1,1,1], tapeColor, 'EdgeColor', 'none' );
		end
		%	Determine the accuracy
		accStr = fGetAccStr(ops.TnglClassifiers.UniClassifier(iRow).Acc);
		hStr = fGetHstr(ops.TnglClassifiers.UniClassifier(iRow).h);
	end
	thx = cell(size(hx));
	for iCat = 1:nCats
		thx{iCat} = interp1(bins, hx{iCat}, xTape, 'cubic');
		plot(xTape, thx{iCat}*0.9, 'Color', ops.Colors{iCat}, 'LineWidth', ops.LineWidth);
	end
	for iCat = nCats:-1:1
		plot(xTape, thx{iCat}*0.9, 'Color', ops.Colors{iCat}, 'LineWidth', 1);
	end
	pixHeight = lbwh(4)*ops.FigureSizeXy(2);
	ax = [ ops.Axes(iRow).T_Lims(1), ops.Axes(iRow).T_Lims(2), 0, min(2,pixHeight/(pixHeight-2*(6+ops.FontSize_Acc))) ];

	axis(ax);
	set(gca, 'XTickMode', 'manual');
	set(gca, 'XTick', ops.Axes(iRow).T_TickVals);
	set(gca, 'XTickLabel', ops.Axes(iRow).TickLabels);
	set(gca, 'YTickMode', 'manual');
	set(gca, 'YTick', []);
	set(gca, 'TickLength', [0.05 0.05]);
	set(gca, 'TickDir', 'out');
	set(gca, 'FontName', ops.FontName_Axis, 'FontSize', ops.FontSize_Axis);
	
	%	Draw labeled ticks inside of box.
	for iTick=1:length(ops.Axes(iRow).TickLabels)
		if ~isempty(ops.Axes(iRow).TickLabels{iTick})
			plot( ops.Axes(iRow).T_TickVals(iTick)+[0 0], ax(3)+[0, 0.1*(ax(4)-ax(3))], 'Color', [0 0 0], 'LineWidth', 0.5);
		end
	end
	
	xPerPix = (ax(2)-ax(1))/(72/ops.Dpi*szPaper(1)*ops.SubplotRects(iRow,iRow,3));
	yPerPix = (ax(4)-ax(3))/(72/ops.Dpi*szPaper(2)*ops.SubplotRects(iRow,iRow,4));
	
	%	Draw black text without stroking.
	xStroke=0;
	yStroke=0;
	col = [0 0 0];
	if ops.DrawAcc && ~isempty(accStr)
		text( ax(1)+2*xPerPix+xStroke*xPerPix*72/ops.Dpi, ax(4)+yStroke*yPerPix*72/ops.Dpi, accStr,'Color', col, ...
			'FontName', ops.FontName_Acc, 'FontWeight', 'Normal', 'FontSize', ops.FontSize_Acc+1, ...
			'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
	end
	%	Draw entropy
	if ops.DrawH && ~isempty(hStr)
		text( ax(1)+2*xPerPix+xStroke*xPerPix*72/ops.Dpi, ax(4)-(2*72/ops.Dpi+ops.FontSize_Acc)*yPerPix+yStroke*yPerPix*72/ops.Dpi, hStr,'Color', col, ...
			'FontName', ops.FontName_Acc, 'FontWeight', 'Normal', 'FontSize', ops.FontSize_Acc+1, ...
			'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
	end
	if ~isempty(ops.Axes(iRow).TextLabel)
		text( ax(2)-ops.FontSize_Label/2*xPerPix, ax(4)-4*yPerPix, ops.Axes(iRow).TextLabel,'Color', col, ...
			'FontName', ops.FontName_Label, 'FontWeight', 'normal', 'FontSize', ops.FontSize_Label, ...
			'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter',  ops.Axes(iRow).TextInterpreter);
	end
	%	Ensure that axes are black and not covered by any density plots.
	plot(ax(1:2), ax(3)+[0 0], 'Color', [0 0 0]);
	plot(ax(1)+[0 0], ax(3:4), 'Color', [0 0 0]);
	%	Flash the drawing queue
	drawnow;
end
%--------------------------------------------------------------------

%--------------------------------------------------------------------
%	Lower-Left corner.  Show distribution of category membership probabilities for the correct category.
%--------------------------------------------------------------------
if ~isempty(ops.FullClassifier) && nDims>1
	%	Initiate the subplot
	lbwh = squeeze(ops.SubplotRects(nDims,1,:))';
	lbwh(3:4) = lbwh(3:4) * max(1,(nDims-1)/2);
	%	Shorten the plot.
	if RAISE_LOWER_LEFT
		lbwh(2) = lbwh(2) + 0.2*lbwh(4);
		lbwh(4) = lbwh(4) * 0.6;
	else
		lbwh(4) = lbwh(4) * 0.8;
	end
	pixHeight = lbwh(4)*ops.FigureSizeXy(2);
	ax = [0 1 0 min(2,pixHeight/(pixHeight-2*(6+ops.FontSize_Acc)))];
	
	xPerPix = (ax(2)-ax(1))/(72/ops.Dpi*szPaper(1)*lbwh(3));
	yPerPix = (ax(4)-ax(3))/(72/ops.Dpi*szPaper(2)*lbwh(4));
	
	subplot('position', lbwh);
	set(gca, 'TickLength', [0.025 0.025]);
	set(gca, 'XTickMode', 'manual');
	set(gca, 'XTick', (0:nCats)/nCats);
	xTickLabels = cell(1,nCats+1);
	xTickLabels{1} = '  0';
	xTickLabels{nCats+1} = '1  ';
	for iCat=1:(nCats-1)
		xTickLabels{iCat+1} = num2str(round(100*iCat/nCats)/100);
	end
	set(gca, 'XTickLabel', xTickLabels);
	set(gca, 'YTickMode', 'manual');
	set(gca, 'YTick', []);
	hold on;
	%	Compute the bin centers for the probability histograms
	npBins = floor(0.5*ops.HistSize_1d);
	binCenters = (0.5:npBins)/npBins;
	tBinCenters = [0, (0.5:(npBins*5))/(npBins*5), 1];
	%	Figure out the bin heights for the probability histograms
	if isequal('Mcl_Exemplar' ,ops.FullClassifier.Class) || isequal('Mcl_Poly', ops.FullClassifier.Class)
		%	The density for individual bins
		dProb = zeros(nCats,length(binCenters));
		%	Working memory for classifier
		for iCat=1:nCats
			%	Assume Xcell is the training set and use training set (probabilities are already given by the classifier).
			dProb(iCat,:) = hist(ops.FullClassifier.P(ops.FullClassifier.Cat==int32(iCat-1),iCat),binCenters) * ops.FullClassifier.CatWeight(iCat);
		end
		dProb = dProb./max(dProb(:));
		if nCats==2
			dProb(2,:) = dProb(2,length(binCenters):-1:1);
		end
		tProb = cell(nCats,1);
		for iCat=1:nCats
			%plot([0, binCenters, 1], 0.8*[dProb(iCat,1), dProb(iCat,:), dProb(iCat,size(dProb,2))],'LineWidth',ops.LineWidth,'Color',ops.Colors{iCat});
			tProb{iCat} = interp1([0, binCenters, 1], [dProb(iCat,1), dProb(iCat,:), dProb(iCat,size(dProb,2))], tBinCenters, 'cubic');
			plot(tBinCenters, tProb{iCat},'LineWidth',ops.LineWidth,'Color',ops.Colors{iCat});
		end
		for iCat=(nCats:-1:1)
			plot(tBinCenters, tProb{iCat},'LineWidth',1,'Color',ops.Colors{iCat});
		end
		%	Draw accuracy
		if ops.DrawAcc
			accStr = fGetAccStr(ops.FullClassifier.Acc);
			text( ax(1)+12*xPerPix*72/ops.Dpi, ax(4)-2*yPerPix*72/ops.Dpi,...
				accStr,'Color', [0 0 0], ...
				'FontName', ops.FontName_Acc, 'FontWeight', 'Normal', 'FontSize', ops.FontSize_Acc+2, ...
				'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
		end
		%	Draw entropy
		if ops.DrawH
			hStr = fGetHstr(ops.FullClassifier.h);
			text( ax(1)+12*xPerPix*72/ops.Dpi, ax(4)-2*yPerPix*72/ops.Dpi-(2*72/ops.Dpi+ops.FontSize_Acc+2)*yPerPix,...
				hStr,'Color', [0 0 0], ...
				'FontName', ops.FontName_Acc, 'FontWeight', 'Normal', 'FontSize', ops.FontSize_Acc+2, ...
				'VerticalAlignment', 'top', 'HorizontalAlignment', 'left' );
		end
	else
		error(['Unrecognized class of Full classifier:  ' ops.FullClassifier.Class]);
	end
	set(gca, 'TickLength', [0.05 0.05]);
	set(gca, 'TickDir', 'out');
	set(gca, 'FontName', ops.FontName_Axis, 'FontSize', ops.FontSize_Axis);
	%	p(c|x) label
	if ~isempty(ops.FullLabel)
		if RAISE_LOWER_LEFT
			xlabel(ops.FullLabel, 'Interpreter', ops.TextInterpreter, 'FontName', ops.FontName_Label, 'FontSize', ops.FontSize_Label);
		else
			fs = ops.FontSize_Label+ADDEND_FONTSIZE_LABEL_LOWERLEFT;
			text( ax(2)-fs*xPerPix, ax(4)-2*yPerPix,...
				ops.FullLabel,'Color', col, ...
				'FontName', ops.FontName_Label, 'FontWeight', 'Normal', 'FontSize', fs, ...
				'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', ops.FullLabel_TextInterpreter );
		end
	end
	%	Ensure that axes are black and not covered by any density plots.
	plot(ax(1:2), ax(3)+[0 0], 'Color', [0 0 0]);
	plot(ax(1)+[0 0], ax(3:4), 'Color', [0 0 0]);
	%	Paint an axis on the top side
	plot(ax(1:2), ax(4)+[0 0], 'Color', [0 0 0]);
	%	Paint an axis on the right side
	plot(ax(2)+[0 0], ax(3:4), 'Color', [0 0 0]);
	%	Paint a gray line at the decision rule
	plot([0 0]+1/nCats, [0 1], 'Color', [0 0 0]+0.7);
	%	Set axis range.
	axis(ax);
end
%--------------------------------------------------------------------

%--------------------------------------------------------------------
%	Triangle
%--------------------------------------------------------------------
for iRow=1:nDims
	for iCol = (iRow+1):nDims
		%	Get the axis range
		ax = [ops.Axes(iCol).T_Lims(1) ops.Axes(iCol).T_Lims(2) ops.Axes(iRow).T_Lims(1) ops.Axes(iRow).T_Lims(2)];
		%	Build the histogram for x, scale to a max of 1.0, and invert rows to work with the image plotter
		hx = cell(nCats,1);
		for iCat = 1:nCats
			hx{iCat} = double(Mcl_TruncHist2(Xcell{iCat}(:,iCol), Xcell{iCat}(:,iRow), ax, ops.HistSize_2d, ops.BlurLevel  ));
			%	Truncates very large pixels
			[q, qi] = Mcl_Quantiles(hx{iCat}(:), prod(ops.HistSize_2d-1)/prod(ops.HistSize_2d)-0.01);
			hx{iCat}(qi{2}) = q;
			hx{iCat} = hx{iCat}/(max(hx{iCat}(:)));
			%	Build rgb table
			if iCat==1
				rgb = 0;
				rgbDark = 0;
			end
			%	Cross-multiply (tall,wide) expands to a 2D (tall*wide) array.
			rgb		= rgb		+ hx{iCat}(:)*(ops.Colors{iCat}/(1-darks(iCat)));
			rgbDark = rgbDark	+ hx{iCat}(:)*darks(iCat);
		end
		%	Calculate max RGB
		maxRgb = max(rgb(:,1)+rgb(:,2));
		maxRgb = max(max(rgb(:,1)+rgb(:,3)),maxRgb);
		maxRgb = max(max(rgb(:,2)+rgb(:,3)),maxRgb);
		%	Scale
		rgb		=	(1/maxRgb)*rgb;
		rgbDark =	(1/maxRgb)*rgbDark; % Add 1 
		%	Create an image, impose limits, and gamma correct.
		img = min(1,max(0,(1-[...
			rgb(:,2)+rgb(:,3), ...
			rgb(:,1)+rgb(:,3), ...
			rgb(:,1)+rgb(:,2) ]))).^ops.Gamma;
		%	Multiply with the darkness (after gamma correction).
		rgbDark = 1-rgbDark;
		img(:,1) = img(:,1) .* rgbDark;
		img(:,2) = img(:,2) .* rgbDark;
		img(:,3) = img(:,3) .* rgbDark;
		%	Reshape as image.
		img = reshape(img, [ops.HistSize_2d, 3]);
		%	Begin the subplot by drawing the image into the plot area
		ops.SubplotHandles(iRow,iCol) = subplot('position', reshape(ops.SubplotRects(iRow,iCol,:),1,4));
		image(ax(1:2), ax(3:4), img);
		hold on;
		axis(ax);
		set(gca, 'XTickMode', 'manual');
		set(gca, 'XTick', []);
		set(gca, 'YTickMode', 'manual');
		set(gca, 'YTick', []);
		%	Get the pairwise classifier.  If it does not exist, then iPair will be empty.
		iPair = [];
		iPairReversed = false;
		if ~isempty(ops.TnglClassifiers)
			iPair = find(ops.TnglClassifiers.DimPairs(:,1)==int32(iRow) & ops.TnglClassifiers.DimPairs(:,2)==int32(iCol));
			if isempty(iPair)
				iPair = find(ops.TnglClassifiers.DimPairs(:,1)==int32(iCol) & ops.TnglClassifiers.DimPairs(:,2)==int32(iRow));
				iPairReversed = true;
			end
		end
		
		%	If a model exists, then plot it
		if ~isempty(iPair)
			%	Draw boundaries
			if ops.DrawBounds
				%	Later, these values will be used to help stroke any text that must be drawn over the scatterplot:
				xPerPix = (ax(2)-ax(1))/(72/ops.Dpi*szPaper(1)*ops.SubplotRects(iRow,iCol,3));
				yPerPix = (ax(4)-ax(3))/(72/ops.Dpi*szPaper(2)*ops.SubplotRects(iRow,iCol,4));
				%	Generate a mesh of coordinates over which the classifier will be evaluated.
				xvec	= (0.5:ops.HistSize_2d(1))/ops.HistSize_2d(1)*(ax(2)-ax(1))+ax(1);
				yvecRev = (0.5:ops.HistSize_2d(2))/ops.HistSize_2d(2)*(ax(4)-ax(3))+ax(3);
				%	Note that yvec is always reversed because the y-axis is upside-down on these image plots.
				yvec = yvecRev(length(yvecRev):-1:1);
				if iPairReversed
					[X1mesh,X2mesh] = meshgrid(xvec,yvec);
				else
					[X2mesh,X1mesh] = meshgrid(xvec,yvec);
				end
				if isequal(ops.TnglClassifiers.Class, 'Mcl_Exemplar_Tngl') || isequal(ops.TnglClassifiers.Class, 'Mcl_Poly_Tngl')
					%----------------------------------------------------------------------------
					%	For classifier:  Compute the predicted probabilities over the mesh and draw.
					%----------------------------------------------------------------------------
					%	Init memory (if necessary).
					if size(clsfrDv,1)~=numel(X1mesh)
						clsfrDv = zeros(numel(X1mesh),nCats);
						clsfrP = zeros(numel(X1mesh),nCats);
					end
					%	Prepare mesh coordinates for classification.
					Xmesh = [X1mesh(:),X2mesh(:)];
					%	Evaluate decision values for entire transformed mesh.
					if isequal('Mcl_Exemplar', ops.TnglClassifiers.PairClassifier(iPair).Class)
						Mcl_Exemplar_CalcDv(ops.TnglClassifiers.PairClassifier(iPair), clsfrDv, Xmesh, false);
					elseif isequal('Mcl_Poly', ops.TnglClassifiers.PairClassifier(iPair).Class)
						Mcl_Poly_CalcDv(ops.TnglClassifiers.PairClassifier(iPair), clsfrDv, ...
							Mcl_Poly_Expand(Xmesh,ops.TnglClassifiers.PairClassifier(iPair).CoeffDims) );
					else
						error(['Unrecognized bivariate classifier type:  ' ops.TnglClassifiers.PairClassifier(iPair).Class]);
					end
					%	Map the decision values to probabilities for each coordinate of the mesh.
					Mcl_MapDv(clsfrP ,clsfrDv, ops.TnglClassifiers.PairClassifier(iPair).Quant);
					%----------------------------------------------------------------------------
				else
					error(['Unrecognized classifier Tngl type:  ' ops.TnglClassifiers.Class]);
				end

				if nCats==2
					%----------------------------------------------------------------------------
					%	Draw gray contours for decision boundaries when nCats==2
					%----------------------------------------------------------------------------
					%	Find the contours
					ctrs = Mcl_SeparateContours(contourc(xvec,yvecRev,reshape(clsfrP(:,1),size(X1mesh)), [0.5 0.5])); % to be drawn white
					%	Get their lengths
					cLengths = zeros(size(ctrs));
					for ctCon=1:length(ctrs)
						cLengths(ctCon) = size(ctrs{ctCon},1);
					end
					%	Clip the very short contours
					cLengthThresh = 0.1*max(cLengths);
					for ctCon=1:length(ctrs)
						if cLengths(ctCon) > cLengthThresh
							ctrs{ctCon} = Mcl_SmoothContour(ctrs{ctCon},[ax(2)-ax(1), ax(4)-ax(3)]/100);
							plot(ctrs{ctCon}(:,1), ctrs{ctCon}(:,2), 'LineWidth', ops.LineWidth/2, 'Color', [.7 .7 .7]);
						end
					end
					%----------------------------------------------------------------------------
				else
					%----------------------------------------------------------------------------
					%	Draw gray contours for decision boundaries when nCats>2
					%----------------------------------------------------------------------------
					pMin = min(clsfrP(:))-1;
					%	Subtract the second-highest p-value from the p-values of each row.  So the highest p-value will be >0.0, 2nd 0.0, others <0.0.
					%	Also scale so that the maximum p-value reaches 1.0.
					[pMax,iMax] = max(clsfrP,[],2); % max of each row
					%	Subscript to index
					iMax = numel(X1mesh)*(iMax-1)+(1:numel(X1mesh))';
					clsfrP(iMax)=pMin;
					[p2Max,i2Max] = max(clsfrP,[],2); % 2nd max of each row
					%	Subscript to index
					i2Max = numel(X1mesh)*(i2Max-1)+(1:numel(X1mesh))';
					clsfrP(i2Max)=pMin;
					p3MaxSubZero = max(clsfrP,[],2)-pMax; % 3rd max of each row
					%	Set each row to p3MaxSubZero
					for iCat=1:nCats
						clsfrP(:,iCat) = p3MaxSubZero;
					end
					%	Set row max above zero
					clsfrP(iMax) = pMax-p2Max;
					%	Set 2nd row max below zero, but have it approach zero as it approaches pMax
					clsfrP(i2Max) = p2Max-pMax;
					%	Scale to maximum 1.0
					clsfrP = clsfrP/max(0.1,max(clsfrP(:)));
					%--------------------------------------
					%	Gray
					%--------------------------------------
					for iCat=1:nCats
						%if mean(clsfrP(:,iCat)>0.05)>0.05 % If more than 10% of the samples from iCat have a clsfrP value of >0.1, then draw contours.
							%	Find the contours
							ctrs = Mcl_SeparateContours(contourc(xvec,yvecRev,reshape(clsfrP(:,iCat),size(X1mesh)), [0 0])); % to be drawn gray
							%	Get their lengths
							cLengths = zeros(size(ctrs));
							for ctCon=1:length(ctrs)
								cLengths(ctCon) = size(ctrs{ctCon},1);
							end
							%	Clip the very short contours
							cLengthThresh = 0.1*max(cLengths);
							for ctCon=1:length(ctrs)
								if cLengths(ctCon) > cLengthThresh
									ctrs{ctCon} = Mcl_SmoothContour(ctrs{ctCon},[ax(2)-ax(1), ax(4)-ax(3)]/25);
									plot(ctrs{ctCon}(:,1), ctrs{ctCon}(:,2), 'LineWidth', ops.LineWidth/2, 'Color', [0 0 0]+0.7);
								end
							end
						%end
					end
					%----------------------------------------------------------------------------
				end
				%----------------------------------------------------------------------------
				%	Add text
				%----------------------------------------------------------------------------
				%	Flush painter
				drawnow;
				%	Draw stroked text (black on white), text must have stroked context to be perceptible.
				for xStroke = STROKE_SET
					for yStroke = STROKE_SET
						if sqrt(xStroke^2+yStroke^2) <= max(abs(STROKE_SET))
							%	Write black text with a white stroke around it
							if xStroke==0 && yStroke==0
								col = [0 0 0];
							else
								col = [1 1 1];
							end
							%	For some reason, apparently Matlab uses a smaller font size with the print command inside imagesc plots.
							%	Bump up the font size by 1.
							
							%	Draw accuracy
							if ops.DrawAcc
								accStr = fGetAccStr(ops.TnglClassifiers.PairClassifier(iPair).Acc);
								xLoc = ax(1)*0.00+ax(2)*1.00 - 2*xPerPix*72/ops.Dpi;
								yLoc = ax(3)*1.00+ax(4)*0.00 + 2*yPerPix*72/ops.Dpi; % For some reason, y goes in opposite direction
								text( xLoc+xStroke*xPerPix, yLoc+yStroke*yPerPix,...
									accStr,'Color', col, ...
									'FontName', ops.FontName_Acc, 'FontWeight', 'Normal', 'FontSize', ops.FontSize_Acc+1, ...
									'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
							end
							%	Draw entropy
							if ops.DrawH
								hStr = fGetHstr(ops.TnglClassifiers.PairClassifier(iPair).h);
								xLoc = ax(1)*0.00+ax(2)*1.00 - 2*xPerPix*72/ops.Dpi;
								yLoc = ax(3)*1.00+ax(4)*0.00 + 2*yPerPix*72/ops.Dpi + (2*72/ops.Dpi+ops.FontSize_Acc)*yPerPix;   % For some reason, y goes in opposite direction
								text( xLoc+xStroke*xPerPix, yLoc+yStroke*yPerPix,...
									hStr,'Color', col, ...
									'FontName', ops.FontName_Acc, 'FontWeight', 'Normal', 'FontSize', ops.FontSize_Acc+1, ...
									'VerticalAlignment', 'top', 'HorizontalAlignment', 'right' );
							end
						end
					end
				end
				%----------------------------------------------------------------------------
			end
		end
	end
end
%--------------------------------------------------------------------
out = ops;
return;

function hStr = fGetHstr(h)
hStr = '';
if h < 0.995
	if h < 0.095
		hStr = ['0.0' num2str(round(h*100)) ];
	elseif h<0.005
		hStr = '0.00';
	else
		hStr = ['0.' num2str(round(h*100)) ];
	end
else
	hStr = num2str(round(h*100)/100);
	%	Pad with zeros (if necessary).
	if numel(hStr)==3 && h<9.95
		hStr = [hStr, '0'];
	elseif numel(hStr)==1
		hStr = [hStr, '.00'];
	end
end
return;

function accStr = fGetAccStr(acc)
accStr = [num2str(round(acc*100)) '%'];
return;