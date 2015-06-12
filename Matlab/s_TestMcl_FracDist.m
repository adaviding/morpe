dbstop if error;

NDIMS = 2;
NFRAC = 1;
DPRIME = 2;
DPRIME_FRAC = 2.8;
NSAMP = 100;

DO_GRAYSCALE = true;
if DO_GRAYSCALE
	cMap = Mcl_ColormapGray(256);
else
	cMap = Mcl_ColormapBlueRed(256);
end

COLORS = {[0 0 1], [1 0 0]};
SYMBOLS = {'d', 'o'};
FONTSIZE = 24;
FONTNAME = 'Times New Roman';
XYFIGSIZE = [420 420];

MARKERSIZE = 12;

NQUANT = 10;

rBasis = Mcl_RandRotate(int32(NDIMS));

SAVETAG = ['_' num2str(NFRAC)];
oDist = [...
	Mcl_FracDist_Ctor(NDIMS, NFRAC, -DPRIME/2*rBasis(1,:), DPRIME_FRAC), ...
	Mcl_FracDist_Ctor(NDIMS, NFRAC, +DPRIME/2*rBasis(1,:), DPRIME_FRAC) ];

nCats = numel(oDist);
Xcell = cell(size(oDist));
Xtable = [];
for iCat=1:nCats
	Xcell{iCat} = Mcl_FracDist_Sample(oDist(iCat), NSAMP);
	%[ iCat, q(x), x ]
	Xtable = [Xtable; repmat(iCat,[NSAMP,1]), zeros(NSAMP,1), Xcell{iCat}];
end

if NDIMS==2
	%---------------------------------------------------------------------------------------------------
	%	Two dimensional workup to illustrate the polynomial classifier
	%---------------------------------------------------------------------------------------------------
	%	Define the plot area and sample the mesh
	if NFRAC==1
		ax = [-1, 1, -1, 1]*4;
	else
		ax = [-1, 1, -1, 1]*5;
	end
	xVec = ax(1):0.1:ax(2);
	[Xmesh, Ymesh] = meshgrid(xVec, xVec);
	Xmesh = [Xmesh(:), Ymesh(:)];
	clear Ymesh;
	
	%-------------------------------------------------------
	%	Plot raw data
	%-------------------------------------------------------
	Figure1 = figure; hold on;
	set(gcf, 'RendererMode' ,'manual'); set(gcf, 'Renderer', 'painters');
	set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperUnits', 'points');
	set(gcf, 'Position', [50, 50, XYFIGSIZE]); set(gcf, 'PaperPosition', [0, 0, XYFIGSIZE]); % [left, bottom, width, height]
	title('Raw Data', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
	xlabel('$x_{2}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	ylabel('$x_{1}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	for iCat=1:2
		plot(Xcell{iCat}(:,1), Xcell{iCat}(:,2), SYMBOLS{iCat}, 'Color', COLORS{iCat}, 'LineWidth', 1, 'MarkerSize', MARKERSIZE);
	end
	axis square; axis(ax); grid on;
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig1.tiff']);
	%-------------------------------------------------------
	
	%-------------------------------------------------------
	%	Measure the gaussian summary stats and plot
	%-------------------------------------------------------
	Mu		= cell(size(oDist));
	Cov		= cell(size(oDist));
	InvCov	= cell(size(oDist));
	covScale = 0;
	for iCat=1:nCats
		Mu{iCat}		= mean(Xcell{iCat});
		Cov{iCat}		= cov(Xcell{iCat});
		covScale		= covScale + 0.5*sum(sum(diag(Cov{iCat})));
	end
	for iCat=1:nCats
		%Cov{iCat}		= Cov{iCat} * covScale / sum(sum(diag(Cov{iCat})));
		InvCov{iCat}	= inv(Cov{iCat});
	end
	%	For drawing circles
	rCir = (0:0.1:(2*pi+0.05))';
	xyCir = [cos(rCir), sin(rCir)];
	%----------------------------
	%	Plot the Gaussian world
	%----------------------------
	Figure2 = figure; hold on;
	set(gcf, 'RendererMode' ,'manual'); set(gcf, 'Renderer', 'painters');
	set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperUnits', 'points');
	set(gcf, 'Position', [50, 50, XYFIGSIZE]); set(gcf, 'PaperPosition', [0, 0, XYFIGSIZE]); % [left, bottom, width, height]
	title('Gaussian Simplification', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
	xlabel('$x_{2}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	ylabel('$x_{1}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	for iCat=1:2
		plot(Mu{iCat}(1), Mu{iCat}(2), SYMBOLS{iCat}, 'MarkerEdgeColor', COLORS{iCat}, 'MarkerFaceColor', COLORS{iCat}, 'MarkerSize', MARKERSIZE+2);
		el = xyCir*chol(Cov{iCat});
		plot(  el(:,1)+Mu{iCat}(1),   el(:,2)+Mu{iCat}(2), '-', 'Color', COLORS{iCat}, 'LineWidth', 2);
		plot(2*el(:,1)+Mu{iCat}(1), 2*el(:,2)+Mu{iCat}(2), ':', 'Color', COLORS{iCat}, 'LineWidth', 1);
	end
	axis(ax); axis square; grid on;
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig2.tiff']);
	%----------------------------
	%-------------------------------------------------------
	
	%-------------------------------------------------------
	%	Compute the Fisher quadratic discriminant and then plot
	%-------------------------------------------------------
	Q = struct(...
		'A', 0.5*(InvCov{1} - InvCov{2}), ...
		'b', Mu{2}*InvCov{2} - Mu{1}*InvCov{1}, ...
		'c', Mu{1}*InvCov{1}*Mu{1}' - Mu{2}*InvCov{2}*Mu{2}' + log(det(Cov{1})) - log(det(Cov{2})) );
	xVec = ax(1):0.1:ax(2);
	yVec = ax(3):0.1:ax(4);
	[xMesh, yMesh] = meshgrid(xVec, yVec);
	xyEval= [xMesh(:), yMesh(:)];
	qMesh = reshape(    xyEval*Q.b' + Q.c,    size(xMesh)    );
	for i=1:numel(qMesh)
		qMesh(i) = qMesh(i) + xyEval(i,:)*Q.A*xyEval(i,:)';
	end
	%	Get contours
	qCons = Mcl_SeparateContours(contourc(xVec,yVec,qMesh, [0 0]));
	for iCon=1:numel(qCons)
		plot( qCons{iCon}(:,1), qCons{iCon}(:,2), 'Color', [1 1 1], 'LineWidth', 4 );
	end
	for iCon=1:numel(qCons)
		plot( qCons{iCon}(:,1), qCons{iCon}(:,2), 'Color', [0 0 0], 'LineWidth', 2 );
	end
	axis(ax); axis square; grid on;
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig3.tiff']);
	%-------------------------------------------------------
	
	%-------------------------------------------------------
	%	Plot color temperature map of the quadratic surface.
	%-------------------------------------------------------	
	Figure4 = figure; hold on;
	set(gcf, 'RendererMode' ,'manual'); set(gcf, 'Renderer', 'painters');
	set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperUnits', 'points');
	set(gcf, 'Position', [50, 50, XYFIGSIZE]); set(gcf, 'PaperPosition', [0, 0, XYFIGSIZE]); % [left, bottom, width, height]
	title('Fisher Quadratic', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
	xlabel('$x_{2}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	ylabel('$x_{1}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	%	Make the range symmetric (so that middle color is equal to decision bound).
	qRng = [min(qMesh(:)), max(qMesh(:))];
	qRng(1) = min(qRng(1), -qRng(2));
	qRng(2) = max(qRng(2), -qRng(1));
	%	Tighten the range (some data will be outside of range) to make plot look nicer.
	qRng = 0.8*qRng;
	imagesc(xVec, yVec, qMesh);
	caxis(qRng);
	if DO_GRAYSCALE
		sgnMap = cMap-0.5;
		gcMap = sqrt(abs(cMap))/sqrt(2).*sign(sgnMap)+0.5;
		colormap(cMap);
	else
		colormap(cMap);
	end
	colorbar('FontSize', FONTSIZE, 'FontName', FONTNAME);
	axis(ax); axis square; grid on;
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig4.tiff']);
	%-------------------------------------------------------
	
	%-------------------------------------------------------
	%	Plot color temperature map of the probability surface.
	%-------------------------------------------------------
	title('Canonical Probability (Fisher)', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
	pMeshFish = 1./(1+exp(-qMesh));
	imagesc(xVec, yVec, pMeshFish);
	caxis([0 1]);
	colormap(cMap);
	colorbar('FontSize', FONTSIZE+4, 'FontName', FONTNAME);
	axis(ax); axis square; grid on;
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig5.tiff']);
	%-------------------------------------------------------
	
	%-------------------------------------------------------
	%	Quantize raw data into bins
	%-------------------------------------------------------
	for iSamp=1:size(Xtable,1)
		x = Xtable(iSamp,3:(NDIMS+2));
		Xtable(iSamp,2) = x*Q.A*x' + x*Q.b' + Q.c;
	end
	%	Sort table by row 2.
	Xtable = sortrows(Xtable, 2);
	%	Get the quantile separation indices for q(x)
	qSepInd = (1:(NQUANT-1))'/NQUANT*NSAMP*2;
	qSep = interp1(Xtable(:,2),qSepInd);
	figure(Figure1); hold on;
	title('Quantize Raw Data', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
	%	For each quantile separator
	for iq = 1:length(qSep)
		%	Get the color
		iqColor = max(0,min(1,(qSep(iq)-qRng(1))/diff(qRng)));
		qColor = cMap(ceil(iqColor*(size(cMap,1)-1.0001)+0.00005),:);
		%	Get contours and plot
		qCons = Mcl_SeparateContours(contourc(xVec,yVec,qMesh, [0 0]+qSep(iq)));
		for iCon=1:numel(qCons)
			plot( qCons{iCon}(:,1), qCons{iCon}(:,2), 'Color', qColor, 'LineWidth', 1 );
		end
	end
	%	Re-plot the data
	for iCat=1:2
		plot(Xcell{iCat}(:,1), Xcell{iCat}(:,2), SYMBOLS{iCat}, 'Color', COLORS{iCat}, 'LineWidth', 1, 'MarkerSize', MARKERSIZE);
	end
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig6.tiff']);
	%-------------------------------------------------------
	
	%-------------------------------------------------------
	%	Plot the quadratic probability mesh from the quantized decision value
	%-------------------------------------------------------
	%	Quantile the decision values
	qDv = zeros(NQUANT,2);
	%	For each quantile
	for iq = 1:NQUANT
		%	Get the upper quantile bound
		if iq==NQUANT
			dvHigh = Inf;
		else
			dvHigh = qSep(iq);
		end
		%	Get the lower quantile bound
		if iq==1
			dvLow = -Inf;
		else
			dvLow = qSep(iq-1);
		end
		%	Index the quantile
		ind = Xtable(:,2)>=dvLow & Xtable(:,2)<dvHigh;
		%	Compute the mean of the decision value in the quantile
		qDv(iq,1) = mean(Xtable(ind,2));
		%	Compute the probability of category 2 in the quantile
		qDv(iq,2) = mean(Xtable(ind,1))-1;
	end
	%	Compute the quadratic probability based on the quantized decision value
	pMeshQuad = reshape(interp1(qDv(:,1), qDv(:,2), max(qDv(1,1),min(qDv(NQUANT,1),qMesh(:)))), size(qMesh));
	%	Activate figure 4.
	figure(Figure4);
	title('Quantized Probability', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
	imagesc(xVec, yVec, pMeshQuad);
	print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_FracDist' SAVETAG '_Fig7.tiff']);
	%-------------------------------------------------------
	
	%---------------------------------------------------------------------------------------------------
	close(Figure1);
	close(Figure2);
	close(Figure4);
else
	%---------------------------------------------------------------------------------------------------
	%	Plot a triangular scatterplot
	%---------------------------------------------------------------------------------------------------
	ops = Mcl_TnglPlot('OPTIONS');
	ops.Axes  = Mcl_Axis_Ctor(repmat([-6, 6], [NDIMS 1]));
	figure; Mcl_TnglPlot(Xcell, ops);
	%---------------------------------------------------------------------------------------------------
end