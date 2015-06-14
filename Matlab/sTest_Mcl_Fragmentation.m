dbstop if error;

NCATS = 2;
NDIMS = 2;
NCYCLES = 2;
NSAMP = 400;

COLORS = {[0 0 1], [1 0 0]};
SYMBOLS = {'d', 'o'};
FONTSIZE = 24;
FONTNAME = 'Times New Roman';
XYFIGSIZE = [400 400];

MARKERSIZE = 12;

NQUANT = 10;

SAVETAG = ['_' num2str(NCYCLES)];

%---------------------------------------------------------------------------------------------------
%	Simulate fragmented data
%---------------------------------------------------------------------------------------------------
% [Cat (0,1), z1, z2, sin(z1*2*pi*NCYCLES), sin(z2*2*pi*NCYCLES)]
Xtable = [zeros(NSAMP,1), rand(NSAMP,2)];
Xtable = [Xtable, sin(Xtable(:,2:3)*2*pi*NCYCLES)];
Xtable = [Xtable, Xtable(:,4).*Xtable(:,5), cos(Xtable(:,2)*2*pi*NCYCLES).*cos(Xtable(:,3)*2*pi*NCYCLES)];

%	The probability that each data point is in category 1
p1 = 1./(1+exp(8*Xtable(:,4).*Xtable(:,5)));
%	The category labels are assigned stochastically by this probability.
Xtable(:,1) = (rand(size(p1))<p1);
%	The indices
ind0 = find(Xtable(:,1)==0);
ind1 = find(Xtable(:,1)==1);
Xcell = {...
	Xtable(ind0,2:7), ...
	Xtable(ind1,2:7) };

%	The range
ax1 = [0, 1, 0, 1]; %[xmin, xmax, ymin, ymax]
ax2 = [-1, 1, -1, 1];
ax3 = [-1, 1, -1, 1];

%-------------------------------------------------------
%	Plot raw data (very fragmented)
%-------------------------------------------------------
Figure1 = figure; hold on;
set(gcf, 'RendererMode' ,'manual'); set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperUnits', 'points');
set(gcf, 'Position', [50, 50, XYFIGSIZE]); set(gcf, 'PaperPosition', [0, 0, XYFIGSIZE]); % [left, bottom, width, height]
title('More Fragmented', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
%xlabel('$z_{2}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
%ylabel('$z_{1}$', 'Interpreter', 'latex', 'FontSize', FONTSIZE+4, 'FontName', FONTNAME);
for iCat=1:2
	plot(Xcell{iCat}(:,1), Xcell{iCat}(:,2), SYMBOLS{iCat}, 'Color', COLORS{iCat}, 'LineWidth', 1, 'MarkerSize', MARKERSIZE);
end
axis square; axis(ax1); grid on;
print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_Fragmentation' SAVETAG '_Fig1.tiff']);
%-------------------------------------------------------

%-------------------------------------------------------
%	Plot fragmented data
%-------------------------------------------------------
Figure2 = figure; hold on;
set(gcf, 'RendererMode' ,'manual'); set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperUnits', 'points');
set(gcf, 'Position', [50, 50, XYFIGSIZE]); set(gcf, 'PaperPosition', [0, 0, XYFIGSIZE]); % [left, bottom, width, height]
title('Fragmented', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
%xlabel('$sin(kz_{2})$', 'Interpreter', 'latex', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
%ylabel('$sin(kz_{1})$', 'Interpreter', 'latex', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
for iCat=1:2
	plot(Xcell{iCat}(:,3), Xcell{iCat}(:,4), SYMBOLS{iCat}, 'Color', COLORS{iCat}, 'LineWidth', 1, 'MarkerSize', MARKERSIZE);
end
axis square; axis(ax2); grid on;
print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_Fragmentation' SAVETAG '_Fig2.tiff']);
%-------------------------------------------------------

%-------------------------------------------------------
%	Plot least fragmented data
%-------------------------------------------------------
Figure3 = figure; hold on;
set(gcf, 'RendererMode' ,'manual'); set(gcf, 'Renderer', 'painters');
set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperUnits', 'points');
set(gcf, 'Position', [50, 50, XYFIGSIZE]); set(gcf, 'PaperPosition', [0, 0, XYFIGSIZE]); % [left, bottom, width, height]
title('Less Fragmented', 'FontSize', FONTSIZE, 'FontName', FONTNAME);
%xlabel('$sin(kz_{1})*sin(kz_{2})$', 'Interpreter', 'latex', 'FontSize', FONTSIZE/2, 'FontName', FONTNAME);
%ylabel('$cos(kz_{1})*cos(kz_{2})$', 'Interpreter', 'latex', 'FontSize', FONTSIZE/2, 'FontName', FONTNAME);
for iCat=1:2
	plot(Xcell{iCat}(:,5), Xcell{iCat}(:,6), SYMBOLS{iCat}, 'Color', COLORS{iCat}, 'LineWidth', 1, 'MarkerSize', MARKERSIZE);
end
axis square; axis(ax3); grid on;
print(gcf, '-dtiff', '-r396', '-painters', ['s_TestMcl_Fragmentation' SAVETAG '_Fig3.tiff']);
%-------------------------------------------------------
	
%---------------------------------------------------------------------------------------------------
close(Figure1);
close(Figure2);
close(Figure3);