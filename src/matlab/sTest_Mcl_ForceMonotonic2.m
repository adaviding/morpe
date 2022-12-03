PATH_OUT = 'c:/users/ing/desktop/';

%	A script for testing Mcl_ForceMonotonic.
%	Run it repeatedly (N_SIMS=1000) to see how the function works.

%-----------------------------------------------------------------------
%	Specify how the plotsappears
%-----------------------------------------------------------------------
STYLE = 2; % 1=cumulative Gaussian, 2=exponential
PLOTNEG = false;
SHOWSIDES = true;
NSAMP = 100;
FONTNAME = 'Cambria';
AXFONTSIZE = 24;
TIFONTSIZE = 22;
SZFIGURE = [250 150];
DPI = 300;
LINEWIDTH = 6;
SZMARKER = 8;
COLMARKER = [0 0 0];
IDMETHOD = int32(2); % Method for monotonic regression: int32(0), int32(1), or int32(2).
GRAY = [0.6 0.6 0.6];

%-----------------------------------------------------------------------
%	Simulate data
%-----------------------------------------------------------------------
ax = [-3, 3, 0, 1];
x = 6*((0.5:NSAMP)'/NSAMP-0.5);
if STYLE==1
	xx = x + randn(NSAMP,1)/2;
elseif STYLE==2
	xx = x + 1*(rand(NSAMP,1)-0.5);
end

if STYLE==1
	y = normcdf(xx);
else
	y = exp(-abs(xx));
end

%	Try monotonic increasing
yMono1 = zeros(size(y));
Mcl_ForceMonotonic(yMono1, y, IDMETHOD);
rho1 = abs(corr(y,yMono1));

%	Try monotonic decreasing
y = -y;
yMono2 = zeros(size(y));
Mcl_ForceMonotonic(yMono2, y, IDMETHOD);
rho2 = abs(corr(y,yMono2));

%	Get better fit
y = -y;
yMono2 = -yMono2;
if abs(rho1) > abs(rho2)
	rho = -rho1;
	yMono = yMono1;
else
	rho = rho2;
	yMono = yMono2;
end

%	Index the positive and negative parts
iNeg = find(x<0);
iPos = find(x>0);
yMonoNeg = zeros(size(iNeg));
yMonoPos = zeros(size(iPos));
yNeg = y(iNeg);
yPos = y(iPos);
%	Perform monotonic regression for both parts
if STYLE==1
	Mcl_ForceMonotonic(yMonoNeg, yNeg, IDMETHOD);
	Mcl_ForceMonotonic(yMonoPos, yPos, IDMETHOD);
elseif STYLE==2
	yPos = -yPos;
	Mcl_ForceMonotonic(yMonoNeg, yNeg, IDMETHOD);
	Mcl_ForceMonotonic(yMonoPos, yPos, IDMETHOD);
	yPos = -yPos;
	yMonoPos = -yMonoPos;
end

rhoNeg = abs(corr(yMonoNeg,yNeg));
rhoPos = abs(corr(yMonoPos,yPos));

%-----------------------------------------------------------------------
%	Determine the string output (coefficient values)
%-----------------------------------------------------------------------
sRho = num2str(round(   rho*100   )/100);
sRhoNeg = num2str(round(rhoNeg*100)/100);
sRhoPos = num2str(round(rhoPos*100)/100);

if PLOTNEG
	y = 1-y;
	yMono = 1-yMono;
	yNeg = 1-yNeg;
	yPos = 1-yPos;
	yMonoNeg = 1-yMonoNeg;
	yMonoPos = 1-yMonoPos;
	if STYLE==1
		sRhoNeg = ['-' sRhoNeg];
	elseif STYLE==2
		sRhoNeg = ['-' sRhoNeg];
		sRhoPos = ['-' sRhoPos];
	end
else
	if STYLE==1
		sRhoPos = ['-' sRhoPos];
	end
end

%-----------------------------------------------------------------------
%	Create a figure
%-----------------------------------------------------------------------
figure;
set(gcf, 'Renderer', 'painters');
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [10, 50, SZFIGURE]);
accFigure = subplot('position', [0.08, 0.1, 0.88, 0.75]); %[left, bottom, width, height]
subplot(accFigure);
hold on;
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
set(gca, 'XTick', []);
set(gca, 'YTick', []);
%set(gca, 'FontName', FONTNAME);
%set(gca, 'FontSize', AXFONTSIZE);
axis(ax);
xlabel('Signed Distance from Contour', 'FontName', FONTNAME, 'FontSize', AXFONTSIZE);
ylabel('Pixel value', 'FontName', FONTNAME, 'FontSize', AXFONTSIZE);
%	Draw correlation coefficient
text( 0, 1.05, ['$I_{l}\:\rm{=' sRho '}$'],	'FontName', FONTNAME, 'FontSize', TIFONTSIZE, 'Interpreter', 'latex', 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
%	Plot midline
plot([0 0], ax(3:4), '-', 'Color', [0.5 1 0.5], 'LineWidth', LINEWIDTH*3);
%	Plot full monotonic fit
%if STYLE==1
	plot(x		,yMono,		'-', 'Color', GRAY,		'LineWidth', LINEWIDTH*2.5);
%else
%	plot(x		,yMono,		'-', 'Color', GRAY,		'LineWidth', LINEWIDTH);
%end
if SHOWSIDES
	%	Plot negative monotonic fit
	plot(x(iNeg),yMonoNeg,	'-', 'Color', [0 0 1],	'LineWidth', LINEWIDTH);
	%	Plot positive monotonic fit
	plot(x(iPos),yMonoPos,	'-', 'Color', [1 0 0],	'LineWidth', LINEWIDTH);
	%	Draw correlation coefficients
	text(-2.2, 1.05, ['$I_{l-}\:\rm{=' sRhoNeg '}$'],	'FontName', FONTNAME, 'FontSize', TIFONTSIZE, 'Interpreter', 'latex', 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
	text(+2.2, 1.05, ['$I_{l+}\:\rm{=' sRhoPos '}$'],	'FontName', FONTNAME, 'FontSize', TIFONTSIZE, 'Interpreter', 'latex', 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
end
%	Plot symbols
plot(x,y, 'o', ...
	'MarkerEdgeColor', COLMARKER, 'LineWidth', 2, ...'MarkerFaceColor', COLMARKER, ...
	'Marker', 'o', 'MarkerSize', SZMARKER );


%-----------------------------------------------------------------------
%	Print figure
%-----------------------------------------------------------------------
drawnow;
if ~exist(PATH_OUT, 'dir')
	mkdir(PATH_OUT);
end
outputName = ['Mcl_Mono_' num2str(STYLE) '_' num2str(PLOTNEG) '_' num2str(SHOWSIDES)];
disp(['Printing ' PATH_OUT outputName '.tiff']);
print(gcf, '-dtiff', ['-r' num2str(DPI)], '-painters', [PATH_OUT outputName '.tiff']);

close gcf;
pause(1);