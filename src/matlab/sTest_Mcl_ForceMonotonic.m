%	A script for testing Mcl_ForceNondecreasing.
%	Run it repeatedly to see how the function works.

PLOT_SCALE = 4;
NOISE_SIZE = 0.2;

yRng = [0, 1];
nBins = 100;
x = (0.5:nBins)'/nBins;
y = x + NOISE_SIZE*randn(size(x));
%y = rand(size(x));
y = max(yRng(1), min(yRng(2), y));

binNumber = 1:length(x);

figure;
hold on;
set(gcf, 'Renderer', 'painters');
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [10, 50, 300*PLOT_SCALE, 190*PLOT_SCALE]);
axFontName = 'Times';
axFontSize = 9*PLOT_SCALE;
plot(binNumber,y, 'b', 'LineWidth', PLOT_SCALE/2);
axis([1, length(x), 0 1]);
grid on;

yForcedNd = zeros(size(y));
yForced = zeros(size(y));
yForcedBlend = zeros(size(y));

nTrips = Mcl_ForceMonotonic(yForcedNd,		y, int32(0));
nTrips = Mcl_ForceMonotonic(yForced,		y, int32(1));
nTrips = Mcl_ForceMonotonic(yForcedBlend,   y, int32(2));

%plot(binNumber,yForcedNd, 'Color', [0 0.8 0], 'LineWidth', 2*PLOT_SCALE);
%plot(binNumber,yForced, 'Color', [1 0 0], 'LineWidth', 1*PLOT_SCALE);
plot(binNumber,yForcedBlend, 'Color', [1 0 0], 'LineWidth', 1*PLOT_SCALE);
hold off;

disp(num2str([nTrips, corr(x,y), mean(y), mean(y)-mean(yForced), min(diff(yForced))]));

set(gca, 'FontName', axFontName);
set(gca, 'FontSize', axFontSize);
%xlabel('Quantile Number');
%ylabel('Probability');
xlabel('Index');
ylabel('Value');