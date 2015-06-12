function out = Mcl_FastExp_Ctor

% function out = Mcl_FastExp_Ctor
%	Tables the exp function to facilitate faster evaluation of it (via linear interpolation).
out = struct(...
	'LookupTable', zeros(6826,1), ...
	'i', [   0, 1001, 1902, 2803, 3413, 4414, 5315, 6216], ...
	'n', [1001,  901,  901,  610, 1001,  901,  901,  610], ...
	'zLow', [0, 1, 10, 100, 0, -1, -10, -100], ...
	'zHigh', [1, 10, 100, 709, -1, -10, -100, -709], ...
	'zStep', [0.001, 0.01, 0.1, 1, 0.001, 0.01, 0.1, 1]);
	
% public double[] LookupPos_0_1 = new double[1001];
% public double[] LookupPos_1_10 = new double[901];
% public double[] LookupPos_10_100 = new double[901];
% public double[] LookupPos_100_709 = new double[611];
% public double[] LookupNeg_0_1 = new double[1001];
% public double[] LookupNeg_1_10 = new double[901];
% public double[] LookupNeg_10_100 = new double[901];
% public double[] LookupNeg_100_709 = new double[610];

z = zeros(3413,1);

i = 0:1000;
z(i+1) = i*0.001;

i = 0:900;
z(i+1001) = (i+100)*0.01;

z(i+1902) = (i+100)*0.1;

i = 0:609;
z(i+2803) = i+100;

out.Table(   1:3413) = exp(z);
out.Table(3414:6826) = 1./out.Table(1:3413);