function out = Mcl_FracDist_Ctor(Ndims, Nfrac, Mu, Dprime)

out = struct(...
	'Class', 'Mcl_FracDist', ...
	'Ndims', Ndims, ...
	'Nfrac', Nfrac, ...
	'Mu', [], ...
	'Dprime', [], ...
	'MuFrac', []);

if ~exist('Mu', 'var') || isempty(Mu)
	out.Mu = zeros(1,Ndims);
elseif numel(Mu) ~= Ndims
	error('The input Mu must have Ndims elements.');
else
	out.Mu = reshape(Mu, [1, Ndims]);
end
if ~exist('Dprime', 'var') || isempty(Dprime)
	out.Dprime = 1;
elseif numel(Dprime)>1 || Dprime<0
	error('Dprime must be a non-negative scalar.');
else
	out.Dprime = Dprime;
end

%	Generate the mean of each Gaussian.
Nfrac = round(Nfrac);
if Nfrac < 0
	error('The number of fractures must be a positive integer.');
elseif Nfrac==1
	%	Not fractured
	out.MuFrac = out.Mu;
else
	%	Fractured
	ob = [];
	while size(ob,1)<Nfrac
		ob = [ob; Mcl_RandRotate(int32(Ndims))];
	end
	muShift = mean(ob(1:Nfrac,:))*out.Dprime;
	out.MuFrac = ob(1:Nfrac,:)*out.Dprime + repmat(out.Mu-muShift, [Nfrac,1]);
end

