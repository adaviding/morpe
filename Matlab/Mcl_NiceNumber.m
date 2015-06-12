function out = Mcl_NiceNumber(target, direction)

%function out = Mcl_NiceNumber(target, direction)
% Finds a nice number close to "target" suitable for use on plotting axes.
% If direction == -1, "out" will be less than "target".
% If direction == 1, "out" will be greater than "target".
% If direction == 0, "out" will be either greater or less than "target".

if length(target)>1
	error('Target can only be a single number.');
end

if target==0
	out = 0;
	return;
end

if nargin<2
	direction = 0;
end

if isinteger(target)
	out = target;
	return;
end


%	Handle the arbitrary direction
if direction==0
	out1 = Mcl_NiceNumber(target, -1);
	out2 = Mcl_NiceNumber(target, 1);
	if abs(out1-target) < abs(out2-target)
		out = out1;
	else
		out = out2;
	end
	return;
end

%	Remember if target was positive.
pos = target>0;
target = abs(target);
if ~pos
	direction = -direction;
end

targetLog = log10(target);
%targetLogDiv = floor(targetLog);
rngMod = 10.^mod(targetLog,1);
if direction ==1
	rngMod = ceil(10*rngMod)/10;
elseif direction ==-1
	rngMod = floor(10*rngMod)/10;
end
out = 10.^floor(targetLog) * rngMod;

if ~pos
	out = -out;
end