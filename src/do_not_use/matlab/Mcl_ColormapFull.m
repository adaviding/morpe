function out = Mcl_ColormapFull( nrows, doWhite )

%	function OUT = Mcl_ColormapFull( NROWS, DOWHITE )
% Returns a vector of length NROWS and 3 colomns suitable as input for
%	matlab's "colormap" funciton.
% The map transitions through the following colors
%	purple (high negative values)
%	blue
%	green
%	white (zero values)  -> only when when DRAWWHITE is true
%	yellow
%	orange
%	red (high positive values)

if nargin < 2
	doWhite = false;
end

if( nrows < 8 )
	error('nrows must be >= 8');
end

highNeg = floor(0.3*nrows);
lowPos = ceil(0.7*nrows);

%	Indices for each zone
ineg = (1:highNeg)';
ipos = (lowPos:nrows)';
imid = ((highNeg+1):(lowPos-1))';

%	Color angles for each zone
aneg = (ineg-0.5)*(1.8*pi/nrows);
apos = (ipos-0.5)*(1.8*pi/nrows);
amid = (imid-0.5)*(1.8*pi/nrows);

c = cos([aneg; amid; apos]);
s = sin([aneg; amid; apos]);
z = zeros(size(c));

R = ...
	(c>0).*(s>0).*(min(1,c+s)) + ...
	(c<0).*(s>0).*(min(1,s)) + ...
	(c>0).*(s<0).*(min(1,c)) + ...
	(c<0).*(s<0).*(min(1,z));
G = ...
	(c>0).*(s>0).*(min(1,s)) + ...
	(c<0).*(s>0).*(min(1,s-c)) + ...
	(c>0).*(s<0).*(min(1,z)) + ...
	(c<0).*(s<0).*(min(1,-c));
B = ...
	(s>0).*(min(1,z)) + ...
	(s<0).*(min(1,-s));

if doWhite
	wmid = (0:(lowPos-highNeg-2))'/(lowPos-highNeg-2);
	wmid = 1-abs(wmid-0.5)*2;

	R(imid) = R(imid).*(1-wmid) + wmid;
	G(imid) = G(imid).*(1-wmid) + wmid;
	B(imid) = B(imid).*(1-wmid) + wmid;
end

out = [R(length(R):-1:1), G(length(R):-1:1), B(length(R):-1:1)];