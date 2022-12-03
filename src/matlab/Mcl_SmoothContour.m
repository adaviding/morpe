function out = f_SmoothContour(Xy, ScaleXy)

%  function out = f_SmoothContour(Xy)
%	Returns a smoothed contour.
%-----------------------------------------------------------
% INPUT
%-----------------------------------------------------------
% Xy (double 2D array, nVerts * 2)
%	A contour.  Each row is a vertex.  Each column is a spatial coordinate.
% ScaleXy (double vector, 2)
%	The range of the plotting surface.
%-----------------------------------------------------------
% OUTPUT
%-----------------------------------------------------------
% out (double 2D array, nVerts * 2)
%	The smoothed contour, same size as input Xy.
%-----------------------------------------------------------

%	Default scale.
if nargin<2 || isempty(ScaleXy)
	ScaleXy = [1 1];
end

%	The number of vertices
nVerts = size(Xy,1);

%	The distance between adjacent vertices
Dist = sqrt(diff(Xy(:,1)/ScaleXy(1)).^2 + diff(Xy(:,2)/ScaleXy(2)).^2);

%	The mean unwrapped distance between vertices.
meanDist = mean(Dist);

%	The distance for the wrapping part
wrapDist = sqrt(  sum( (Xy(1,:)-Xy(nVerts,:))./ScaleXy ).^2  );

%	Do wrapping
doWrap = false; %wrapDist*2<meanDist;

%	Span of blur kernel from center to tail.
halfSpan = ceil(max(ScaleXy)*4/meanDist);

%	Cumulative distance
Dist = [wrapDist; cumsum(Dist)];

%	The output
out = zeros(nVerts,2);

if doWrap
	error('This needs to be tested.');
	iHighLim = floor((nVerts-1)/2);
	iLowLim = -iHighLim;
	for iMid = 1:nVerts
		iLow = max(iLowLim,iMid-halfSpan);
		iHigh = min(iHighLim,iMid+halfSpan);
		if iLow<1
			%	Wrap low
			indLow = (nVerts+iLow):nVerts;
		else
			%	Do not wrap low
			indLow = [];
		end
		if iHigh>nVerts
			%	Wrap high
			indHigh = 1:(iHigh-nVerts);
		else
			%	Do not wrap high
			indHigh = [];
		end
		%	Unwrapped
		ind = max(1,iLow):min(nVerts,iHigh);
		w = [...
			exp(-0.5*abs(+Dist(indLow )-Dist(iMid)+Dist(nVerts))); ...
			exp(-0.5*abs(+Dist(ind    )-Dist(iMid))); ...
			exp(-0.5*abs(+Dist(indHigh)-Dist(iMid)-Dist(nVerts))) ];
		w = w / sum(w);
		out(iMid,1) = sum(Xy([indLow, ind, indHigh],1).*w);
		out(iMid,2) = sum(Xy([indLow, ind, indHigh],2).*w);
	end
else
	for iMid = 1:nVerts
		iLow = max(1,iMid-halfSpan);
		iHigh = min(nVerts,iMid+halfSpan);
		ind = iLow:iHigh;
		w = exp(-0.5*abs(Dist(ind)-Dist(iMid)));
		w = w / sum(w);
		out(iMid,1) = sum(Xy(ind,1).*w);
		out(iMid,2) = sum(Xy(ind,2).*w);
	end
end