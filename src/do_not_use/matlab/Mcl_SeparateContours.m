function out = Mcl_SeparateContours(ccOutput)

% function out = Mcl_SeparateContours(ccOutput)
%	Given output ("contour matrix") from the matlab function contourc, 
%	this returns a cell array where each element is a unique contour (defined by x,y vertices).

%	The contour matrix has two rows.
%	   Contour 1				Subsequent contours
%	[level  x(1), ..., x(N) | ... 					|  ;
%	 numPts y(1), ..., y(N) | ...					|  ] 

out = {};
sz = size(ccOutput);

ic = 1;
ctCon=1;
while ic<sz(2)
	nv = ccOutput(2,ic);
	out{ctCon} = ccOutput(:,(ic+1):(ic+nv))';
	ic = ic+nv+1;
	ctCon = ctCon+1;
end