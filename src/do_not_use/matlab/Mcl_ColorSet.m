function out = Mcl_ColorSet(n)

%	Returns a list of random colors as elements of a cell array.  Each random color is an rgb triplet in [0 1]

if n==0
	out = {};
elseif n==1
	out = {[0 0 1]};
elseif n==2
	out = {[1 0 0], [0 0 1]};
elseif n==3
	out = {[1 0 0], [0 0.75 0], [0 0 1]};
elseif n==4
	out = {[1 0 1]*0.7, [1 0 0], [0 0.75 0], [0 0 1]};
elseif n==5
	out = {[1 0 1]*0.7, [1 0 0], [0.5 0.5 0], [0 0.75 0], [0 0 1]};
elseif n==6
	out = {[1 0 1]*0.7, [1 0 0], [0.5 0.5 0], [0 0.75 0], [0 0.5 0.5], [0 0 1]};
elseif n==7
	out = {[1 0 1]*0.7, [1 0 0], [0.8 0.3 0], [0.5 0.7 0], [0 0.75 0], [0 0.6 0.6], [0 0 1]};
elseif n>=8
	cmap = Mcl_ColormapFull(n,false);
	out = cell(1,n);
	for i=1:n
		out{i} = cmap(i,:);
	end
end
%--------------------------------------------------------------------