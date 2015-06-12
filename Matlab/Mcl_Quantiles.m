%function [out, indices] = Mcl_Quantiles(x, qVec)
%   Computes quantiles in each column of your data and returns the quantile
%   bin locations and indices for the samples falling within each
%   quantiled range.
%   INPUT
%   x is a column vector of data (or many columns of data)
%   qVec is a vector of desired quantiles, each element between 0 and 1.
%   For each column of x, a column of quantile values is returned.
%   OUTPUT
%   'out' has as many columns as 'x' and as many rows as 'qVec'.
%   'indices' is a cell array with columns corresponding to columns of 'x' 
%   and rows corresponding to a quantile range, so the number of rows
%   is equal to the length of 'qVec'+1.  Each member of that cell array
%   is a vector of indices indexing the rows of 'x', thus, can be used
%	to sample 'x' one quantile at a time.
function varargout = Mcl_Quantiles(x,qVec)

[rows,cols] = size(x);
out = zeros(length(qVec),cols);
indices = cell( length(qVec)+1, cols );

%drow = 1/(rows-1);

for j=1:cols

    [xx, ui] = sort(x(:,j));

    %ind = 1;
    lastind = 1;

    for i=1:length(qVec)
        ind = qVec(i)*rows;
        f = floor(ind);
        c = ceil(ind);
        if(f==c)
            out(i,j) = xx(ind);
        else
            df = ind-f;
            dc = c-ind;
			if f==0
				out(i,j) = xx(c);
			else
				out(i,j) = (xx(f)*dc+xx(c)*df);
			end
        end
        indices{i,j} = ui( lastind:f );
        lastind = f+1;
    end
    
    indices{length(qVec)+1,j} = ui(lastind:rows);

end

varargout = {out, indices};