function out = Mcl_RandVec(setSize, subSetSize)

%function out = Mcl_RandVec(setSize, subSetSize)
% If nargin <= 1
%	This function returns a vector of length setSize with elements equal to
%    a random ordering of the set {1, ..., setSize}.
% If nargin>1
%   This function returns a vector of length subSetSize containing the
%   first subSetSize integers sampled from the randomly ordered set {1, ..., setSize}.

sorted = sortrows([(1:setSize)',rand(setSize,1)], 2);
out = sorted(:,1);

if nargin > 1
    if subSetSize==0
        out = [];
    elseif subSetSize <= setSize && subSetSize>=0
        out = out(1:subSetSize);
    else
        error('Incorrect arguments.  The value of subSetSize must be non-negative and not greater than setSize.');
    end
end