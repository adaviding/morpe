function Xcell = Mcl_CleanSamples(Xcell, dispWarning)

% function Xcell = Mcl_CleanSamples(Xcell, dispWarning)
% Ensures that all samples are free of Nan and Inf values.
%
% Xcell{1}:  Training samples from category 1.  Each row is a sample, each column is a spatial coordinate.
%	...
% Xcell{nCats}:  Training samples from category nCats.

if nargin<2
	dispWarning = true;
end

if iscell(Xcell) 
    nCats = numel(Xcell);
    for iCat=1:nCats
        nSamp = size(Xcell{iCat},1);
        goodInd = find(     ~(  any(isinf(Xcell{iCat}),2) |  any(isnan(Xcell{iCat}),2)  )     );
        if length(goodInd)~=nSamp
            if dispWarning
                disp(['WARNING:  only ', num2str(length(goodInd)/nSamp*100), '% of rows from category ' num2str(iCat) ' did not contain NaN or Inf values.']);
            end
            Xcell{iCat} = Xcell{iCat}(goodInd,:);
        end
    end
else
    nSamp = size(Xcell,1);
    goodInd = find(     ~(  any(isinf(Xcell),2) |  any(isnan(Xcell),2)  )     );
    if length(goodInd)~=nSamp
        if dispWarning
            disp(['WARNING:  only ', num2str(length(goodInd)/nSamp*100), '% of rows did not contain NaN or Inf values.']);
        end
        Xcell = Xcell(goodInd,:);
    end
end