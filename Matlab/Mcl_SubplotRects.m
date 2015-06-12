function out = Mcl_SubplotRects(nTall, nWide, waxSpacer, taxSpacer, frameLtrb)

%----------------------------------------------------------------------------
% function OUT = Mcl_SubplotRects(NTALL,NWIDE,[WAXSPACER],[TAXSPACER],[FRAMELTRB])
%----------------------------------------------------------------------------
% This function is handy when generating a grid of plots using matlab's "subplot" command.
% Specifically, this function computes the position of each sub-plot within the grid.
% So, if you wanted to generate a 3 (wide) by 4 (tall) grid of plots, you might do something like:
%
% rects = Mcl_SubplotRects(3,4);
% for iWide = 1:3
%     for iTall = 1:4
%			subplot( 'position', squeeze(rects(iTall,iWide,:)) );
%			%  TO DO:  DRAW YOUR PLOT
%     end
% end
%
% Note that the upper-left corner is indexed (iTall=1,iWide=1).
% The lower-right corner is indexed (iTall=NTALL,iWide=NWIDE).
%----------------------------------------------------------------------------
% OUT:  A 3-d double array.
%	OUT(iTall,iWide,1) --> The left position of the subplot rectangle at grid location (iTall,iWide) where the possible range is [0, 1].
%	OUT(iTall,iWide,2) --> The bottom position of the subplot rectangle at grid location (iTall,iWide) where the possible range is [0, 1].
%   OUT(iTall,iWide,3) --> The width of the subplot rectangle at grid location (iTall,iWide) where the possible range is [0, 1-OUT(iTall,iWide,1)].
%	OUT(iTall,iWide,4) --> The height of the subplot rectangle at grid location (iTall,iWide) where the possible range is [0, 1].
%----------------------------------------------------------------------------
% NTALL: The height of the grid
% NWIDE: The width of the grid
% WAXSPACER: (default 0.1) The proportion of each square given underneath the horizontal axis (for spacing)
% TAXSPACER: (default 0.1) The proportion of each square given left of the vertical axis (for spacing)
% FRAMELTRB: (default [0,0.02,0.02,0.04]) The spacing around the entire grid as [left, top, right, bottom] 
%----------------------------------------------------------------------------

if nargin<5
    frameLtrb = [0.00 0.02 0.02 0.04];
end
if isempty(frameLtrb)
	frameLtrb = [0.00 0.02 0.02 0.04];
end
if nargin < 4
    taxSpacer = 0.1;
end
if isempty(taxSpacer)
	taxSpacer = 0.1;
end
if nargin < 3
    waxSpacer = 0.1;
end
if isempty(waxSpacer)
	waxSpacer = 0.1;
end

spacer = 0;

tallness = (1-sum(frameLtrb([2 4])))/nTall;
tallVec = frameLtrb(4):tallness:(1-frameLtrb(2));
tspace = spacer*tallness;
htspace = tspace/2;
taxspace = taxSpacer*tallness;

wideness = (1-sum(frameLtrb([1 3])))/nWide;
wideVec = frameLtrb(1):wideness:(1-frameLtrb(3));
wspace = spacer*wideness;
hwspace = wspace/2;
waxspace = waxSpacer*wideness;

out = zeros(nTall, nWide, 4);

for i=1:nTall
    for j=1:nWide
        out(i,j,:) = [wideVec(j)+hwspace+waxspace, tallVec(nTall-i+1)+htspace+taxspace,...
            wideness-wspace-waxspace, tallness-tspace-taxspace];
    end
end