function out = Mcl_Axis_Ctor(Lims, TickVals, Transform, labeledTix)

% function out = Mcl_Axis_Ctor(Lims, TickVals, Transform)
%	Constructs one or many Mcl_Axis structures and outputs them.
%------------------------------------------------------------------------------------------------------
% NOMENCLATURE
%------------------------------------------------------------------------------------------------------
% nAxes (int scalar)
%	The number of axes.  Equal to numel(out), size(Lims,1), size(TickVals,1), and numel(Transforms).
% mTickVals (int scalar)
%	The maximum number of tix on an axis.  Equal to size(TickVals,2).
% "original units"
%	These units are used for display and are the units of data without any transform applied.
% "transformed units"
%	These are the units after the specified transform was applied to the original units.  The function
%	Mcl_TnglPlot will display histograms and scatterplot using the transformed axis.
%------------------------------------------------------------------------------------------------------
% INPUT
%------------------------------------------------------------------------------------------------------
% Lims (double 2D array: nAxes * 2)
%	Inputs the minimum (column 1) and maximum (column 2) of each axis.
% TickVals (double 2D array: nAxes * mTickVals)
%	Inputs tick values for each axis in original units.  To specify fewer than mTickVals ticks for a given
%	axis, pack the associated row with NaN values.  Each NaN value indicates the absence of a tick.
% Transform (Mcl_Trfm structure: nAxes)
%	Inputs the definition of the transform applied to each axis. 
%------------------------------------------------------------------------------------------------------
% OUTPUT
%------------------------------------------------------------------------------------------------------
% For each axis iAxis in {1,...,nAxes}:
%
% out(iAxis).Lims (double vector: 2)
%	Same as Lims input.
% out(iAxis).T_Lims (double vector: 2)
%	The Lims in transformed units.  If a transform is not specified, original units are used.
% out(iAxis).T_TickVals (double vector: nTickVals)
%	The TickVals in transformed units.  If a transform is not specified, original units are used.
% out(iAxis).TextInterpreter
%	'latex'
% out(iAxis).TextLabel (string)
%	The label for the axis (to be displayed on a plot).
% out(iAxis).TickLabels (cell 2D array of strings:  nAxes * nTickVals)
%	Each iTick entry of the cell array contains a string representation of each tick (ideal for labeling an axis).
%	If empty, default num2str() strings are created for the first (1) and last (nTickVals) tick values.  The remaining
%	ticks may be left as entry strings, but some of them will be filled in depending on the Transform provided.
% out(iAxis).TickVals (double vector: nTickVals)
%	Same as TickVals input.  If input is not provided, this variable is created automatically.
%------------------------------------------------------------------------------------------------------

%	Matlab doesn't do floating point equality very well.  Especially for the number 0.4.
EPSILON = 16*eps;

nAxes = size(Lims,1);
if size(Lims,2)~=2
	error('The input Lims must have two columns.');
end
if any(Lims(:,2)<=Lims(:,1))
	error('All Lims(:,2) must be greater or equal to Lims(:,1).');
end
if nargin<2
	TickVals = [];
end
if ~isempty(TickVals) && size(TickVals,1)~=nAxes
	error('The input TickVals must have nAxes rows, also equal to size(Lims,1)).');
end
if nargin<3
	Transform = [];
end
if ~isempty(Transform) && all(numel(Transform)~=[1 nAxes])
	error('The input argument Transform must have nAxes elements.');
end

out = struct(...
	'Class', '', ...
	'Lims', [], ...
	'T_Lims', [], ...
	'T_TickVals', [], ...
	'TextLabel', '', ...
	'TextInterpreter', 'latex', ...
	'TickLabels', [], ...
	'TickVals', [], ...
	'Transform', [] );

if ~exist('labeledTix', 'var')
	if isempty(TickVals)
		labeledTix = [Lims(:,1), 0.5*(Lims(:,1)+Lims(:,2)) ,Lims(:,2)];
	else
		if mod(size(TickVals,2),2)==1
			labeledTix = [TickVals(:,1), TickVals(:,ceil(size(TickVals,2)/2)), TickVals(:,size(TickVals,2))];
		else
			labeledTix = [TickVals(:,1), TickVals(:,size(TickVals,2))];
		end
	end
end
for iAxis = 1:nAxes
	%	Axis structure class name
	out(iAxis).Class = 'Mcl_Axis';
	%	Axis name
	out(iAxis).TextLabel = ['${\it{x}\rm{_{' num2str(iAxis) '}}}$'];
	out(iAxis).TextInterpreter = 'latex';
	%	Axis limits
	out(iAxis).Lims = Lims(iAxis,:);
	%	Set the transform
	if isempty(Transform)
		out(iAxis).Transform = Mcl_Trfm_Ctor('None');
	elseif numel(Transform)==1
		out(iAxis).Transform = Transform;
	else
		out(iAxis).Transform = Transform(iAxis);
	end
	%	Get the transformed Lims.
	out(iAxis).T_Lims = Mcl_Trfm_Vec(out(iAxis).Transform,true,Lims(iAxis,:));
	%	Axis tick values
	if isempty(TickVals)
		%	Default 5 tick values.
		out(iAxis).T_TickVals = out(iAxis).T_Lims(1) + [0, 0.25, 0.5, 0.75, 1] * (out(iAxis).T_Lims(2)-out(iAxis).T_Lims(1));
		out(iAxis).TickVals = Mcl_Trfm_Vec(out(iAxis).Transform,false,out(iAxis).T_TickVals);
	else
		%	Copy input TickVals
		out(iAxis).TickVals = TickVals(iAxis,:);
		out(iAxis).TickVals = out(iAxis).TickVals(~isnan(out(iAxis).TickVals));
		out(iAxis).T_TickVals = Mcl_Trfm_Vec(out(iAxis).Transform,true,out(iAxis).TickVals);
	end
	%	Calculate the tick labels.
	out(iAxis).TickLabels = cell(1,numel(out(iAxis).TickVals));
	for iTick=1:numel(out(iAxis).TickLabels)
		if any(    abs( labeledTix(iAxis,:)-out(iAxis).TickVals(iTick) )  <  EPSILON    )
			%	Get the string representation
			str = num2str(out(iAxis).TickVals(iTick));
			if abs(out(iAxis).TickVals(iTick)-out(iAxis).Lims(1))  <  EPSILON
				%	Pad left
				%out(iAxis).TickLabels{iTick} = [repmat(' ', [1, numel(str)*2-1]), str];
				%	No padding
				out(iAxis).TickLabels{iTick} = str;
			elseif abs(out(iAxis).TickVals(iTick)-out(iAxis).Lims(2))  <  EPSILON
				%	Pad right
				out(iAxis).TickLabels{iTick} = [str, repmat(' ', [1, numel(str)*2-1])];
			else
				%	No padding
				out(iAxis).TickLabels{iTick} = str;
			end
		end
	end
end