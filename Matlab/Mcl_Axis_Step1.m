function out = Mcl_Axis_Step1(ax, nSteps, leaveSpaceAtEnds)

% function out = Mcl_Axis_Step1(ax, nSteps, leaveSpaceAtEnds)
%	Creates a uniform sampling of the axis range in transformed units.  Sampling is linear in transformed units.
%-----------------------------------------------------------------------------------------------------
% INPUT
%-----------------------------------------------------------------------------------------------------
% ax (Mcl_Axis structure scalar)
%	The axis.  An Mcl_Axis structure.
% nSteps (int32 scalar)
%	The number of steps.  This will be numel(out).
% leaveSpaceAtEnds (logical scalar)
%	If false, the uniform sampling spans the full axis limits such that out(1)==ax.T_Lims(1) and out(nSteps)==ax.T_Lims(2).
%	If true, uniform sampling does not extend the full span as a little space is left between ax.T_Lims and the most
%		extreme samples.
%-----------------------------------------------------------------------------------------------------
% OUTPUT
%-----------------------------------------------------------------------------------------------------
% out (double vector: nSteps)
%	A vector of uniform steps along the axis in transformed units.  Sampling is linear in transformed units.
%-----------------------------------------------------------------------------------------------------

if nargin<3 || isempty(leaveSpaceAtEnds)
	leaveSpaceAtEnds = false;
end
if nargin<2 || isempty(nSteps)
	nSteps=10;
end
if nargin<1 || isempty(ax) || ~isstruct(ax) || isfield(ax, 'Class') || ~isequal(ax.Class, 'Mcl_Axis')
	error('The first input argument (ax) must be an Mcl_Axis structure.');
end

if leaveSpaceAtEnds
	out = (0.5:nSteps)/nSteps*(ax.T_Lims(2)-ax.T_Lims(1))+ax.T_Lims(1);
else
	out = (0.5:(nSteps-1))/(nSteps-1)*(ax.T_Lims(2)-ax.T_Lims(1))+ax.T_Lims(1);
end