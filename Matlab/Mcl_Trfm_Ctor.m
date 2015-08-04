function out = Mcl_Trfm_Ctor(Name, Vals)

% funciton out = Mcl_Trfm_Ctor(Name, [Vals])
% This creates a structure that can identify a transform applied to data.  Use with other Mcl_Trfm methods
%	to easily transform and inverse-transform the columns of a dataset.
%
% The Mcl_Trfm structure is used by the Mcl classifiers to label and transform spatial dimensions that require
%	transformation between two domains.  One domain is usually based on the data's raw unit, but such data is sometimes
%	skewed or distributed in a way that is troublesome for Mcl classifiers.  These classifiers work best when the data
%	is distributed as Gaussian as possible.  Thus, the Mcl_Trfm object specifies a univariate and usually isotonic
%	transform that can be used to map data from one domain into another.  Thesecond domain is defined by the transform
%	(specified by an Mcl_Trfm structure) and a reverse transform (also specified by the Mcl_Trfm structure).
%-----------------------------------------------------------------------------------
% INPUT
%-----------------------------------------------------------------------------------
% Name (string)
%	'None'		An empty placeholder, no transform is applied.
%	'FisherZ'	A safe version of the FisherZ transform is applied.
%	'LogAbs'	A safe version of log(abs(x)) is applied.
% [Vals] (double vector: nParams)
%	Any optional parameter values associated with the transform.
%-----------------------------------------------------------------------------------
% OUTUPT
%-----------------------------------------------------------------------------------
% out (struct)
%	A structure that encapsulates information to identify the transform.
%-----------------------------------------------------------------------------------


if	~isequal(Name, 'None') &&...
	~isequal(Name, 'FisherZ') &&...
	~isequal(Name, 'LogAbs') &&...
	~isequal(Name, 'SgnPow') &&...
    ~isequal(Name, 'Affine')
	error('Unrecognized transform Name.');
end

out = struct(...
	'Class', 'Mcl_Trfm', ...
	'Name', Name, ...
	'Vals', []);
if exist('Vals', 'var')
	out.Vals = Vals;
	if numel(out.Vals)==1		&& isequal(Name, 'FisherZ')
		out.Vals = [out.Vals, 1];
	elseif numel(out.Vals)==1	&& isequal(Name, 'LogAbs')
		out.Vals = [out.Vals, 0];
    elseif numel(out.Vals)==1	&& isequal(Name, 'Affine')
		out.Vals = [out.Vals, 1];
	end
else
	if	isequal(Name, 'None')
	elseif isequal(Name, 'SgnPow')
		%	Signed square-root
		out.Vals = 0.5;
	elseif isequal(Name, 'FisherZ')
		% Limits the extent of the transform
		out.Vals = [1.0001, 1];  %1.00000001;
	elseif isequal(Name, 'LogAbs')
		% Bounds transform domain to > -4.6
		out.Vals = [0.00000001, 0]; %0.00000001;
    elseif isequal(Name, 'Affine')
		% Default subtract 0 divide by 1
		out.Vals = [0 1];
	else
		error('Unrecognized transform Name.');
	end
end