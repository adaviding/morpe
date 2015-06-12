%=================================================================================
% WHAT DOES THIS SCRIPT FILE DO?
%	After the classifier is trained, this script loads the classifier and exports it to XML.
%---------------------------------------------------------------------------------

%	This is the name of the input MAT file and output XML file.
FILENAME = 'TnglPlot';

%	Input and output files.
fileIn = [cd '\' FILENAME  '.mat'];
fileOut = [cd '\' FILENAME  '.xml'];

%	Load the trained classifier.
load(fileIn);

%	Get the trained classifier.
data = struct('Morpe', '');
data.Morpe = rmfield(TnglPlot.FullClassifier, {'Cat', 'X', 'Dv', 'H', 'P', 'Class', 'Xmeans', 'Xvars', 'Xmean', 'Xvar', 'Xscale', 'Xacc', 'wScale', 'wInit', 'DisplayModulus', 'SolverOptions', 'SolverOutput'});

%	Save the file
Struct2Xml(data, fileOut);