function out = f_ConditionalEntropy_UnidimensionalQuadratic(a, stim)

%  function out = f_ConditionalEntropy_UnidimensionalQuadratic(a, )
% -------------------------------------------------------
%	A function that computes the conditional entropy
% -------------------------------------------------------
% -------------------------------------------------------
% -------------------------------------------------------

pPred = 1./(1+ exp(a(1)*stim{1}(:,2).^2 + a(2)*stim{1}(:,2) + a(3)));
out = mean(log(pPred));

pPred = 1-1./(1+ exp(a(1)*stim{2}(:,2).^2 + a(2)*stim{2}(:,2) + a(3)));
out = -0.5* (out + mean(log(pPred)));