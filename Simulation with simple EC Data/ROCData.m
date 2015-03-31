function [TPR, FPR] = ROCData(Theta, Theta_hat, rho)
% The function to compute true positive rates and false positive rates
%
% Input Parameters: 
% Theta:: true precision matrix
% Theta_hat: estimated precision matrix
% rho:: threshold to avoid inperfect convergency

% Output Parameters:
% TPR:: True Positive Rates
% FNR:: False Positive Rates

Theta(abs(Theta) <rho) = 0;
Theta_hat(abs(Theta_hat) < rho) = 0;

ZeroInd = Theta == 0;
NonzeroInd = Theta ~= 0;

TPR = sum(Theta_hat(NonzeroInd) ~= 0)/ sum(NonzeroInd(:));
FPR = sum(Theta_hat(ZeroInd) ~= 0)/ sum(ZeroInd(:));




