function [TPR, FPR] = ROC_Data(Theta, Theta_hat, rho)
% Input Parameters: 
% Theta:: true precision matrix
% Theta_hat: estimated precision matrix
% rho:: threshold to avoid inperfect convergency

% Output Parameters:
% TPR:: True Positive Rates
% FNR:: False Negative Rates

Theta(abs(Theta) <rho) = 0;
Theta_hat(abs(Theta_hat) < rho) = 0;

NumZeros = sum(Theta(~~tril(Theta,1)) == 0);
NumNonzeros = sum(Theta(~~tril(Theta,1)) ~= 0);

% NumZeros_hat = sum(Theta_hat(~~tril(Theta_hat,1)) == 0);
NumNonzeros_hat = sum(Theta_hat(~~tril(Theta_hat,1)) ~= 0);

TPR = NumNonzeros_hat/NumNonzeros;
FPR = NumNonzeros_hat/NumZeros;




