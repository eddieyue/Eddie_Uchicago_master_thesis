function Theta_hat = Truncate_Theta(Theta, eta)
%
% Truncate the precision matrix with due to imperfect estimation, bascially
% I treat the entry with absolute value less or equal to eta as zero
%
% Input:
%       Theat: the estimated precision matrix
%       eta: the threshold for imperfect estimation, default setting is
%       1e-4
%
% Outputs:
%
%       Theta_hat: the truncated version of precision matrix 
if nargin < 2; eta = 1e-4;end

Theta(abs(Theta) <= eta) = 0;
Theta_hat = Theta;
end