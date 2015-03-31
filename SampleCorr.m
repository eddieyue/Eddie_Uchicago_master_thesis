function Thetahat = SampleCorr(X,Lambda)
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       Lambda: the tuning parameters for CLIME method
% Output:
%       Thetahat: estimated precision matrices based on different lambda
if nargin < 5, plotkey = 0; end


Tau = corr(X);
N = length(Lambda);

Thetahat = zeros([size(Tau), N]);

for n = 1:N
    Thetahat(:,:,n) = CLIME(Tau,Lambda(n));
end

end
