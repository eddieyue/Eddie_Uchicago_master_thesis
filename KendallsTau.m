function Thetahat = KendallsTau(X,Lambda)
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       Lambda: the tuning parameters for CLIME method
% Output:
%       Thetahat: estimated precision matrices based on different lambda

Tau = sin(pi/2*corr(X,X,'type','Kendall'));
N = length(Lambda);

Thetahat = zeros([size(Tau), N]);


for n = 1:N
    Thetahat(:,:,n) = CLIME(Tau,Lambda(n));
end

end
