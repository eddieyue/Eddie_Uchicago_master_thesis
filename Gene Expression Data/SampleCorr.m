function Thetahat = SampleCorr(X,L,lambdaMax)
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       L: Number of bins in lambda sequence
%       lambdaMax: the maximum value of Lambda, default setting is 1
% Output:
%       Thetahat: estimated precision matrices based on different lambda
if nargin < 3; lambdaMax = 1; end
if nargin < 2; L = 10; end

Tau = corr(X);
lambdaMin = ADMM_minlambda(Tau);
sprintf('Minimum lambda in Sample Correlation is found')
lambda = lambda_range(lambdaMin, L);

Thetahat = ADMM(Tau,lambda,2);
sprintf('Whole Process in Sample Correlation is completed')

end
