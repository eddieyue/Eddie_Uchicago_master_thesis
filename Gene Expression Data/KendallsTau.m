function Thetahat = KendallsTau(X,L,lambdaMax)
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       L: Number of dice in lambda sequence
%       lambdaMax: the maximum value of Lambda, default setting is 1
% Output:
%       Thetahat: estimated precision matrices based on different lambda

% Tau = sin(pi/2*corr(X,X,'type','Kendall'));
if nargin < 3; lambdaMax = 1; end
if nargin < 2; L = 10; end
Tau = sin(pi/2*kendalltau_fast(X));
lambdaMin = ADMM_minlambda(Tau);
sprintf('Minimum lambda in Kendalls Tau is found ')
lambda = lambda_range(lambdaMin, L);

Thetahat = ADMM(Tau,lambda,1);
sprintf('Whole Process in Kendalls Tau is completed')

end
