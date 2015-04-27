function Thetahat = Nonparanormal(X,L,lambdaMax)
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       L: Number of dice in lambda sequence
%       lambdaMax: the maximum value of Lambda, default setting is 1
% Output:
%       Thetahat: estimated precision matrices based on different lambda
if nargin < 3; lambdaMax = 1; end
if nargin < 2; L = 10; end

[Nsim, Ndim] = size(X); 
Zhat = zeros(Nsim, Ndim);

quant_normal = norminv(T_Delta((1:Nsim)/(Nsim+1),Nsim), 0, 1);


for i = 1:Ndim
    [~, ~, Xind] = unique(X(:,i));
    Zhat(:,i) = quant_normal(Xind);
end 

Tau = corr(Zhat,Zhat);
lambdaMin = ADMM_minlambda(Tau);
sprintf('Minimum lambda in Nonparanormal is found')
lambda = lambda_range(lambdaMin, L);


Thetahat = ADMM(Tau,lambda,3);
sprintf('Whole Process in Nonparanormal is completed')

end
