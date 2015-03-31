function Thetahat = Nonparanormal(X,Lambda )
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       Lambda: the tuning parameters for CLIME method
% Output:
%       Thetahat: estimated precision matrices based on different lambda


[Nsim, Ndim] = size(X); 
Zhat = zeros(Nsim, Ndim);

quant_normal = norminv(T_Delta((1:Nsim)/(Nsim+1),Nsim), 0, 1);


for i = 1:Ndim
    [~, ~, Xind] = unique(X(:,i));
    Zhat(:,i) = quant_normal(Xind);
end 

Tau = corr(Zhat,Zhat);
N = length(Lambda);

Thetahat = zeros([size(Tau), N]);


for n = 1:N
    Thetahat(:,:,n) = CLIME(Tau,Lambda(n));
end

end
