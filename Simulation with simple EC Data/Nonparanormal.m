function Thetahat = Nonparanormal(X,Theta,Lambda,eta, plotkey )
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       mu, sigma: the parameter for gaussian wanted to restore
%       plotkey: the key to decide to plot the correlation graph or not
% Output:
%       Zhat: the data we have them marginally retore to gaussian
%             distribution


if nargin < 5, plotkey = 0; end

[Nsim, Ndim] = size(X); 
Zhat = zeros(Nsim, Ndim);

quant_normal = norminv(T_Delta((1:Nsim)/(Nsim+1),Nsim), 0, 1);


for i = 1:Ndim
    [~, ~, Xind] = unique(X(:,i));
    Zhat(:,i) = quant_normal(Xind);
end 

Tau = corr(Zhat,Zhat);
N = length(Lambda);

Thetahat = zeros([size(Theta), N]);
TPR = zeros(1,N);
FPR = zeros(1,N);

for n = 1:N
    Thetahat(:,:,n) = CLIME(Tau,Lambda(n));
    [TPR(n), FPR(n)] = ROCData(Theta, Thetahat(:,:,n), eta);
end

if plotkey == 1
    figure()
    plot(FPR,TPR)
    title('ROC Curve of Nonparanormal Method')
    xlabel('FPR')
    ylabel('TPR')
end
end
