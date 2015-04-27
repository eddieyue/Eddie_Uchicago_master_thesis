function Thetahat = KendallsTau(X,Theta,Lambda,eta, plotkey )
% Input: 
%       X: Input Data, Nsim*Ndim matrix
%       Theta: the original precision matrix 
%       Lambda: the tuning parameters for CLIME method
%       eta: imperfect parameter 
%       plotkey: the key to decide to plot the correlation graph or not
% Output:
%       Thetahat: estimated precision matrices based on different lambda
if nargin < 5, plotkey = 0; end

Tau = sin(pi/2*corr(X,X,'type','Kendall'));
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
    title('ROC Curve of Kendalls Tau Method')
    xlabel('FPR')
    ylabel('TPR')
end

end
