function [X, Theta] = GenerateData(Nsim,Ndim, CorrMatrix, TransVector, mu, sigma)
% Input Parameters: 
% Nsim:: Int, number of simulations
% Ndim:: Int, number of variates for simulation
% CorrMatrix:: Ndim* Ndim+1 Lower Triangular Matrix, correlation function for each Z
% TransVector:: Ndim Cell Array, mapping function from each Z to each X

% Output Parameters:
% OutputData:: Nsim*Ndim Matrix, simulated data we wanted

if nargin < 6, sigma =1; end
if nargin < 5, mu = 0; end

X = zeros(Nsim, Ndim);
Z = zeros(Nsim, Ndim);
Z(:,1) = CorrMatrix(1,1)*normrnd(mu, sigma, [Nsim 1]);
X(:,1) = TransVector{1}(Z(:,1));
for i = 2:Ndim
    Z0 = normrnd(mu, sigma, [Nsim 1]); 
    Corri = CorrMatrix(i,1:i);
    Z(:,i) = [Z0 Z(:,1:i-1)]*diag(Corri)*ones(i,1);
    X(:,i) = TransVector{i}(Z(:,i));
end
Theta = cov(Z)\eye(Ndim);
end
