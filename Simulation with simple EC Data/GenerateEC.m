function X = GenerateEC(Nsim, Ndim,nu,Mu,Theta)
% function to generate the Elliptical Data, I use Student T distribution to
% simulation the value of Xi
%
% Input: 
%       Nsim: Number of simulation
%       Ndim: Number of dimension
%       nu: Parameter for T-distribution
%       Mu: Parameter for elliptical distribution
%       Theta: precision matrix s.t. Theta = inv(Sigma)
% Output:
%       X: Generated EC Data

if nargin < 5, Theta = eye(Ndim); end
if nargin < 4, Mu = zeros(Ndim,1); end


Sigma = Theta\eye(Ndim);
A = chol(Sigma,'lower'); % Compute deterministic matrix 

U = normr(randn([Ndim Nsim])); % random vector uniformly distributed on unit sphere

Xi = abs(trnd(nu, [1, Nsim])); 

if Mu == zeros(Ndim,1)
    X = bsxfun(@times,A*U, Xi);
else
    X = bsxfun(@times,A*U, Xi);
    X = bsxfun(@plus,X, Mu);
end

X = X';
end


