function X = GenerateEC(Nsim, Ndim, Sigma, nu,Mu)
% function to generate the Elliptical Data, I use Student T distribution to
% simulation the value of Xi
%
% Inputs: 
%       Nsim: Number of simulation
%       Ndim: Number of dimension
%       nu: Degree of freedom for T-distribution
%       Mu: Parameter for elliptical distribution
%       Sigma: the inverse of precision matrix 
% Output:
%       X: Generated EC Data

if nargin < 5, Mu = zeros(Ndim,1); end
if nargin < 4, nu = 1; end

% Compute deterministic matrix
A = chol(Sigma,'lower');  

% Random vector uniformly distributed on unit sphere
U = normr(randn([Ndim Nsim])); 

Xi = abs(trnd(nu, [1, Nsim])); 

if Mu == zeros(Ndim,1)
    X = bsxfun(@times,A*U, Xi);
else
    X = bsxfun(@times,A*U, Xi);
    X = bsxfun(@plus,X, Mu);
end

X = X';
end


