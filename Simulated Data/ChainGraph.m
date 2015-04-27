function [Theta, Sigma] = ChainGraph(Ndim, Pvec, Vvec)
%
% function to generate tridiagonal precision matrix, it's diagnoal should
% be all one and it's sub-diagnoal and sup-diagoal absolute value should
% less than 1
%
% Inputs:
%       Ndim: number of dimensions 
%       Pvec: the position vector
%       Vvec: the value vector corresponding to position vector 
% Output:
%       Theta: the generated precision matrix
%       Sigma: the variance-covariance matrix to generate data and it
%              equals to inverse of Theta

Theta = eye(Ndim);

for i = 1:length(Pvec)
    P = Pvec(i);
    V = Vvec(i);
    Theta = Theta + diag(V*ones(Ndim-P,1),P) + diag(V*ones(Ndim-P,1),-P);
end

Sigma = Theta\eye(Ndim);
end 