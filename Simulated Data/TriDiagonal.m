function [Theta, Sigma] = TriDiagonal(Ndim, V)
%
% function to generate tridiagonal precision matrix, it's diagnoal should
% be all one and it's sub-diagnoal and sup-diagoal absolute value should
% less than 1
%
% Inputs:
%       Ndim: number of dimensions 
%       V: the value on sub-diagnoal and sup-diagnoal
% Output:
%       Theta: the generated precision matrix
%       Sigma: the variance-covariance matrix to generate data and it
%              equals to inverse of Theta


Theta = eye(Ndim);
Theta = Theta + diag(V*ones(Ndim-1,1),1) + ...
    diag(V*ones(Ndim-1,1),-1);

Sigma = Theta\eye(Ndim);
end 