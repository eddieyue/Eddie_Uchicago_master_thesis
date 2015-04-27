function R = lambda_range(lambdamin, L, lower, upper)
%
% The function to generate range for lambda 
%
% Inputs:
%       lambdamin: the minimum lambda could solve the optimization problem.
%       L: the number of bins in the range, default setting is L = 100;
%       lower: the lower bound for the range, default setting is lower_bound = 1.01 * lambdamin
%       upper: the upper bound for the range, default setting is upper_bound = 1
%
% Outputs:
%       R: a sequnce of lambdas for optimization problem

if nargin < 4; upper = 1; end 
if nargin < 3; lower = 1.01; end
if nargin < 2; L = 100; end 

lower = lambdamin*lower;
R = lower +(upper - lower)/(L-1)*((L-1):-1:0);
end