
%   minimize   ||x||_1
%       s.t.   ||Ax-b||_infty <= c
%   input: matrix A, vector b, scalar c



% problem dimensions (m inequalities in n-dimensional space)
m = 100;
n = 50;

% construct a feasible set of inequalities
% (this system is feasible for the x0 point)
A  = randn(m,n);
x0 = randn(n,1);
b  = A*x0 + rand(m,1);
% feasibility: ||A*x0 - b||_infty <= c
c = max(abs(A*x0-b))*2;

cvx_begin
  variable x(n) % defines variable
  minimize( norm( x, 1 ) ) % minimize   ||x||_1
  subject to
    norm((A*x-b), Inf) <= c; % s.t.   ||Ax-b||_infty <= c
cvx_end
