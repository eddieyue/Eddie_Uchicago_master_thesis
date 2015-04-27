function yk1 = Truncate_B(b,mu)
%
% function to solve the optimization problem:
%           min mu/2*||yk-b||^2_2 + ||yk||_{infty}
%
% Inputs: 
%        b: the constant part in optimization problem 
%        mu: the coefficient of optimization problem
% Output:
%        yk1: updated y

V = zeros(1,length(b)); 
c = abs(b); 
for i = 1:length(b)
    y = b;
    T = c(i);
    y( c>T ) = T * sign( b(c>T) );
    V(i) = mu/2*norm(y-b,2)^2 + norm(y,'inf');
end
[~,Ind] = min(V);
Tm = c(Ind);
yk = b;
yk( c>Tm ) = Tm * sign( b(c>Tm) );
yk1 = yk;
end
    