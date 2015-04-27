function yk1 = Truncate_B_alt(b,mu)
%
% function to solve the optimization problem:
%           min mu/2 * ||yk-b||^2_2 + ||yk||_{infty}
%
% Inputs: 
%        yk: old y
%        b: the constant part in optimization problem 
% Output:
%        yk1: updated y

n = length(b);
V = zeros(1,n); 
c = abs(b); 

ymat=min(c*ones(1,n),ones(n,1)*c');
objfun=sum((ymat-c*ones(1,n)).^2)*mu/2+c';
[~,Ind]=min(objfun);
yk1=sign(b).*ymat(:,min(Ind));
end
    