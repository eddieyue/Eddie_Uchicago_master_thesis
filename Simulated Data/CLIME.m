function Theta_hat = CLIME(Theta,lambda)
%
% I use the CVX package to achieve CLIME method 
%
% Input: 
%       Theta: the latent generalized concentration matrix need to
%       estimated
%       lambda: the constrained parameter, need to be greater than zero
% Output:
%       Theta_hat: the estimated Theta by CLIME method

n = size(Theta,1);
Theta_hat = zeros(n);
B = eye(n);

for i = 1:n
    cvx_begin quiet
      variable x(d) 
      minimize( norm( x, 1 ) ) 
      subject to
        norm( (Theta*x-B(:,i)), Inf ) <= lambda; 
    cvx_end
    Theta_hat(:,i) = x;
end 