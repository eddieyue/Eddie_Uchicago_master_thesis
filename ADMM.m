function Theta_hat = ADMM(Sigma,lambda, mu )
%
% ADMM method to estimate the Theta_hat
% 
% Input: 
%       Sigma: estimated correlation matrix
%       lambda: the constrained parameter, need to be greater than zero
%       mu:  the constrained parameter, need to be greater than zero
% Output:
%       Theta_hat: the estimated Theta by ADMM method

if nargin < 3; mu = 0.1; end

n = size(Sigma,1);
err = 1e-4; alpha = eigs(Sigma, 1);

% S =  (alpha^2 * eye(n) - Sigma^2)/alpha;
S =  (alpha^2 * eye(n) - Sigma^2);
R = Sigma/alpha; L = ones(n,1); 
A = chol(S);
% A = chol(S,'lower');
Theta_hat = zeros(n,n);
B = eye(n);

soft_thresh = @(b, rambda) sign(b).*max(abs(b) - 1/rambda,0);
for i = 1:n
    bi = B(:,i);
    Ri = R(:,i);
    xk = L; uk = L; 
%     mu = 1; 
    cond = true;
    while cond 
        yk = min(max(Sigma*xk -bi + uk/mu, -lambda*L),lambda*L);
%         tempb = Ri + R * (yk - uk/mu) + S/alpha * xk;
        tempb = Ri + R * (yk - uk/mu) + A/alpha * xk;
        xk1 = soft_thresh(tempb, mu*alpha);
        uk = uk + mu * (Sigma*xk1 - bi - yk);
        cond = norm( xk1 - xk) > err ;
        %         cond = norm( xk1 - xk) > err && norm( yk1 - yk) > err && norm(uk1 - uk) > err;
        xk = xk1;
    end
    Theta_hat(:,i) = xk;
end

