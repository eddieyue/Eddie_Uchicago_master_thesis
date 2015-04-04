function Theta_hat = ADMM(Sigma,lambda, err, mu )
%
% ADMM method to estimate the precision matrix
% 
% Input: 
%       Sigma: estimated correlation matrix
%       lambda: the constrained parameter, need to be greater than zero
%       mu:  the constrained parameter, need to be greater than zero
% Output:
%       Theta_hat: the estimated Theta by ADMM method

if nargin < 4; mu = 1; end
if nargin < 3; err = 1;end

n = size(Sigma,1);  I = eye(n);L = diag(I); eta = 1e-4;
soft_threshS = @(b,lambda) sign(b).*max(abs(b) - lambda/2,0);

% Precomputed 
alpha = eigs(Sigma, 1)+ 0.01;
S =  (alpha^2 * eye(n) - Sigma'*Sigma)/alpha^2; R = Sigma/alpha^2;  
Theta_hat = zeros(n);

% Optimization with each column once
for i = 1:n
    ei = I(:,i);
    xk = L; uk = L; 
%     yk = L; 
    
    % Start with iteration
    cond = true;
    while cond 

        yk = min(max(Sigma*xk -ei + uk/mu, (-lambda*L)),lambda*L);
        v = R * (ei + yk - uk/mu);
        b = v + S * xk;
        xk1 = soft_threshS(b, 2/(mu*alpha^2));
        uk1 = uk + mu * (Sigma*xk1 - ei - yk);
        cond = norm( xk1 - xk) > err && norm(uk1 - uk) > err;
%         cond = norm( xk1 - xk) > err;
        xk = xk1; uk = uk1;
    end
    Theta_hat(:,i) = xk;

end


% Theta_hat(abs(Theta_hat) <= eta) = 0;
% Theta_hat