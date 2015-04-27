function lambdaMin = ADMM_minlambda(Sigma, err, mu )
%
% ADMM method to estimate the precision matrix
% 
% Input: 
%       Sigma: estimated correlation matrix
%       err: the constrained parameter for convergence
%       mu:  the constrained parameter, need to be greater than zero
% Output:
%       lambdaMin: the minimun lambda allows the optimization problem has
%       solution 
if nargin < 3; mu = 1; end
if nargin < 2; err = 1e-3;end

n = size(Sigma,1);  I = eye(n);
% soft_threshS = @(b,lambda) sign(b).*max(abs(b) - lambda/2,0);

% Precomputed 
alpha = eigs(Sigma, 1) * 1.01;
S =  (alpha^2 * eye(n) - Sigma'*Sigma)/alpha^2; R = Sigma/alpha^2;  
lambda = zeros(1,n);

% Optimization with each column once
for i = 1:n
    ei = I(:,i);
    xk = zeros(n,1); uk = zeros(n,1);yk = ones(n,1); 
    
    % Start with iteration
    cond = true;iter=0;
%     while cond && iter<=3000
    while cond
        iter=iter+1;
        b = Sigma*xk - ei + uk/mu;
%         b
        yk1 = Truncate_B_alt(b,mu);
%         yk1
        v = R * (ei + yk1 - uk/mu); c = v + S * xk;
        xk1 = c;
        uk1 = uk + mu * (Sigma*xk1 - ei - yk1);
        cond = norm( xk1 - xk) > err || norm(uk1 - uk) > err || norm(yk1 - yk) > err;
        xk = xk1; uk = uk1; yk = yk1;
%         if(mod(iter,50)==0), 
%             [max(abs(Sigma*xk-ei)) toc]
%             tic
%         end
    end
    lambda(i) = norm(Sigma*xk-ei,'inf');   
end
lambdaMin = max(lambda);
% lambdaMin = lambda;