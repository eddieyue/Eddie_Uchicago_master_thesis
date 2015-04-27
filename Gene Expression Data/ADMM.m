function Theta_hat = ADMM(Sigma,lambda,method_ind, quiet_ind, err, mu )
%
% ADMM method to estimate the precision matrix with a sequence of constrain
% parameters. 
% 
% Input: 
%       Sigma: N*N symmetric matrix, estimated correlation matrix
%       lambda: 1*nR vector, the constrain parameter vector from large to less, need to be greater than zero
%       method_ind: [1,2,3], 1 for Kendalls Tau, 2 for Sample Correlation and 3 for Nonparanormal
%       quiet_ind: 'quiet' switch into quiet mode
%       err: the constrained parameter for convergence, defualt setting is
%       1e-6
%       mu:  the constrained parameter, need to be greater than zero, defualt setting is
%       1
% Output:
%       Theta_hat: N*N*nR matrix, the estimated Theta by ADMM method for
%       all lambdas

if nargin < 6; mu = 1; end
if nargin < 5; err = 1e-6;end
if nargin < 4; quiet_ind = 'not quiet';end
if nargin < 3; method_ind = 4; end;
if method_ind == 1; method = 'by Kendalls Tau'; ...
elseif method_ind == 2; method = 'by Sample Correlation'; ...
elseif method_ind == 3; method = 'by Nonparanormal'; ...
else method = '';
end
if strcmp('quiet', quiet_ind); quiet_cond = 0; else quiet_cond = 1; end


nR = length(lambda);n = size(Sigma,1);  I = eye(n);L = diag(I);
soft_threshS = @(b,lambda) sign(b).*max(abs(b) - lambda/2,0);

% Precomputed 
alpha = eigs(Sigma, 1)*1.01;
S =  (alpha^2 * eye(n) - Sigma'*Sigma)/alpha^2; R = Sigma/alpha^2;  
Theta_hat = zeros(n,n,nR);

% Optimization with each column 
for i = 1:n
    ei = I(:,i);
    xk = zeros(n,1); uk = zeros(n,1);    yk = zeros(n,1); 
    
    % Start with iteration for each lambda and save the convergence value
    % of xk, yk and uk 
    for j = 1:nR
        lambdaj = lambda(j);
        cond = true;
        tic;
        while cond 
            yk1 = min(max(Sigma*xk -ei + uk/mu, (-lambdaj*L)),lambdaj*L);
            v = R * (ei + yk1 - uk/mu); b = v + S * xk;
            xk1 = soft_threshS(b, 2/(mu*alpha^2));
            uk1 = uk + mu * (Sigma*xk1 - ei - yk1);
            cond = norm( xk1 - xk) > err || norm(uk1 - uk) > err || norm(yk1 - yk) > err;
            xk = xk1; uk = uk1; yk = yk1;
        end
        t = toc; 
        if quiet_cond
            sprintf('Spend %.02f seconds on the %d th lambda for %d th colunm %s',t, j, i, method)
        end
        Theta_hat(:,i,j) = xk;   
    end
end