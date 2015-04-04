clear
clc
% Generate Elliptical Data
% Generate precision matrix 
Ndim = 10; Nsim = 100; 
Theta = eye(Ndim);
Theta = Theta + diag(0.2*ones(Ndim-1,1),1) + diag(0.2*ones(Ndim-1,1),-1);

% Generate with the T-distribution with parameter = Nu and the
% randomized vector Mu 
Mu = rand(Ndim,1); Nu = 10;
% Theta = Theta' * Theta; % Make theta semi-positive definite
DataEC = GenerateEC(Nsim,Ndim,Nu, Mu, Theta);

% Estimate variance-covariance matrix by correlation matrix
Sigma = corr(DataEC, DataEC);

% Initial Setting for ADMM
lambda = 0.05; 
err = 1e-4;
Theta_hat = ADMM(Sigma, lambda, err);



% mu = 1; n = size(Sigma,1); err = 1; I = eye(n);L = diag(I);
% soft_threshS = @(b,lambda) sign(b).*max(abs(b) - lambda/2,0);
% 
% % Precomputed 
% alpha = eigs(Sigma, 1)+ 0.01;
% S =  (alpha^2 * eye(n) - Sigma'*Sigma)/alpha; R = Sigma/alpha;  
% Theta_hat = zeros(n);
% 
% 
% for i = 1:n
%     ei = I(:,i);
%     %     Ri = R(:,i);
%     xk = L; uk = L; yk = L; 
%     cond = true;
%     while cond 
% 
%         yk = min(max(Sigma*xk -ei + uk/mu, (-lambda*L)),lambda*L);
%         v = R * (ei + yk - uk/mu);
%         b = v + S * xk;
%         xk1 = soft_threshS(b, 2/mu*alpha);
%         uk1 = uk + mu * (Sigma*xk1 - ei - yk);
%         cond = norm( xk1 - xk) > err && norm(uk1 - uk) > err;
% %         cond = norm( xk1 - xk) > err;
%         xk = xk1; uk = uk1;
%     end
% 
%     xk
%     Theta_hat(:,i) = xk;
% 
% end
% Theta_hat
% 
