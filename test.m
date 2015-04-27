clear
clc
% Generate Elliptical Data
% Generate precision matrix 
Ndim = 50; Nsim = 1000; 
Theta = eye(Ndim);
Theta = Theta + diag(0.2*ones(Ndim-1,1),1) + diag(0.2*ones(Ndim-1,1),-1);

% Generate with the T-distribution with parameter = Nu and the
% randomized vector Mu 
Mu = rand(Ndim,1); Nu = 10;
% Theta = Theta' * Theta; % Make theta semi-positive definite
DataEC = GenerateEC(Nsim,Ndim,Nu, Mu, Theta);

% Estimate variance-covariance matrix by correlation matrix
Sigma = corr(DataEC, DataEC);
tic;lambdamin = ADMM_minlambda(Sigma, 1e-4);toc;
lambdaMin = lambdamin*1.01; 
L =60; lambdaMax =1;
lambda = lambdaMin +(lambdaMax - lambdaMin)/(L-1)*((L-1):-1:0);
% lambda = lambdaMin * (1.30:-0.01:1.00);

% Initial Setting for ADMM
err = 1e-6;
tic;Theta = ADMM(Sigma, lambda, err);toc;
Theta_hat = Truncate_Theta(Theta, 1e-4);

for i = 1:L 
    lambda(i)
    figure()
    imagesc(Theta_hat2(:,:,i));
end

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
