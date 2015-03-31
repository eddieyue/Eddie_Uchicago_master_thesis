%% Repeat for elliptical data 
clear
clc
% R = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]; 
lambda = 1;
N = 10; Nu = 1; eta = 1e-4; plotkey = 1;
% TPR = zeros(7, N);FPR = zeros(7, N);
Ndim = 10;
Nsim = 100;

Theta = eye(Ndim);
Theta = Theta + diag(0.2*ones(Ndim-1,1),1) + diag(0.2*ones(Ndim-1,1),-1);
Mu = rand(Ndim,1);
% Theta = Theta' * Theta; % Make theta semi-positive definite
DataEC = GenerateEC(Nsim,Ndim,Nu, Mu, Theta);
Sigma = corr(DataEC, DataEC);
% Thetahat = ADMM(Tau,lambda);

n = size(Sigma,1);
err = 1e-4; alpha = eigs(Sigma, 1);

S =  (alpha^2 * eye(n) - Sigma^2);R = Sigma/alpha; L = ones(n,1); 
A = chol(S);
Theta_hat = zeros(n);
E = eye(n);

soft_thresh = @(b, rambda) sign(b).*max(abs(b) - 1/rambda,0);
% for i = 1:n
i = 1;
    ei = E(:,i);
    Ri = R(:,i);
    xk = L; uk = L; yk = 0; 
    mu = 0.01; 
    cond = true;
    while cond 
        xk
        yk1 = min(max(Sigma*xk -ei + uk/mu, -lambda*L),lambda*L);
        tempb = R * (ei + yk1 - uk/mu) + A * xk/alpha;
        xk1 = soft_thresh(tempb, mu*alpha);
        uk1 = uk + mu * (Sigma*xk1 - ei - yk1);
%         cond = norm( xk1 - xk) > err && norm( yk1 - yk) > err && norm(uk1 - uk) > err;
        cond = norm( xk1 - xk) > err;
        xk = xk1; yk = yk1; uk = uk1;
    end
    Theta_hat(:,i) = xk;

% end

