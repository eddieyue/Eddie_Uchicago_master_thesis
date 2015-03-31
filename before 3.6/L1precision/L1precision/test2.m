%% Experiment of ROC 

Nsim = 1000; Ndim = 50;
Sigma = eye(Ndim);
Sigma = Sigma + diag(0.2*ones(Ndim-1,1),1) + diag(0.2*ones(Ndim-1,1),-1);
Theta = (Sigma\eye(Ndim));
X = mvnrnd(zeros(Ndim,1),Sigma,Nsim);

tau = sin(pi/2*corr(X,X,'type','Kendall'));
Theta_hat = L1precisionBCD(tau,0.01);
[TPR, FPR] = ROC_Data(Theta, Theta_hat, 1e-3);


