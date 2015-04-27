%% Experiment of ROC 
clear 
clc
addpath('L1precision')
Nsim = 100; Ndim = 50;
Mu = zeros(Ndim,1);
Nu = 1;  %% heavy tailed
Theta = eye(Ndim);
Theta = Theta + diag(0.2*ones(Ndim-1,1),1) + diag(0.2*ones(Ndim-1,1),-1);
Sigma = Theta\eye(Ndim);
A = chol(Sigma,'lower');
DataGauss = mvnrnd(Mu,Sigma,Nsim);
DataEC = GenerateEC(Nsim,Ndim,Nu,Mu, A);

%%




TauGauss = sin(pi/2*corr(DataGauss,DataGauss,'type','Kendall'));
TauEC = sin(pi/2*corr(DataEC,DataEC,'type','Kendall'));
RhoGauss = corr(DataGauss,DataGauss);
RhoEC = corr(DataEC,DataEC);
R = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]; % do this over a range of penalty parameters (replace 0.01)



% TPRGauss = zeros(7,1);
% TPRSGauss = zeros(7,1);
% FPRGauss = zeros(7,1); %% false positives
% FPRSGauss = zeros(7,1);
% TPREC = zeros(7,1);
% TPRSEC = zeros(7,1);
% FPREC = zeros(7,1);
% FPRSEC = zeros(7,1);

TPRCGauss = zeros(7,1);
TPRSCGauss = zeros(7,1);
FPRCGauss = zeros(7,1); 
FPRSCGauss = zeros(7,1);
TPRCEC = zeros(7,1);
TPRSCEC = zeros(7,1);
FPRCEC = zeros(7,1);
FPRSCEC = zeros(7,1);


n = 1 
for r = R
%     tic;
% 	Theta_hatGuass = L1precisionBCD(TauGauss,r);
% 	Theta_hatEC = L1precisionBCD(TauEC,r);
% 	Theta_hatSGuass = L1precisionBCD(RhoGauss,r);
% 	Theta_hatSEC = L1precisionBCD(RhoEC,r);
%     t1= toc
% %     fprintf('time for Grahical Lasso %d', t1);
%     tic;
    Theta_hatCGuass = CLIME(TauGauss,r);
	Theta_hatCEC = CLIME(TauEC,r);
	Theta_hatSCGuass = CLIME(RhoGauss,r);
	Theta_hatSCEC = CLIME(RhoEC,r);
%     t2= toc
% %     fprintf('time for CLIME %d', t2);
% 	[TPRGauss(n), FPRGauss(n)] = ROCData(Theta, Theta_hatGuass, 1e-4);
% 	[TPREC(n), FPREC(n)] = ROCData(Theta, Theta_hatEC, 1e-4);
% 	[TPRSGauss(n), FPRSGauss(n)] = ROCData(Theta, Theta_hatSGuass, 1e-4);
% 	[TPRSEC(n), FPRSEC(n)] = ROCData(Theta, Theta_hatSEC, 1e-4);
	[TPRCGauss(n), FPRCGauss(n)] = ROCData(Theta, Theta_hatCGuass, 1e-4);
	[TPRCEC(n), FPRCEC(n)] = ROCData(Theta, Theta_hatCEC, 1e-4);
	[TPRSCGauss(n), FPRSCGauss(n)] = ROCData(Theta, Theta_hatSCGuass, 1e-4);
    [TPRSCEC(n), FPRSCEC(n)] = ROCData(Theta, Theta_hatSCEC, 1e-4);
    n = n +1
end	

%   & plot curve

% compare to GLasso using sample correlation rather than tau

% figure()
% subplot(2,2,1)
% plot(FPRSGauss,TPRSGauss)
% title('sample correlation gaussian')
% xlabel('FPR')
% ylabel('TPR')
% subplot(2,2,2)
% plot(FPRGauss,TPRGauss)
% title('Tau elliptical gaussian')
% xlabel('FPR')
% ylabel('TPR')
% subplot(2,2,3)
% plot(FPRSEC,TPRSEC)
% title('sample correlation elliptical')
% xlabel('FPR')
% ylabel('TPR')
% subplot(2,2,4)
% plot(FPREC,TPREC)
% title('tau elliptical')
% xlabel('FPR')
% ylabel('TPR')

figure()
subplot(2,2,1)
plot(FPRSCGauss,TPRSCGauss)
title('sample correlation gaussian')
xlabel('FPR')
ylabel('TPR')
subplot(2,2,2)
plot(FPRCGauss,TPRCGauss)
title('Tau elliptical gaussian')
xlabel('FPR')
ylabel('TPR')
subplot(2,2,3)
plot(FPRSCEC,TPRSCEC)
title('sample correlation elliptical')
xlabel('FPR')
ylabel('TPR')
subplot(2,2,4)
plot(FPRCEC,TPRCEC)
title('tau elliptical')
xlabel('FPR')
ylabel('TPR')


%% Repeat for nongaussian data elliptical
R = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]; 
N = 10;
TPR = zeros(7, N);
FPR = zeros(7, N);

for n = 1:N
    i = 1;
    for r = R
        Theta = rand(Ndim);
        Theta = Theta * Theta % Make theta semi-positive definite
        A = chol(Sigma,'lower');
        DataEC = GenerateEC(Nsim,Ndim,Nu,Mu, A);
        TauEC = sin(pi/2*corr(DataEC,DataEC,'type','Kendall'));
        Theta_hat= CLIME(TauEC,r);
        [TPR(i,n), FPR(i,n)] = ROCData(Theta, Theta_hat, 1e-4);
        i = i + 1; 
    end
end




