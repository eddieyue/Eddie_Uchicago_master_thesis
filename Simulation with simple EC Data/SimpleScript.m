%% Repeat for elliptical data 
clear
clc
R = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]; 
N = 10; Nu = 1; eta = 1e-4; plotkey = 1;
% TPR = zeros(7, N);FPR = zeros(7, N);
Ndim = 10;
Nsim = 100;

Theta = eye(Ndim);
Theta = Theta + diag(0.2*ones(Ndim-1,1),1) + diag(0.2*ones(Ndim-1,1),-1);
Mu = rand(Ndim,1);
% Theta = Theta' * Theta; % Make theta semi-positive definite
DataEC = GenerateEC(Nsim,Ndim,Nu, Mu, Theta);
ThetaHatK = KendallsTau(DataEC,Theta,R,eta,plotkey);
ThetaHatS = SampleCorr(DataEC,Theta,R,eta,plotkey);
ThetaHatN = Nonparanormal(DataEC,Theta,R,eta,plotkey);
