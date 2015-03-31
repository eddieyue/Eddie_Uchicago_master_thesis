%% Nonparanormal Model Test for 2 and 3 vairbles 

%% For 2 Variables

f1 = @(x) exp(x);
f2 = @(x) x.^3;
TransVector = {f1, f2};

Nsim = 1000; Ndim = 2;
CorrMatrix = [1, 0, 0; sqrt(1-.3^2), 0.3, 0];
X = GenerateData(Nsim, Ndim, CorrMatrix, TransVector); 
Zhat = RestoreData(X);
figure()
% plotmatrix(Zhat)
scatter(Zhat(:,1), Zhat(:,2))

XX = randn(Nsim,Ndim)*[1,.5;.5,1];
XX = diag(abs(trnd(5,[Nsim,1])))*XX;
%XX(:,1)= exprnd(1, [Nsim, 1]);
%XX(:,2) = XX(:,1)+ trnd(2, [Nsim, 1]);
ZZhat = RestoreData(XX);
figure()
% plotmatrix(ZZhat)
scatter(ZZhat(:,1), ZZhat(:,2))
%% For 3 Variables

f1 = @(x) exp(x);
f2 = @(x) x.^3;
f3 = @(x) x.^5;
TransVector = {f1, f2, f3};

Nsim = 2000; Ndim = 3;
CorrMatrix = [1, 0, 0, 0; sqrt(1-.3^2), 0.3, 0, 0; sqrt(1-.2^3), 0.2, 0.1, 0];
X = GenerateData(Nsim, Ndim, CorrMatrix, TransVector); 
Zhat = RestoreData(X);
figure()
plotmatrix(Zhat)
% scatter(Zhat(:,1), Zhat(:,2))

XX = zeros(Nsim, 3);
XX(:,1)= exprnd(1, [Nsim, 1]);
XX(:,2) = XX(:,1) + trnd(2, [Nsim, 1]);
XX(:,3) = 0.2 * XX(:,1) + 0.7 * XX(:,2) +  trnd(4, [Nsim, 1]);
ZZhat = RestoreData(XX);
figure()
plotmatrix(ZZhat)
% scatter(ZZhat(:,1), ZZhat(:,2))