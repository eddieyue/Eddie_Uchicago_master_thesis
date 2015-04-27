clear; clc;

% load Gene Expression Dataset
genedata = get_arabidopsis_gene_data();
% [genedata, genenames, genepathways, is_isoprenoid] = ...
%     get_arabidopsis_gene_data();
genedata = genedata';

% Select the Feature with one pathway
% unipathway = unique(genepathways);
% label = unipathway(36); 
% In this case, pathway is Inositolphosphatemetabolism
% genesubdata = PathWay_subset(label, gecd nepathways, genedata);
genesubdata = genedata(:,1:39);
[Nsim, Ndim] = size(genesubdata);

% test part
X = genesubdata; 
Tau_test = sin(pi/2*kendalltau_fast(X));
tic;lambdamin_test = ADMM_minlambda(Tau_test);toc;
sprintf('Minimum lambda in Kendalls Tau is found ')
L_test = 500;
lambda_test = lambda_range(lambdamin_test, L_test);
lambda_test2 = lambda_test(751:800);

start_test = tic;Thetahat_test = ADMM(Tau_test,lambda_test,1);elapsed_test = toc(start_test);
% start_test = tic;Thetahat_test = ADMM(Tau_test,lambda_test,1,'quiet');elapsed_test = toc(start_test);
sprintf('Whole Process in Kendalls Tau is completed')
% end of test part

% StARS sampling
% set Number of sample and Siza of sample
Nsam = 20; Ssam = 20;
% Inds = stars_split(Nsim, Ssam, Nsam);
Inds = bootstrip_split(Nsim, Ssam, Nsam);

% for eveyr subsample, I implement three kind of method to estimate the
% precison matrix
L = 20;  N = size(Inds,1);

% Start for ADMM
ThetaK = zeros(Ndim, Ndim, L, N);
ThetaS = zeros(Ndim, Ndim, L, N);
ThetaN = zeros(Ndim, Ndim, L, N);
for n = 1:N
    SampleIndex = Inds(n,:);
    InputData = genesubdata(SampleIndex,:);
    sprintf('Working on the %d th subsample for Kendalls Tau', n)
    ThetaK(:,:,:,n) = KendallsTau(InputData, L);
    sprintf('Working on the %d th subsample for Sample Correlation', n)
    ThetaS(:,:,:,n) = SampleCorr(InputData, L);
    sprintf('Working on the %d th subsample for Nonparanormal', n)
    ThetaN(:,:,:,n) = Nonparanormal(InputData,L);
end

eta =1e-6;
ThetaHatK = Truncate_Theta(ThetaK,eta);
ThetaHatS = Truncate_Theta(ThetaS,eta);
ThetaHatN = Truncate_Theta(ThetaN,eta);

% for each lambda
PowerHatK = zeros(L,1);
PowerHatS = zeros(L,1);
PowerHatN = zeros(L,1);
InstabilityHatK = zeros(L,1);
InstabilityHatS = zeros(L,1);
InstabilityHatN = zeros(L,1);
for j = 1:L
    PsiK = zeros(Ndim);PsiS = zeros(Ndim);PsiN = zeros(Ndim);
    % for each subsample
    for n = 1:N
        PsiK  = PsiK + (ThetaHatK(:,:,j,n) ~= 0)/N; 
        PsiS  = PsiS + (ThetaHatS(:,:,j,n) ~= 0)/N; 
        PsiN  = PsiN + (ThetaHatN(:,:,j,n) ~= 0)/N;
    end
    PowerHatK(j) = MyPower(PsiK);
    PowerHatS(j) = MyPower(PsiS);
    PowerHatN(j) = MyPower(PsiN);
    InstabilityHatK(j) = MyStablility(PsiK);
    InstabilityHatS(j) = MyStablility(PsiS);
    InstabilityHatN(j) = MyStablility(PsiN);
end

figure()
hold on
plot(PowerHatK, InstabilityHatK,'-ro')
plot(PowerHatS, InstabilityHatS, '-bo')
plot(PowerHatN, InstabilityHatN,'-go')
hold off
xlabel('Power');
ylabel('Instability')
legend('Kendalls Tau','Sample Correlation','Nonparanormal')


% PowerHatK = zeros(length(lambda),1);
% PowerHatS = zeros(length(lambda),1);
% PowerHatN = zeros(length(lambda),1);
% StabilityHatK = zeros(length(lambda),1);
% StabilityHatS = zeros(length(lambda),1);
% StabilityHatN = zeros(length(lambda),1);
% for i = 1:length(lambda)
%     r = lambda(i);
%     PsiK = zeros(Ndim);PsiS = zeros(Ndim);PsiN = zeros(Ndim);
%     for n = 1:size(Inds,1)
%         SampleIndex = Inds(n,:);
%         InputData = genesubdata(SampleIndex,:); % correction needed
%         ThetaHatK = KendallsTau(InputData, r);
%         ThetaHatS = SampleCorr(InputData, r);
%         ThetaHatN = Nonparanormal(InputData,r);
%         
%         % 
%         PsiK  = PsiK + (ThetaHatK ~= 0)/size(Inds,1); 
%         PsiS  = PsiS + (ThetaHatS ~= 0)/size(Inds,1); 
%         PsiN  = PsiN + (ThetaHatN ~= 0)/size(Inds,1); 
%     end
%     PowerHatK(i) = MyPower(PsiK);
%     PowerHatS(i) = MyPower(PsiS);
%     PowerHatN(i) = MyPower(PsiN);
%     StabilityHatK(i) = MyStablility(PsiK);
%     StabilityHatS(i) = MyStablility(PsiS);
%     StabilityHatN(i) = MyStablility(PsiN);
% end
% Measure the stability