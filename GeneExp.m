clear; clc;
lambda = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]; 

% load Gene Expression Dataset
GeneExpData = importdata('Processed Gene expression Data');
[Nsim, Ndim] = size(GeneExpData);

% StARS sampling
Nsam = 100; Ssam = 500;
IndS = stars_split(Nsim, Ssam, Nsam);

% for eveyr subsample, I implement three kind of method to estimate the
% precison matrix
for n = 1:size(Inds,1)
    InputData = GeneExpData(IndS(n,:),:);
    ThetaHatK = KendallsTau(InputData, lambda);
    ThetaHatS = SampleCorr(InputData, lambda);
    ThetaHatN = Nonparanormal(InputData,lambda);
end

% Measure the stability