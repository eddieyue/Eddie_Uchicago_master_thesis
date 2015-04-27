function Zhat = RestoreData(X, mu, sigma)
% Input Parameter: 
% X:: Input Data

if nargin < 3, sigma =1; end
if nargin < 2, mu = 0; end

[Nsim, Ndim] = size(X); 
Zhat = zeros(Nsim, Ndim);
quant_normal = norminv((1:Nsim)/(Nsim+1), mu, sigma);
%[~, Xind] = sort(X);
for i = 1:Ndim
    [~, ~, Xind] = unique(X(:,i));
    Zhat(:,i) = quant_normal(Xind);
end 
end