function R = stars_split(Nsim, b, N)

% 

V = 1:Nsim ;
[~, x] = sort(rand(N,Nsim),2) ;
x = sort(x(:,1:b),2) ;
R = V(x); % N random combinations of K unique elements from V
R = unique(R, 'rows');

end
