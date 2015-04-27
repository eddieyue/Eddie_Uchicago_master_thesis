function D = MyInstablility(theta)

N = size(theta,1);
xi = 2 * (theta .* (1-theta));
D = sum(sum(triu(xi,1)))/((N^2-N)/2); %stability or instability??
end