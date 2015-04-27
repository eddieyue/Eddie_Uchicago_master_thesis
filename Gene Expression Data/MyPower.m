function P = MyPower(theta)

N = size(theta,1);
P = sum(sum(triu(theta,1)))/((N^2-N)/2);
end