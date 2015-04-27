function Y = T_Delta(X, N)

% delta = 1/(N+1);
delta = 1/(4*N^0.25*(pi*log(N)));

X(X < delta) = delta;
X(X > (1-delta)) = 1-delta;

Y = X;
end