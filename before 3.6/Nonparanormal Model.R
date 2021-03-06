#############################
### nonparanormal model
#############################


n=1000
Z1=rnorm(n)
Z2=.3*Z1+sqrt(1-.3^2)*rnorm(n)
plot(Z1,Z2,pch=20)



f1=function(z){exp(z)}
f2=function(z){z^3}


X1=f1(Z1)
X2=f2(Z2)


plot(X1,X2,pch=20) # doesn't look multivariate normal


hist(X1) # doesn't look normal
hist(X2) # doesn't look normal


## tranform X1 & X2 back to normal
quant_normal=qnorm(1:n/(n+1))
hist(quant_normal)


Z1hat=quant_normal[rank(X1)]
Z2hat=quant_normal[rank(X2)]


# sanity check
plot(X1,Z1hat,pch=20)
hist(Z1hat)


# check if data is actually nonparanormal
plot(Z1hat,Z2hat,pch=20)















#############################
### does not fit nonparanormal model
#############################



XX1=rexp(n)
XX2=XX1+rt(n,2)





plot(XX1,XX2,pch=20) # doesn't look multivariate normal


hist(XX1) # doesn't look normal
hist(XX2) # doesn't look normal



ZZ1hat=quant_normal[rank(XX1)]
ZZ2hat=quant_normal[rank(XX2)]


# sanity check
plot(XX1,ZZ1hat,pch=20)
hist(ZZ1hat)


# check if data is actually nonparanormal
plot(ZZ1hat,ZZ2hat,pch=20)