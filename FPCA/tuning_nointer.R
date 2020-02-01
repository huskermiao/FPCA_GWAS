# Choose the optimal tuning parameter when there is no fixed effects

tuning_nointer = function(lower, upper, Omega, Xmat, Y.vec)
{
lam.list=exp(seq(lower,upper,1))
gcv=rep(0,length(lam.list))
for(ii in 1:length(lam.list))
{
A <- solve(t(Xmat)%*%Xmat+adiag(Omega*lam.list[ii]))
Y.vec.hat <- (Xmat%*%A) %*% (t(Xmat)%*%Y.vec)
diag.mean <- sum(diag(t(Xmat)%*%Xmat%*%A))/(dim(Xmat)[1]) # the original mean(diag(Hmat))
gcv[ii] <- mean((Y.vec-Y.vec.hat)^2)/(1-diag.mean)^2
}
ind=which(gcv==min(gcv))
lam.list[ind]
}