### basic setup ###
rm(list=ls())	# remove ALL objects
library(fda)
library(magic)

# read data
data = read.csv("example_unfit.csv") 
head(data)

# choose the first genotype
subdata = data[data$Sample=="284_CM003",]
subdata

K.int = 6	# number of interior knots
tt = subdata$day	# evaluation points
J = length(tt)		# number of evaluation points
t.min = min(tt)
t.max = max(tt)

### basis functions ###
order = 4 # cubic splines, for second derivative, using order 6 rather than 4
knots = t.min + (t.max-t.min)*(1:K.int)/(1+K.int)
K = length(knots) + order # number of basis functions
basis = create.bspline.basis(c(t.min,t.max),K,norder=order)

# evaluate original, the first, and second derivatives of basis functions
BS = eval.basis(tt,basis,0)
BS1 = eval.basis(tt,basis,1)
BS2 = eval.basis(tt,basis,2)

# penalty matrix
Omega = inprod(basis,basis,2,2)	# for second derivative, using 4 rather than 2

# function for tuning parameter selection using GCV
tuning_nointer = function(lower, upper, Omega, Xmat, Y.vec, xlen)
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


################################ Mean function estimation ##################################
### Design matrix ###
Y.vec = as.vector(subdata$Height)
N = length(Y.vec)
# design matrix
Xmat = BS

### Penalized least squares estimates ###
lam = tuning_nointer(-10,15,Omega,Xmat,Y.vec,xlen)	# tunning parameter
bhat = solve(t(Xmat)%*%Xmat+adiag(Omega*lam))%*%t(Xmat)%*%Y.vec


############# Figures #############
# evaluation time points
n.eval = range(tt)[2]-range(tt)[1]	
t.eval = seq(t.min,t.max,by=(t.max-t.min)/n.eval)
mu.eval = eval.basis(t.eval,basis,0)%*%bhat
y_lim = range(Y.vec)
png(filename = 'sm_fitting.png')
plot(tt,Y.vec,type="p",ylim=c(y_lim[1]-20,y_lim[2]+25),xlab="DAP",ylab="Plant Height (CM)",main='test', cex.axis=1.5,cex.lab=1.5,col="blue")
lines(t.eval,mu.eval,lty=1,lwd=2,col="red")
dev.off()
