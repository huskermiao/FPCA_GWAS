## install the dependent libraries
#install.packages('fda')
#install.packages('magic')

library(fda)
library(magic)

### read data ###
data = readRDS('data.Rdata')
# read Y and T, both are lists
Y = data[[1]]	# a list of observations/measurements
T = data[[2]]	# a list of time points

# load R functions
source("pca_fun.R")
source("pca_score.R")
source("tuning_nointer.R")

# range of the data
Y.vec = Reduce(c,Y)	# vectorize the observation values
T.vec = Reduce(c,T)	# vectorize the time points
ylimit = range(Y.vec)	# range of observation values
tlimit = range(T.vec)	# range of measurement time points
n = length(Y)		# sample size n

### draw a picture for the raw data ###
# empirical mean curve
t.all = unique(sort(T.vec))
mu.emp = rep(0,length(t.all))
for (i in 1:length(t.all))
{
	mu.emp[i] = mean(Y.vec[T.vec==t.all[i]])
}

# raw data with the empirical mean function
par(mar = c(4.2, 4.2, 0.1, 0.1))
plot(tlimit,ylimit,type="n",xlab="Day",ylab="Measurements",cex.axis=1.5,cex.lab=1.5)
for (i in 1:n)
{
	lines(T[[i]],Y[[i]],lty=1,lwd=2)
}
lines(t.all,mu.emp,col="green",lwd=2)


############################################################################################
# Stage 0: Setup before the analysis
############################################################################################
tmin = tlimit[1]	# minimum observation time
tmax = tlimit[2]	# maximum observation time

### basis functions ###
K.int = 12	# K.int is the number of interior knots
order = 4 	# cubic splines
knots = tmin + (tmax-tmin)*(1:K.int)/(1+K.int)	# time points of interios knots
K = length(knots) + order 				# number of basis functions
basis = create.bspline.basis(tlimit,K,norder=order)

### penalty matrix ###
Omega = inprod(basis,basis,2,2)			# penalty matrix
inte1 = kronecker(inprod(basis,basis,2,2),inprod(basis,basis))
inte2 = kronecker(inprod(basis,basis,1,1),inprod(basis,basis,1,1))
inte3 = kronecker(inprod(basis,basis),inprod(basis,basis,2,2))
Omega2 = inte1+2*inte2+inte3


############################################################################################
# Stage 1: Estimation of the mean function
############################################################################################
################################ Mean function estimation ##################################
# design matrix
N = length(Y.vec)
Xmat = matrix(0,N,K)
start.temp = 1
for (i in 1:n)
{
	n.i = length(T[[i]])
	Xmat[(start.temp:(start.temp+n.i-1)),] = eval.basis(T[[i]],basis,0)
	start.temp = start.temp+n.i
}

### Penalized least squares estimates ###
lam = tuning_nointer(-10,15,Omega,Xmat,Y.vec)	# tunning parameter
# print(log(lam))
bhat = solve(t(Xmat)%*%Xmat+adiag(Omega*lam))%*%t(Xmat)%*%Y.vec

### draw a picture for the raw data with the estimated mean function ###
# J = tmax-tmin+1
# tt = tmin:tmax
J = 500
tt = seq(tmin, tmax, length.out=J)	# evaluation time points
BS = eval.basis(tt,basis,0)		# evaluation of the basis functions over the evaluation time points
mu = BS%*%bhat	# estimate mean function

par(mar = c(4.2, 4.2, 0.1, 0.1))
plot(tlimit,ylimit,type="n",xlab="Day",ylab="Measurements",cex.axis=1.5,cex.lab=1.5)
for (i in 1:n)
{
	lines(T[[i]],Y[[i]],lty=1,lwd=2)
}
lines(tt,mu,col="red",lwd=2)
lines(t.all,mu.emp,col="green",lwd=2)
legend("topleft",col=c("red","green"),legend=c("Estimated","Empirical"),cex=1.5,lwd=2)
abline(v = 0, lwd=2, lty=2)

# output the empirial and estimated mean curves
write.table(data.frame(t.all,mu.emp), file = "empirical.csv",row.names=FALSE,sep=",")
write.table(data.frame(tt,mu), file = "estimated.csv",row.names=FALSE,sep=",")


############################################################################################
# Stage 2: Estimation of eigenvalues and eigenfunctions by FPCA
############################################################################################
temp = pca_fun(data,Y.vec,Xmat,bhat,K,n,BS,basis,Omega,Omega2)
v1 = temp[[1]]
V1 = temp[[2]]
sigma.e2.hat = temp[[3]]

# The eigen-function matrix
Phi1 = list(n)
for (i in 1:n)
{
	Phi1[[i]] = eval.basis(T[[i]],basis,0)%*%V1
}

# Choose the number of PCs by PVE method
prop = c(0.85,0.9,0.95,0.99)
v1.trim = v1[v1>0]
n1 = length(v1.trim)
s1 = sum(v1.trim)
m = length(prop)
K.prop = rep(0,m)
for (i in 1:m)
{
	for (j in 1:n1) 
	{
		if((sum(v1.trim[1:j])/s1)>prop[i])
		{
			K.prop[i] = j
			break
		}
	}
}
K.prop
K.est = 2
sum(v1.trim[1:K.est])/s1	# first two 97.42%

### draw a picture for the eigenfunctions ###
phi.fun = BS%*%V1
par(mfrow=c(2,2))
par(mar = c(4.2, 4.4, 1.2, 0.2))
plot(tt,phi.fun[,1],type='l',main=expression(psi[1]('80.44%')),xlab='Day',ylab=expression(psi[1]*(t)),
ylim=range(phi.fun[,1]),lwd=2,cex.axis=1.5,cex.lab=1.5)
plot(tt,phi.fun[,2],type='l',main=expression(psi[2]('16.98%')),xlab='Day',ylab=expression(psi[2]*(t)),
ylim=range(phi.fun[,2]),lwd=2,cex.axis=1.5,cex.lab=1.5)
plot(tt,phi.fun[,3],type='l',main=expression(psi[3]('2.07%')),xlab='Day',ylab=expression(psi[3]*(t)),
ylim=range(phi.fun[,3]),lwd=2,cex.axis=1.5,cex.lab=1.5)
plot(tt,phi.fun[,4],type='l',main=expression(psi[4]('0.43%')),xlab='Day',ylab=expression(psi[4]*(t)),
ylim=range(phi.fun[,4]),lwd=2,cex.axis=1.5,cex.lab=1.5)

# output for the esimated eigen-functions
write.table(data.frame(tt,phi.fun[,1:4]), file = "eigen_fun.csv",row.names=FALSE,sep=",")


############################################################################################
# Stage 3: Prediction of PC scores by BLUP
############################################################################################
### Principal component scores ###
PC = pca_score(data,Xmat,bhat,K.est,v1,Phi1,n,sigma.e2.hat)

# save the PC scores
write.table(PC, file = "pc_scores.csv",row.names=FALSE,sep=",")
