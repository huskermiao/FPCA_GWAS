# Uni-level FPCA with irregular time points

pca_fun = function(data,Y.vec,Xmat,bhat,K,n,BS,basis,Omega,Omega2)
{
##################################################################################
# Estiamtion of G1
##################################################################################
# residuals
e = Y.vec-Xmat%*%bhat
# initialize the response and covariates in the least squares
ee.vec = numeric()
Xmat.tensor = numeric()
start.temp = 1
for (i in 1:n)
{
n.i = length(Y[[i]])
e.i = e[(start.temp:(start.temp+n.i-1))]
start.temp = start.temp+n.i

# the response part
ee.i = kronecker(e.i,e.i)
ee.i.j = ee.i[-c((0:(n.i-1))*(n.i+1)+1)]	# remove the case j = j'
ee.vec = c(ee.vec,ee.i.j)

# the covariates part
BS.i = eval.basis(T[[i]],basis,0)
BS.i.tensor = kronecker(BS.i,BS.i)
BS.i.tensor.j = BS.i.tensor[-c((0:(n.i-1))*(n.i+1)+1),]	# remove the case j = j'
Xmat.tensor = rbind(Xmat.tensor,BS.i.tensor.j)
}
lam.G1 = tuning_nointer(-10,10,Omega2,Xmat.tensor,ee.vec)	# tunning parameter
bhat.G1 = solve(t(Xmat.tensor)%*%Xmat.tensor+lam.G1*Omega2) %*%t(Xmat.tensor)%*%ee.vec

############ variance function ############
# \hat{\sigma}_Y^2
e2.vec = e^2
lam = tuning_nointer(-10,15,Omega,Xmat,e2.vec)	# tunning parameter
bhat.e2 = solve(t(Xmat)%*%Xmat+adiag(Omega*lam))%*%t(Xmat)%*%e2.vec


#########################################################################
# Estiamtion of eigen-values and eigen-functions
#########################################################################
Jmat=inprod(basis,basis)
Jmat.sqrt=chol(Jmat)

### PCA for level 1 ###
A1 = matrix(bhat.G1,K,K)
fpca1 = eigen(Jmat.sqrt%*%A1%*%t(Jmat.sqrt))
v1 = fpca1$values
V1 = solve(Jmat.sqrt)%*%fpca1$vectors
# estimation of eigen-functions corresponding eigen-values, evaluated at tt (50 points)

# the estimated covariance functions
BS.tensor = kronecker(BS,BS)
K.hat.vec = BS.tensor%*%bhat.G1
K.hat = matrix(K.hat.vec,J,J)

# estimates for \hat{\sigma}_e^2
# If J is too small like 20 or 50, the estimate of sigma_e^2 will be negative
sigma.e2.hat = mean(BS%*%bhat.e2-diag(K.hat))


##################### return a list of results ##################
final = list(v1,V1,sigma.e2.hat)
return(final)
}