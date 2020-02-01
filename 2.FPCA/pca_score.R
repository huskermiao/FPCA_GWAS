# Principal component scores prediction

pca_score = function(data,Xmat,bhat,K.est,v1,Phi1,n,sigma.e2.hat)
{
#########################################################################
# Prediction of the principle component scores
#########################################################################
K.pc = K.est	# number of PCs
PC = as.data.frame(matrix(0,n,K.pc))	# initialization the PC matrix

# The eigenvalues matrix: Lambda
v1.pc = v1[1:K.pc]
G = diag(rep(v1.pc,n))

# The eigenfunction matrix: Phi
Z = matrix(0,N,n*K.pc)
start.temp = 1
for (i in 1:n)
{
	n.i = length(T[[i]])
	Z[(start.temp:(start.temp+n.i-1)),(((i-1)*K.pc+1):(i*K.pc))] = Phi1[[i]][,1:K.pc]
	start.temp = start.temp+n.i
}

# The matrix: Sigma^{-1}
S.inv = (1/sigma.e2.hat)*diag(N) - (sigma.e2.hat)^(-2)*Z%*%solve(solve(G)+t(Z)%*%Z/sigma.e2.hat,t(Z))

# BLUP
e = Y.vec-Xmat%*%bhat
xi = G%*% t(Z)%*%(S.inv%*%e)
PC = t(matrix(xi,K.pc,n))

##################### return the results  ##################
return(PC)
}
