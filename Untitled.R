library(GHS)

f_GHS_est <- function(S,n,burnin,nmc,alpha ){
  p <- nrow(S)
  omega_save <- array(0, dim=c(p,p,nmc))
  lambda_sq_save <- array(0,dim=c(p*(p-1)/2,nmc))
  tau_sq_save <- array(0,dim=c(1,nmc))
  
  ind_all = array(0,dim=c(p-1,p))
  for (i in seq(1,p,1))
  {
    if (i==1)
    {
      ind <- t(array(seq(2,p,1),dim=c(1,p-1)))
    }else if (i==p)
    {
      ind <- t(array(seq(1,p-1,1),dim=c(1,p-1)))
    }else
    {
      ind <- t(array(c(seq(1,i-1,1),seq(i+1,p,1)),dim=c(1,p-1)))
    }
    ind_all[,i] <- ind
  }
  
  # set initial values
  Omega <- diag(p)
  Sigma <- diag(p)
  tau_sq <- 1
  xi <- 1
  Lambda_sq <- array(1,dim=c(p,p))
  Nu <- array(1,dim=c(p,p))
  
  for (iter in seq(1,burnin+nmc,1))
  {
    for (i in seq(1,p,1))
    {
      ind <- ind_all[,i]
      Sigma_11 <- Sigma[ind,ind]
      sigma_12 <- Sigma[ind,i]
      sigma_22 <- Sigma[i,i]
      s_21 <- S[ind,i]*alpha
      s_22 <- S[i,i]*alpha
      lambda_sq_12 <- Lambda_sq[ind,i]
      nu_12 <- Nu[ind,i];
      
      # sample gamma and beta, gamma in the code is gamma_sample here and beta in code is beta_def here
      gamma_sample <- rgamma(1,alpha*(n/2)+1,s_22/2) # random gamma with shape=n/2+1, rate=2/s_22
      inv_Omega_11 <- Sigma_11 - (sigma_12%*%t(sigma_12))/sigma_22
      inv_C <- s_22*inv_Omega_11 + diag(1/(as.vector(lambda_sq_12)*tau_sq))
      inv_C_chol <- chol(inv_C)
      mu_i <- matrix(solve(-inv_C,s_21,tol = 1e-25))
      beta_def <- mu_i+ matrix(solve(inv_C_chol,rnorm(p-1,0,1)))
      omega_12 <- beta_def
      omega_22 <- gamma_sample + (t(beta_def)%*%inv_Omega_11%*%beta_def)
      rate <- omega_12^2/(2*tau_sq) + 1/nu_12
      lambda_sq_12 <- matrix(1/rgamma(rep(1,dim(rate)[1]),rep(1,dim(rate)[1]),rate=rate),dim(rate)[1],1)  #' random inv gamma with shape=1, rate=rate
      nu_12 <- matrix(1/rgamma(rep(1,dim(lambda_sq_12)[1]),rep(1,dim(lambda_sq_12)[1]),rate=1+1/lambda_sq_12),dim(lambda_sq_12)[1],1)#' random inv gamma with shape=1, rate=1+1/lambda_sq_12
      
      # update Omega, Sigma, Lambda_sq, Nu
      Omega[i,ind] <- omega_12
      Omega[ind,i] <- omega_12
      Omega[i,i] <- omega_22
      temp <- inv_Omega_11%*%beta_def
      Sigma_11 <- inv_Omega_11 + temp%*%t(temp)/gamma_sample
      sigma_12 <- -temp/gamma_sample
      sigma_22 <- 1/gamma_sample
      Sigma[ind,ind] <- Sigma_11
      Sigma[i,i] <- sigma_22
      Sigma[i,ind] <- sigma_12
      Sigma[ind,i] <- sigma_12
      Lambda_sq[i,ind] <- lambda_sq_12
      Lambda_sq[ind,i] <- lambda_sq_12
      Nu[i,ind] <- nu_12
      Nu[ind,i] <- nu_12
    }
    # Sample tau_sq and xi
    omega_vector <- matrix(Omega[lower.tri(Omega)])
    lambda_sq_vector <- matrix(Lambda_sq[lower.tri(Omega)])
    rate <- (1/xi) + sum((omega_vector^2)/(2*lambda_sq_vector));
    tau_sq = 1/rgamma(1,((p*(p-1)/2)+1)/2,rate)
    xi = 1/rgamma(1,1,1+1/tau_sq);    #' inv gamma w/ shape=1, rate=1+1/tau_sq
    
    # save Omega, lambda_sq, tau_sq
    if (iter >burnin)
    {
      omega_save[,,iter-burnin] <- Omega
      lambda_sq_save[,iter-burnin] <- lambda_sq_vector
      tau_sq_save[iter-burnin] <- tau_sq
    }
  }
  output <- list(omega_save,lambda_sq_save,tau_sq_save)
  return(output)
}
Posdef <- function (n,ev){
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}


#library(GGMselect)
#p = 60
#eta = 0.5 # sparsity
#Gr <- simulateGraph(p,eta)
# plot the graph
#library(network)
#par(mfrow=c(1,1))
#plot(network(Gr$G),jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)

