
source('Untitled.R')

p = 50
n = 30

Sigmmma = matrix(0.2, ncol = p, nrow = p)
diag(Sigmmma) = 1
Omega0 = solve(Sigmmma)

  Y = mvrnorm(n = n, mu=rep(0,p), Sigma = Sigmmma)
  S <- t(Y)%*%Y
  
  out <- GHS_est(S,n = n,burnin = 100,nmc = 1000); est_matrix <- apply(out[[1]],c(1,2),mean)
  
  out_09 <- f_GHS_est(S,n = n,burnin = 100,nmc = 1000, alpha = .9) ; 
  est_mat_09 <- apply(out_09[[1]],c(1,2),mean)
  out_05 <- f_GHS_est(S,n = n,burnin = 100,nmc = 1000, alpha = .5) ; 
  est_mat_05 <- apply(out_05[[1]],c(1,2),mean)
  out_01 <- f_GHS_est(S,n = n,burnin = 100,nmc = 1000, alpha = .1) ; 
  est_mat_01 <- apply(out_01[[1]],c(1,2),mean)
  
 mean((est_matrix - Omega0)^2); 
  mean((est_mat_09 - Omega0)^2)
  mean((est_mat_05 - Omega0)^2)
 mean((est_mat_01 - Omega0)^2)

