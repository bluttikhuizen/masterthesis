coupling_fixed <- function(data, n, m, k, incidence_MCMC, iterations, step_save, a_1, T_0, a_2, V){
  # We initialize a list to store the sampled incidence matrices.
  MCMC <- list()
  # We fix a number of hyperparameters.
  v <- 1
  mu <- numeric(n)
  dof <- a_2
  prec <- V
  step_counter <- 1
  # We start our MCMC sampling scheme.
  for (i in 1:iterations) {
    # We run single MCMC iterations here, because we need to update T_0 after each iteration.
    MCMC_iteration <- list()
    # T_0 is scaled like discussed in the paper.
    scaling <- a_1
    for (j in 1:k) {
      MCMC_iteration[[j]] <- strMCMC(data[[j]],incidence_MCMC[[j]],1,1,scaling*T_0, a=a_1)
    }
    # We now sample a precision matrix from the posterior for each data set
    T_m <- list()
    W <- list()
    for (j in 1:k) {
      T_m[[j]] <- scaling*T_0 + (m-1)* cov(t(data[[j]])) + ((v*m)/(v+m))* (mu - rowMeans(data[[j]]))%*%t(mu - rowMeans(data[[j]]))
      W[[j]] <- rWishart(1, a_1+m, ginv(T_m[[j]]))[,,1]
    }
    
    # The next step is to make the precision matrices graph specific.
    covariance <- list()
    covariance_top <- list()
    incidence_top <- list()
    for (j in 1:k) {
      covariance[[j]] <- ginv(W[[j]])
      incidence_MCMC[[j]] <- MCMC_iteration[[j]][[1]][[2]]
      
      # The topological order can be computed from the incidence matrix.
      top <- top_order(incidence_MCMC[[j]])
      
      # We now reorder the relevant matrices according to the topological order
      covariance_top[[j]] <- covariance[[j]][top,top]
      incidence_top[[j]] <- incidence_MCMC[[j]][top,top]
    }
    
    # We now run a loop that computes all the conditional variances and network specific precision matrices.
    conditionals <- list()
    W_network <- list()
    for (j in 1:k) {
      conditionals[[j]] <- array(data = NA, dim = n)
      W_network[[j]] <- array(data=NA,dim = c(n,n))
      for (c in 0:(n-2)) {
        conditionals[[j]][c+1] <- conddist(n-c, covariance_top[[j]], incidence_top[[j]])
      }
      # The final element of the conditional distributions is not conditioned on anything.
      conditionals[[j]][n] <- covariance_top[[j]][1,1]
      
      # We now reverse the order of the conditional variances to make it more intuitive.
      conditionals[[j]] <- rev(conditionals[[j]])
      
      # We are now going to use these conditional distributions to find the precision matrix.
      W_network[[j]] <- precision(covariance_top[[j]],conditionals[[j]],incidence_top[[j]])
    }
    
    # We now store the incidence matrices according to the burnin and thinning.
    if(i>burnin & i%%step_save==0){
      MCMC[[step_counter]] <- incidence_MCMC
      step_counter <- step_counter+1
      print(step_counter)
    }
    
    # We can now resample T_0 for the next iteration.
    dof <- dof+k*a_1
    prec <- prec+scaling*Reduce('+', W_network)
    T_0 <- rWishart(1, dof, ginv(prec))[,,1]
  }
  return(MCMC)
}
