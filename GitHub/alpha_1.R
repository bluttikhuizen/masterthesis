# We fix the interaction coefficients.
edges <- 20
a1 <- runif(edges,0.5,2)
a2 <- sample(c(-1,1),edges, replace=TRUE)
a <- a1*a2    # vector with regression coefficients for the 20 edges

# We first sample a large data set to derive a representative T_0.
m <- 100
n <- 11
data_1 <- make_test_Data_2(a,10^6,1)

# We now create the incidence matrix suggested for the RAF signaling pathway in literature.
incidence <- make_true_Net()

# Computation of the empirical covariance matrix.
covariance <- cov(t(data_1))

# The topological order can be computed from the incidence matrix.
top <- top_order(incidence)

# Reorder the relevant matrices according to the topological order.
covariance_top <- covariance[top,top]
incidence_top <- incidence[top,top]

# Run a loop that computes all the conditional variances.
conditionals <- c()
for (k in 0:(n-2)) {
  conditionals[k+1] <- conddist(n-k, covariance_top, incidence_top)
}
# The final element of the conditional distributions is equal to the variance of the first node.
conditionals[n] <- covariance_top[1,1]

# We now reverse the order of the conditional variances to make it more intuitive.
conditionals <- rev(conditionals)

# We are now going to use these conditional distributions to make W graph specific.
W <- precision(covariance_top,conditionals,incidence_top)

# We use the fact that E[W]=a_1*T_0^(-1)
T_0 <- ginv(W)
T_0 <- diag(n)

# We initialize a_1 equal to the minimal allowed value of n+2.
a_1 <- n+2
v <- 1

# For the optimization of a_1, we need a second data set.
data_2 <- make_test_Data_2(a,m,1)

# We now compute the BGE scores
mu <- numeric(n)
a_stored <- c()
bge_old_stored <- c()
bge_prop_stored <- c()
sampleinterval <- 1
# We now run a MCMC sampling scheme.
iterations <- 100000
for (i in 1:iterations) {
  # We now propose a new a_1 with a uniform step centered around 0.
  a_1_prop <- a_1+runif(1, min = -sampleinterval, max = sampleinterval)
  
  # computation of the (logarithmizid) BGe Score of the current alpha
  a <- a_1
  bge_old <- BGE_alpha(a,mu,data_2,incidence, T_0)
  bge_old_stored[i] <- bge_old
  
  # computation of the (logarithmizid) BGe Score of the PROPOSED alpha
  a <- a_1_prop
  bge_prop <- BGE_alpha(a,mu,data_2,incidence, T_0)
  bge_prop_stored[i] <- bge_prop
  
  # We can now compute the acceptance ratio
  acceptance <- min(1, exp(bge_prop - bge_old))
  rand <- runif(1)
  
  if(acceptance > rand & a_1_prop>n+2){
    a_1<-a_1_prop
  }
  a_stored[i]<-a_1
}
plot(a_stored, main= "Markov Trajectory of a_1",xlab="Iterations",ylab="a_1")
