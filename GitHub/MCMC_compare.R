# We first sample the data needed for the different structure MCMC algorithms.
k <- 5    # We set the number of datasets.
m <- 15  # We set the sample size for the datasets.
n <- 11   # We set the number of nodes in the graph.

# We also fix the regression coefficients, so that they remain the same over the data sets.
edges <- 20
a1 <- runif(edges,0.5,2)
a2 <- sample(c(-1,1),edges, replace=TRUE)
a <- a1*a2    # vector with regression coefficients for the 20 edges
data <- list()
for (j in 1:k) {
  data[[j]] <- make_test_Data_2(a,m,1)
}

# We also initialize an empty incidence matrix.
incidence <- matrix(data=numeric(n*n), nrow = n, ncol = n)
incidence_MCMC <- list()
for (j in 1:k) {
  incidence_MCMC[[j]] <- incidence
}

# Here we set the number of iterations, step save and burnin period.
iterations <- 100000
step_save <- 1000
burnin <- 10000

# We sample an initial T_0.
a_2 <- n+2   # Set the degrees of freedom for the prior distribution of T_0.
V <- diag(a_2,n) # Set the precision matrix for the prior distribution.
T_0 <- rWishart(1,a_2, ginv(V))[,,1] # Sample T_0 from the prior.
a_1 <- a_2 # We initialize a_1 at the minimal value.

# We call the functions that run the different sampling schemes.
update <- coupling_1_update(data, n, m, k, incidence_MCMC, iterations, step_save, a_1, T_0, a_2, V)
fixed <- coupling_fixed(data, n, m, k, incidence_MCMC, iterations, step_save, a_1, T_0, a_2, V)
structureMarco <- list()
for (i in 1:5) {
  structureMarco[[i]] <- strMCMC(Data=data[[i]],incidence_MCMC[[i]],iterations,step_save, T_0, a=a_1)
}


# We extract the ROC curves for the k data sets.
CPDAG_update <- list()
postP_update <- list()
CPDAG_fixed <- list()
postP_fixed <- list()
CPDAG_structure <- list()
postP_structure <- list()

# We now create the incidence matrix suggested for the RAF signaling pathway in literature.
true_incidence <- make_true_Net()


# The CPDAG is calculated for the true DAG
trueEdges <- extract_cpdag_of_dag(true_incidence)
for (j in 1:k) {
  CPDAG_update[[j]] <- cpdag_list(sapply(update[[1]], `[`, j),1)
  CPDAG_fixed[[j]] <- cpdag_list(sapply(fixed, `[`, j),1)
  CPDAG_structure[[j]] <- cpdag_list(structureMarco[[j]][[1]],250)
  postP_update[[j]] <- CPDAG_update[[j]][[3]]
  postP_fixed[[j]] <- CPDAG_fixed[[j]][[3]]
  postP_structure[[j]] <- CPDAG_structure[[j]][[3]]
  draw_ROC(postP_update[[j]],trueEdges)
  compute_AUROC(postP_update[[j]], trueEdges)
  draw_ROC(postP_fixed[[j]],trueEdges)
  compute_AUROC(postP_fixed[[j]], trueEdges)
  draw_ROC(postP_structure[[j]],trueEdges)
  compute_AUROC(postP_structure[[j]], trueEdges)
}

auroc_update <- c()
auroc_fixed <- c()
auroc_structure <- c()

for (r in 1:5) {
  auroc_update[r]<- compute_AUROC(postP_update[[r]], trueEdges)
  auroc_fixed[r]<- compute_AUROC(postP_fixed[[r]], trueEdges)
  auroc_structure[r]<- compute_AUROC(postP_structure[[r]], trueEdges)
}

plot(update[[2]], main="Scaling of T_0", xlab="Iterations", ylab="Scaling")