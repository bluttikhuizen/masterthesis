library(MASS)
# Function that returns the conditional covariance and mean based on the full one
conddist <- function(data, mean, covariance, incidence){
  # Determine the row/column index of the child node
  C <- nrow(covariance)
  # We select the submatrices using the parent sets
  Sigma_11 <- covariance[which(incidence[,C]==1),which(incidence[,C]==1)]
  Sigma_12 <- covariance[C,which(incidence[,C]==1)]
  Sigma_21 <- t(Sigma_12)
  Sigma_22 <- covariance[C,C]
  # We can now compute the conditional mean and covariance
  Sigma_cond <- Sigma_22+Sigma_21%*%ginv(Sigma_11)%*%Sigma_12
  Mean_cond <- mean[C]+Sigma_21%*%ginv(Sigma_11)%*%(data[which(incidence[,C]==1),
                                    dim(data)[2]]-mean[which(incidence[,C]==1)])
  return(Mean_cond,Sigma_cond)
}