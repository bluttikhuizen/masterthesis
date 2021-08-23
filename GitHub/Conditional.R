library(MASS)
# Function that returns the conditional covariance and mean based on the full one
conddist <- function(child_index, covariance, incidence){
  # Determine the row/column index of the child node and those of the parents
  C <- child_index
  P <- which(incidence[,C]==2)
  # We select the submatrices using the parent sets
  Sigma_11 <- covariance[P,P]
  Sigma_12 <- covariance[P,C]
  Sigma_21 <- t(Sigma_12)
  Sigma_22 <- covariance[C,C]
  
  # We can now compute the conditional variance
  if (length(P)>0) {
    Sigma_cond <- Sigma_22-Sigma_21%*%ginv(Sigma_11)%*%Sigma_12}
    else{
      Sigma_cond <- Sigma_22
      }
  return(Sigma_cond)
}