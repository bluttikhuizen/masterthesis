# m_obs: number of samples
# var_noise: variance of Gaussian distributed noise terms

make_test_Data <- function(m_obs, var_noise){
edges <- 20
a1 <- runif(edges,0.5,2)
a2 <- sample(c(-1,1),edges, replace=TRUE)
a <- a1*a2    # vector with regression coefficients for the 20 edges
  

# 1. pip3
x_pip3 <- rnorm(m_obs, sd=1)
pip3 <- (x_pip3 - mean(x_pip3))/sd(x_pip3)

# 2. plcg
x_plcg <- a[1]* pip3 + rnorm(m_obs, sd=sqrt(var_noise))
plcg <- (x_plcg - mean(x_plcg))/sd(x_plcg)

# 3. pip2
x_pip2 <- a[2]* pip3 + a[3]*plcg + rnorm(m_obs, sd=sqrt(var_noise))
pip2 <- (x_pip2 - mean(x_pip2))/sd(x_pip2)

# 4. pkc
x_pkc <- a[4]* pip2 + a[5]*plcg + rnorm(m_obs, sd=sqrt(var_noise))
pkc  <- (x_pkc - mean(x_pkc))/sd(x_pkc)

# 5. pka
x_pka <- a[6]* pkc + rnorm(m_obs, sd=sqrt(var_noise))
pka  <- (x_pka - mean(x_pka))/sd(x_pka)

# 6. jnk
x_jnk <- a[7]* pkc + a[8]* pka + rnorm(m_obs, sd=sqrt(var_noise))
jnk  <- (x_jnk - mean(x_jnk))/sd(x_jnk)

# 7. p38
x_p38 <- a[9]* pkc + a[10]* pka + rnorm(m_obs, sd=sqrt(var_noise))
p38  <- (x_p38 - mean(x_p38))/sd(x_p38)

# 8. raf
x_raf <- a[11]* pkc + a[12]* pka + rnorm(m_obs, sd=sqrt(var_noise))
raf  <- (x_raf - mean(x_raf))/sd(x_raf)

# 9. mek
x_mek <- a[13]* pkc + a[14]* pka + a[15]* raf + rnorm(m_obs, sd=sqrt(var_noise))        
mek  <- (x_mek - mean(x_mek))/sd(x_mek)

# 10. erk
x_erk <- a[16]* pka + a[17]* mek + rnorm(m_obs, sd=sqrt(var_noise))
erk  <- (x_erk - mean(x_erk))/sd(x_erk)

# 11. akt
x_akt <- a[18]* pip3 + a[19]* pka + a[20]* erk + rnorm(m_obs, sd=sqrt(var_noise))
akt  <- (x_akt - mean(x_akt))/sd(x_akt)    

daten <- cbind(pip3, plcg, pip2, pkc, pka, jnk, p38, raf, mek, erk, akt)

return(t(daten))
}

# m_obs: number of samples
# var_noise: variance of Gaussian distributed noise terms

make_test_Data_2 <- function(a,m_obs, var_noise){
  # 1. pip3
  x_pip3 <- rnorm(m_obs, sd=1)
  pip3 <- (x_pip3 - mean(x_pip3))/sd(x_pip3)
  
  # 2. plcg
  x_plcg <- a[1]* pip3 + rnorm(m_obs, sd=sqrt(var_noise))
  plcg <- (x_plcg - mean(x_plcg))/sd(x_plcg)
  
  # 3. pip2
  x_pip2 <- a[2]* pip3 + a[3]*plcg + rnorm(m_obs, sd=sqrt(var_noise))
  pip2 <- (x_pip2 - mean(x_pip2))/sd(x_pip2)
  
  # 4. pkc
  x_pkc <- a[4]* pip2 + a[5]*plcg + rnorm(m_obs, sd=sqrt(var_noise))
  pkc  <- (x_pkc - mean(x_pkc))/sd(x_pkc)
  
  # 5. pka
  x_pka <- a[6]* pkc + rnorm(m_obs, sd=sqrt(var_noise))
  pka  <- (x_pka - mean(x_pka))/sd(x_pka)
  
  # 6. jnk
  x_jnk <- a[7]* pkc + a[8]* pka + rnorm(m_obs, sd=sqrt(var_noise))
  jnk  <- (x_jnk - mean(x_jnk))/sd(x_jnk)
  
  # 7. p38
  x_p38 <- a[9]* pkc + a[10]* pka + rnorm(m_obs, sd=sqrt(var_noise))
  p38  <- (x_p38 - mean(x_p38))/sd(x_p38)
  
  # 8. raf
  x_raf <- a[11]* pkc + a[12]* pka + rnorm(m_obs, sd=sqrt(var_noise))
  raf  <- (x_raf - mean(x_raf))/sd(x_raf)
  
  # 9. mek
  x_mek <- a[13]* pkc + a[14]* pka + a[15]* raf + rnorm(m_obs, sd=sqrt(var_noise))        
  mek  <- (x_mek - mean(x_mek))/sd(x_mek)
  
  # 10. erk
  x_erk <- a[16]* pka + a[17]* mek + rnorm(m_obs, sd=sqrt(var_noise))
  erk  <- (x_erk - mean(x_erk))/sd(x_erk)
  
  # 11. akt
  x_akt <- a[18]* pip3 + a[19]* pka + a[20]* erk + rnorm(m_obs, sd=sqrt(var_noise))
  akt  <- (x_akt - mean(x_akt))/sd(x_akt)    
  
  daten <- cbind(pip3, plcg, pip2, pkc, pka, jnk, p38, raf, mek, erk, akt)
  
  return(t(daten))
}


################################################################################
##################### function for plotting the ROC curve

library(limma)
library(ROC)

draw_ROC <- function(postP, trueEdges, steps=0.01){


n1 <- numeric(0)
for(i in 1:dim(postP)[1]){
        for(j in 1:dim(postP)[2]){
                                  if (abs(i-j)>0){
                                  n1 <- c(n1, postP[i,j])
                                  }
                                  }
                          }

n2 <- numeric(0)
for(i in 1:dim(trueEdges)[1]){
       for(j in 1:dim(trueEdges)[2]){
                              if (abs(i-j)>0){     
                              n2 <- c(n2, trueEdges[i,j])
                                        }
                                  }
                          }

sp <- sort(n1)         # posterior probabilities
st <- n2[order(n1)]    # true ordered by posterior probabilities

rocc.obj <- rocdemo.sca(n2, n1, rule=NULL, cutpts=NA)#seq(0,1,steps))

return(list(plot((1-rocc.obj@"spec"), rocc.obj@"sens", type="l",las=1, xlab="1 - specificity", ylab="sensitivity", main='ROC curve'), abline(0,1, lty=2)))


}



################################################################################

make_true_Net <- function(){

NETWORK <- matrix(numeric(11*11),11,11)

# 1. pip3

# 2. plcg
NETWORK[1,2] <- 1

# 3. pip2
NETWORK[1,3] <- 1
NETWORK[2,3] <- 1

# 4. pkc
NETWORK[2,4] <- 1
NETWORK[3,4] <- 1

# 5. pka
NETWORK[4,5] <- 1

# 6. jnk
NETWORK[4,6]  <- 1
NETWORK[5,6]  <- 1

# 7. p3B
NETWORK[4,7]  <- 1
NETWORK[5,7]  <- 1

# 8. raf
NETWORK[4,8]  <- 1
NETWORK[5,8]  <- 1

# 9. mek
NETWORK[4,9]  <- 1
NETWORK[5,9]  <- 1
NETWORK[8,9]  <- 1

# 10. erk
NETWORK[5,10]  <- 1
NETWORK[9,10]  <- 1

# 11. akt
NETWORK[1,11]  <- 1
NETWORK[5,11]  <- 1
NETWORK[10,11]  <- 1

return(NETWORK)
}





################################################################################
##################### function for computing the AUROC value and plotting the ROC

library(limma)
library(ROC)

compute_AUROC <- function(postP, trueEdges, steps=0.01){

n1 <- numeric(0)
for(i in 1:dim(postP)[1]){
        for(j in 1:dim(postP)[2]){
                                  if (abs(i-j)>0){
                                  n1 <- c(n1, postP[i,j])
                                  }
                                  }
                          }

n2 <- numeric(0)
for(i in 1:dim(trueEdges)[1]){
       for(j in 1:dim(trueEdges)[2]){
                              if (abs(i-j)>0){     
                              n2 <- c(n2, trueEdges[i,j])
                                        }
                                  }
                          }
                          
                          
sp <- sort(n1)         # posterior probabilities
st <- n2[order(n1)]    # true ordered by posterior probabilities

rocc.obj <- rocdemo.sca(n2, n1, rule=NULL, cutpts=NA)#seq(0,1,steps))

return(auROC(st,sp))

}

################################################################################
# Extrahiere das CPDAG des wahren Netwzerks

extract_cpdag_of_dag <- function(true_incidence){    

L <- list()

nodes <- dim(true_incidence)[1]
    
k <- cpdag(true_incidence)

    
dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
    
if(length(nrow(k))!=0){
    dummy[k[,1]] <- k[,5]
    L <- dummy
}
   
if(length(nrow(k))==0 && length(k)>0){
     dummy[k[1]] <- k[5]
     L <- dummy
}

   mat.com <- matrix(numeric(nodes*nodes),nrow=nodes)
   mat.re  <- matrix(numeric(nodes*nodes),nrow=nodes)

   com <- which(L>0)
   re  <- which(L<0)

   mat.com[com] <- 1
   mat.re[re] <- 1
   mat <- mat.com + mat.re + t(mat.re)

 return(mat)
}

### function for the computation of c(n, alpha)
c_function <- function(N,A){
  fact <- numeric(N)
  for (i in 1:N){
    fact[i] <- -lgamma((A+1-i)/2)
  }
  product <- sum(fact) -(A*N/2)*log(2)- (N*(N-1)/4)*log(pi)
  return(product)
}

#########################################
# Boris added 3 more functions.

library(MASS)
# Function that returns the conditional covariance and mean based on the full one
conddist <- function(child_index, covariance, incidence){
  # Determine the row/column index of the child node and those of the parents
  C <- child_index
  P <- which(incidence[,C]==1)
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

precision <- function(covariance, Sigma_cond, incidence){
  # We initialize an empty precision matrix of the correct dimensions.
  dim <- dim(covariance)
  W<-matrix(data = rep(NA,times=length(covariance)),nrow = dim[1],ncol = dim[2])
  # We precompute all of the b_ij's
  b<-W
  # We set all the b_ij's that correspond to a non-existent parent-child connection to 0.
  b[which(t(incidence)==0)]<-0
  for (j in 2:dim[1]) {
    # We now compute the matrix based on the parent sets from which to select the b_ij's
    P<-which(incidence[,j]==1)
    Sigma_11 <- covariance[P,P]
    Sigma_12 <- covariance[P,j]
    if(length(P)>0){
      Interaction <- t(Sigma_12)%*%ginv(Sigma_11)
    }
    else{
      Interaction <- 0
    }
    count<-1
    for (k in P) {
      b[j,k] <- Interaction[count]
      count <- count+1
    }
  }
  # We now compute the elements of W iteratively
  W[1,1] <- 1/Sigma_cond[1] # 1/v_1
  for (i in 1:(dim[1]-1)) {
    b_temp<-b[i+1,1:i]
    v_temp<-1/Sigma_cond[i+1] # 1/v_i
    W[1:i,1:i] <- W[1:i,1:i]+b_temp%*%t(b_temp)*v_temp
    W[1:i,i+1] <- -b_temp*v_temp
    W[i+1,1:i] <- -t(b_temp)*v_temp
    W[i+1,i+1] <- v_temp
  }
  return(W)
}

BGE_alpha <- function(a,mu,data,incidence,T_0,v=1){
  n <- nrow(data)
  m <- ncol(data)
  P_local_num <- numeric(n)   ###  numerator of the factors
  P_local_den <- numeric(n)   ### denumerator of the factors
  T_0_iteration <- a*T_0
  T_m <- T_0_iteration + (m-1)* cov(t(data)) + ((v*m)/(v+m))* (mu - rowMeans(data))%*%t(mu - rowMeans(data))
  
  for (j in 1:n)  {
    n_nodes <- which(incidence[,j]==1)         # parents of j
    l <- length(n_nodes)+1
    P_local_num[j] <- (-(length(n_nodes)+1)*m/2)*log(2*pi) + ((length(n_nodes)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes)+1),a-n+l)-c_function((length(n_nodes)+1),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0_iteration[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))
    if(sum(incidence[,j])>0){          # if j has at least one parent
      l <- length(n_nodes)
      P_local_den[j] <- (-(length(n_nodes))*m/2)*log(2*pi) + (length(n_nodes)/2)*log(v/(v+m)) + c_function(length(n_nodes),a-n+l)- c_function(length(n_nodes),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0_iteration[n_nodes,n_nodes])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[n_nodes,n_nodes])))
    }
    else{                              # if j has no parents
      P_local_den[j] <- 0
    }
  }
  bge <- (sum(P_local_num))-(sum(P_local_den))
  return(bge)
}
