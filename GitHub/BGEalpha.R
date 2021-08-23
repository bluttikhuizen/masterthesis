BGE_alpha <- function(a,mu,data,incidence){
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
