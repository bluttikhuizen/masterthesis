# Peregrine Script

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

strMCMC <- function(Data,incidence,iterations,step_save, T_0, a, fan.in=nrow(Data)-1, v=1, mu=numeric(nrow(Data))){
  n <- nrow(Data)                                # number of nodes
  m <- ncol(Data)                                # number of observations
  
  T_m <- T_0 + (m-1)*cov(t(Data)) + ((v*m)/(v+m))*(mu - rowMeans(Data))%*%t(mu - rowMeans(Data))
  
  L1 <- list()    # incidence matrix
  L2 <- list()    # log BGe score
  ################################################################################
  ##### functions we need in the algorithm
  
  ### calculation of the first ancestor matrix:
  ancestor <- function(incidence){
    incidence1 <- incidence
    incidence2 <- incidence
    k <- 1
    while (k < nrow(incidence)){
      incidence1 <- incidence1%*%incidence
      incidence2 <- incidence2 + incidence1
      k <-k+1
    }
    incidence2[which(incidence2[,]>0)] <- 1
    return(t(incidence2))}
  
  ### function for the computation of c(n, alpha)
  c_function <- function(N,A){
    fact <- numeric(N)
    for (i in 1:N){
      fact[i] <- -lgamma((A+1-i)/2)
    }
    product <- sum(fact) -(A*N/2)*log(2)- (N*(N-1)/4)*log(pi)
    return(product)}
  
  
  top_order <- function(incidence){
    n<-nrow(incidence)
    Order <- numeric(n)
    fan_in <- numeric(n)
    no_fan_in <- numeric(0)
    m <- 1
    for (p in 1:n){                                       # number of parent nodes at the beginning
      fan_in[p] <- sum(incidence[,p])
    }
    no_fan_in <- which(fan_in==0)
    while (length(which(Order==0))>0){                    # as long as there is a node without an order
      fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
      no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
      Order[m] <- no_fan_in[1]
      no_fan_in <- no_fan_in[-1]
      m <- m+1
    }
    return(Order)
  }
  
  
  ### assign the topological order of the descendants of the child
  des_top_order <- function(incidence, ancest1,child){
    top <- top_order(incidence)
    position_child <- which(top==child)
    top_all_after <- top[position_child:n]                # top. order without the "first" nodes
    desc <- which(ancest1[,child]==1)                     # descendants of the child
    inter_step <- c(child,desc,top_all_after)
    des_top <- inter_step[which(duplicated(inter_step))]
    return(des_top)
  }
  
  ################################################################################
  ### computation of the (logarithmizid) BGe Score of the FIRST graph
  P_local_num <- numeric(n)   ###  numerator of the factors
  P_local_den <- numeric(n)   ### denumerator of the factors
  
  for (j in 1:n)  {
    n_nodes <- which(incidence[,j]==1)         # parents of j
    l <- length(n_nodes)+1
    P_local_num[j] <- (-(length(n_nodes)+1)*m/2)*log(2*pi) + ((length(n_nodes)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes)+1),a-n+l)-c_function((length(n_nodes)+1),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))
    if(sum(incidence[,j])>0){          # if j has at least one parent
      l <- length(n_nodes)
      P_local_den[j] <- (-(length(n_nodes))*m/2)*log(2*pi) + (length(n_nodes)/2)*log(v/(v+m)) + c_function(length(n_nodes),a-n+l)- c_function(length(n_nodes),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[n_nodes,n_nodes])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[n_nodes,n_nodes])))
    }
    else{                              # if j has no parents
      P_local_den[j] <- 0
    }
  }
  bge_old <- (sum(P_local_num))-(sum(P_local_den))
  
  # first ancestor matrix
  ancest1 <- ancestor(incidence)
  
  ####### ... the number of neighbour graphs/proposal probability for the FIRST graph
  ### 1.) number of neighbour graphs obtained by edge deletions
  num_deletion <- sum(incidence)
  
  ### 2.) number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
  inter_add <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence - ancest1 >0)
  add <- matrix(numeric(n*n),nrow=n)
  add[inter_add] <- 1
  add[,which(colSums(incidence)>fan.in-1)] <- 0
  num_addition <- sum(add)
  
  ### 3.) number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
  inter_rev <- which(incidence - t(t(incidence)%*% ancest1)==1)
  re <- matrix(numeric(n*n),nrow=n)
  re[inter_rev] <- 1
  re[which(colSums(incidence)>fan.in-1),] <- 0 # CORRECTED!!!???!!!
  num_reversal <- sum(re)
  
  ##### total number of neighbour graphs:
  total <- sum(num_deletion,num_addition,num_reversal)
  
  ### proposal probability:
  proposal <- 1/total
  
  ############## sampling a new graph (or rather sampling an edge to shift)
  ### sample one of the three single edge operations
  random <- sample(1:total,1)
  
  operation <- 0                           # memorise, if the single edge operation is (will be) an edge reversal
  if (random > total - num_reversal){
    operation <- 1}
  
  #### shifting of the incidence matrix
  incidence_new <- incidence
  
  if (random <= num_deletion){             # if edge deletion was sampled
    if(length(which(incidence>0))>1){
      new_edge <- sample(which(incidence>0),1)} # sample one of the existing edges
    else
    {new_edge <- which(incidence>0)}
    incidence_new[new_edge] <- 0}            # and delete it
  
  if (random > (total - num_reversal)){      # if edge reversal was sampled
    if(num_reversal>1){
      new_edge <- sample(which(re==1),1)}  # sample one of the existing edges where a reversal leads to a valid graph
    else{
      new_edge <- which(re==1)}
    incidence_new[new_edge] <- 0             # delete it
    junk <- matrix(numeric(n*n),nrow=n)      # creating a matrix with all entries zero
    junk[new_edge] <- 1                      # an only a "1" at the entry of the new (reversed) edge
    incidence_new <- incidence_new + t(junk)}# sum the deleted matrix and the "junk-matrix"
  
  if (random <= (total - num_reversal) & random > num_deletion){     # if edge addition was sampled
    if(num_addition>1){
      new_edge <- sample(which(add==1),1)} # sample one of the existing edges where a addition leads to a valid graph
    else{
      new_edge <- which(add==1)}
    incidence_new[new_edge] <- 1             # and add it
  }
  
  
  #################### Updating the ancestor matrix
  
  # creating a matrix with dimensions of the incidence matrix and all entries zero except for the entry of the chosen edge
  help_matrix <- matrix(numeric(n*n),nrow=n)
  help_matrix[new_edge] <- 1
  
  # numbers of the nodes that belong to the shifted egde
  parent <- which(rowSums(help_matrix)==1)
  child <- which(colSums(help_matrix)==1)
  
  ### updating the ancestor matrix (after edge reversal)
  ## edge deletion
  ancestor_new <- ancest1
  if (operation==1){
    ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0           # delete all ancestors of the child and its descendants                                           #
    top_name <- des_top_order(incidence_new, ancest1, child)
    for (d in top_name){
      for(g in which(incidence_new[,d]==1)) {
        ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
      }
    }
    ## edge addition
    anc_parent <- which(ancestor_new[child,]==1)                     # ancestors of the new parent
    des_child <- which(ancestor_new[,parent]==1)                     # descendants of the child
    ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
  }
  
  ### updating the ancestor matrix (after edge deletion)
  if (random <= num_deletion){
    ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0           # delete all ancestors of the child and its descendants                                           #
    top_name <- des_top_order(incidence_new, ancest1, child)
    for (d in top_name){
      for(g in which(incidence_new[,d]==1)) {
        ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
      }
    }
  }
  
  # updating the ancestor matrix (after edge addition)
  if (random <= total - num_reversal & random > num_deletion){
    anc_parent <- which(ancest1[parent,]==1)        # ancestors of the new parent
    des_child <- which(ancest1[,child]==1)          # descendants of the child
    ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
  }
  
  ####### ... the number of neighbour graphs/proposal probability for the proposed graph
  ### 1.) number of neighbour graphs obtained by edge deletions
  num_deletion_new <- sum(incidence_new)
  
  ### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
  inter_add.new <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence_new - ancestor_new >0)
  add.new <- matrix(numeric(n*n),nrow=n)
  add.new[inter_add.new] <- 1
  add.new[,which(colSums(incidence_new)>fan.in-1)] <- 0
  num_addition_new <- sum(add.new)
  
  ### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
  inter_rev.new <- which(incidence_new - t(t(incidence_new)%*% ancestor_new)==1)
  re.new <- matrix(numeric(n*n),nrow=n)
  re.new[inter_rev.new] <- 1
  re.new[which(colSums(incidence_new)>fan.in-1),] <- 0  # CORRECTED!!!???!!!
  num_reversal_new <- sum(re.new)
  
  ##### total number of neighbour graphs:
  total_new <- sum(num_deletion_new,num_addition_new,num_reversal_new)
  
  ### proposal probability:
  proposal_new <- 1/total_new
  
  ### BGe Score for the new graph
  P_local_num_new <- P_local_num
  P_local_den_new <- P_local_den
  n_nodes_new <- which(incidence_new[,child]==1)
  l <- length(n_nodes_new)+1
  
  P_local_num_new[child] <- (-(length(n_nodes_new)+1)*m/2)*log(2*pi) + ((length(n_nodes_new)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes_new)+1),a-n+l)-c_function((length(n_nodes_new)+1),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))
  
  if(sum(incidence_new[,child])>0){           # if child at least one parent
    l <- length(n_nodes_new)
    P_local_den_new[child] <- (-(length(n_nodes_new))*m/2)*log(2*pi) + (length(n_nodes_new)/2)*log(v/(v+m)) + c_function(length(n_nodes_new),a-n+l)- c_function(length(n_nodes_new),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[n_nodes_new,n_nodes_new])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[n_nodes_new,n_nodes_new])))
  }
  else{                                       # if child has no parents
    P_local_den_new[child] <- 0
  }
  
  if (operation==1){                          # if single edge operation was an edge reversal
    n_nodesP <- which(incidence_new[,parent]==1)
    l <- length(n_nodesP)+1
    P_local_num_new[parent] <- (-(length(n_nodesP)+1)*m/2)*log(2*pi) + ((length(n_nodesP)+1)/2)*log(v/(v+m)) + c_function((length(n_nodesP)+1),a-n+l)-c_function((length(n_nodesP)+1),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))
    if(sum(incidence_new[,parent])>0){          # if parent at least one parent
      l <- length(n_nodesP)
      P_local_den_new[parent] <- (-(length(n_nodesP))*m/2)*log(2*pi) + (length(n_nodesP)/2)*log(v/(v+m)) + c_function(length(n_nodesP),a-n+l)- c_function(length(n_nodesP),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[n_nodesP,n_nodesP])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[n_nodesP,n_nodesP])))
    }
    else{                                       # if parent has no parents
      P_local_den_new[parent] <- 0
    }
  }
  bge_new <- (sum(P_local_num_new))-(sum(P_local_den_new))
  
  
  L1[[1]] <- incidence                        # initial graph
  L2[[1]] <- bge_old                          # and it`s BGe score
  
  acceptance <- min(1, exp((bge_new + log(proposal_new)) - (bge_old  + log(proposal))))
  rand <- runif(1)
  
  if(acceptance > rand){
    incidence <- incidence_new
    bge_old <- bge_new
    P_local_num <- P_local_num_new
    P_local_den <- P_local_den_new
    proposal <- proposal_new
    ancest1 <- ancestor_new
    total <- total_new
    num_deletion <- num_deletion_new
    num_addition <- num_addition_new
    num_reversal <- num_reversal_new
    add <- add.new
    re <- re.new
  }
  
  ####################################################################################################################################
  #################################################################################
  
  for (z in 2:((iterations/step_save)+1)){
    for (count in 1:step_save){
      
      ############## sampling a new graph (or rather sampling an edge to shift)
      ### sample one of the three single edge operations
      random <- sample(1:total,1)
      
      operation <- 0                            # memorise, if the single edge operation is (will be) an edge reversal
      if (random > total - num_reversal){
        operation <- 1}
      
      #### shifting of the incidence matrix
      incidence_new <- incidence
      
      if (random <= num_deletion){              # if edge deletion was sampled
        if(length(which(incidence>0))>1){
          new_edge <- sample(which(incidence>0),1)} # sample one of the existing edges
        else
        {new_edge <- which(incidence>0)}
        incidence_new[new_edge] <- 0}            # and delete it
      
      if (random > (total - num_reversal)){    # if edge reversal was sampled
        if(num_reversal>1){
          new_edge <- sample(which(re==1),1)}      # sample one of the existing edges where a reversal leads to a valid graph
        else{
          new_edge <- which(re==1)}
        incidence_new[new_edge] <- 0             # delete it
        junk <- matrix(numeric(n*n),nrow=n)      # creating a matrix with all entries zero
        junk[new_edge] <- 1                      # an only a "1" at the entry of the new (reversed) edge
        incidence_new <- incidence_new + t(junk)}# sum the deleted matrix and the "junk-matrix"
      
      if (random <= (total - num_reversal) & random > num_deletion){     # if edge addition was sampled
        if(num_addition>1){
          new_edge <- sample(which(add==1),1)} # sample one of the existing edges where a addition leads to a valid graph
        else{
          new_edge <- which(add==1)}
        incidence_new[new_edge] <- 1             # and add it
      }
      
      
      ### Updating the ancestor matrix
      
      # creating a matrix with dimensions of the incidence matrix and all entries zero except for the entry of the chosen edge
      help_matrix <- matrix(numeric(n*n),nrow=n)
      help_matrix[new_edge] <- 1
      
      # numbers of the nodes that belong to the shifted egde
      parent <- which(rowSums(help_matrix)==1)
      child <- which(colSums(help_matrix)==1)
      
      ### updating the ancestor matrix (after edge reversal)
      ## edge deletion
      ancestor_new <- ancest1
      if (operation==1){
        ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants                                           #
        
        top_name <- des_top_order(incidence_new, ancest1, child)
        for (d in top_name){
          for(g in which(incidence_new[,d]==1)) {
            ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
          }
        }
        
        anc_parent <- which(ancestor_new[child,]==1)          # ancestors of the new parent
        des_child <- which(ancestor_new[,parent]==1)          # descendants of the child
        ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
      }
      
      ### updating the ancestor matrix (after edge deletion)
      if (random <= num_deletion){
        ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants                                           #
        top_name <- des_top_order(incidence_new, ancest1, child)
        for (d in top_name){
          for(g in which(incidence_new[,d]==1)) {
            ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
          }
        }
      }
      
      # updating the ancestor matrix (after edge addition)
      if (random <= total - num_reversal & random > num_deletion){
        anc_parent <- which(ancest1[parent,]==1)             # ancestors of the new parent
        des_child <- which(ancest1[,child]==1)               # descendants of the child
        ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
      }
      
      ####### ... the number of neighbour graphs/proposal probability for the proposed graph
      ### 1.) number of neighbour graphs obtained by edge deletions
      num_deletion_new <- sum(incidence_new)
      
      ### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
      inter_add.new <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence_new - ancestor_new >0)
      add.new <- matrix(numeric(n*n),nrow=n)
      add.new[inter_add.new] <- 1
      add.new[,which(colSums(incidence_new)>fan.in-1)] <- 0
      num_addition_new <- sum(add.new)
      
      ### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
      inter_rev.new<- which(incidence_new - t(t(incidence_new)%*% ancestor_new)==1)
      re.new <- matrix(numeric(n*n),nrow=n)
      re.new[inter_rev.new] <- 1
      re.new[,which(colSums(incidence_new)>fan.in-1)] <- 0
      num_reversal_new <- sum(re.new)
      
      ##### total number of neighbour graphs:
      total_new <- sum(num_deletion_new, num_addition_new, num_reversal_new)
      
      ### proposal probability:
      proposal_new <- 1/total_new
      
      ### BGe Score for the new graph
      P_local_num_new <- P_local_num
      P_local_den_new <- P_local_den
      n_nodes_new <- which(incidence_new[,child]==1)
      l <- length(n_nodes_new)+1
      
      P_local_num_new[child] <- (-(length(n_nodes_new)+1)*m/2)*log(2*pi) + ((length(n_nodes_new)+1)/2)*log(v/(v+m)) + c_function((length(n_nodes_new)+1),a-n+l)-c_function((length(n_nodes_new)+1),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))
      
      if(sum(incidence_new[,child])>0){       # if child at least one parent
        l <- length(n_nodes_new)
        P_local_den_new[child] <- (-(length(n_nodes_new))*m/2)*log(2*pi) + (length(n_nodes_new)/2)*log(v/(v+m)) + c_function(length(n_nodes_new),a-n+l)- c_function(length(n_nodes_new),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[n_nodes_new,n_nodes_new])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[n_nodes_new,n_nodes_new])))
      }
      else{                                   # if child has no parents
        P_local_den_new[child] <- 0
      }
      
      if (operation==1){                      # if single edge operation was an edge reversal
        n_nodesP <- which(incidence_new[,parent]==1)
        l <- length(n_nodesP)+1
        P_local_num_new[parent] <- (-(length(n_nodesP)+1)*m/2)*log(2*pi) + ((length(n_nodesP)+1)/2)*log(v/(v+m)) + c_function((length(n_nodesP)+1),a-n+l)-c_function((length(n_nodesP)+1),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))
        if(sum(incidence_new[,parent])>0){      # if parent at least one parent
          l <- length(n_nodesP)
          P_local_den_new[parent] <- (-(length(n_nodesP))*m/2)*log(2*pi) + (length(n_nodesP)/2)*log(v/(v+m)) + c_function(length(n_nodesP),a-n+l)- c_function(length(n_nodesP),a+m-n+l)+ ((a-n+l)/2)*log(det(as.matrix(T_0[n_nodesP,n_nodesP])))+ (-(a+m-n+l)/2)*log(det(as.matrix(T_m[n_nodesP,n_nodesP])))
        }
        else{                                   # if parent has no parents
          P_local_den_new[parent] <- 0
        }
      }
      bge_new <- (sum(P_local_num_new))-(sum(P_local_den_new))
      
      acceptance <- min(1, exp((bge_new +log(proposal_new)) - (bge_old +log(proposal))))
      rand <- runif(1)
      
      if(acceptance > rand){
        incidence <- incidence_new
        bge_old <- bge_new
        P_local_num <- P_local_num_new
        P_local_den <- P_local_den_new
        proposal <- proposal_new
        ancest1 <- ancestor_new
        total <- total_new
        num_deletion <- num_deletion_new
        num_addition <- num_addition_new
        num_reversal <- num_reversal_new
        add <- add.new
        re <- re.new
      }
    }
    
    L1[[z]] <- incidence
    L2[[z]] <- bge_old
  }
  return(list(L1,L2))
}

################################################################################

child <- function(edges,n){         # input: the numbers of the edges in the incidence matrix and the number of nodes
  p <- ceiling(edges/n)
  return(p)
}

parent <- function(edges,n){
  ch <- edges + n - child(edges,n)*n
  return(ch)
}


top_order <- function(incidence){
  n <- nrow(incidence)
  Order <- numeric(n)
  fan_in <- numeric(n)
  no_fan_in <- numeric(0)
  m <- 1
  for (p in 1:n){                                       # number of parent nodes at the beginning
    fan_in[p] <- sum(incidence[,p])
  }
  no_fan_in <- which(fan_in==0)
  while (length(which(Order==0))>0){                    # as long as there is a node without an order
    fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
    no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
    Order[m] <- no_fan_in[1]
    no_fan_in <- no_fan_in[-1]
    m <- m+1
  }
  return(Order)
}

################################################################################
order.edges <- function(incidence){
  top.order <- top_order(incidence)
  n <- length(top.order)
  edges <- which(incidence!=0)
  children <- child(edges,n)
  parents <- parent(edges,n)
  m <- length(edges)
  ordered_edges  <- numeric(m)
  incidence_n <- incidence
  tog <- matrix(c(edges,parents,children,ordered_edges),ncol=4, byrow=FALSE)
  k <- 1
  while(any(tog[,4]==0)){
    node1 <- top.order[which(colSums(incidence_n[,top.order])>0)][1]    # first node in top. order that has at least one parent
    par1<- tog[which(tog[,3]==node1),2]                # find the parents of  first child in the top. order that has an unordered edge incident into it
    g <- par1[which(par1>0)]
    f1 <- numeric(length(g))
    for (i in 1:length(g)){
      f1[i] <- which(top.order==g[i])
    }
    par2 <- g[which.max(f1)]                           # find the highest ordered node that has an edge leading into node1
    tog[which(tog[,2]==par2 & tog[,3]==node1),4] <- k
    k <- k + 1
    incidence_n[tog[which(tog[,2]==par2 & tog[,3]==node1),1]] <- 0     # delete the edge in the "incidence" matrix
    tog[which(tog[,2]==par2 & tog[,3]==node1),2] <- 0
  }
  to <- matrix(c(edges,parents,children,tog[,4]),ncol=4,byrow=FALSE)
  return(to)                                          # return the whole matrix, the order is the fourth column
}

#################################################################################
### DAG-to-CPDAG algorithm
# +1 if the edge is "compelled"
#   -1 if the edge is "reversible"
#######################################################################
cpdag <- function(incidence){
  z <- order.edges(incidence)
  new_mat <- cbind(z,numeric(nrow(z)))    # edges, parents, children, order, zeros
  n_mat <- new_mat[order(new_mat[,4]),]   # sort the edges by its order
  vec <- numeric(nrow(z))
  while(any(vec==0)){                                  # while there are unlabeled edges            l.3
    if (length(vec)>1){                                  # if there are at least 2 edges
      first <- which(n_mat[,5]==0)[1]                    # first EDGE that ist labeled "unknown" (0)  l.4
      parent1 <- n_mat[first,2]                          # x   parent NODE
      child1 <- n_mat[first,3]                           # y   child NODE
      comp1 <- n_mat[which(n_mat[,3]==parent1 & n_mat[,5]==1),2]      # w NODES that have an edge incident into the parent labeled compelled)
    }
    if (length(vec)==1){
      first <- which(n_mat[5]==0)                      # first edge that ist labeled "unknown" (0)
      parent1 <- n_mat[2]                             # x   parent
      child1 <- n_mat[3]                              # y   child
      comp1 <- numeric(0)
    }
    for (j in comp1){                                   #                                            l.5
      if (incidence[j,child1]==0){                     # if w is not a parent of the child          l.6
        n_mat[first,5] <- 1                             # label x -> y compelled                     l.7
        n_mat[which(n_mat[,3]==child1),5] <- 1          # label every edge incident into y compelled l.7
        vec[first] <- 1
        vec[which(n_mat[,3]==child1)] <- 1
        break
      }
      if (incidence[j,child1]!=0)    {
        n_mat[which(n_mat[,2]==j & n_mat[,3]==child1),5] <- 1  # label w -> y compelled                l.10
        vec[which(n_mat[,2]==j & n_mat[,3]==child1)] <- 1
      }
    }
    if (length(vec)>1){
      if(n_mat[first,5]==0){
        
        moep <- n_mat[which(n_mat[,3]==child1 & n_mat[,2]!=parent1),2]      # other parents of the child
        if(length(moep)>0){                              #                     l.11
          for(o in moep){
            if(incidence[o,parent1]==0){
              vec[first] <- 1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- 1
              n_mat[first,5] <- 1                                     # label x -> y compelled
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- 1   # label all "unknown" edges incident into y compelled
              break
            }
            if(all(incidence[moep,parent1]!=0)){
              vec[first] <- -1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
              n_mat[first,5] <- -1                                    # label x -> y reversible
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1  # label all "unknown" edges incident into y reversible
            }
          }
        }
        if(length(moep)==0){
          vec[first] <- -1
          vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
          n_mat[first,5] <- -1                                    # label x -> y reversible
          n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1  # label all "unknown" edges incident into y reversible
        }
      }
    }
    if (length(vec)==1){
      n_mat[5] <- -1                                    # label x -> y reversible
      vec <- -1
    }
  }
  return(n_mat)
}


################################################################################

cpdag_list <- function(list.inc,E){    # E: end of burnIn phase
  L <- list()
  G <- list()
  nodes <- dim(list.inc[[1]])[1]
  mat.sum <- matrix(numeric(nodes*nodes),nrow=nodes)
  for (i in E:length(list.inc)){
    k <- cpdag(list.inc[[i]])
    dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
    if(length(nrow(k))!=0){
      dummy[k[,1]] <- k[,5]
      L[[i]] <- dummy
    }
    if(length(nrow(k))==0 && length(k)>0){
      dummy[k[1]] <- k[5]
      L[[i]] <- dummy
    }
    mat.com <-matrix(numeric(nodes*nodes),nrow=nodes)
    mat.re <- matrix(numeric(nodes*nodes),nrow=nodes)
    com <- which(L[[i]]>0)
    re <- which(L[[i]]<0)
    mat.com[com] <- 1
    mat.re[re] <- 1
    mat <- mat.com + mat.re + t(mat.re)
    G[[i]] <- mat
    mat.sum <- mat.sum + mat
  }
  return(list(L,G, (mat.sum/(length(list.inc)- E+1))))
}

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

coupling_1_update <- function(data, n, m, k, incidence_MCMC, iterations, step_save, a_1, T_0, a_2, V){
  # We initialize a list to store the sampled incidence matrices.
  MCMC <- list()
  scaling_stored <- c()
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
    scaling_stored[i] <- scaling
    for (j in 1:k) {
      MCMC_iteration[[j]] <- strMCMC(data[[j]],incidence_MCMC[[j]],1,1,scaling*T_0,a=a_1)
    }
    # We now sample a precision matrix from the posterior for each data set
    T_m <- list()
    W <- list()
    for (j in 1:k) {
      T_m[[j]] <- scaling*T_0 + (m-1)*cov(t(data[[j]])) + ((v*m)/(v+m))*(mu - rowMeans(data[[j]]))%*%t(mu - rowMeans(data[[j]]))
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
    
    # We also resample a_1 to determine the coupling strength.
    a_1_prop <- a_1+runif(1, min = -1, max = 1)
    
    # computation of the (logarithmizid) BGe Score of the current alpha
    a <- a_1
    bge_old <- BGE_alpha(a,mu,data[[1]],incidence_MCMC[[1]],T_0)+BGE_alpha(a,mu,data[[2]],incidence_MCMC[[2]],T_0)+BGE_alpha(a,mu,data[[3]],incidence_MCMC[[3]],T_0)+BGE_alpha(a,mu,data[[4]],incidence_MCMC[[4]],T_0)+BGE_alpha(a,mu,data[[5]],incidence_MCMC[[5]],T_0)
    
    # computation of the (logarithmizid) BGe Score of the PROPOSED alpha
    a <- a_1_prop
    bge_prop <- BGE_alpha(a,mu,data[[1]],incidence_MCMC[[1]],T_0)+BGE_alpha(a,mu,data[[2]],incidence_MCMC[[2]],T_0)+BGE_alpha(a,mu,data[[3]],incidence_MCMC[[3]],T_0)+BGE_alpha(a,mu,data[[4]],incidence_MCMC[[4]],T_0)+BGE_alpha(a,mu,data[[5]],incidence_MCMC[[5]],T_0)
    
    # We can now compute the acceptance ratio
    acceptance <- min(1, exp(bge_prop - bge_old))
    rand <- runif(1)
    
    if(acceptance > rand & a_1_prop>n+2){
      a_1<-a_1_prop
    }
  }
  return(list(MCMC,scaling_stored))
}

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
V <- diag(n) # Set the precision matrix for the prior distribution.
T_0 <- rWishart(1,a_2, ginv(V))[,,1] # Sample T_0 from the prior.
a_1 <- a_2 # We initialize a_1 at the minimal value.

# We call the functions that run the different sampling schemes.
update <- coupling_1_update(data, n, m, k, incidence_MCMC, iterations, step_save, a_1, T_0, a_2, V)
fixed <- coupling_fixed(data, n, m, k, incidence_MCMC, iterations, step_save, a_1, T_0, a_2, V)
structureMarco <- list()
for (i in 1:5) {
  structureMarco[[i]] <- strMCMC(Data=data[[i]],incidence_MCMC[[i]],iterations,step_save, a_1*T_0, a=a_1)
}

# Save the output generated by peregrine
save.image(file = "peregrine_output7.RData")