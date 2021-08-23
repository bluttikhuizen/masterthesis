# We first create test data to produce a covariance matrix.
data <- make_data_simple(100,1)

# We now compute the covariance matrix and the mean.
covariance <- cov(t(data))

# We model consider the following simple linear regression model for node_2.
regression_model_1 <- lm(data[2,]~data[1,])
regression_model_2 <- lm(data[3,]~data[2,])

# We can now extract the interaction coefficient with node 1.
beta_1 <- regression_model_1$coefficients[2]
beta_2 <- regression_model_2$coefficients[2]

# We now compute the coefficients b_ij calculated in the iterative computation of W.
b_21 <- t(covariance[1,2])%*%ginv(covariance[1,1])
b_32 <- t(covariance[2,3])%*%ginv(covariance[2,2])