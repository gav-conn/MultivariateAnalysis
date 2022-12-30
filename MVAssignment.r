# Load in required packages
library('e1071')
library('pls')
library("matlib")

# Q1
data = read.csv("Milk_MIR_Traits_data.csv")
# Set seed to student number
set.seed(18308483)
n = nrow(data)
#generate random number between 1 & n
rdm_number = sample(1:n, 1)
data <- data[-rdm_number,] # delete row number from dataset
n = n-1

#Q2
ncol(data)
spectra <- data[,52:582] # last 531 columns correspond to spectral readings
# plot column means of spectral values to get idea of their distribution
plot(colMeans(spectra), pch = 20, col = 4,
    xlab = 'Spectral Values', ylab = "Column Means of Spectral Values")
# plot column standard deviations of spectral values
plot(apply(spectra, 2, sd),  pch = 20, col = 8,
     xlab = 'Spectral Values', ylab = "Column sds of Spectral Values")

# 9th column corresponds to measurements of alphas1-casein
alphas1 <- data[,9]
plot(alphas1, pch = 20, main = 'alphas1-casein values') # large portion of missing data apparent in plot
# remove observations which have no value for alphas1-casein
cdata <- subset(data, !is.na(alphas1))

# calculate standard deviation of alphas1-casein observations
sd_as1 = sqrt(var(cdata[,9]))
# set upper/lower bounds for data to be used in our analysis
sd_up <- mean(cdata[,9]) + 3*sd_as1
sd_down <- mean(cdata[,9]) - 3*sd_as1

which(cdata[,9] > sd_up) # observations 71, 73 & 81 are 3 sds higher than mean
which(cdata[,9] < sd_down) # no observations are 3 sds lower than mean
cdata <- cdata[-which(cdata[,9] > sd_up), ] #remove observations which are over 3 sds from eman
# plot alphas1-casein data without the removed values
plot(cdata$alpha_s1_casein, pch = 20, main = 'alphas1-casein values')

#Q3
spectra <- as.matrix(cdata[,52:582])
sc_spectra <- scale(spectra)

# Hierarchical Clustering
# compare average linkage & complete linkage solutions
cl.average = hclust(dist(sc_spectra, method = 'euclidean'), method="average")
plot(cl.average, labels = FALSE) # Looks to be roughly 6 clusters
hcl.average = cutree(cl.average, k = 6)
table(hcl.average) # average linkage creates 2 large clusters as well as many smaller clusters

plot(cdata[,9], col = hcl.average, pch = 20,
     xlab = 'Spectral Values', ylab = 'alphas1 values')
plot(colMeans(spectra), col = hcl.average, pch = 20,
     xlab = 'Spectral Values', ylab = 'Column means of spectral values')
plot(apply(spectra, 2, sd), col = hcl.average, pch = 20,
     xlab = 'Spectral Values', ylab = 'Column sds of spectral values')
# average linkage gives us one very big cluster with lots of smaller ones which is not ideal

cl.complete = hclust(dist(sc_spectra, method = 'euclidean'), method="complete")
plot(cl.complete, labels = FALSE) 
# dendrogram for complete linkage looks much better
hcl.complete = cutree(cl.complete, k = 5) # seems to roughly  main clusters in data
table(hcl.complete)

plot(cdata[,9], col = hcl.complete, pch = 20,
     xlab = 'Spectral Values', ylab = 'alphaS1 values')
plot(colMeans(spectra), col = hcl.complete, pch = 20,
     xlab = 'Spectral Values', ylab = 'Column means of spectral values')
plot(apply(spectra, 2, sd), col = hcl.complete, pch = 20,
     xlab = 'Spectral Values', ylab = 'Column sds of spectral values')

table(hcl.average, hcl.complete)# clustering methods seem to be giving similar solutions
classAgreement(table(hcl.average, hcl.complete))

# K-means Clustering
WGSS = rep(0, 20)
n = nrow(sc_spectra)
WGSS[1] = (n-1)*sum(apply(spectra, 2, var))
for (k in 2:20) {
  WGSS[k] = sum(kmeans(spectra, centers = k)$withinss)
}
plot(1:20, WGSS, pch = 20, type = 'b', xlab = 'k', ylab = 'within group sum of squares', 
     main = 'WGSS vs number of clusters')
# Choosing k = 4 for our number of centres, seems optimal as it provides us with a low within group sum of squares without letting number of clusters grow too large 
k = 4
cl.kmeans = kmeans(spectra, centers = k)
plot(cdata[,9], pch = 20, col = cl.kmeans$cluster)
plot(colMeans(spectra), pch = 20, col = cl.kmeans$cluster)
plot(apply(spectra,2,sd), pch = 20, col = cl.kmeans$cluster)
table(cl.kmeans$cluster)

tab <- table(hcl.complete, cl.kmeans$cluster)
classAgreement(tab) # seems to be some agreement between solutions, but agreement not particularly strong

# Q4
fit <- prcomp(spectra, scale = TRUE)
CPV <- summary(fit)$importance[3,1:10]
plot(CPV, type = "b", pch = 20, col = 4,
     xlab = "Principal Components", ylab = "Cumulative Proportion of Variance",
     main = 'Plot of Cum. Prop. of Var. explained by No. of Components')

# Pretty much all of the variance explained by the first 4 principal components (over 99%) 
# with only marginally more variance explained by additional PCs
pairs(fit$x[,1:4], col = cl.kmeans$cluster)

# Q5
sc_spectra <- scale(spectra)
cov_spectra <- cov(sc_spectra)
eigen_spectra <- eigen(cov_spectra)
evect_spectra <- eigen_spectra$vectors
eval_spectra <- eigen_spectra$values
# Check Principal Components are same as those calculated in Q4
CPV_1 <- matrix(data = NA, nrow = 1, ncol = 10)
for (i in 1:10) {
  CPV_1[i] <- sum(eval_spectra[1:i]/sum(eval_spectra))
}  
rbind(CPV, CPV_1)
# Check sum of eigenvalues equal to sum of variances
sum(apply(sc_spectra, 2, var))
sum(eval_spectra)

# Calculate PC scores for each of the first 4 Principal Components
Y_PC <- matrix(data = NA, nrow = nrow(spectra), ncol = 4)
for(i in 1:nrow(spectra)){
  for (j in 1:4) {
    Y_PC[i, j] = evect_spectra[,j]%*%sc_spectra[i,] 
  }
}

colMeans((Y_PC-fit$x[,1:4])) # calcuated PC scores are same as outputted from prcomp
pairs(Y_PC, col = cl.kmeans$cluster, main = 'Pairs of 1st 4 PC Scores coloured by k-means cluster')

# Q7
alphas1 <- cdata[,9]
plsr_data <- cbind(alphas1, spectra)
plsr_data <- as.data.frame(plsr_data)
train <- sample(1:nrow(plsr_data), 2/3*nrow(plsr_data), replace = FALSE)

plsr.fit <- plsr(alphas1~., ncomp = 10, data = plsr_data, subset = train, scale = FALSE, validation = 'CV')
summary(plsr.fit) # Lowest adjusted Cross-Validation error for 5 component model
# Model explains 59.5% of the variance in Y & 98.6% of variation in X
# Predict alphas1-casein levels using the 5 component model
pred_test <- predict(plsr.fit, newdata = plsr_data[-train,], ncomp = 5)
# Examine performance of model
#plot(pred_test-plsr_data[-train,]$alphas1, pch = 20, ylab = 'Difference between Fitted & Actual Values') # errors seem to be distributed with mean 0
mean((pred_test-plsr_data[-train,]$alphas1)^2) # relatively low mean squared error on our test set

plot(alphas1[-train], pch = 20, main = 'Plot of Fitted vs Actual Values', 
     ylab = 'alphas1-casein values')
points(pred_test, col = 2, pch = 20) # model seems to do a decent job in estimating alphas1-casein levels

# Q8
plsr_train_X <- as.matrix(spectra[train,])
plsr_test_X <- as.matrix(spectra[-train,])

X <- plsr_train_X
Y <- alphas1[train]
m = 6

W <- matrix(data = NA, nrow = ncol(spectra), ncol = m)
T <- matrix(data = NA, nrow = nrow(plsr_train_X), ncol = m)
P <- matrix(data = NA, nrow = ncol(spectra), ncol = m)
Q <- matrix(data = NA, nrow = 1, ncol = m)

E <- X
F <- Y
for (i in 1:m) {
  #print(i)
  # Get SVD of S
  S = t(E)%*%F
  svd_S <- svd(S)
  # Define w/q as the left/right singular vectors
  w <- svd_S$u
  # Define t & normalise it
  t <- E%*%w
  t <- t/as.vector(sqrt(t(t)%*%t))
  # define p.q
  p <- t(E)%*%t
  q <- t(F)%*%t
  # Deflate data to remove info relating to previous latent variable
  E <- E-t%*%t(p)
  F <- F-t%*%t(q)
  # Store each of w,t,p & q as a column of relevant matrix
  W[,i] <- w
  T[,i] <- t
  P[,i] <- p
  Q[,i] <- q
}

#Calculate Regression Coefficients
R <- W%*%inv(t(P)%*%W)
B <- R%*%t(Q)

# Predict the alphas1-casein values of the test set
# Multiply scaled data by regression coefficients & add on mean of training data
pred2_test <- plsr_test_X%*%B
# Plot fitted values vs actual test alphas1-values
plot(alphas1[-train], pch = 20, main = 'Fitted vs Actual Values', ylab = 'alphas1-casein values')
points(pred_test, pch = 15, col = '#a60a0a')
points(pred2_test, pch = 20, col = '#fa4b4b')
# We see both predictions give the same estimates for alphas1-casein level