# =================================
# Application Script -
# Standard Errors for IID Returns
# =================================

# Required libraries
library(RPESE)
library(metRology)

# ==============
# SHARPE RATIO
# ==============

# Number of replications
M <- 30000

# Sample sizes
sample.sizes <- c(60, 120, 240)

# Distribution parameters
mean <- 0.01
true.val <- c(0.2, 0.5)[1]
N.scale <- c(0.05, 0.02)[1]
t.scale <- c(0.039, 0.0155)[1]
data.dist <- c("N", "t")[1]

# Matrix for output
n.out <- matrix(nrow=3, ncol=5)

for(n in sample.sizes){
  
  # Setting seed
  set.seed(0)
  
  # Simulating data
  if(data.dist=="N")
    full.data <- sapply(1:M, function(x) rnorm(n, mean, N.scale), simplify=TRUE) else if(data.dist=="t")
      full.data <- sapply(1:M, function(x) rt.scaled(n, df=5, mean=mean, sd=t.scale), simplify=TRUE)
  # Casting as TS
  full.data <- as.ts(full.data)
  
  # Vectors to store output
  point.est <- SE.est <- numeric(M)
  
  # Computation over replications
  for(rep in 1:M){

    # Computation of SE output
    out <- SharpeRatio.SE(full.data[,rep, drop=FALSE], se.method=c("IFiid"))
    point.est[rep] <- as.numeric(out$SR)
    SE.est[rep] <- out$IFiid$se
  }
  
  # Computation of SEs output
  SE.mc <- sd(point.est)
  SE.IF <- mean(SE.est)
  rej.prob <- 1-mean((point.est-qnorm(0.975)*SE.IF<=true.val)*(point.est+qnorm(0.975)*SE.IF>=true.val))
  SE.bias <- (SE.IF - SE.mc)/SE.mc
  
  # Simulation output
  n.out[n==sample.sizes,] <- c(round(mean(point.est),4), round((mean(point.est)-true.val)/true.val*100, 1), round(SE.mc, 4), round(SE.IF, 4), round(SE.bias*100, 1), round(rej.prob, 3)*100)
}

# Sample size output
n.out

# ========================
# DOWNSIDE SHARPE RATIO
# ========================

# Number of replications
M <- 30000

# Sample sizes
sample.sizes <- c(60, 120, 240)

# Distribution parameters
mean <- 0.01
true.val <- c(0.2, 0.5)[1]
N.scale <- c(0.05, 0.02)[1]
t.scale <- c(0.039, 0.0155)[1]
data.dist <- c("N", "t")[1]

# Matrix for output
n.out <- matrix(nrow=3, ncol=5)

for(n in sample.sizes){
  
  # Setting seed
  set.seed(0)
  
  # Simulating data
  if(data.dist=="N")
    full.data <- sapply(1:M, function(x) rnorm(n, mean, N.scale), simplify=TRUE) else if(data.dist=="t")
      full.data <- sapply(1:M, function(x) rt.scaled(n, df=5, mean=mean, sd=t.scale), simplify=TRUE)
  # Casting as TS
  full.data <- as.ts(full.data)
    
  # Vectors to store output
  point.est <- SE.est <- numeric(M)
  
  # Computation over replications
  for(rep in 1:M){

    # Computation of SE output
    out <- SortinoRatio.SE(full.data[,rep, drop=FALSE], se.method=c("IFiid"), threshold=c("mean","const")[1])
    point.est[rep] <- as.numeric(out$SoR)/sqrt(2)
    SE.est[rep] <- out$IFiid$se/sqrt(2)
  }
  
  # Computation of SEs output
  SE.mc <- sd(point.est)
  SE.IF <- mean(SE.est)
  rej.prob <- 1-mean((point.est-qnorm(0.975)*SE.IF<=true.val)*(point.est+qnorm(0.975)*SE.IF>=true.val))
  SE.bias <- (SE.IF - SE.mc)/SE.mc
  
  # Simulation output
  n.out[n==sample.sizes,] <- c(round((mean(point.est)-true.val)/true.val*100, 1), round(SE.mc, 4), round(SE.IF, 4), round(SE.bias*100, 1), round(rej.prob, 3)*100)
}

# Sample size output
n.out


# ========================
# STANDARD DEVIATION
# ========================

# Number of replications
M <- 30000

# Sample sizes
sample.sizes <- c(60, 120, 240)

# Distribution parameters
mean <- 0.01
true.val <- c(0.02, 0.05)[1]
N.scale <- c(0.02, 0.05)[1]
t.scale <- c(0.0155, 0.039)[1]
data.dist <- c("N", "t")[1]

# Matrix for output
n.out <- matrix(nrow=3, ncol=5)

for(n in sample.sizes){
  
  # Setting seed
  set.seed(0)
  
  # Simulating data
  if(data.dist=="N")
    full.data <- sapply(1:M, function(x) rnorm(n, mean, N.scale), simplify=TRUE) else if(data.dist=="t")
      full.data <- sapply(1:M, function(x) rt.scaled(n, df=5, mean=mean, sd=t.scale), simplify=TRUE)
  # Casting as TS
  full.data <- as.ts(full.data)
  
  # Vectors to store output
  point.est <- SE.est <- numeric(M)
  
  # Computation over replications
  for(rep in 1:M){

    # Computation of SE output
    out <- StdDev.SE(full.data[,rep, drop=FALSE], se.method=c("IFiid"))
    point.est[rep] <- as.numeric(out$SD)
    SE.est[rep] <- out$IFiid$se
  }
  
  # Computation of SEs output
  SE.mc <- sd(point.est)
  SE.IF <- mean(SE.est)
  rej.prob <- 1-mean((point.est-qnorm(0.975)*SE.IF<=true.val)*(point.est+qnorm(0.975)*SE.IF>=true.val))
  SE.bias <- (SE.IF - SE.mc)/SE.mc
  
  # Simulation output
  n.out[n==sample.sizes,] <- c(round((mean(point.est)-true.val)/true.val*100, 1), round(SE.mc, 4), round(SE.IF, 4), round(SE.bias*100, 1), round(rej.prob, 3)*100)
}

# Sample size output
n.out


# ==========================
# SEMI-STANDARD DEVIATION
# ==========================

# Number of replications
M <- 30000

# Sample sizes
sample.sizes <- c(60, 120, 240)

# Distribution parameters
mean <- 0.01
true.val <- c(0.02, 0.05)[2]
N.scale <- c(0.02855, 0.0715)[1]
t.scale <- c(0.022, 0.055)[2]
data.dist <- c("N", "t")[2]

# Matrix for output
n.out <- matrix(nrow=3, ncol=5)

for(n in sample.sizes){
  
  # Setting seed
  set.seed(0)
  
  # Simulating data
  if(data.dist=="N")
    full.data <- sapply(1:M, function(x) rnorm(n, mean, N.scale), simplify=TRUE) else if(data.dist=="t")
      full.data <- sapply(1:M, function(x) rt.scaled(n, df=5, mean=mean, sd=t.scale), simplify=TRUE)
  # Casting as TS
  full.data <- as.ts(full.data)
  
  # Vectors to store output
  point.est <- SE.est <- numeric(M)
  
  # Computation over replications
  for(rep in 1:M){
    
    cat(rep, "\n")
    
    # Computation of SE output
    out <- SemiSD.SE(full.data[,rep, drop=FALSE], se.method=c("IFiid"))
    point.est[rep] <- as.numeric(out$SSD)
    SE.est[rep] <- out$IFiid$se
  }
  
  # Computation of SEs output
  SE.mc <- sd(point.est)
  SE.IF <- mean(SE.est)
  rej.prob <- 1-mean((point.est-qnorm(0.975)*SE.IF<=true.val)*(point.est+qnorm(0.975)*SE.IF>=true.val))
  SE.bias <- (SE.IF - SE.mc)/SE.mc
  
  # Simulation output
  n.out[n==sample.sizes,] <- c(round((mean(point.est)-true.val)/true.val*100, 1), round(SE.mc, 4), round(SE.IF, 4), round(SE.bias*100, 1), round(rej.prob, 3)*100)
}

# Sample size output
n.out


