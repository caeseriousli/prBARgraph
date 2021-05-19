poisson_BAR_StARS <- function(dt, sequen = sequen, eps = .01) {
  library(dplyr)
  library(doParallel)
  library(glmnet)
  library(prBARgraph)
  
  
  # Define col and row numbers to reduce repeated computations in loops
  n = nrow(dt)
  p = ncol(dt)
  # Define subsample size
  # used to be set to 2.5
  bn = floor(2.5*sqrt(n))
  # Calculate how many subsamples we can have at each iteration
  m = floor(n/bn)
  # How many combinations (potential edges) do we have (p choose 2)
  combination = p*(p-1)/2
  # number of lambda's
  nlambda = length(sequen)
  
  # Subsampling for partitioning
  # Subsampling, partitioning into m many subsamples
  partitions <- sample(rep(1:m, bn), m*bn, replace = FALSE)
  # if n is not divisible by subsample size bn, fill in the last few assignments at random
  partitions <- c(partitions, rep(m+1, n-m*bn))
  # use the partition factor vector to split data
  dt.temp = split(as.data.frame(dt), partitions)
  
  # Multithreading loops
  # Make sure the core count does not exceed your virtual core count  (or physical cores x 2, 
  # given that CPU is Intel Iris series) or R (and the machine) would crash
  # If not sure use detectCores(). It is recommended to use this number-1 to aviod crashes
  #cores = detectCores()
  #cl <- makeCluster((cores-1), type = "PSOCK")
  #registerDoParallel(cl)
  #registerDoParallel(4)
  
  #iterate through features
  instability.stacked = foreach(j = 1:p, .combine = cbind, .packages = c('doParallel')) %:%
    
    # iterate through subsamples
    foreach(i = 1:m, .combine = '+') %dopar% {
      
      subsample = as.matrix(dt.temp[[i]])
      
      # func = function(lamb) {
      #   fit.temp = prBARgraph::poisBAR(subsample[, j], subsample[, -j], lambda = lamb)
      #   return(fit.temp$coef)
      # }
      
      # fit = prBARgraph::poisBAR(subsample[, j], subsample[, -j], lambda = sequen)
      fit = prBARgraph::poisBAR(dt[, j], dt[, -j], lambda = sequen, xi = log(n), eps = 1E-6, xitype = "original", max.iter = 1000)
      
      #coeffs = sapply(sequen, func, simplify = "array")
      coeffs = fit$coef
      
      # if(sum(is.na(coeffs)) > 0) {
      # 
      #   print(paste("No good:", j, sequen[1], coeffs))
      #   write.csv(subsample, paste("./output/subsample_j",
      #                              j, ".csv", collapse = ""), row.names = F)
      #   
      #   # testing phase, setting to 0 for now
      #   coeffs[is.na(coeffs)] = 0
      # 
      # }
      
      #coeffs = as.matrix(coeffs)
      coeffs = ifelse(coeffs == 0, 0, 1)
      
      coeffs = coeffs/m
      
      #Add the diagnal element (vector) in the vector before binding into a matrix
      # Using an if statement because the indicator j would go out of bound at last iteration
      if (j > nrow(coeffs)) {
        return(rbind(coeffs, rep(0, nlambda)) )
      } else {
        return(rbind(coeffs[0:(j-1), ], rep(0, nlambda), coeffs[j:nrow(coeffs), ]))
      }
      
      
    }
    
    # disagreement <- disagreement/m
    # 
    # # Add the diagnal element (vector) in the vector before binding into a matrix
    # # Using an if statement because the indicator j would go out of bound at last iteration
    # if (j > nrow(disagreement)) {
    #   return(rbind(disagreement, rep(0, nlambda)) )
    # } else {
    #   return(rbind(disagreement[0:(j-1), ], rep(0, nlambda), disagreement[j:nrow(disagreement), ]))
    # }
    
  
  instability.results <- foreach(j = 1:nlambda, .combine = c) %dopar% {
    
    # we've stacked p (number of features) matrices each of which has lambda sequen many vectors
    # for each lambda, pick out the vector at that position from all these matrices
    A = instability.stacked[ ,seq(j, j+nlambda*(p-1), by = nlambda)]
    # We use the union of asymmetric adjacency matrix parameters, Poisson GGM equation 5
    instability = (t(A)+A)/2
    
    # STARS page 6, parameter xi
    instability = 2*instability*(1-instability)
    
    # Calculate and output total instability, STARS page 6, D
    # Note that the matrix is symmetric now with each edge added twice (I made diagnal)
    # elements zero. Thus it must be divided by 2 first and then (p choose 2)
    return(sum(instability)/2/combination)
  }
  
  #stopCluster(cl)
  #closeAllConnections()
  
  ####### Plot Results ##############
  stability.graph = rbind(sequen, instability.results)
  stability.graph = t(stability.graph)
  names(stability.graph) = c("lambda", "instability")
  
  #plot(stability.graph, xlab = expression(lambda), ylab = "Instability")
  
  return(instability.results)
}