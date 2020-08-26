#' Fit BAR Poisson Graphical Model.
#' 
#' @param dt nxp dataset.
#' @param instabilities if not NA, shoose to use StARS instability for selecting a lambda. Run poisson_BAR_StARS() first to get StARS
#' @param sequen to be used with instability a grid of l0/l1 penalization parameters to be chosen from.
#' @param beta stability 
#' @param regularization default to l0 (BAR). The function also supports l1 (LASSO Poisson)
#' @param lchosen if not NA, do not use StARS to select Lambda. Instead, run all penalizations and return a list of adjacency matrices
#' @return either an adjacenty matrix (selected by StARS), or a list of adjacency matrices, corresponding to each lambda supplied by 'lchose' option.
#' @export 
#' @examples
#' require(XMRF)
#' # Simulate RNA data from XMRF package (Wan, Y. W., et. al, 2016)
#' Xsim <- XMRF.Sim(n=100, p=200, model="LPGM", graph.type="hub")
#' 
#' # Specify a grid of l0 penalization parameter
#' sequenBAR = 2^seq(logb(18, base = 2), logb(1e-4, base = 2), length = 15)
#' 
#' # Run StARS to calculate stability (Liu, H. et. al, 2010)
#' stability_results = poisson_BAR_StARS(t(Xsim$X), sequen = sequenBAR)
#' 
#' # Use stability results to select a lambda and fit Poisson BAR graphical model
#' # beta grid specifies stability thresholds by user. Function returns as many 
#' # estimated adjacency matrices as the beta's supplied
#' fitModel(t(Xsim200sf$X), stability_results, sequenBAR, beta=c(1,2), 
#'          regularization = 'l0', lchosen = NA, xitype = "original")
#'          
#' # Alternatively, skip selecting tuning parameter by StARS, fit models on all lambdas
#' # when  lchosen is supplied, the function will ignore stability results and beta grid
#' fitl0 = fitModel(t(Xsim200sf$X), lchosen = sequenBAR, xitype = "original")
#' 
#' @references 
#' Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. In Advances in neural information processing systems (pp. 1432-1440).
#' Wan, Y. W., Allen, G. I., Baker, Y., Yang, E., Ravikumar, P., Anderson, M., & Liu, Z. (2016). XMRF: an R package to fit Markov Networks to high-throughput genetics data. BMC systems biology, 10(3), 69.

fitModel <- function(dt, instabilities, sequen, beta = .05, regularization = "l0", 
                     lchosen = NA, th = NA, first.timer = TRUE, xitype = "original") {
  library(doParallel)
  
  if(is.na(lchosen[1])) {
    Dbar = data.frame(1/sequen, instabilities)
  
    for (i in 2:nrow(Dbar)) {
      if(Dbar[i, 2] < Dbar[(i-1), 2]) {
        Dbar[i,2] = Dbar[(i-1), 2]
      }
    }
  
    lchosen = rep(NA, length(beta))
    for(i in 1:length(beta)) {
      candidates = Dbar[Dbar[, 2] <= beta[i], ]
      if(nrow(candidates) == 0) {
        Lambda.chosen = 1/sequen[1]
      } else {
        Lambda.chosen = max(candidates[, 1])
      }
      lchosen[i] = 1/Lambda.chosen
      if(is.na(Lambda.chosen)) print(paste(i, "th is no good", sep = ""))
    }
    lchosen = lchosen[!is.na(lchosen)]
  }
  
  nlambda = length(lchosen)
  p = ncol(dt)
  n = nrow(dt)
  print(log(n)/2)
  
  if(first.timer == TRUE) {
    # cores = detectCores()
    # cl <- makeCluster((5), type = "PSOCK")
    # registerDoParallel(cl)
    registerDoParallel(4)
  }
  #iterate through features
  stackedA = foreach(j = 1:ncol(dt), .combine = cbind, .packages = c('doParallel', 'glmnet')) %dopar% {
    
    if(regularization == 'l1') {  
      fit = glmnet::glmnet(dt[, -j], dt[, j], family = "poisson", 
                           alpha = 1, lambda = lchosen, intercept = F, standardize = F)
      coeffs = as.matrix(fit$beta)
    } else if(regularization == 'l0') {
      fit = prBARgraph::poisBAR(dt[, j], dt[, -j], lambda = lchosen, xi = log(n)/2, eps = 1E-6, xitype = xitype, max.iter = 3000)
      coeffs = fit$coef
      
      # func = function(lamb) {
      #   fit.temp = prBARgraph::poisBAR(dt[, j], dt[, -j], lambda = lamb)
      #   return(fit.temp$coef)
      # }
      # coeffs = sapply(sequen, func, simplify = "array")
    }
    
    if(is.na(th)) {
      coeffs = ifelse(coeffs == 0, 0, 1)
    } else {
      coeffs = ifelse(coeffs <= th, 0, 1)
    }
    coeffs = as.matrix(coeffs)
    
    # Add the diagnal element (vector) in the vector before binding into a matrix
    # Using an if statement because the indicator j would go out of bound at last iteration
    if(nlambda > 1) {
      if (j > nrow(coeffs)) {
        return(rbind(coeffs, rep(0, nlambda)) )
      } else if(j > 1) {
        return(rbind(coeffs[0:(j-1), ], rep(0, nlambda), coeffs[j:nrow(coeffs), ] ))
      } else if(j == 1) {
        return(rbind(rep(0, nlambda), coeffs))
      }
    } else {
      coeffs = as.numeric(coeffs)
      if (j > length(coeffs)) {
        return(c(coeffs, 0))
      } else if(j > 1) {
        return(c(coeffs[0:(j-1)], 0, coeffs[j:length(coeffs)] ))
      } else if(j == 1) {
        return(c(0, coeffs))
      }
    }
    
  }
  
  A.results <- foreach(j = 1:length(lchosen)) %dopar% {
    
    A = stackedA[ ,seq(j, j+nlambda*(p-1), by = nlambda)]
    if(regularization == 'l1') {
      A = A+t(A)
      #A = A*t(A)
    } else {
      A = A+t(A)
      #A = A*t(A)
    }
    A[lower.tri(A)] = 0
    A = ifelse(!A == 0, 1, 0)
    return(A)
  }

  #stopCluster()
  closeAllConnections()
  
  return(A.results)
  
}

# Make an adjacency matrix upper triangular
UpperTriangulrize <- function(A) {
  A = A + t(A)
  A[lower.tri(A)] = 0
  A = ifelse(!A == 0, 1, 0)
  return(A)
}

ROC <- function(A, estimatedA) {
  # To upper triangular matrix
  A = UpperTriangulrize(A)
  estimatedA = lapply(estimatedA, UpperTriangulrize)
  p = nrow(A)
  
  roc = data.frame(falsePos = rep(0,length(estimatedA)), truePos = rep(0,length(estimatedA)))
  for(i in 1:length(estimatedA)) {
    totalest = p*(p-1)/2
    totalA = sum(A)
    truePos = sum(A*estimatedA[[i]])/totalA
    falsePos = (sum(estimatedA[[i]]) - sum(A*estimatedA[[i]]))/totalest
    #falsePos = (sum(estimatedA[[i]]) - sum(A*estimatedA[[i]]))/sum(estimatedA[[i]])
    #print(paste(truePos, falsePos, sep = " / "))
    roc[i, ] = c(falsePos, truePos)
  }
  return(roc)
}
