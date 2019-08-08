#' @export
#' @useDynLib prBARgraph ccd_ridge
poisBAR <- function(y, x, lambda = 0, xi = 5, delta = 0,
                   eps = 1E-6, tol = 1E-6,
                   max.iter = 1000, xitype = "fixed"){
  # Caesar added, idea from elastic net
  #if(xi == 0) xi = 10*lambda
  #if(xi == 0) xi = ifelse(lambda > 1, lambda, 
  if(xitype == "fixed") {
    ridgename = "ccd_ridge0"
    barname = "ccd_bar0"
    #barname = "ccd_bar_hybrid"
  } else if (xitype == "timesLam") {
    #if(!xi == 0) xi = ifelse(lambda > 5, 50*lambda, 200)
    xi = xi/lambda
    #xi = lambda
    ridgename = "ccd_ridge"
    barname = "ccd_bar"
  } else if(xitype == "original") {
    ridgename = "ccd_ridge0"
    barname = "ccd_hybrid_ridge"
  } else stop("Error: xitype has to be either fixed or timesLam.")
  #####################################
  
  ## Error checking
  #if(xi < 0) stop("xi must be a non-negative number.")
  if(max.iter < 1) stop("max.iter must be positive integer.")
  if(eps <= 0) stop("eps must be a positive number.")
  if(tol <= 0) stop("tol must be a positive number.")
  p <- dim(x)[2]

  # Fit the PSH Ridge Model here w/ tuning parameter xi
  ridgeFit   <- .Call(ridgename, x, as.numeric(y),
                      xi, eps, as.integer(1000),
                      penalty.factor = rep(1, p), PACKAGE = "prBARgraph")
  ridgeCoef  <- ridgeFit[[1]] #Divide coeff estimates by sdev
  ridgeIter  <- ridgeFit[[3]]

  # ridgeFit <- glmnet::glmnet(x, as.numeric(y), family = "poisson",
  #                            alpha = 0, lambda = xi, intercept = F)
  # ridgeCoef <- as.double(ridgeFit$beta)
  # ridgeIter <- 5

  #Results to store:
  btmp <- ridgeFit[[1]]
  #btmp <- ridgefit2$beta[,1]
  #ridgeCoef <- btmp
  # btmp = ridgeCoef
  barFit <- .Call(barname,  x, as.numeric(y),
                  as.vector(lambda), as.double(eps), as.integer(max.iter),
                  as.vector(btmp), PACKAGE = "prBARgraph")
  
  # turn long beta veto into matrix
  betas = barFit[[1]]
  if(length(lambda) > 1) betas = matrix(betas, ncol = length(lambda))

  ## Output
  val <- structure(list(coef = betas,
                        logLik = barFit[[2]] / -2,
                        iter = barFit[[3]],
                        ridgeCoef = ridgeCoef,
                        ridgeIter = ridgeIter,
                        xi = xi,
                        tst = barFit[5],
                        call = sys.call()),
                   class = "prBARgraph")

  val
}
