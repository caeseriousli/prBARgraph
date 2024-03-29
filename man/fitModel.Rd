% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitModel.R
\name{fitModel}
\alias{fitModel}
\title{Fit BAR Poisson Graphical Model.}
\usage{
fitModel(dt, instabilities, lams, beta = 0.05, regularization = "l0",
  lchosen = NA, th = NA)
}
\arguments{
\item{dt}{any nxp matrix, or data frame. Design matrix to be fitted.}

\item{instabilities}{a numeric vector, if not NA, shoose to use StARS instability for selecting a lambda. Run poisson_BAR_StARS() first to get StARS}

\item{lams}{a numeric vector, to be used with instability a grid of l0/l1 penalization parameters to be chosen from.}

\item{beta}{a number, specifying stability threshold}

\item{regularization}{a character, either default to "l0" (BAR). The function also supports "l1" (LASSO Poisson)}

\item{lchosen}{a number or a vector of penalization (either for l0 or l1) parameters, if not NA, function will ignore StARS (instabilities, lams, and beta) to select Lambda. Instead, run all penalizations and return a list of adjacency matrices}
}
\value{
Either an adjacenty matrix (selected by StARS), or a list of adjacency matrices, corresponding to each lambda supplied by 'lchose' option.
}
\description{
Fit BAR Poisson Graphical Model.
}
\examples{
require(XMRF)
require(doParallel)
# Simulate RNA data from XMRF package (Wan, Y. W., et. al, 2016)
Xsim <- XMRF.Sim(n=100, p=200, model="LPGM", graph.type="hub")

# Specify a grid of l0 penalization parameter
sequenBAR = 2^seq(logb(18, base = 2), logb(1e-4, base = 2), length = 15)

# Run StARS to calculate stability (Liu, H. et. al, 2010)
stability_results = poisson_BAR_StARS(t(Xsim$X), sequen = sequenBAR)

# Use stability results to select a lambda and 
# fit Poisson BAR graphical model
# beta grid specifies stability thresholds by user. 
# Function returns as many estimated adjacency
# matrices as the beta's supplied
fitModel(t(Xsim200sf$X), stability_results, sequenBAR, beta=c(1,2), 
         regularization = 'l0', lchosen = NA, xitype = "original")
         
# Alternatively, skip selecting tuning parameter by StARS, 
# fit models on all lambdas
# when  lchosen is supplied, the function will ignore 
# stability results and beta grid
fitl0 = fitModel(t(Xsim200sf$X), lchosen = sequenBAR, xitype = "original")

}
\references{
Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to 
regularization selection (stars) for high dimensional graphical models. 
In Advances in neural information processing systems (pp. 1432-1440).
Wan, Y. W., Allen, G. I., Baker, Y., Yang, E., Ravikumar, P., Anderson, 
M., & Liu, Z. (2016). XMRF: an R package to fit Markov Networks 
to high-throughput genetics data. BMC systems biology, 10(3), 69.
"Reference to our paper"
}
