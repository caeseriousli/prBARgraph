#include <math.h>
#include <Rmath.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#define LEN sizeof(double)


// get sum of log factorials for loglik
double getSumLogFactorial(int y)
{
  double val = 0;
  for(int i = 1; i < y; i++)
  {
    val += log(i + 1);
  }
  return val;
}

// get log-likelihood for Poisson model
double getLogLikelihood(double *y, double *eta, int nin)
{
  const int n = nin;
  double likli = 0;
  for (int i = 0; i < n; i++)
  {
    //likli += (y[i] * eta[i] - exp(eta[i]) - getSumLogFactorial(y[i]));
    likli += (y[i] * eta[i] - exp(eta[i]));
  }
  return likli / n;
}

double sgn(double z) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  return(s);
}
//Standardize design matrix
SEXP standardize(SEXP X_) {
  // Declarations
  int n = nrows(X_);
  int p = ncols(X_);
  SEXP XX_, c_, s_;
  PROTECT(XX_ = allocMatrix(REALSXP, n, p));
  PROTECT(c_ = allocVector(REALSXP, p));
  PROTECT(s_ = allocVector(REALSXP, p));
  double *X = REAL(X_);
  double *XX = REAL(XX_);
  double *c = REAL(c_);
  double *s = REAL(s_);

  for (int j = 0; j < p; j++) {

    // Center (Calculate mean and subtract)
    c[j] = 0;
    for (int i = 0; i < n; i++) {
      c[j] += X[j * n + i];
    }
    c[j] = c[j] / n;
    for (int i = 0; i < n; i++) XX[j * n + i] = X[j * n + i] - c[j];

    // Scale (Calculate sdev and divide)
    s[j] = 0;
    for (int i = 0; i < n; i++) {
      s[j] += pow(XX[j * n + i], 2);
    }
    s[j] = sqrt(s[j] / n);
    for (int i = 0; i < n; i++) XX[j * n + i] = XX[j * n + i] / s[j];
  }

  // Return list
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, XX_); // Standardized design matrix
  SET_VECTOR_ELT(res, 1, c_); // Mean
  SET_VECTOR_ELT(res, 2, s_); // Standard deviations
  UNPROTECT(4);
  return(res);
}

// Weighted cross product of y with jth column of x
double getGradient(double *X, double *y, double *eta, int n, int j) {
  int nn = n * j;
  double val = 0;
  for (int i = 0; i < n; i++) val += X[nn + i] * (y[i] - exp(eta[i]));
  return(val);
}

// Weighted sum of squares of jth column of X
double getHessian(double *X, double *eta, int n, int j) {
  int nn = n * j;
  double val = 0;
  for (int i = 0; i < n; i++) val += exp(eta[i]) * pow(X[nn + i], 2);
  return(val);
}


// Criterion for convergence: All coefficients must pass the following |(b_new - b_old) / b_old| < eps
int checkConvergence(double *beta, double *beta_old, double eps, int p, int l) {
  int converged = 1;
  double allsum = 0;
  for (int j = 0; j < p; j++) {
    if (fabs((beta[l * p + j] - beta_old[j]) / beta_old[j]) > eps) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

///////////////////////////////////////////////////////////////////////////////////////
SEXP cleanupCRR(double *a, double *eta, double *diffBeta,
                SEXP beta, SEXP Dev, SEXP iter, SEXP linpred, SEXP tstp) {
  Free(a);
  Free(eta);
  Free(diffBeta);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta); //coefficient estimates
  SET_VECTOR_ELT(res, 1, Dev); //deviance = -2*loglik
  SET_VECTOR_ELT(res, 2, iter); //iterations until convergence
  SET_VECTOR_ELT(res, 3, linpred); //hessian
  SET_VECTOR_ELT(res, 4, tstp);
  UNPROTECT(6);
  return(res);
}

//////////////////////////////////////////////////////////////////////////////////
//start cordinate descent

//x_ = design matrix
//y_ = outcome
//lambda = tuning parameter
//esp_ = epsilon (thershold)
//max_iter_ = max iterations

SEXP ccd_ridge(SEXP x_, SEXP y_, SEXP lambda_,
               SEXP esp_, SEXP max_iter_, SEXP multiplier) {

  //Declaration
  int n = length(y_);
  int p = length(x_) / n;
  double nullDev;

  //Output
  SEXP res, beta, Dev, iter, linpred, tst;
  PROTECT(beta = allocVector(REALSXP, p));
  double *b = REAL(beta);
  for (int j = 0; j < p; j++) b[j] = 0;
  PROTECT(Dev = allocVector(REALSXP, 1));
  for (int i = 0; i < 1; i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, 1));
  for (int i = 0; i < 1; i++) INTEGER(iter)[i] = 0;
  PROTECT(linpred = allocVector(REALSXP, n));
  double *lp = REAL(linpred);
  for (int i = 0; i <  n; i++) lp[i] = 0;
  PROTECT(tst = allocVector(REALSXP, p));
  double *tstp = REAL(tst);
  for (int i = 0; i <  p; i++) tstp[i] = 0;

  //Intermediate quantities for internal use
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0;
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *diffBeta = Calloc(p, double);
  for (int j = 0; j < p; j++) diffBeta[j] = 1;

  double grad, hess, l1, shift, si, delta;
  int converged;

  //Pointers
  double *x = REAL(x_);
  double *y = REAL(y_);
  double lam = REAL(lambda_)[0];
  double esp = REAL(esp_)[0];
  double *m = REAL(multiplier);
  int max_iter = INTEGER(max_iter_)[0];


  //end of declaration;

  //initialization
  nullDev = -2 * getLogLikelihood(y, eta, n); // Calculate null deviance at beta = 0

  for (int j = 0; j < p; j++) a[j] = b[j];

  while (INTEGER(iter)[0] < 1000) {
    if (REAL(Dev)[0] - nullDev > 0.99 * nullDev) break;

    INTEGER(iter)[0]++;


    // calculate gradient & hessian andupdate beta_j
    for (int j = 0; j < p; j++) {
      grad = -getGradient(x, y, eta, n, j); // jth component of gradient [l'(b)]
      hess = getHessian(x, eta, n, j); // jth component of hessian [l''(b)]
      l1 = m[j] * lam ; //divide by n since we are minimizing the following: -(1/n)l(beta) + lambda * p(beta)

      delta =  -(grad / n + a[j] * l1 / n) / (hess / n + l1);
      //u   = grad / n + (hess / n) * a[j]; // z in paper
      //v   = hess / n;


      // Do one dimensional ridge update.
      // Employ trust region as in Genkin et al. (2007) for quadratic approximation.
      b[j] = a[j] + sgn(delta) * fmin(fabs(delta), diffBeta[j]);
      diffBeta[j] = fmax(2 * fabs(delta), diffBeta[j] / 2);
      //b[j] = a[j] + (xwr / n - a[j] * l1) / (xwx / n + l1);

      // Update r
      shift = b[j] - a[j];
      if (shift != 0) {
        for (int i = 0; i < n; i++) {
          si = shift * x[j * n + i];
          eta[i] += si;
        }
      } //end shift

    } //for j = 0 to (p - 1)
    // Check for convergence
    converged = checkConvergence(b, a, esp, p, 0);
    for (int i = 0; i < p; i++)
      a[i] = b[i];

    //Calculate deviance
    REAL(Dev)[0] = -2 * getLogLikelihood(y, eta, n);

    for (int i = 0; i < n; i++){
      lp[i] = eta[i];
    }
    if (converged)  break;
    //for converge
  } //for while loop

  res = cleanupCRR(a, eta, diffBeta, beta, Dev, iter, linpred, tst);
  return(res);
}

// BAR:
//////////////////////////////////////////////////////////////////////////////////
double newBarL0(double h, double g, double b, double l) {
  double tmp;
  double s = 0;
  tmp = h * b + g;
  if (tmp > 0) s = 1;
  else if (tmp < 0) s = -1;
  if (fabs(tmp) < 2 * sqrt(h * l)) return(0);
  else return((tmp + s * sqrt(pow(tmp, 2) - 4 * l * h)) / (2 * h));
}

//start cordinate descent

//x_ = design matrix
//y_ = outcome
//lambda = tuning parameter
//esp_ = epsilon (thershold)
//max_iter_ = max iterations

// Right now, changes are: disabling deviance check
SEXP ccd_bar(SEXP x_, SEXP y_, SEXP lambda_,
               SEXP esp_, SEXP max_iter_, SEXP beta0_) {

  //Declaration
  int n = length(y_);
  int p = length(x_) / n;
  int L = length(lambda_);
  double nullDev;

  //Output
  SEXP res, beta, Dev, iter, linpred, tst;
  PROTECT(beta = allocVector(REALSXP, L * p));
  double *b = REAL(beta);
  for (int j = 0; j < (L * p); j++) b[j] = 0;
  PROTECT(Dev = allocVector(REALSXP, 1));
  for (int i = 0; i < 1; i++) REAL(Dev)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i = 0; i < L; i++) INTEGER(iter)[i] = 0;
  PROTECT(linpred = allocVector(REALSXP, n));
  double *lp = REAL(linpred);
  for (int i = 0; i <  n; i++) lp[i] = 0;
  PROTECT(tst = allocVector(REALSXP, p));
  double *tstp = REAL(tst);
  for (int i = 0; i <  p; i++) tstp[i] = 0;

  //Intermediate quantities for internal use
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j = 0; j < p; j++) a[j] = 0;
  double *eta = Calloc(n, double);
  for (int i = 0; i < n; i++) eta[i] = 0;
  double *diffBeta = Calloc(p, double);
  for (int j = 0; j < p; j++) diffBeta[j] = 0;

  double grad, hess, l1, u, v, shift, si;
  double *m = REAL(beta0_);

  int converged;

  //Pointers
  double *x = REAL(x_);
  double *y = REAL(y_);
  double *lam = REAL(lambda_);
  //double lam = REAL(lambda_)[0];
  double esp = REAL(esp_)[0];
  int max_iter = INTEGER(max_iter_)[0];


  //end of declaration;
  
  //Outer loop for each lambda
  for(int l = 0; l < L; l++) {
  
    for (int i = 0; i < n; i++) eta[i] = 0;
    //initialization
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < p; j++) {
        eta[i] += m[j] * x[j * n + i];
      } 
    }
    
    // CAESAR Initialize eta and nullDev
    nullDev = -2 * getLogLikelihood(y, eta, n); // Calculate null deviance at beta = 0
    for (int j = 0; j < p; j++) a[j] = m[j];
  
  
    while (INTEGER(iter)[l] < max_iter) {
      //if (REAL(Dev)[0] - nullDev > 0.99 * nullDev) break;
  
      INTEGER(iter)[l]++;
  
  
      // calculate gradient & hessian andupdate beta_j
      for (int j = 0; j < p; j++) {
        grad = getGradient(x, y, eta, n, j); // jth component of gradient [l'(b)]
        hess = getHessian(x, eta, n, j); // jth component of hessian [l''(b)]
        //u   = grad / n + (hess / n) * a[j]; // z in paper
        //v   = hess / n;
  
        // Update b_j
        l1 = lam[l] ; //divide by n since we are minimizing the following: -(1/n)l(beta) + lambda * p(beta)
  
        //Do one dimensional ridge update.
        b[l * p + j] = newBarL0(hess / n, grad / n, a[j], l1);
        
        //tstp[j] = 0;
        //for (int i = 0; i < n; i++) tstp[j] += eta[i];
        //tstp[j] = eta[i];
        //if( newBarL0(hess / n, grad / n, a[j], l1) > 1E+05 ) {
        //  tstp[j] = grad;
        //  tstp[j+1] = hess;
        //}
  
        // Update r
        shift = b[l * p + j] - a[j];
        if (shift != 0) {
          for (int i = 0; i < n; i++) {
            si = shift * x[j * n + i];
            eta[i] += si;
          }
        } //end shift
        
        tstp[j] = hess;
  
      } //for j = 0 to (p - 1)
      // Check for convergence
      //CAESAR added iteration number l
      converged = checkConvergence(b, a, esp, p, l);
  
      //CAESAR update
      //for (int i = 0; i < n; i++) eta[i] = 0;
      //for (int i = 0; i < n; i++) {
      //  for (int j = 0; j < p; j++) {
      //    eta[i] += b[l * p + j] * x[j * n + i];
      //  }
      //}
  
      // CAESAR check
      //for (int j = 0; j < p; j++) {
      //  tstp[j] = fabs( (b[j] - a[j]) / a[j]) > esp;
      //}
  
      for (int i = 0; i < p; i++)
        a[i] = b[l * p + i];
  
      //Calculate deviance
      REAL(Dev)[0] = -2 * getLogLikelihood(y, eta, n);
  
      //if (INTEGER(iter)[0] <= p)
        //tstp[INTEGER(iter)[0]] = ( REAL(Dev)[0] - nullDev ) / nullDev;
        //tstp[INTEGER(iter)[0]] = REAL(Dev)[0];
  
      for (int i = 0; i < n; i++){
        lp[i] = eta[i];
      }
      
      if (converged) break;
      
      //for converge
    } //for while loop
    
  } //lambda loop
  
  res = cleanupCRR(a, eta, diffBeta, beta, Dev, iter, linpred, tst);
  return(res);
}


