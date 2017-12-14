#ifndef __EEGPAL_H__
#define __EEGPAL_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <fstream>

using namespace std;
using namespace arma;
using namespace Rcpp;


/************************************************************************************/
/* Data Structure                                                                   */
/************************************************************************************/
struct dta{
  icube miss; // Indicator cube for missing observations
  int ns;     // Number of subjects
  int nt;     // Number of time segments
  int np;     // Number of covariates (groups)
  int nr;     // Number of brain regions
  vec nmiss;  // Number of missing observations
};


/************************************************************************************/
/* Parameters Structure                                                             */
/************************************************************************************/
struct pars{
    // Model parameters -------------------------------------------------------------
    cube bi;        // Individual mean coefficients
    cube mbi;       // prior mean of bi
    cube beta;      // Population means coefficients
    mat S;          // Spatial Covariance
    mat He;         // Error Variance
    double lbeta;   // Population smoothing
    double lbi;     // Subject-level smoothing
    // Static Summaries -------------------------------------------------------------
    int nbase;      // Number of basis functions
    mat Df;         // Smoothing Matrix
    mat BtB;        // Quadratic basis product
    mat XtX;        // Qadratic regression design
    mat XtXi;       // Inverse of Quadratic design matrix
    cube fit;       // Current subject-level fit
    mat cholXX;     // Cholesky of X'X
    mat cholDf;     // Cholesky of Df
};
struct parsV{
   // Model parameters for vectorized model -----------------------------------------
   cube bi;         // Individual coefficient
   cube mbi;        // Individual mean coefficient (fixed)
   cube beta;       // Coefficients mean;
   mat Ld;          // Latent Factor loading
   mat H;           // Latent Factors
   mat Sigp;        // Regional variance, diagonal
   mat Sigq;        // Temporal variance, diagonal
   cube Omega;      // Spatial covariance, free
   vec pen;         // Column covariance penalty on Loadings
   double sige;     // System error variances
   mat HtH;         // H'H
   // Static Summaries --------------------------------------------------------------
   int nfac;        // Number of factors;
   int nbase;       // Number of basis functions;
   cube fit;        // Current sub-level fit
};

/************************************************************************************/
/* Prior Structure                                                                  */
/************************************************************************************/
struct prior{
  // Priors for hierarchical models
  // Error variance
  double ae;
  double be;

  // subject level smoothing
  double ab;
  double bb;

  // spatial smoothing
  double aS;
  mat W;

};

struct priorV{
  // Prior for vectorized models
  // Error variance
  double ae;
  double be;
  
  // Coefficients variances (Sigp and Sigq)
  double athe;
  double bthe;
  
  // Penalties shapes, the first and the rest
  double ap1;
  double ap2;
  
  // spatial covs (IW)
  double aOm;
  mat Om0;
};

/************************************************************************************/
/* Estimate Structure                                                               */
/************************************************************************************/
struct ergodics{
   cube fit;    // Posterior mean fit
   cube beta;   // Posterior mean population regression functions
   cube beta2;  // Second posterior moment of population regression functions
   double n;    // number of samples used in mc estimators
};


// Exposed utility functions -----------------------------------------------------
//
// Bayesian Regression and Triangular Systems ------------------------------------
vec BayesReg(vec const &b, mat const &Q);
vec BayesRegL(vec const &b, mat const &L);
void trisolve(vec &solution, mat const &L, vec const &x, bool const &lower);
// Multivariate Samplers
mat rWish(mat const &S, double v);
mat rMN(mat const &m, mat const &S, mat const& V);
// Other utilities
mat cov2cor(mat S);
double meanNA(vec &y);
mat vec2mat(vec const &v, int const& nr, int const& nc);
//
// slicing of Armadillo cubes ---------------------------
mat longslice(cube &A, int r);
mat flatslice(cube &A, int r);
//
// shared functions among different implementations
void completeY(cube &y, mat const& Bs);

#endif
