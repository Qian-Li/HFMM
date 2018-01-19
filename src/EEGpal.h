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
struct modelpars{
  // Model parameters shared by both V&S implementations
  cube bi;         // Individual coefficient
  cube mbi;        // Individual mean coefficient (fixed)
  cube beta;       // Coefficients mean;
  double taue;     // System error variances;
  mat taup;        // coef residual precision, diagonal
  mat tauq;        // coef residual precision, diagonal
  // Static Summaries --------------------------------------------------------------
  int nbase;       // Number of basis functions;
  cube fit;        // Current sub-level fit
  mat BtB;         // Bs'Bs;
  mat iBtB;        // pinv(Bs'Bs);
  mat XtXi;        // pinv(X'X)
  mat X;           // design matrix X
};
struct parsV{
  // Parameters specific to vectorized implementation
  cube Ld;         // Latent Factor loading
  mat H;           // Latent Factors
  cube Omega;      // Spatial covariance, free
  vec pen;         // Column covariance penalty on Loadings
  mat HtH;         // H'H
  // static
  int nfac;        // Number of factors;
};
struct parsS{
  // Parameters specific to sandwich implementation
  mat LL;          // Latent Factor loading, Left
  mat LR;          // Latent Factor loading, Right
  cube H;          // Latent Factors
  mat tauL;        // Precision, Left
  mat tauR;        // Precision, Right
  vec penL;        // Penalty, Left
  vec penR;        // Penalty, Right
  // mat HtH;         // H'H
  // Static Summaries --------------------------------------------------------------
  int nfacL;       // Number of factors;
  int nfacR;       // Number of factors;  
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
  
  // Latent factors;
  vec eta0;
};
struct priorS{
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
  
  // Loading precision
  double nuL;
  double nuR;
  
  // Latent factors;
  vec eta0;
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

struct ergodicsV{
  cube fit;     // Posterior mean fit, (mbi + bi) * Bs'
  cube beta;    // Group mean fit, (mbi) * Bs'
  cube beta2;   // Spline coefficients
  cube betai;   // Individual mean fit (mbi + Lambda eta) * Bs'
  // cube L;       // E(lf)
  double n;     // sampel included for ergodics sum
};

// Exposed utility functions -----------------------------------------------------
//
// Bayesian Regression and Triangular Systems ------------------------------------
vec BayesReg(vec const &b, mat const &Q);
vec BayesRegL(vec const &b, mat const &L);
void trisolve(vec &solution, mat const &L, vec const &x, bool const &lower);
//
// Multivariate Samplers
mat rWish(mat const &S, double v);
mat rMN(mat const &m, mat const &S, mat const& V);
// Other utilities
mat cov2cor(mat S);
double meanNA(vec &y);
mat vec2mat(vec &v, int const& nr, int const& nc);
cube mat2cube(mat const &m, int const& nr, int const& nc);
void vec2long(cube &a, vec const &v, int const& j);
cube cubexmat(cube &c, mat const&m);
cube matxcube(mat const& L, cube const&C, mat const&R);
mat cube2mat(cube const &c);
mat safechol(mat const& S);
//
// slicing of Armadillo cubes ---------------------------
mat longslice(cube &A, int r);
mat flatslice(cube &A, int r);
//
#endif
