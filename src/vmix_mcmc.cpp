#include "EEGpal.h"

/***********************************************************************************/
/* FILES & CONSTANTS ***************************************************************/
/***********************************************************************************/
#define fileOutput1 "post_v/bi.txt"               //coeffs per subject
#define fileOutput2 "post_v/beta.txt"             //betas
#define fileOutput3 "post_v/precision.txt"        //tau's
#define fileOutput4 "post_v/factors.txt"          //latent factors


/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
dta datv;                   // Data summaries
parsV thv;                  // Parameters
// ergodics hat;              // Ergodic summaries
priorV prv;                  // Prior parameters
vec yvec, bvec;             // Global subj/region data/pars
//
// // Utility constants
// double w1;
// double w2;
// Output files ------------------------------------------------------------------
ofstream outv1;              // B-spline coeffs bi's
ofstream outv2;              // betas (group mean coeffs)
ofstream outv3;              // precisions
ofstream outv4;              // latent factors

/***********************************************************************************/
/* completeY()                                                                     */
/***********************************************************************************/
void completeY2(cube &y, mat const& Bs)
{
  int i, r, j;
  vec yn;
  //
  for(r=0; r<datv.nr; r++){
    for(i=0; i<datv.ns; i++){
      yvec = y.tube(i,r);
      for(j=0; j<datv.nt; j++){
        if(datv.miss(i,r,j) == 1){
          yn = randn(1);
          y(i,r,j) = thv.fit(i,r,j) + yn(0) * std::pow(thv.sige, -0.5);
        }
      }
    }
  }
}

/***********************************************************************************/
/* sampleLoadings()                                                                */
/***********************************************************************************/
void sampleLoadings()
{// Sample Lambdam Omega, delta | Theta, M, SigP, SigQ, Fac
  int k, l, i, j;
  vec lfvec, tau, invp, b, thdif, samp;
  mat S, Sp, sampM, W;
  double tt;
  //
  //
  double nu = prv.aOm + thv.nfac;
  // ---- Global* vars
  cube cdif = thv.bi - thv.mbi;
  // invp = Sigma_p^-1
  invp = 1.0/diagvec(thv.Sigp);
  // tau = prod pen's
  tau = 1.0 / cumprod(thv.pen);
  mat T = diagmat(tau);
  for(k=0; k<thv.nbase; k++){
    thdif = vectorise(longslice(cdif, k));
    Sp = diagmat(1.0 / thv.Sigq(k,k) * invp);
    S = kron(T,inv(thv.Omega.slice(k))) + kron(Sp, thv.HtH);
    b = kron(diagmat(1.0 / thv.Sigq(k,k) * invp), thv.H.t()) * thdif;
    // Sample loadings;
    samp = BayesReg(b, S);
    // Matricize;
    sampM = vec2mat(samp, datv.nr, thv.nfac);
    // Update Loading matrix --------------------------------------------------------
    thv.Ld.submat(k*datv.nr, (k+1)*datv.nr-1, 0, thv.nfac) = sampM;
    // Update Omega -----------------------------------------------------------------
    W = prv.Om0 + sampM * T * sampM.t();
    thv.Omega.slice(k) = rWish(W, nu);
  }
  // shrinkage update
  // delta0;
  {
    double a1 = prv.ap1 + datv.nr * thv.nbase * thv.nfac * 0.5;
    double a2 = 1.0;
    for(i=0; i<thv.nfac; i++){
      vec copy = thv.pen; copy[0] = 1.0;
      tau = cumprod(copy);
      for(j=0; j<thv.nbase; j++){
        mat tempvec = thv.Ld.submat(j*datv.nr, (j+1)*datv.nr-1, i,i);
        mat tempsum = tempvec.t() * thv.Omega.slice(j) * tempvec;
        a2 += 0.5*tau[j]* tempsum[0];
      }
    }
    thv.pen[0] = Rf_rgamma(a1, a2);
  }
  // delta1+;
  for(l=1; l<thv.nfac; l++){
    double a1 = prv.ap2 + datv.nr * thv.nbase * (thv.nfac-l) * 0.5;
    double a2 = 1.0;
    for(i=0; i<thv.nfac; i++){
      vec copy = thv.pen; copy[l] = 1.0;
      tau = cumprod(copy);
      for(j=0; j<thv.nbase; j++){
        mat tempvec = thv.Ld.submat(j*datv.nr, (j+1)*datv.nr-1, i,i);
        mat tempsum = tempvec.t() * thv.Omega.slice(j) * tempvec;
        a2 += 0.5*tau[j]* tempsum[0];
      }
    }
    thv.pen[l] = Rf_rgamma(a1, a2);
  }
}

/***********************************************************************************/
/* sampleFacCovs()                                                                 */
/***********************************************************************************/
void sampleFacCovs()
{// Sample SigQ, Sigp | Theta, Lambda H
  // Sample row covs SigQ;

  // Sample col covs SigP;
}

/***********************************************************************************/
/* sampleErrCovs()                                                                 */
/***********************************************************************************/
void sampleErrCovs()
{// Sample sige | Yi, Theta
  //sample residual covariances;
}

/***********************************************************************************/
/* sampleBetas()                                                                   */
/***********************************************************************************/
void sampleBetas()
{// Sample M = B'xi; B| Sigma, Theta, X

}

/***********************************************************************************/
/* sampleFactors()                                                                 */
/***********************************************************************************/
void sampleFactors()
{// Sample factors | everything else;
  //sample latent factors;

  //sample coefficient mats;
}


// /***********************************************************************************/
// /* ErgodicAverages()                                                               */
// /***********************************************************************************/
// void ErgodicAverages(mat const& Bs)
// {
//   int r;
//   mat beta, betat;
//
//   w1 = 1.0/(1.0 + hat.n);
//   w2 = 1.0 - w1;
//   hat.n++;
//   // Average fit ------------------------------------------
//   hat.fit  = w1*th.fit + w2*hat.fit;
//   // Average regression coefficients ----------------------
//   for(r=0; r<dat.nr; r++){
//     beta  = th.beta.slice(r);                                          // np x nbase
//     betat = Bs * trans(beta);                                          // nt x np
//     hat.beta.slice(r)  = w1*betat + w2*hat.beta.slice(r);              // nt x np
//     hat.beta2.slice(r) = w1*betat%betat + w2*hat.beta2.slice(r);       // nt x np
//   }
// }
//
// /***********************************************************************************/
// /* OutputSample()                                                                  */
// /***********************************************************************************/
// void OutputSample()
// {
//   int j, r, p;
//   mat invS;
//   { // Population Splines --------------------------------
//     std::stringstream ss;
//     for(p=0; p<dat.np; p++){
//       for(r=0; r<dat.nr; r++){
//         for(j=0; j<th.nbase; j++){
//           ss << th.beta(p,j,r) << " ";
//         }
//       }
//     }
//     ss << "\n"; out1 << ss.str();
//   }
//   { // Precisions ----------------------------------------
//     std::stringstream ss;
//     for(r=0; r<dat.nr; r++){
//       out2 << th.He(r,r) << " "; // Error precision
//     }
//     ss << th.lbi << " ";         // Subject smoothing
//     ss << "\n"; out2 << ss.str();
//   }
//   { // Spatial -------------------------------------------
//     std::stringstream ss;
//     invS = (inv(th.S+0.0));
//     for(r=0; r<dat.nr; r++){
//       for(j=r; j<dat.nr; j++){
//         ss << invS(r,j) << " "; // Spatial matrix
//       }
//     }
//     ss << "\n"; out3 << ss.str();
//   }
// }

// Vectorized HFMM
/***********************************************************************************/
/* vmix_mcmc()                                                                     */
/***********************************************************************************/

// [[Rcpp::export]]
void vmix_mcmc(arma::cube y, 
               arma::mat const& X, 
               arma::mat const& Bs,
               int const& burnin, int const& nsim, int const& thin
               )
{
//   int i, j, r;                // Indices
//   mat W, D, D2;               // AR matrices
//   vec rho;                    // regularization eigenvalues
//   // Get initialization summaries ---------------------------------------------------
//   dat.ns = y.n_rows;           // number of subjects
//   dat.nr = y.n_cols;           // number of regions
//   dat.nt = y.n_slices;         // number of segments
//   dat.np = X.n_cols;           // number of covariates (groups)
//   //
//   th.nbase = Bs.n_cols;        // number of basis functions
//   //
//
//   /* ----------------------------------------------------------------------------- */
//   /* Initialize Chain and static summaries                                         */
//   /* ----------------------------------------------------------------------------- */
//
//   // Prior parameters -----------------------------------------------
//   pr.ae = 0.00001;         pr.be = 0.00001;             // Error variance
//   pr.ab = 0.9*10.0 + 1.0;  pr.bb = 0.9;                 // Subject-level smoothing
//   pr.aS = 1.0;
//   pr.W = eye<mat>(dat.nr, dat.nr) * 1.0;                // spatial smoothing
//
//   // Smoothing Matrix (AR-1) ----------------------------------------
//   W  = eye<mat>(th.nbase, th.nbase)*1.0;
//   D  = eye<mat>(th.nbase, th.nbase)*3.0;
//   D2 = eye<mat>(th.nbase, th.nbase)*1.0;
//   // AR1 adjecency -----------------------------------
//   for(i=1; i<th.nbase; i++){ W(i,(i-1)) = 1.0; }
//   for(i=0; i<(th.nbase-1); i++){ W(i,(i+1)) = 1.0; }
//   // AR1 no. neighbors -------------------------------
//   D(0,0) = D((th.nbase-1),(th.nbase-1)) = 2.0;
//   // Regularizing eigenvalue -------------------------
//   for(i=0; i<th.nbase; i++) D2(i,i) = std::pow(D(i,i), -0.5);
//   rho = randu<vec>(th.nbase);
//   eig_sym(rho, D2*W*D2);
//   // Final penalization matrix -----------------------
//   th.Df = (D - rho[(th.nbase-2)]*W)*10.0;
//
//   // Initialize missing data ----------------------------------------
//   vec rn1;
//   dat.miss = icube(dat.ns, dat.nr, dat.nt);
//   dat.nmiss= vec(dat.ns);
//   yir = randu<vec>(dat.ns);
//   for(i=0; i<dat.ns; i++){
//     for(r=0; r<dat.nr; r++){
//       yir = y.tube(i,r);
//       for(j=0; j<dat.nt; j++){
//         dat.miss(i,r,j) = 0;
//         if( yir(j)==12345.12345 ){
//           dat.miss(i,r,j) = 1;
//           rn1 = randn(1)*0.5;
//           y(i,r,j) = meanNA(yir) + rn1(0);
//         }
//       }
//     }
//     dat.nmiss[i] = accu(dat.miss.tube(0,0));
//   }
//   // Individual functional means and smoothing ----------------------
//   th.lbi = 10.0;                                  // precision
//   th.bi  = randu<cube>(dat.ns, dat.nr, th.nbase); // basis functions
//   th.fit = randu<cube>(dat.ns, dat.nr, dat.nt);   // current fit
//   th.BtB = trans(Bs)*Bs;
//   bir    = randu<vec>(th.nbase);
//   for(i=0; i<dat.ns; i++){
//     for(r=0; r<dat.nr; r++){
//       yir = y.tube(i,r);
//       bir = solve((th.BtB + th.lbi*th.Df), trans(Bs)*yir);
//       th.bi.tube(i,r) = bir;
//       th.fit.tube(i,r) = Bs*bir;
//     }
//   }
//   // Error Precision and Spatial Covariance -------------------------
//   double err;
//   int nn1 = dat.ns*dat.nt;
//   th.He   = zeros<mat>(dat.nr, dat.nr);
//   for(r=0; r<dat.nr; r++){
//     err = 0.0;
//     for(i=0; i<dat.ns; i++){
//       for(j=0; j<dat.nt; j++){
//         err += pow(y(i,r,j) - th.fit(i,r,j), 2.0)/nn1;
//       } // j
//     } //i
//     th.He(r,r) = 5.0/err;
//   }
//   th.S = eye<mat>(dat.nr, dat.nr)*1.0;
//   // Population means means and smoothing ---------------------------
//   th.lbeta = 0.01/(dat.ns*dat.nr);
//   th.beta  = randu<cube>(dat.np, th.nbase, dat.nr);
//   th.XtX   = trans(X)*X;
//   mat br   = randu<mat>(dat.ns, th.nbase);
//   for(r=0; r<dat.nr; r++){
//     br               = longslice(th.bi, r);
//     th.beta.slice(r) = solve(th.XtX , trans(X)*br);
//   };
//   th.mbi   = randu<cube>(dat.ns, dat.nr, th.nbase);
//
//   // Initialize and decompose regression matrices -------------------
//   th.XtXi   = pinv(th.XtX);           // (X'X)^-1
//   th.cholXX = chol(th.XtXi);          // Cholesky of (X'X)^-1 (upper)
//   th.cholDf = chol(pinv(th.Df));      // Cholesky of Df^-1 (upper)
//
//   // Initialize ergodics --------------------------------------------
//   hat.n     = 0.0;
//   hat.beta  = randu<cube>(dat.nt, dat.np, dat.nr);
//   hat.beta2 = randu<cube>(dat.nt, dat.np, dat.nr);
//   hat.fit   = randu<cube>(dat.ns, dat.nr, dat.nt);
//
//   // Open output files ----------------------------------------------
//   out1.open(fileOutput1);
//   out2.open(fileOutput2);
//   out3.open(fileOutput3);
//
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */

  for(int rep=0; rep<(nsim+burnin+1); rep++){
    // mcmc samples --------------------
    completeY2(y, Bs);
    sampleLoadings();
    sampleFacCovs();
    sampleErrCovs();
    sampleBetas();
    sampleFactors();
    // Store mcmc samples --------------
    // if( (rep > burnin) && (rep%thin == 0) ){
    //   OutputSample();           // Write mcmc samples to file
    //   ErgodicAverages(Bs);      // Update ergodic averages
    // }
  }
//   // Close output files ---------------------------------------------
//   out1.close(); out2.close(); out3.close();
//
//   // Output to R --------------------------------------------------------------------
//   return List::create(
//     Named("fit")      = hat.fit,
//     Named("coef")     = hat.beta,
//     Named("coef2")    = hat.beta2);
}

// // END vmix_mcmc --------------------------------------------------------------------