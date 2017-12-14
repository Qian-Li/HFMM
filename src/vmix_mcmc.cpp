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
    S = kron(T,pinv(thv.Omega.slice(k))) + kron(Sp, thv.HtH);
    b = kron(diagmat(1.0 / thv.Sigq(k,k) * invp), trans(thv.H)) * thdif;
    // Sample loadings;
    samp = BayesReg(b, S);
    // Matricize;
    sampM = vec2mat(samp, datv.nr, thv.nfac);
    // Update Loading matrix ---------------------------------------------------- /
    thv.Ld.submat(k*datv.nr, (k+1)*datv.nr-1, 0, thv.nfac) = sampM;
    // Update Omega --------------------------------------------------------------/
    W = prv.Om0 + sampM * T * trans(sampM);
    thv.Omega.slice(k) = rWish(W, nu);
  }
  // shrinkage update
  // delta0;----------------------------------------------------------------------/
  {
    double a1 = prv.ap1 + datv.nr * thv.nbase * thv.nfac * 0.5;
    double a2 = 1.0;
    for(i=0; i<thv.nfac; i++){
      vec copy = thv.pen; copy[0] = 1.0;
      tau = cumprod(copy);
      for(j=0; j<thv.nbase; j++){
        mat tempvec = thv.Ld.submat(j*datv.nr, (j+1)*datv.nr-1, i,i);
        mat tempsum = trans(tempvec) * thv.Omega.slice(j) * tempvec;
        a2 += 0.5*tau[j]* tempsum[0];
      }
    }
    thv.pen[0] = Rf_rgamma(a1, 1.0/a2);
  }
  // delta1+;---------------------------------------------------------------------/
  for(l=1; l<thv.nfac; l++){
    double a1 = prv.ap2 + datv.nr * thv.nbase * (thv.nfac-l) * 0.5;
    double a2 = 1.0;
    for(i=0; i<thv.nfac; i++){
      vec copy = thv.pen; copy[l] = 1.0;
      tau = cumprod(copy);
      for(j=0; j<thv.nbase; j++){
        mat tempvec = thv.Ld.submat(j*datv.nr, (j+1)*datv.nr-1, i,i);
        mat tempsum = trans(tempvec) * thv.Omega.slice(j) * tempvec;
        a2 += 0.5*tau[j]* tempsum[0];
      }
    }
    thv.pen[l] = Rf_rgamma(a1, 1.0/a2);
  }
}

/***********************************************************************************/
/* sampleFacCovs()                                                                 */
/***********************************************************************************/
void sampleFacCovs()
{// Sample SigQ, Sigp | Theta, Lambda H
  double  ap, aq;
  mat     lat;
  vec     bp(datv.nr), bq(thv.nbase);
  int     i, j, k;
  lat = thv.Ld * trans(thv.H);                  //pq * n latent RE's
  {
    ap = prv.athe + 0.5 * datv.ns * thv.nbase;
    aq = prv.athe + 0.5 * datv.ns * datv.nr;
    bp.fill(prv.bthe);
    bq.fill(prv.bthe);
    mat sum;
    for(i=0; i<datv.ns; i++){
      vec latvi=conv_to<vec>::from(lat.col(i));
      mat lati = vec2mat(latvi, datv.nr, thv.nbase);
      mat diff = thv.bi.slice(i) - latvi;       // p*q diff
      bp += 0.5 * diagvec(diff * pinv(thv.Sigq) * trans(diff));
      bq += 0.5 * diagvec(trans(diff) * pinv(thv.Sigp) * diff);
    }
  }
  // Sample col covs SigP;---------------------------------------------------------/
  for(j=0; j<datv.nr; j++){thv.Sigp(j,j) = 1.0/Rf_rgamma(ap, 1.0/bp[j]);}
  // Sample row covs SigQ;---------------------------------------------------------/
  for(k=0; k<thv.nbase; k++){thv.Sigq(k,k) = 1.0/Rf_rgamma(aq, 1.0/bq[k]);}
}

/***********************************************************************************/
/* sampleErrCovs()                                                                 */
/***********************************************************************************/
void sampleErrCovs(cube &y)
{// Sample sige | Yi, Theta
  int i, r, t;
  double a, b, v, resb;
  vec  h1;
  mat  resE, bi, mi, betar, S;
  
  b    = prv.be;
  a = (datv.ns*datv.nt-accu(datv.nmiss))/2.0*datv.nr + prv.ae;
  
  // Error precision ---------------------
  for(r=0; r<datv.nr; r++){
    // a
    resE = longslice(y, r) - longslice(thv.fit, r);
    // b
    for(i=0; i<datv.ns; i++){
      for(t=0; t<datv.nt; t++){
        b += (datv.miss(i,r,t)==0)? std::pow(resE(i,t),2.0) : 0.0;
      }
    }
  }
  //sample residual covariances;----------------------------------------------------/
  thv.sige = 1.0/Rf_rgamma(a, 1.0/b);
}

/***********************************************************************************/
/* sampleBetas()                                                                   */
/***********************************************************************************/
void sampleBetas() // needs to be figured out
{// Sample M = B'xi; B| Sigma, Theta, X--------------------------------------------------------*****************************
  
}

/***********************************************************************************/
/* sampleFactors()                                                                 */
/***********************************************************************************/
void sampleFactors(cube &y, mat const& Bs)
{// Sample factors | everything else;
  int i;
  mat res, Q, iSig, Qth;
  vec b, bi;
  
  iSig = pinv(kron(thv.Sigp, thv.Sigq) + kron(eye(datv.nr, datv.nr), thv.iBtB));
  for(i=0; i<datv.ns; i++){
    Q = trans(thv.Ld) * iSig * thv.Ld + eye(thv.nfac, thv.nfac);
    //residual of pq
    res = (y.slice(i) - thv.mbi.slice(i) * trans(Bs)) * Bs * thv.iBtB;
    b = prv.eta0 + trans(thv.Ld) * iSig * conv_to<colvec>::from(vectorise(res));
    // Sample latent Factors;--------------------------------------------------------/
    thv.H.row(i) = conv_to<rowvec>::from(BayesReg(b,Q));
    // Update thv.mbi!!!!! --------------------------------------------------------*****************************
    
    // Sample coefficient matrices --------------------------------------------------/
    Qth = (1.0/thv.sige * kron(eye(datv.nr, datv.nr), thv.BtB))
      + pinv(kron(thv.Sigp, thv.Sigq));
    thv.bi.slice(i) = vec2mat(BayesReg(vectorise(thv.mbi.slice(i)), Qth), 
                 datv.nr, thv.nbase);
    // Update thv.fit!!!!--------------------------------------------------------*****************************
    
  }
  // Update H'H;
  thv.HtH = trans(thv.H) * thv.H;
}

// Thinking on what to record
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

/***********************************************************************************/
/* OutputSample()                                                                  */
/***********************************************************************************/
void OutputSampleV()
{
  int j, r, p;
  // mat invS;
  // { // Population Splines --------------------------------
  //   std::stringstream ss;
  //   for(p=0; p<datv.np; p++){
  //     for(r=0; r<dat.nr; r++){
  //       for(j=0; j<th.nbase; j++){
  //         ss << th.beta(p,j,r) << " ";
  //       }
  //     }
  //   }
  //   ss << "\n"; out1 << ss.str();
  // }
  { // Precisions ----------------------------------------
    std::stringstream ss;
    ss << thv.sige << " ";          // Error covs;
    for(r=0; r<datv.nr; r++){
      ss << thv.Sigp(r,r) << " ";   // Spatial res covs;
    }
    for(j=0; j<thv.nbase; j++){
      ss << thv.Sigq(j,j) << " ";   // Temp res covs;
    }
    for(p=0; p<thv.nfac; p++){
      ss << thv.pen(p) << " ";      // penalty;
    }
    ss << "\n"; outv3 << ss.str();
  }
  { // Factors -------------------------------------------
    std::stringstream ss;
    for(r=0; r<datv.ns; r++){
      for(j=0; j<thv.nfac; j++){
        ss << thv.H(r,j) << " "; // Spatial matrix
      }
    }
    ss << "\n"; outv4 << ss.str();
  }
}

// Vectorized HFMM
/***********************************************************************************/
/* vmix_mcmc()                                                                     */
/***********************************************************************************/

// [[Rcpp::export]]
void vmix_mcmc(arma::cube y, 
               arma::mat const& X, 
               arma::mat const& Bs,
               int const& nfac,
               int const& burnin, int const& nsim, int const& thin
               )
{
  // vmix var declaration:
  int i, j, r;                // Indices
//   mat W, D, D2;               // AR matrices
//   vec rho;                    // regularization eigenvalues
  // Get initialization summaries ---------------------------------------------------
  datv.ns = y.n_rows;           // number of subjects
  datv.nr = y.n_cols;           // number of regions
  datv.nt = y.n_slices;         // number of segments
  datv.np = X.n_cols;           // number of covariates (groups)
  //
  thv.nbase = Bs.n_cols;        // number of basis functions
  thv.nfac  = nfac;             // number of latent factors
  //
  /* ----------------------------------------------------------------------------- */
  /* Initialize Chain and static summaries                                         */
  /* ----------------------------------------------------------------------------- */
  //
  // Prior parameters -----------------------------------------------
  prv.ae = 0.00001;         prv.be = 0.00001;             // Error variance
  prv.athe = 0.001;         prv.bthe = 0.001;             // coeff res row/col var
  prv.ap1= 1.5;             prv.ap2= 1.5;                 // penalties
  prv.aOm= 1.0;                                           // --------------could be less informative
  prv.Om0= eye<mat>(datv.nr, datv.nr) * 1.0;              // spatial smoothing
  prv.eta0 = zeros<vec>(thv.nfac);                        // Latent factor priors

//   // Smoothing Matrix (AR-1) ----------------------------------------//needed for Beta's
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

  // Initialize missing data ----------------------------------------
  {
    vec rn1, yir;
    datv.miss = icube(datv.ns, datv.nr, datv.nt);
    datv.nmiss= vec(datv.ns);
    yir = randu<vec>(datv.ns);
    for(i=0; i<datv.ns; i++){
      for(r=0; r<datv.nr; r++){
        yir = y.tube(i,r);
        for(j=0; j<datv.nt; j++){
          datv.miss(i,r,j) = 0;
          if( yir(j)==12345.12345 ){
            datv.miss(i,r,j) = 1;
            rn1 = randn(1)*0.5;
            y(i,r,j) = meanNA(yir) + rn1(0);
          }
        }
      }
      datv.nmiss[i] = accu(datv.miss.tube(0,0));
    }
  }
  // RANDOM initialization                            // ----------------------
  thv.H   = randu<mat>(datv.ns, nfac);                // Latent factors
  thv.Ld  = randu<mat>(datv.nr*thv.nbase, nfac);      // Loading matrix
  thv.HtH = trans(thv.H) * thv.H;
  thv.bi  = randu<cube>(datv.ns, datv.nr, thv.nbase); // basis functions
  thv.mbi = randu<cube>(datv.ns, datv.nr, thv.nbase); // mean bases funcs
  thv.fit = randu<cube>(datv.ns, datv.nr, datv.nt);   // current fit
  thv.beta= randu<cube>(datv.np, thv.nbase, datv.nr); // Group means
  thv.Omega=randu<cube>(datv.nr, datv.nr, thv.nbase); // Spatial correlations
  thv.BtB = trans(Bs)*Bs;
  thv.iBtB= pinv(thv.BtB);
  thv.sige= 1.0/Rf_rgamma(prv.ae, 1.0/prv.be);        // Residual covs
  thv.Sigp= eye(datv.ns, datv.ns);
  thv.Sigq= eye(thv.nbase, thv.nbase);
  // Meaningful initialization:                       // ----------------------
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
  // Open output files ----------------------------------------------
  outv1.open(fileOutput1);
  outv2.open(fileOutput2);
  outv3.open(fileOutput3);
  outv4.open(fileOutput4);
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */

  for(int rep=0; rep<(nsim+burnin+1); rep++){
    // mcmc samples --------------------
    completeY2(y, Bs);
    sampleLoadings();
    sampleFacCovs();
    sampleErrCovs(y);
    sampleBetas();
    sampleFactors(y, Bs);
    // Store mcmc samples --------------
    if( (rep > burnin) && (rep%thin == 0) ){
      OutputSampleV();           // Write mcmc samples to file
      // ErgodicAverages(Bs);      // Update ergodic averages
    }
  }
  // Close output files ---------------------------------------------
  outv1.close(); outv2.close(); outv3.close(); outv4.close();

//   // Output to R --------------------------------------------------------------------
//   return List::create(
//     Named("fit")      = hat.fit,
//     Named("coef")     = hat.beta,
//     Named("coef2")    = hat.beta2);
}

// // END vmix_mcmc --------------------------------------------------------------------