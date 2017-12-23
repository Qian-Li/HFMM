#include "EEGpal.h"
//
/***********************************************************************************/
/* FILES & CONSTANTS ***************************************************************/
/***********************************************************************************/
#define fileOutput1 "post_v/beta.txt"             //fixed-effects related coeffs
#define fileOutput2 "post_v/precision.txt"        //sigP,sigq and sige^-1
// #define fileOutput3 "post_v/spatial"              //Omega's
#define fileOutput4 "post_v/factors.txt"          //factors and penalties
//
/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
dta       datv;                   // Data summaries
parsV     thv;                    // Parameters
ergodicsV hatv;                   // Ergodic summaries
priorV    prv;                    // Prior parameters
vec       vbi;                    // vbi as vec(bi)
mat       Sinv;                   // Sigma^{-1} of size pq*pq
double    g = 100.0;              // Mean Prior constant, set large
//
// Utility constants
//
vec bit, yit;
double ww1, ww2;                    // For Ergodic updates
// Output files ------------------------------------------------------------------
ofstream outv1;
ofstream outv2;
// ofstream outv3;
ofstream outv4;
//
/***********************************************************************************/
/* completeY2()                                                                    */
/***********************************************************************************/
void completeY2(cube &y)
{
  int i, r, j;
  vec yn(1);
  //
  for(r=0; r<datv.nr; r++){
    for(i=0; i<datv.ns; i++){
      for(j=0; j<datv.nt; j++){
        if(datv.miss(i,r,j) == 1){
          yn = randn(1);
          y(i,r,j) = thv.fit(i,r,j) + yn[0] * std::pow(thv.taue, -0.5);
        }
      }
    }
  }
}
//
/***********************************************************************************/
/* sampleLoadings()                                                                */
/***********************************************************************************/
void sampleLoadings()
{// Sample Lambda, Omega, delta | Theta, M, SigP, SigQ, Fac
  int k, l;
  vec samp, tau, b, copy;
  mat Sp, sampM, W, T, sch;
  bool sinv;
  // cumsum for penalty update
  vec bg = 1.0 * ones<vec>(thv.nfac);
  //
  double nu = prv.aOm + thv.nfac;
  // tau = prod pen's
  tau = cumprod(thv.pen);
  T = diagmat(tau);
  for(k=0; k<thv.nbase; k++){
    vbi = vectorise(longslice(thv.bi, k));
    Sp  = thv.taup * thv.tauq(k,k);
    Sinv= kron(T, thv.Omega.slice(k));
    Sinv= Sinv + kron(thv.HtH, Sp);
    b   = kron(trans(thv.H), Sp) * vbi;
    // Sample loadings, pm * 1
    sinv= chol(sch, Sinv, "lower");
    Rcout << "ok 2.1" << std::endl;
    if(!sinv){
      thv.pen.print("penalty:");
      T.print("T:");
      thv.Omega.slice(k).print("Omega:");
      Sp.print("Sp:");
    }
    samp = BayesRegL(b, sch);
    Rcout << "ok 2.2" << std::endl;
    // Matricize; p*m
    sampM = vec2mat(samp, datv.nr, thv.nfac);
    Rcout << "ok 2.3" << std::endl;
    // Update Loading matrix ---------------------------------------------------- /
    vec2long(thv.Ld, samp, k);
    Rcout << "ok 2.4" << std::endl;
    // Update Omega --------------------------------------------------------------/
    W = prv.Om0 + sampM * T * trans(sampM);
    thv.Omega.slice(k) = rWish(W, nu);
    Rcout << "ok 2.5" << std::endl;
    // Update cumsum for penalties
    for(l=0; l<thv.nfac; l++){
      copy = thv.pen;     copy[l] = 1.0;
      tau = cumprod(copy);    T = diagmat(tau);
      bg[l] += 0.5 * trace(T.submat(l,l,thv.nfac-1,thv.nfac-1) * 
          trans(sampM.cols(l, thv.nfac-1)) * 
          thv.Omega.slice(k) * sampM.cols(l, thv.nfac-1));
    }
  }
  // shrinkage(penalty) update
  // delta0;----------------------------------------------------------------------/
  {
    double a1 = prv.ap1 + datv.nr * thv.nbase * thv.nfac * 0.5;
    thv.pen[0] = Rf_rgamma(a1, 1.0/bg[0]);
  }
  // delta1+;---------------------------------------------------------------------/
  for(l=1; l<thv.nfac; l++){
    double a1 = prv.ap2 + datv.nr * thv.nbase * (thv.nfac-l) * 0.5;
    thv.pen[l] = Rf_rgamma(a1, 1.0/bg[l]);
  }
}
//
/***********************************************************************************/
/* sampleFacCovs()                                                                 */
/***********************************************************************************/
void sampleFacCovs()
{// Sample SigQ, Sigp | Theta, Lambda H
  double  ap, aq;
  cube    res;
  vec     bp(datv.nr), bq(thv.nbase);
  int     i, j, k;
  //
  res = thv.bi - cubexmat(thv.Ld, thv.H);               //Theta_i-M_i-mat(Lambda*H)
  {
    ap = prv.athe + 0.5 * datv.ns * thv.nbase;
    aq = prv.athe + 0.5 * datv.ns * datv.nr;
    bp.fill(prv.bthe);
    bq.fill(prv.bthe);
    mat sump, sumq;
    for(i=0; i<datv.ns; i++){
      sump = res.slice(i) * thv.tauq * trans(res.slice(i));
      sumq = trans(res.slice(i)) * thv.taup * res.slice(i);
      bp += 0.5 * diagvec(sump);
      bq += 0.5 * diagvec(sumq);
    }
  }
  // Sample col covs SigP;---------------------------------------------------------/
  for(j=0; j<datv.nr; j++){thv.taup(j,j) = Rf_rgamma(ap, 1.0/bp[j]);}
  // Sample row covs SigQ;---------------------------------------------------------/
  for(k=0; k<thv.nbase; k++){thv.tauq(k,k) = Rf_rgamma(aq, 1.0/bq[k]);}
}
//
/***********************************************************************************/
/* sampleErrCovs()                                                                 */
/***********************************************************************************/
void sampleErrCovs(cube &y)
{// Sample sige | Yi, Theta
  int i, r, t;
  double a, b;
  vec  h1;
  mat  resE, bi, mi, betar, S;
  
  b    = prv.be;
  a    = prv.ae + (datv.ns*datv.nt-accu(datv.nmiss)) * 0.5 * datv.nr;
  
  // Error sum of squares ---------------------
  for(r=0; r<datv.nr; r++){
    // a
    resE = longslice(y, r) - longslice(thv.fit, r);
    // b
    for(i=0; i<datv.ns; i++){
      for(t=0; t<datv.nt; t++){
        b += 0.5 * ((datv.miss(i,r,t)==0)? std::pow(resE(i,t),2.0) : 0.0);
      }
    }
  }
  //sample residual covariances;----------------------------------------------------/
  thv.taue = Rf_rgamma(a, 1.0/b);
}
//
/***********************************************************************************/
/* sampleBetas()                                                                   */
/***********************************************************************************/
void sampleBetas() // needs to be figured out
{// Sample M = B'xi; B| Sigma, Theta, X
  mat S, V, Theta, Bn, Bsamp;
  mat Lmat = cube2mat(thv.Ld);        //pq*m
  Theta = cube2mat(thv.mbi+thv.bi);   //pq*n
  //
  S = g / (g+1.0) * thv.XtXi;
  V = kron(pinv(thv.tauq), pinv(thv.taup));
  V = V + Lmat * trans(Lmat);
  //
  Bn = S * trans(thv.X) * trans(Theta); 
  //sample population means: beta;-------------------------------------------------/
  Bsamp = rMN(Bn, S, V);                                 // pq*d
  thv.beta = mat2cube(trans(Bsamp), datv.nr, thv.nbase); // p*q*d
  // Update mbi based on beta------------------------------------------------------/
  thv.mbi = cubexmat(thv.beta, thv.X);                   // p*q*n
}

/***********************************************************************************/
/* sampleFactors()                                                                 */
/***********************************************************************************/
void sampleFactors(cube &y, mat const& Bs)
{// Sample factors | everything else;
  int i, r, t;
  mat lmat, res, Q, Stinv, Sipq, fit, qch, sch;
  cube LFcube;
  vec b, bi, thi, bsamp;
  bool decom1, decom2;
  //
  lmat = cube2mat(thv.Ld);
  LFcube= cubexmat(thv.Ld, thv.H);
  Sinv = kron(pinv(thv.tauq), pinv(thv.taup));
  Sinv = Sinv + kron((thv.iBtB), 1.0/thv.taue*eye<mat>(datv.nr, datv.nr));
  Q = trans(lmat) * pinv(Sinv) * lmat + 1.0 * eye<mat>(thv.nfac, thv.nfac); ////could be improved by performing CHOL here
  decom1 = chol(qch, Q, "lower");
  Sipq = kron(thv.taup, thv.tauq);
  Stinv= Sipq + kron(thv.taue * eye<mat>(datv.nr, datv.nr), thv.BtB);
  decom2 = chol(sch, Stinv, "lower");
  if(!decom1){Q.print();}
  if(!decom2){Stinv.print();}
  //
  for(i=0; i<datv.ns; i++){
    //residual of pq
    res = (flatslice(y,i) * Bs * thv.iBtB - thv.mbi.slice(i));
    b   = prv.eta0 + trans(lmat) * pinv(Sinv) * conv_to<colvec>::from(vectorise(res));
    // Sample latent Factors;--------------------------------------------------------/
    bsamp   = BayesRegL(b,qch);
    thv.H.row(i) = conv_to<rowvec>::from(bsamp);
    // residual of Theta_i
    bi  = vectorise(thv.taue * trans(Bs) * trans(flatslice(y,i)));
    bi  = bi + Sipq * vectorise(trans(thv.mbi.slice(i) + LFcube.slice(i)));
    // Sample coefficient matrices --------------------------------------------------/
    thi = BayesRegL(bi, sch);
    // Update bi = Theta_i - mbi
    thv.bi.slice(i) = trans(vec2mat(thi, thv.nbase, datv.nr)) - thv.mbi.slice(i);
    fit = (thv.bi.slice(i)+thv.mbi.slice(i)) * trans(Bs);       //p*nt      
    // Update fit_i = Theta_i * Bs'
    for(r=0; r<datv.nr; r++){
      for(t=0; t<datv.nt; t++){
        thv.fit(i,r,t) = fit(r,t);
      }
    }
  }
  // Update H'H;
  thv.HtH = trans(thv.H) * thv.H;
}
//
/***********************************************************************************/
/* ErgodicAverages()                                                               */
/***********************************************************************************/
void ErgodicAveragesV(mat const& Bs)
{
  int i, r, t;
  mat beta1, beta2;
  vec tj, tji;
  //
  ww1 = 1.0/(1.0 + hatv.n);
  ww2 = 1.0 - ww1;
  hatv.n++;
  //
  // Average fit ------------------------------------------
  hatv.fit  = ww1*thv.fit + ww2*hatv.fit;
  hatv.beta2= ww1*thv.beta+ ww2*hatv.beta2;
  hatv.L    = ww1*thv.Ld + ww2*hatv.L;
  cube lcube= cubexmat(thv.Ld, thv.H);
  //
  for(i=0; i<datv.ns; i++){
    beta1 = thv.mbi.slice(i) * trans(Bs);                       //p*nt
    beta2 = (thv.mbi.slice(i) + lcube.slice(i)) * trans(Bs);    //p*nt
    for(r=0; r<datv.nr; r++){
      for(t=0; t<datv.nt; t++){
        hatv.beta(i,r,t) = ww1*beta1(r,t) + ww2*hatv.beta(i,r,t);
        // hatv.beta2(i,r,t)= w1*std::pow(beta1(r,t), 2.0) + w2*hatv.beta2(i,r,t);
        hatv.betai(i,r,t)= ww1*beta2(r,t) + ww2*hatv.betai(i,r,t);
      }
    }
  }
}

/***********************************************************************************/
/* OutputSample()                                                                  */
/***********************************************************************************/
void OutputSampleV()
{
  int j, r, p;
  { // Population Splines --------------------------------
    std::stringstream ss;
    for(p=0; p<datv.np; p++){
      for(r=0; r<datv.nr; r++){
        for(j=0; j<thv.nbase; j++){
          ss << thv.beta(r,j,p) << " ";
        }
      }
    }
    ss << "\n"; outv1 << ss.str();
  }
  { // Precisions ----------------------------------------
    std::stringstream ss;
    ss << thv.taue << " ";          // Error covs;
    for(r=0; r<datv.nr; r++){
      ss << thv.taup(r,r) << " ";   // Spatial res covs;
    }
    for(j=0; j<thv.nbase; j++){
      ss << thv.tauq(j,j) << " ";   // Temp res covs;
    }
    ss << "\n"; outv2 << ss.str();
  }
  { // Factors -------------------------------------------
    std::stringstream ss;
    for(r=0; r<datv.ns; r++){
      for(j=0; j<thv.nfac; j++){
        ss << thv.H(r,j) << " "; // Spatial matrix
      }
    }
    for(p=0; p<thv.nfac; p++){
      ss << thv.pen(p) << " ";      // penalty;
    }
    ss << "\n"; outv4 << ss.str();
  }
}
//
// Vectorized HFMM
/***********************************************************************************/
/* vmix_mcmc()                                                                     */
/***********************************************************************************/
// [[Rcpp::export]]
List vmix_mcmc(arma::cube y, 
               arma::mat const& X, 
               arma::mat const& Bs,
               int const& nfac,
               int const& burnin, int const& nsim, int const& thin
               )
{
  // vmix var declaration:
  int i, j, r;               // Indices
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
  prv.ae = 0.0001;         prv.be = 0.0001;               // Error variance
  prv.athe = 0.001;         prv.bthe = 0.001;             // coeff res row/col var
  prv.ap1= 1.5;             prv.ap2= 1.5;                 // penalties
  prv.aOm= 1.0;                                           // --------------could be less informative
  prv.Om0= eye<mat>(datv.nr, datv.nr) * 1.0;              // spatial smoothing
  prv.eta0 = zeros<vec>(thv.nfac);                        // Latent factor priors
  //
  // Initialize missing data ----------------------------------------
  {
    vec rn1, yit;
    datv.miss = icube(datv.ns, datv.nr, datv.nt);
    datv.nmiss= vec(datv.ns);
    yit = randu<vec>(datv.ns);
    for(i=0; i<datv.ns; i++){
      for(r=0; r<datv.nr; r++){
        yit = y.tube(i,r);
        for(j=0; j<datv.nt; j++){
          datv.miss(i,r,j) = 0;
          if( yit(j)==12345.12345 ){
            datv.miss(i,r,j) = 1;
            rn1 = randn(1)*0.5;
            y(i,r,j) = meanNA(yit) + rn1(0);
          }
        }
      }
      datv.nmiss[i] = accu(datv.miss.tube(0,0));
    }
  }
  // RANDOM initialization                            // ----------------------
  thv.H   = randu<mat>(datv.ns, nfac);                // Latent factors
  thv.Ld  = randu<cube>(datv.nr, thv.nbase, nfac);    // Loading matrix
  thv.HtH = trans(thv.H) * thv.H;
  thv.bi  = randu<cube>(datv.nr, thv.nbase, datv.ns); // basis functions
  thv.mbi = randu<cube>(datv.nr, thv.nbase, datv.ns); // mean bases funcs
  thv.fit = randu<cube>(datv.ns, datv.nr, datv.nt);   // current fit
  thv.beta= randu<cube>(datv.nr, thv.nbase, datv.np); // Group means
  thv.Omega=randu<cube>(datv.nr, datv.nr, thv.nbase); // Spatial correlations
  thv.BtB = trans(Bs)*Bs;
  thv.iBtB= pinv(thv.BtB);
  thv.taue= Rf_rgamma(prv.ae, 1.0/prv.be);        // Residual covs
  thv.taup= 0.001*eye(datv.nr, datv.nr);
  thv.tauq= 0.001*eye(thv.nbase, thv.nbase);
  thv.pen = randu<vec>(nfac);
  for(i=0; i<nfac; i++){thv.pen[i] = Rf_rgamma(prv.ap1, 1.0/prv.ap2);}
  for(j=0; j<thv.nbase; j++){thv.Omega.slice(j) = eye<mat>(datv.nr, datv.nr);}
  // Meaningful initialization:                       // ----------------------
  bit    = randu<vec>(thv.nbase);
  for(i=0; i<datv.ns; i++){
    for(r=0; r<datv.nr; r++){
      yit = y.tube(i,r);
      bit = solve((thv.BtB + thv.taue * eye<mat>(thv.nbase, thv.nbase)), trans(Bs)*yit);
      for(j=0; j<thv.nbase; j++){
        thv.bi(r,j,i) = bit[j];
      }
      thv.fit.tube(i,r) = vectorise(Bs*bit);
    }
  }
  //
  // Initialize and decompose regression matrices -------------------
  thv.X     = X;                      // X
  thv.XtXi  = pinv(trans(X)*X);       // pinv(X'X)
  // Initialize ergodics --------------------------------------------
  hatv.n     = 0.0;
  hatv.beta  = randu<cube>(datv.ns, datv.nr, datv.nt);
  hatv.beta2 = randu<cube>(datv.nr, thv.nbase, datv.np);
  hatv.betai = randu<cube>(datv.ns, datv.nr, datv.nt);
  hatv.L     = randu<cube>(datv.nr, thv.nbase, thv.nfac);
  hatv.fit   = randu<cube>(datv.ns, datv.nr, datv.nt);

  // Open output files ----------------------------------------------
  outv1.open(fileOutput1);
  outv2.open(fileOutput2);
  // outv3.open(fileOutput3);
  outv4.open(fileOutput4);
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */

  for(int rep=0; rep<(nsim+burnin+1); rep++){
    // mcmc samples --------------------
    Rcout << "ok 1" << std::endl;
    completeY2(y);
    Rcout << "ok 2" << std::endl;
    sampleLoadings();
    Rcout << "ok 3" << std::endl;
    sampleFacCovs();
    Rcout << "ok 4" << std::endl;
    sampleErrCovs(y);
    Rcout << "ok 5" << std::endl;
    sampleBetas();
    Rcout << "ok 6" << std::endl;
    sampleFactors(y, Bs);
    Rcout << "ok 7" << std::endl;
    // Store mcmc samples --------------
    if( (rep > burnin) && (rep%thin == 0) ){
      OutputSampleV();           // Write mcmc samples to file
      ErgodicAveragesV(Bs);      // Update ergodic averages
    }
  }
  // Close output files ---------------------------------------------
  outv1.close(); outv2.close(); outv4.close(); // outv3.close();

  // Output to R --------------------------------------------------------------------
  return List::create(
    Named("fit")      = hatv.fit,
    Named("coef")     = hatv.beta2,
    Named("gmean")    = hatv.beta,
    Named("imean")    = hatv.betai);
}

// // END vmix_mcmc --------------------------------------------------------------------