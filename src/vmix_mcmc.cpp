#include "EEGpal.h"
//
/***********************************************************************************/
/* FILES & CONSTANTS ***************************************************************/
/***********************************************************************************/
#define fileOutput1 "post_v/beta.txt"             //fixed-effects related coeffs  -V
#define fileOutput2 "post_v/precision.txt"        //sigP,sigq and sige^-1         -V
#define fileOutput3 "post_v/factors.txt"          //factors and penalties         -V
#define fileOutput4 "post_s/beta.txt"             //fixed-effects related coeffs  -S
#define fileOutput5 "post_s/precision.txt"        //sigP,sigq and sige^-1         -S
#define fileOutput6 "post_s/factors.txt"          //factors and penalties         -S
//
/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
dta       datn;                   // Data summaries                             -V -S
modelpars pars;                   // Common model parameters                    -V -S
parsV     thv;                    // Parameters                                 -V
parsS     ths;                    // Parameters                                 -S
ergodicsV hatn;                   // Ergodic summaries                          -V -S
priorV    prv;                    // Prior parameters                           -V
priorS    prs;                    // Prior parameters                           -S
//
// Utility constants ---------------------------------------------------------------/
//
vec       vbi;                    // vbi as vec(bi)
mat       Sinv;                   // Sigma^{-1} of size pq*pq
double    g = 100.0;              // Mean Prior constant, set large
vec       bit, yit;
double    ww1, ww2;               // For Ergodic updates
int       i,r,j,k,l,t,s;
// Rcout << "ok 2.3" << std::endl; //Debugger
// Output files ------------------------------------------------------------------
ofstream outv1;
ofstream outv2;
ofstream outv3;
ofstream f("log.txt");
//
/***********************************************************************************/
/* completeY2()                                                                    */
/***********************************************************************************/
void completeY2(cube &y)
{
  vec yn(1);
  //
  for(r=0; r<datn.nr; r++){
    for(i=0; i<datn.ns; i++){
      for(j=0; j<datn.nt; j++){
        if(datn.miss(i,r,j) == 1){
          yn.randn();
          y(i,r,j) = pars.fit(i,r,j) + yn[0] * std::pow(pars.taue, -0.5);
        }
      }
    }
  }
}
//
/***********************************************************************************/
/* sampleLoadingsV()                                                                */
/***********************************************************************************/
void sampleLoadingsV()
{// Sample Lambda, Omega, delta | Theta, M, SigP, SigQ, Fac
  vec     samp, tau, b, copy;
  mat     Sp, sampM, W, T, sch;
  // cumsum for penalty update
  vec     bg = 1.0 * ones<vec>(thv.nfac);
  //
  double  nu = prv.aOm + thv.nfac;
  // tau = prod pen's
  tau     = cumprod(thv.pen);
  T       = diagmat(tau);
  for(k=0; k<pars.nbase; k++){
    vbi   = vectorise(longslice(pars.bi, k));
    Sp    = pars.taup * pars.tauq(k,k);
    Sinv  = kron(T, thv.Omega.slice(k));
    Sinv  +=kron(thv.HtH, Sp);
    b     = kron(trans(thv.H), Sp) * vbi;
    // Sample loadings, pm * 1
    //---------------------Safe Cholesky-------------------------//
    sch   = safechol(Sinv);
    samp  = BayesRegL(b, sch);
    // Matricize; p*m
    sampM = vec2mat(samp, datn.nr, thv.nfac);
    // Update Loading matrix ---------------------------------------------------- /
    vec2long(thv.Ld, samp, k);
    // Update Omega --------------------------------------------------------------/
    W     = prv.Om0 + sampM * T * trans(sampM);
    thv.Omega.slice(k) = rWish(W, nu);
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
    double a1   = prv.ap1 + datn.nr * pars.nbase * thv.nfac * 0.5;
    thv.pen[0]  = Rf_rgamma(a1, 1.0/bg[0]);
  }
  // delta1+;---------------------------------------------------------------------/
  for(l=1; l<thv.nfac; l++){
    double a1   = prv.ap2 + datn.nr * pars.nbase * (thv.nfac-l) * 0.5;
    thv.pen[l]  = Rf_rgamma(a1, 1.0/bg[l]);
  }
}
//
/***********************************************************************************/
/* sampleLoadingsS()                                                                */
/***********************************************************************************/
void sampleLoadingsS()
{ // Sample LL LR, tauL tauR, penL penR | bi, taup tauq 
  // Variable declaration
  vec     pL(ths.nfacL), pR(ths.nfacR), sampL, sampR, ssL(ths.nfacL), ssR(ths.nfacR);
  mat     TL, TR, schL, schR;
  double  ga, gb;
  // Tau as prod of pen's
  pL    = cumprod(ths.penL); TL = diagmat(pL);
  pR    = cumprod(ths.penR); TR = diagmat(pR);
  // ------ Left Loading ------
  for(j=0; j<datn.nr; j++){
    vec temp = pL % conv_to<vec>::from(ths.tauL.row(j));
    vec b = zeros<vec>(ths.nfacL);
    mat S = zeros<mat>(ths.nfacL, ths.nfacL);
    for(i=0; i<datn.ns; i++){
      b += ths.H.slice(i) * ths.LR.t() * pars.tauq * trans(pars.bi.slice(i).row(j));
      S += ths.H.slice(i) * ths.LR.t() * pars.tauq * ths.LR * ths.H.slice(i).t();
    }
    b = b * pars.taup(j,j);
    S = S * pars.taup(j,j) + diagmat(temp);
    //---------------------Safe Cholesky-------------------------//
    schL   = safechol(S);
    sampL = BayesRegL(b, schL);
    // Update loading matrix by rows;
    ths.LL.row(j) = conv_to<rowvec>::from(sampL);
    // Update loading precisions;
    for(r=0; r<ths.nfacL; r++){
      ga = 0.5 * (prs.nuL + 1.0);
      gb = 0.5 * (prs.nuL + pL[r] * std::pow(ths.LL(j,r), 2.0));
      ths.tauL(j,r) = Rf_rgamma(ga, 1.0/gb);
    }
  }
  // ------ Left Penalties ------
  ssL = conv_to<vec>::from(sum((ths.tauL % ths.tauL),0));
  for(r=0; r<ths.nfacL; r++){
    ga = (r==0) ? prs.ap1 : prs.ap2;
    ga += 0.5 * datn.nr * (ths.nfacL - r);
    vec tp = ths.penL;  tp[r] = 1.0;
    pL     = cumprod(tp);
    gb = 1.0;
    for(k=r; k<ths.nfacL; k++){
      gb += 0.5 * pL[k] * ssL[k];
    }
    ths.penL[r] = Rf_rgamma(ga, 1.0/gb);
  }
  //
  // ------ Right Loading ------
  for(j=0; j<pars.nbase; j++){
    vec temp = pR % conv_to<vec>::from(ths.tauR.row(j));
    vec b = zeros<vec>(ths.nfacR);
    mat S = zeros<mat>(ths.nfacR, ths.nfacR);
    for(i=0; i<datn.ns; i++){
      b += ths.H.slice(i).t() * ths.LL.t() * pars.taup * pars.bi.slice(i).col(j);
      S += ths.H.slice(i).t() * ths.LL.t() * pars.taup * ths.LL * ths.H.slice(i);
    }
    b = b * pars.tauq(j,j);
    S = S * pars.tauq(j,j) + diagmat(temp);
    //---------------------Safe Cholesky-------------------------//
    schR   = safechol(S);
    sampR  = BayesRegL(b, schR);
    // Update loading matrix by rows;
    ths.LR.row(j) = conv_to<rowvec>::from(sampR);
    // Update loading precisions;
    for(r=0; r<ths.nfacR; r++){
      ga = 0.5 * (prs.nuR + 1.0);
      gb = 0.5 * (prs.nuR + pR[r] * std::pow(ths.LR(j,r), 2.0));
      ths.tauR(j,r) = Rf_rgamma(ga, 1.0/gb);
    }
  }
  // ------ Reft Penalties ------
  ssR = conv_to<vec>::from(sum((ths.tauR % ths.tauR),0));
  for(r=0; r<ths.nfacR; r++){
    ga = (r==0) ? prs.ap1 : prs.ap2;
    ga += 0.5 * datn.nr * (ths.nfacR - r);
    vec tp = ths.penR;  tp[r] = 1.0;
    pR     = cumprod(tp);
    gb = 1.0;
    for(k=r; k<ths.nfacR; k++){
      gb += 0.5 * pR[k] * ssR[k];
    }
    ths.penR[r] = Rf_rgamma(ga, 1.0/gb);
  }
}
//
/***********************************************************************************/
/* sampleFacCovs()                                                                 */
/***********************************************************************************/
void sampleFacCovs(bool const& v)
{// Sample SigQ, Sigp | Theta, Lambda H
  double  ap, aq;
  cube    res;
  vec     bp(datn.nr), bq(pars.nbase);
  //
  if(v){
    res = pars.bi - cubexmat(thv.Ld, thv.H);               //Theta_i-M_i-mat(Lambda*H)
  } else {
    res = pars.bi - matxcube(ths.LL, ths.H, ths.LR);
  }
  {
    ap = prv.athe + 0.5 * datn.ns * pars.nbase;
    aq = prv.athe + 0.5 * datn.ns * datn.nr;
    bp.fill(prv.bthe);
    bq.fill(prv.bthe);
    mat sump, sumq;
    for(i=0; i<datn.ns; i++){
      sump = res.slice(i) * pars.tauq * trans(res.slice(i));
      sumq = trans(res.slice(i)) * pars.taup * res.slice(i);
      bp += 0.5 * diagvec(sump);
      bq += 0.5 * diagvec(sumq);
    }
  }
  // Sample col covs SigP;---------------------------------------------------------/
  for(j=0; j<datn.nr; j++){pars.taup(j,j) = Rf_rgamma(ap, 1.0/bp[j]);}
  // Sample row covs SigQ;---------------------------------------------------------/
  for(k=0; k<pars.nbase; k++){pars.tauq(k,k) = Rf_rgamma(aq, 1.0/bq[k]);}
}
//
/***********************************************************************************/
/* sampleErrCovs()                                                                 */
/***********************************************************************************/
void sampleErrCovs(cube &y)
{// Sample sige | Yi, Theta
  double a, b;
  vec  h1;
  mat  resE, bi, mi, betar, S;
  
  b    = prv.be;
  a    = prv.ae + (datn.ns*datn.nt-accu(datn.nmiss)) * 0.5 * datn.nr;
  
  // Error sum of squares ---------------------
  for(r=0; r<datn.nr; r++){
    // a
    resE = longslice(y, r) - longslice(pars.fit, r);
    // b
    for(i=0; i<datn.ns; i++){
      for(t=0; t<datn.nt; t++){
        b += 0.5 * ((datn.miss(i,r,t)==0)? std::pow(resE(i,t),2.0) : 0.0);
      }
    }
  }
  //sample residual covariances;----------------------------------------------------/
  pars.taue = Rf_rgamma(a, 1.0/b);
}
//
/***********************************************************************************/
/* sampleBetas()                                                                   */
/***********************************************************************************/
void sampleBetas() // needs to be figured out
{// Sample M = B'xi; B| Sigma, Theta, X
  mat S, V, Theta, Bn, Bsamp;
  mat Lmat = cube2mat(thv.Ld);        //pq*m
  Theta = cube2mat(pars.mbi+pars.bi);   //pq*n
  //
  S = g / (g+1.0) * pars.XtXi;
  V = kron(pinv(pars.tauq), pinv(pars.taup));
  V = V + Lmat * trans(Lmat);
  //
  Bn = S * trans(pars.X) * trans(Theta); 
  //sample population means: beta;-------------------------------------------------/
  Bsamp = rMN(Bn, S, V);                                   // pq*d
  pars.beta = mat2cube(trans(Bsamp), datn.nr, pars.nbase); // p*q*d
  // Update mbi based on beta------------------------------------------------------/
  pars.mbi = cubexmat(pars.beta, pars.X);                   // p*q*n
}
//
/***********************************************************************************/
/* sampleFactors()                                                                 */
/***********************************************************************************/
void sampleFactors(cube &y, mat const& Bs, bool const& v)
{// Sample factors | everything else;
  mat     lmat, res, Q, Stinv, Sipq, fit, qch, sch;
  cube    LFcube;
  vec     b, bi, thi;
  int     ind;
  // -- common variables:
  Sinv =  kron(pinv(pars.tauq), pinv(pars.taup));
  Sinv += kron((pars.iBtB), 1.0/pars.taue*eye<mat>(datn.nr, datn.nr));
  mat pSinv = pinv(Sinv);
  if(v){
    lmat    = cube2mat(thv.Ld);               //loading matrix
    LFcube  = cubexmat(thv.Ld, thv.H);        //latent factor
    Q = trans(lmat) * pSinv * lmat + 1.0 * eye<mat>(thv.nfac, thv.nfac);
  } else{
    LFcube  = matxcube(ths.LL, ths.H, ths.LR);//latent factor
    lmat    = zeros<mat>(datn.nr*pars.nbase, ths.nfacL*ths.nfacR);
    for(r=0; r<ths.nfacL;r++){
      for(s=0; s<ths.nfacR; s++){
        ind = r*ths.nfacR + s;
        lmat.col(ind) = vectorise(ths.LL.col(r)*ths.LR.col(s).t());
      }
    }
    Q = trans(lmat) * pSinv * lmat + 1.0 * eye<mat>(ths.nfacL*ths.nfacR, 
              ths.nfacL*ths.nfacR);
  }
  qch     = safechol(Q);
  Sipq    = kron(pars.taup, pars.tauq);
  Stinv   = Sipq + kron(pars.taue * eye<mat>(datn.nr, datn.nr), pars.BtB);
  sch     = safechol(Stinv);
  //
  for(i=0; i<datn.ns; i++){
    //residual of pq
    res = (flatslice(y,i) * Bs * pars.iBtB - pars.mbi.slice(i));
    b   = (v) ? prv.eta0 : prs.eta0;
    b   += trans(lmat) * pSinv * conv_to<colvec>::from(vectorise(res));
    // Sample latent Factors;--------------------------------------------------------/
    vec bsamp   = BayesRegL(b,qch);
    if(v){
      thv.H.row(i) = conv_to<rowvec>::from(bsamp);
    } else {
      ths.H.slice(i) = vec2mat(bsamp, ths.nfacR, ths.nfacL).t();
    }
    // residual of Theta_i
    bi  = vectorise(pars.taue * trans(Bs) * trans(flatslice(y,i)));
    bi  = bi + Sipq * vectorise(trans(pars.mbi.slice(i) + LFcube.slice(i)));
    // Sample coefficient matrices --------------------------------------------------/
    thi = BayesRegL(bi, sch);
    // Update bi = Theta_i - mbi
    pars.bi.slice(i) = trans(vec2mat(thi, pars.nbase, datn.nr)) - pars.mbi.slice(i);
    fit = (pars.bi.slice(i) + pars.mbi.slice(i)) * trans(Bs);       //p*nt      
    // Update fit_i = Theta_i * Bs'
    for(r=0; r<datn.nr; r++){
      for(t=0; t<datn.nt; t++){
        pars.fit(i,r,t) = fit(r,t);
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
void ErgodicAveragesV(mat const& Bs, bool const& v)
{
  mat beta1, beta2;
  vec tj, tji;
  cube lcube;
  //
  ww1 = 1.0/(1.0 + hatn.n);
  ww2 = 1.0 - ww1;
  hatn.n++;
  //
  // Average fit ------------------------------------------
  hatn.fit  = ww1 * pars.fit + ww2 * hatn.fit;
  hatn.beta2= ww1 * pars.beta + ww2 * hatn.beta2;
  // hatn.L    = ww1 * thv.Ld + ww2*hatv.L;
  if(v){
    lcube = cubexmat(thv.Ld, thv.H);
  } else {
    lcube = matxcube(ths.LL, ths.H, ths.LR);
  }
  //
  for(i=0; i<datn.ns; i++){
    beta1 = pars.mbi.slice(i) * trans(Bs);                       //p*nt
    beta2 = (pars.mbi.slice(i) + lcube.slice(i)) * trans(Bs);    //p*nt
    for(r=0; r<datn.nr; r++){
      for(t=0; t<datn.nt; t++){
        hatn.beta(i,r,t) = ww1*beta1(r,t) + ww2*hatn.beta(i,r,t);
        // hatv.beta2(i,r,t)= w1*std::pow(beta1(r,t), 2.0) + w2*hatv.beta2(i,r,t);
        hatn.betai(i,r,t)= ww1*beta2(r,t) + ww2*hatn.betai(i,r,t);
      }
    }
  }
}

/***********************************************************************************/
/* OutputSample()                                                                  */
/***********************************************************************************/
void OutputSampleV(bool const& v)
{
  { // Population Splines --------------------------------
    std::stringstream ss;
    for(s=0; s<datn.np; s++){
      for(r=0; r<datn.nr; r++){
        for(j=0; j<pars.nbase; j++){
          ss << pars.beta(r,j,s) << " ";
        }
      }
    }
    ss << "\n"; 
    outv1 << ss.str();
  }
  { // Precisions ----------------------------------------
    std::stringstream ss;
    ss << pars.taue << " ";          // Error covs;
    for(r=0; r<datn.nr; r++){
      ss << pars.taup(r,r) << " ";   // Spatial res covs;
    }
    for(j=0; j<pars.nbase; j++){
      ss << pars.tauq(j,j) << " ";   // Temp res covs;
    }
    ss << "\n"; 
    outv2 << ss.str();
  }
  { // Factors -------------------------------------------
    std::stringstream ss;
    for(r=0; r<datn.ns; r++){
      if(v){
        for(j=0; j<thv.nfac; j++){
          ss << thv.H(r,j) << " "; // LF scores
        }
      } else {
        for(j=0; j<ths.nfacR; j++){
          for(k=0; k<ths.nfacL; k++){
            ss << ths.H(k,j,r) << " "; // LF-s scores
          }
        }
      }
    }
    if(v){
      for(t=0; t<thv.nfac; t++){
        ss << thv.pen(t) << " ";      // penalty;
      }
    } else {
      for(t=0; t<ths.nfacL; t++){
        ss << ths.penL(t) << " ";
      }
      for(s=0; s<ths.nfacR; s++){
        ss << ths.penR(s) << " ";
      }
    }
    ss << "\n"; outv3 << ss.str();
  }
}
//
/***********************************************************************************/
/* init_HFMM()                                                                     */
/***********************************************************************************/
//
void init_HFMM(arma::cube y,
               arma::mat const& X,
               arma::mat const& Bs,
               bool const& v)
{
  // Get initialization summaries ---------------------------------------------------
  datn.ns = y.n_rows;           // number of subjects
  datn.nr = y.n_cols;           // number of regions
  datn.nt = y.n_slices;         // number of segments
  datn.np = X.n_cols;           // number of covariates (groups)
  //
  pars.nbase = Bs.n_cols;        // number of basis functions
  //
  double ae = 0.0001; double be = 0.0001; // error init
  // Initialize missing data ----------------------------------------
  {
    vec rn1, yit;
    datn.miss = icube(datn.ns, datn.nr, datn.nt);
    datn.nmiss= vec(datn.ns);
    yit = randu<vec>(datn.ns);
    for(i=0; i<datn.ns; i++){
      for(r=0; r<datn.nr; r++){
        yit = y.tube(i,r);
        for(j=0; j<datn.nt; j++){
          datn.miss(i,r,j) = 0;
          if( yit(j)==12345.12345 ){
            datn.miss(i,r,j) = 1;
            rn1 = randn(1)*0.5;
            y(i,r,j) = meanNA(yit) + rn1(0);
          }
        }
      }
      datn.nmiss[i] = accu(datn.miss.tube(0,0));
    }
  }
  // Initialize and decompose regression matrices -------------------
  pars.X     = X;                      // X
  pars.XtXi  = pinv(trans(X)*X);       // pinv(X'X)
  // Random initialization;
  pars.bi   = randu<cube>(datn.nr, pars.nbase, datn.ns);    // basis functions
  pars.mbi  = randu<cube>(datn.nr, pars.nbase, datn.ns);    // mean bases funcs
  pars.fit  = randu<cube>(datn.ns, datn.nr, datn.nt);       // current fit
  pars.beta = randu<cube>(datn.nr, pars.nbase, datn.np);    // Group means
  thv.Omega =randu<cube>(datn.nr, datn.nr, pars.nbase);     // Spatial correlations
  pars.BtB = trans(Bs)*Bs;
  pars.iBtB= pinv(pars.BtB);
  pars.taup= 0.001 * eye(datn.nr, datn.nr);
  pars.tauq= 0.001 * eye(pars.nbase, pars.nbase);
  pars.taue= Rf_rgamma(ae, 1.0/be);                         // Residual vars
  // Meaningful initialization:                       // ----------------------
  bit    = randu<vec>(pars.nbase);
  for(i=0; i<datn.ns; i++){
    for(r=0; r<datn.nr; r++){
      yit = y.tube(i,r);
      bit = solve((pars.BtB + pars.taue * eye<mat>(pars.nbase, pars.nbase)), trans(Bs) * yit);
      for(j=0; j<pars.nbase; j++){
        pars.bi(r,j,i) = bit[j];
      }
      pars.fit.tube(i,r) = vectorise(Bs*bit);
    }
  }
  // Initialize ergodics --------------------------------------------
  hatn.n     = 0.0;
  hatn.beta  = randu<cube>(datn.ns, datn.nr, datn.nt);
  hatn.beta2 = randu<cube>(datn.nr, pars.nbase, datn.np);
  hatn.betai = randu<cube>(datn.ns, datn.nr, datn.nt);
  // hatv.L     = randu<cube>(datv.nr, thv.nbase, thv.nfac);
  hatn.fit   = randu<cube>(datn.ns, datn.nr, datn.nt);
}
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
  set_cerr_stream(f);
  // number of factors in v;
  thv.nfac  = nfac;             // number of latent factors
  /* ----------------------------------------------------------------------------- */
  /* Initialize Chain and static summaries                                         */
  /* ----------------------------------------------------------------------------- */
  init_HFMM(y, X, Bs, true);
  // Prior parameters -----------------------------------------------
  prv.ae = 0.0001;          prv.be = 0.0001;              // Error variance
  prv.athe = 0.001;         prv.bthe = 0.001;             // coeff res row/col var
  prv.ap1 = 1.5;            prv.ap2 = 1.5;                // penalties
  prv.aOm = 1.0;                                          // --------------could be less informative
  prv.Om0 = eye<mat>(datn.nr, datn.nr) * 1.0;             // spatial smoothing
  prv.eta0 = zeros<vec>(thv.nfac);                        // Latent factor priors
  //
  // RANDOM initialization                                // ----------------------
  thv.H   = randu<mat>(datn.ns, nfac);                    // Latent factors
  thv.Ld  = randu<cube>(datn.nr, pars.nbase, nfac);       // Loading matrix
  thv.HtH = trans(thv.H) * thv.H;
  thv.pen = randu<vec>(nfac);
  for(i=0; i<nfac; i++){thv.pen[i] = Rf_rgamma(prv.ap1, 1.0/prv.ap2);}
  for(j=0; j<pars.nbase; j++){thv.Omega.slice(j) = eye<mat>(datn.nr, datn.nr);}
  //
  // Open output files ----------------------------------------------
  outv1.open(fileOutput1);
  outv2.open(fileOutput2);
  outv3.open(fileOutput3);
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */
  for(int rep=0; rep<(nsim+burnin+1); rep++){
    // mcmc samples --------------------
    completeY2(y);
    sampleLoadingsV();
    sampleFacCovs(true);
    sampleErrCovs(y);
    sampleBetas();
    sampleFactors(y, Bs, true);
    // Store mcmc samples --------------
    if( (rep > burnin) && (rep%thin == 0) ){
      OutputSampleV(true);             // Write mcmc samples to file
      ErgodicAveragesV(Bs, true);      // Update ergodic averages
    }
  }
  // Close output files ---------------------------------------------
  outv1.close(); outv2.close(); outv3.close(); // outv3.close();

  // Output to R --------------------------------------------------------------------
  return List::create(
    Named("fit")      = hatn.fit,
    Named("coef")     = hatn.beta2,
    Named("gmean")    = hatn.beta,
    Named("imean")    = hatn.betai);
}
// Sandwich HFMM
/***********************************************************************************/
/* smix_mcmc()                                                                     */
/***********************************************************************************/
// [[Rcpp::export]]
List smix_mcmc(arma::cube y, 
               arma::mat const& X, 
               arma::mat const& Bs,
               int const& nfacL, int const& nfacR,
               int const& burnin, int const& nsim, int const& thin
)
{
  set_cerr_stream(f);
  // number of factors in v;
  ths.nfacL  = nfacL; ths.nfacR  = nfacR;                // number of latent factors
  /* ----------------------------------------------------------------------------- */
  /* Initialize Chain and static summaries                                         */
  /* ----------------------------------------------------------------------------- */
  init_HFMM(y, X, Bs, false);
  // Prior parameters -----------------------------------------------
  prs.ae = 0.0001;          prs.be = 0.0001;              // Error variance
  prs.athe = 0.001;         prs.bthe = 0.001;             // coeff res row/col var
  prs.ap1 = 1.5;            prs.ap2 = 1.5;                // penalties
  prs.nuL = 5.0;            prs.nuR = 5.0;                // LF precision
  prs.eta0 = zeros<vec>(nfacL*nfacR);                     // Latent factor priors
  //
  // RANDOM initialization                                // ----------------------
  ths.H   = randu<cube>(nfacL, nfacR, datn.ns);           // Latent factors
  ths.LL  = randu<mat>(datn.nr, nfacL);                   // Loading matrix (left)
  ths.LR  = randu<mat>(pars.nbase, nfacR);                // Loading matrix (right)
  ths.penL = randu<vec>(nfacL);
  for(i=0; i<nfacL; i++){ths.penL[i] = Rf_rgamma(prs.ap1, 1.0/prs.ap2);}
  ths.penR = randu<vec>(nfacR);
  for(r=0; r<nfacR; r++){ths.penR[r] = Rf_rgamma(prs.ap1, 1.0/prs.ap2);}
  ths.tauL= randu<mat>(datn.nr, nfacL);
  ths.tauR= randu<mat>(pars.nbase, nfacR);
  for(i=0; i<ths.tauL.n_elem; i++){ths.tauL[i] = Rf_rgamma(prs.nuL, 1.0/prs.nuL);}
  for(j=0; j<ths.tauR.n_elem; j++){ths.tauR[j] = Rf_rgamma(prs.nuR, 1.0/prs.nuR);}
  //
  // Open output files ----------------------------------------------
  outv1.open(fileOutput1);
  outv2.open(fileOutput2);
  outv3.open(fileOutput3);
  /* ----------------------------------------------------------------------------- */
  /* Gibbs Sampling                                                                */
  /* ----------------------------------------------------------------------------- */
  for(int rep=0; rep<(nsim+burnin+1); rep++){
    // mcmc samples --------------------
    completeY2(y);
    sampleLoadingsV();
    sampleFacCovs(false);
    sampleErrCovs(y);
    sampleBetas();
    sampleFactors(y, Bs, false);
    // Store mcmc samples --------------
    if( (rep > burnin) && (rep%thin == 0) ){
      OutputSampleV(false);             // Write mcmc samples to file
      ErgodicAveragesV(Bs, false);      // Update ergodic averages
    }
  }
  // Close output files ---------------------------------------------
  outv1.close(); outv2.close(); outv3.close(); // outv3.close();
  
  // Output to R --------------------------------------------------------------------
  return List::create(
    Named("fit")      = hatn.fit,
    Named("coef")     = hatn.beta2,
    Named("gmean")    = hatn.beta,
    Named("imean")    = hatn.betai);
}

// // END vmix_mcmc --------------------------------------------------------------------