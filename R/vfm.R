#' Regionally Referenced Hierarchical Functional Model (Vectorized Version)
#'
#' This function is a wrapper for a Cpp Gibbs Sampler. Region-specific regression is modeled through a varying coefficient functions.
#'
#' @param y A matrix with vectorized region measurements in columns.
#' @param t A vector with vectorized observed time points or segments.
#' @param X A regression design matrix (subject by predictor).
#' @param nsegs An integer vector with the number of observed segments per subject.
#' @param nLF Number of latent factors. Default is (nkonts * nregions)/2
#' @param knots Number of spline knots or a vector of spline knots.
#' @param mcmc A list of mcmc specifications. Default is (burnin = 1000, nsim = 2000, thin = 1)
#' @return
#' A list of with model fit summaries list(fit, coef, y, t, Bs).
#' \itemize{
#' \item{fit:}{ Posterior mean fit to y in a 3d array [subject, region, segment]}
#' \item{coef:}{ Array containing mean posterior regression coefficients [region, spline, predictor]}
#' \item{gfit:}{ Posterior Fixed effects fit in a 3d array [subject, region, segment]}
#' \item{ifit:}{ Posterior Latent Factor model fit in a 3d array [subject, region, segment]}
#' \item{y:}{ Reashaped data array [subject, region, segment]}
#' \item{t:}{ Matrix of observed segments [subject, segments]}
#' \item{Bs:}{ Matrix of spline polynomials used in the model}
#' \item{X:}{ Covariate matrix}
#' \item{nsim:}{ Number of MCMC samples}
#' }
#'
#' @details
#' For region r, subject i, let y_i(t) ...
#'
#' @examples
#'
#'
#' @export
vfm <- function(y, t, X, nsegs, nLF = NULL, knots=NULL, mcmc=NULL){
  #
  # Normalize longitudinal time scale to be between 0 and 1 -------------------
  #
  tmin <- min(t, na.rm=TRUE)
  tmax <- max(t, na.rm=TRUE)
  t    <- (t - tmin)/(tmax - tmin)
  if(!is.null(knots)) knots <- (knots - tmin)/(tmax - tmin)
  #
  ## Get basic data summaries -------------------------------------------------
  #
  nreg <- ncol(y)            # number of regions
  seg1 <- sort(unique(t))    # unique segments or time-points
  ns1  <- length(seg1)       # number of unique segments
  nsub <- nrow(X)            # number of subjects
  #
  # Put Observed Segments in square Matrix and pad NA for unobserved -------------
  #
  i1 <- cumsum(c(1, nsegs))[1:nsub]
  i2 <- cumsum(c(1, nsegs))[2:(nsub+1)] - 1
  #
  x1 <- matrix(NA, nrow=nsub, ncol=ns1)
  y3 <- array(NA, c(nsub, nreg, ns1))
  #
  for(i in 1:nsub){
    index <- (seg1 %in% t[i1[i]:i2[i]])
    x1[i,index] <- seg1[index]
    y3[i,,index] <- t(y[i1[i]:i2[i],])
  }
  #
  # To be edited out -------------
  #
  # Remove margins if too many NA's
  yy <- y3[,1,]
  nn <- rep(0, ncol(yy))
  for(j in 1:ncol(yy)){
    yj    <- yy[,j]
    nn[j] <- length(yj[is.na(yj)])/ncol(yy)
  }
  # Leave aside for now.
  # Trim segments with too many NA's
  # y3   <- y3[ , , nn < 0.20]
  # x1   <- x1[ , nn < 0.20]
  # seg1 <- seg1[nn < 0.20]
  # ns1  <- length(seg1)       # number of unique segments
  #
  #
  y3_b <- y3
  y3[is.na(y3)] <- 12345.12345
  
  ## Get B-Spline Matrix ------------------------------------------------------
  if(is.null(knots)){
    ns1 <- round(length(seg1)/5)
    knots1 <- seq(0,seg1[ns1],length=(ns1+2))[2:(ns1+1)]
  } else if(length(knots)==1){
    knots1 = seq(0,seg1[ns1],length=(knots+2))[2:(knots+1)]
  } else if(length(knots)>1) {
    knots1 <- (knots - min(seg1))/(max(seg1) - min(seg1))  ## renormalize knots
  } 
  BSX    <- splines::bs(seg1, knots=knots1, intercept=T)
  nbase  <- ncol(BSX)
  #
  # Create Output directory ---------------------------------------------------
  files <- list.files();
  dir   <- files == "post_v";
  if(length(files[dir]) == 0) dir.create("post_v");
  #
  ## Get MCMC summaries -------------------------------------------------------
  if(is.null(mcmc)){
    mcmc$burnin <- 1000
    mcmc$nsim   <- 2000
    mcmc$thin   <- 1
  }
  if(length(nLF) == 2){
    nLF = ceiling(nLF)
    ## sandwich HFMM- MCMC---
    cppfit = smix_mcmc(y3, X, BSX, nLF[1], nLF[2], mcmc$burnin, mcmc$nsim, mcmc$thin)
  } else if(length(nLF == 1)){
    nLF = ceiling(nLF)
    ## vectorized HFMM- MCMC---
    cppfit <- vmix_mcmc(y3, X, BSX, nLF,mcmc$burnin, mcmc$nsim, mcmc$thin)
  } else {
    nLF = ceiling(nreg * nbase / 2.0);
    warning('Number of Latent Factors not given, default nLF = ',nLF)
    ## vectorized HFMM- MCMC---
    cppfit <- vmix_mcmc(y3, X, BSX, nLF,mcmc$burnin, mcmc$nsim, mcmc$thin)
  }
  ## END MCMC -----------------------------------------------------------------
  y3[is.na(y3_b)] <- NA;
  return(list(fit=cppfit$fit, coef=cppfit$coef, gfit=cppfit$gmean, ifit=cppfit$imean, y=y3, t=x1, Bs=BSX, nsim=mcmc$nsim))
  # return(nLF)
}

#' @useDynLib HFMM
#' @importFrom Rcpp sourceCpp
NULL
