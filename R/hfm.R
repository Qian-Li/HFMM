#' Regionally Referenced Hierarchical Functional Model
#'
#' This function is a wrapper for a Cpp Gibbs Sampler. Region-specific regression is modeled through a varying coefficient functions.
#'
#' @param y A matrix with vectorized region measurements in columns.
#' @param t A vector with vectorized observed time points or segments.
#' @param X A regression design matrix (subject by predictor).
#' @param nsegs An integer vector with the number of observed segments per subject.
#' @param knots Number of spline knots or a vector of spline knots
#' @param mcmc A list of mcmc specifications. Default is (burnin = 1000, nsim = 2000, thin = 1)
#' @return
#' A list of with model fit summaries list(fit, coef, y, t, Bs).
#' \itemize{
#' \item{fit:}{ Posterior mean fit to y in a 3d array [subject, region, segmenr]}
#' \item{coef:}{ Array containing mean posterior regression coefficients [time, predictor, region]}
#' \item{coef2:}{ Array containing second moment of posterior regression coefficients [time, predictor, region]}
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
hfm <- function(y, t, X, nsegs, knots=NULL, mcmc=NULL, spatial = FALSE){
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
  dir   <- files == "post";
  if(length(files[dir]) == 0) dir.create("post");
  #
  ## Get MCMC summaries -------------------------------------------------------
  if(is.null(mcmc)){
    mcmc$burnin <- 1000
    mcmc$nsim   <- 2000
    mcmc$thin   <- 1
  }
  ## Run MCMC -----------------------------------------------------------------
  #
  fit1 <- hmix_mcmc(y3, X, BSX, mcmc$burnin, mcmc$nsim, mcmc$thin, spatial)
  #
  ## END MCMC -----------------------------------------------------------------
  y3[is.na(y3_b)] <- NA;
  return(list(fit=fit1$fit, coef=fit1$coef, coef2=fit1$coef2, y=y3, t=x1, Bs=BSX, nsim=mcmc$nsim))
}

#' @useDynLib HFMM
#' @importFrom Rcpp sourceCpp
NULL
