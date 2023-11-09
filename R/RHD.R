#' @title Regularized halfspace depth for functional data
#' 
#' @description Compute depths and the associated rankings of the functional data from the regularizaed haflspace depth (RHD) and find outliers therein based on the outlier detection method linked to the RHD.
#'
#'
#' @param x An n by p matrix of curves where the depths are evaluated. Each row represents one curve. 
#' @param X An N by p matrix of observed curves that are used for constructing the depth function. Each row represents one observed curve. 
#' @param tGrid A vector of p densely equi-spaced time grid points where the curve are evaluated.
#' @param ld_mat A matrix of regularization parameters. The number of rows should match the length of J_vec.
#' @param prob_vec A vector of numeric values in (0,1) that determines the regularization parameters for each truncation level in J_vec based on the RKHS norms from the randomly drawn directions. For each truncation level J in J_vec, (1000) potential random directions are drawn and its p-th quantile is determined as a regularization parameter, where p is a probability in prob_vec.
#' @param J_vec A vector of truncation levels used for estimating the covariance operator of the observed curves. 
#' @param f_iqr_vec A vector of adjustment factors that are used for outlier detection methods linked to the RHD.
#' @param M_vec A vector of integers. Each integer represents the number of randomly drawn directions for rejection sampling approach to compute the RHD. The lengths of J_vec and M_vec should be equal.
#' @param M_type An integer 0 or 1. If M_type=0, the integers in M_vec mean the total number of drawn directions; if M_type=1, the integers in M_vec mean the number of accepted directions.
#' @param ties A character string specifying how ties are treated as rank function in R base.
#' @return A list containing the following fields:
#' \item{depth}{Resulting depth values and depth rankings for each truncation level J and regularization parameter lambda.} 
#' \item{out}{The indices of outliers from the outlier detection method linked to the RHD for each truncation level J, regularization parameter lambda, and adjustment factor and f_iqr.}   
#' \item{ld_mat}{The matrix of regularization parameters that are used for computing the RHD.}
#' \item{RKHSnorm}{The RKHS norms of all randomly generated directions.}
#' \item{acpt_rate}{The rates of random directions whose RKHS norms are less than or equal to each regularization parameter.}
#' \item{M_vec}{The numbers of all randomly generated directions for each truncation level J.}
#' 
#' 
#' 
#' 
#' @examples
#' 
#' # example
#' 
#' library(RHD)
#' 
#' # generating inlying smooth curves
#' J=15
#' c=2; a=5; dt = c*(1:J)^(-a)
#' g1 = c*VGAM::zeta(a)
#' ga = c(g1, sapply(1:(J-1), function(j){g1 - sum(dt[1:j])}))
#' 
#' # eigenfunctions
#' tt = 50; tGrid = seq(0, 1, len=tt); J=15
#' phi= t(fda::fourier(tGrid, J))
#' 
#' n = 50; s = sqrt(3)
#' xi = matrix(runif(n*J, -s, s), n, J) * runif(n, -s, s)
#' Xin = xi %*% (phi * sqrt(ga))
#' 
#' 
#' # outlying curves
#' 
#' eta = 3*sqrt(sum(ga))
#' 
#' fND = function(x){
#'  (x>=-0.05 & x<0.05) * (0.05*x) +
#'  (x>=0.05 & x<0.15) * (-0.05*(x-0.1)) +
#'  (x>=0.15 & x<0.25) * (0.05*(x-0.2)) +
#'  (x>=0.25 & x<0.35) * (-0.05*(x-0.3)) +
#'  (x>=0.35 & x<0.45) * (0.05*(x-0.4)) +
#'  (x>=0.45 & x<0.55) * (-0.05*(x-0.5)) +
#'  (x>=0.55 & x<0.65) * (0.05*(x-0.6)) +
#'  (x>=0.65 & x<0.75) * (-0.05*(x-0.7)) +
#'  (x>=0.75 & x<0.85) * (0.05*(x-0.8)) +
#'  (x>=0.85 & x<0.95) * (-0.05*(x-0.9)) +
#'  (x>=0.95 & x<1.05) * (0.05*(x-1))
#' }
#' 
#' XoutND = fND(tGrid)*eta*400*0.8   # non-differentiable curve
#' 
#' XoutLin = (2*tGrid-1)*eta*0.6     # linear curve
#' 
#' 
#' # apply the RHD
#' 
#' X = rbind(Xin[1:(n-2),], XoutND, XoutLin)
#' 
#' RHD(
#'   X, X, tGrid, 
#'   prob_vec = c(0.95), J_vec=c(6,10), f_iqr_vec=c(1.5, 3),
#'   M_vec=rep(1000,2)
#'   )
#' 
#' 
#' @export
RHD = function(
    x, X, tGrid, ld_mat=0, prob_vec=0.95, J_vec = round(nrow(X)^(1/(7+0.1))), f_iqr_vec = 1.5, 
    M_vec = 1000L, M_type = 0L, ties = "min"
){
  
  if(length(J_vec)==1){J_vec = c(J_vec)}
  if(length(M_vec)==1){M_vec = c(M_vec)}
  if(length(prob_vec)==1){prob_vec = c(prob_vec)}
  
  resFPCA = FPCA(X, tGrid)
  gaHat = resFPCA$ev
  phiHat = resFPCA$ef
  
  if(any(ld_mat==0)){
    ld_size = length(prob_vec)
    res = RHD_inner_prob(
      x, X, tGrid, gaHat, phiHat, 
      prob_vec, J_vec, f_iqr_vec, M_vec,  
      ties, out_type=1, num_M = 100, stt=0
      )
    ld_mat = res$misc$ld_mat
  }else{
    ld_size = ncol(ld_mat)
    res = RHD_inner(
      x, X, tGrid, gaHat, phiHat, 
      ld_mat, J_vec, f_iqr_vec, M_vec, M_type, 
      ties, out_type=1, num_M = 100, stt=0
      )
  }
  
  # outliers
  for(ii_J in 1:length(J_vec)){
    for(ii_ld in 1:ld_size){
      for(ii_f in 1:length(f_iqr_vec)){
        res$out[[ii_J,ii_ld,ii_f]] = 
          intersect(
            c(res$out[[ii_J,ii_ld,ii_f]]),
            c(res$misc$idx_depth_min[[ii_J, ii_ld, 1]])
          )
      }
    }
  }
  dimnames(res$out) = list(
    paste0("J=", J_vec),
    paste0("lambda", 1:ld_size),
    paste0("f_iqr=", f_iqr_vec)
  )
  
  # depth and ranking
  temp = abind::abind(res$depth, res$rank, along=4)
  dimnames(temp) = list(
    paste0("J=", J_vec),
    paste0("lambda", 1:ld_size),
    1:n,
    c("depth", "ranking")
  )
  
  # regularization parameters
  dimnames(ld_mat) = list(
    paste0("J=", J_vec),
    paste0("lambda", 1:ld_size)
  )
  
  # RKHS norms
  dimnames(res$misc$normRKHS) = list(
    paste0("J=", J_vec), NULL, NULL
  )
  # acceptance rates
  dimnames(res$misc$acpt_rate) = list(
    paste0("J=", J_vec), paste0("lambda", 1:ld_size)
  )
  
  # return
  return(
    list(
      depth=temp, out=res$out, ld_mat=ld_mat, 
      RKHSnorms = res$misc$normRKHS, acpt_rate = res$misc$acpt_rate, M_vec=M_vec
      )
  )
}

