## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----genDataIn----------------------------------------------------------------
# generating inlying smooth curves

set.seed(20220704)

J=15; c=2; a=5; dt = c*(1:J)^(-a)
g1 = c*VGAM::zeta(a)
ga = c(g1, sapply(1:(J-1), function(j){g1 - sum(dt[1:j])}))

# eigenfunctions
tt = 50; tGrid = seq(0, 1, len=tt); 
phi= t(fda::fourier(tGrid, J))

n = 50; s = sqrt(3)
xi = matrix(runif(n*J, -s, s), n, J) * runif(n, -s, s)
Xin = xi %*% (phi * sqrt(ga))

## ----genDataOut---------------------------------------------------------------
# outlying curves

eta = 3*sqrt(sum(ga))

fND = function(x){
 (x>=-0.05 & x<0.05) * (0.05*x) +
 (x>=0.05 & x<0.15) * (-0.05*(x-0.1)) +
 (x>=0.15 & x<0.25) * (0.05*(x-0.2)) +
 (x>=0.25 & x<0.35) * (-0.05*(x-0.3)) +
 (x>=0.35 & x<0.45) * (0.05*(x-0.4)) +
 (x>=0.45 & x<0.55) * (-0.05*(x-0.5)) +
 (x>=0.55 & x<0.65) * (0.05*(x-0.6)) +
 (x>=0.65 & x<0.75) * (-0.05*(x-0.7)) +
 (x>=0.75 & x<0.85) * (0.05*(x-0.8)) +
 (x>=0.85 & x<0.95) * (-0.05*(x-0.9)) +
 (x>=0.95 & x<1.05) * (0.05*(x-1))
}

XoutND = fND(tGrid)*eta*400*0.8   # non-differentiable curve

XoutLin = (2*tGrid-1)*eta*0.6     # linear curve

## ----genData------------------------------------------------------------------
X = rbind(Xin[1:(n-2),], XoutND, XoutLin)

## ---- echo=FALSE, out.width='60%', fig.align='center', fig.width=7, fig.height=7----
matplot(tGrid, t(Xin[1:(n-2),]), type='l', col="grey70", lty=1,
        xlab="Time", ylab="",
        main="Smooth inlying curves (grey) with two outliying curves (black)")
lines(tGrid, XoutND, col=1, lwd=2)
lines(tGrid, XoutLin, col=1, lwd=2)

## -----------------------------------------------------------------------------
library(RHD)
resRHD = RHD(
  X, X, tGrid,
  prob_vec = c(0.95), J_vec=c(6,10), f_iqr_vec=c(1.5, 3), M_vec = rep(1000,2)
)

## -----------------------------------------------------------------------------
resRHD$depth["J=6", ,c(10,20,30,49,50),]

## -----------------------------------------------------------------------------
resRHD$out[["J=6", 1, "f_iqr=1.5"]]
resRHD$out[["J=6", 1, "f_iqr=3"]]

## -----------------------------------------------------------------------------
resRHD$ld_mat

## ----adjF---------------------------------------------------------------------
M_adj = 100; nSim = 100; covMat = cov(X)
ld_mat = matrix(1:4, nrow=2); J_vec = c(5,8); f_iqr_vec = c(1.5, 2, 2.5, 3, 3.5);
M_vec = c(1000, 1000); M_type=0; num_M = 20; sst = 1

resAF = AFSG(
 M_adj, nSim, covMat, tGrid, ld_mat, J_vec, f_iqr_vec,
 M_vec, M_type, num_M, sst
 )

## -----------------------------------------------------------------------------
resAF$rateIn

## -----------------------------------------------------------------------------
resAF$adjF

