#' Parameter estimation in a Mh robust design
#' 
#' This function compute demographic parameters and capture probabilities per PSP 
#' using a Mh robust design model in capture-recapture experiments.
#' 
#' @param betaui A vector of size 2 containing the starting values for the optimization problem. 
#'             The vector is obtained by running the function \code{\link{closedp.0}} of \pkg{Rcapture}. 
#' @param Ci A vector containing the total number of captures per PSP.
#' @param Cipsi A vector for the heterogeneity model computed by multiplying  \code{Ci} to the heterogeneity function.
#' @param vi A vector containing the number of units seen for the last time per PSP.
#' @param ui A vector containing the number of first captures per PSP.
#' @param wi A vector containing the number of units captured before a defining PSP i, that will be seen at least
#'           once more, either at PSP i or later.
#'                    
#' @return \item{betai}{A vector of length I containing the estimates for the first log-linear parameter per PSP.}
#' @return \item{taui}{A vector of length I containing the estimates for the second log-linear parameter per PSP.}
#' 
#' @author Mamadou  Yauck and Louis-Paul Rivest
#'
#' @export
MaxPMh_2 <- function(betaui,Ci,Cipsi,vi,ui,wi)
{
  pistar<-1-(1/sum(C*exp(betaui[1]*A+betaui[2]*B)))
  tilden<-(Ci*sum(C*exp(betaui[1]*A+betaui[2]*B))*pistar)/sum(A*C*exp(betaui[1]*A+betaui[2]*B))
  y <- numeric(2)
  y[1] <- (pistar-max((tilden-ui),0)*max((tilden-vi),0)/(max((tilden-ui),0)*max((tilden-vi),0)+tilden*max((wi-tilden+ui),0)))
  y[2] <- ((Ci/Cipsi)-(sum(A*C*exp(betaui[1]*A+betaui[2]*B))/sum(B*C*exp(betaui[1]*A+betaui[2]*B))))
  y
}