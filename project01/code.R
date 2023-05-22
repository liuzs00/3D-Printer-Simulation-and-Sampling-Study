#' s2097920, Zongsheng Liu
#' Add your own function definitions on this file.

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}


#1
library(ggplot2)
library(ggExtra)
library(StatCompLab)
library(shiny)
library(mvtnorm)



#2
#load the required data set
data("filament1", package = "StatCompLab")

#2.1
#compute the joint prior density 
log_prior_density <- function(theta, params){ 
  #Product operator changed into sigma since logarithm was used
  l_prior_d <- (sum(lpd_theta_1 <- dnorm(theta[1], mean = 0, sd = sqrt(params[1]), log = TRUE),
                    lpd_theta_2 <- dnorm(theta[2], mean = 1, sd = sqrt(params[2]), log = TRUE),
                    lpd_theta_3 <- dlogexp(theta[3], rate = params[3], log = TRUE),
                    lpd_theta_4 <- dlogexp(theta[4], rate = params[4], log = TRUE))) 
  l_prior_d
}

#2.2
log_like <- function(theta, x, y){
  #sum up all the log likelihood 
  LogLikelihood <- (sum(dnorm(y, 
                              mean = theta[1]+x*theta[2], 
                              sd = sqrt(exp(theta[3])+x^2*exp(theta[4])), 
                              log = TRUE)))
  LogLikelihood
}


#2.3
log_posterior_density <- function(theta, x, y, params){
  l_posterior_d <- (log_prior_density(theta, params)+log_like(theta, x, y))
  l_posterior_d
}

#2.4
posterior_mode <- function(theta_start, x, y, params){
  opt <- optim(par = theta_start,  #compute the mode which max the log posterior density
               fn = log_posterior_density, 
               x=x, 
               y=y, 
               params=params,
               hessian = TRUE, # compute Hessian at mode
               control=list(fnscale=-1))  #additional argument to change default minimization process to maximization  
  
  mode <- opt$par # mode
  hessian <- opt$hessian # Hessian
  S <- solve(-hessian) #negated inverse of the Hessian
  return(list(mode=mode, Hessian=hessian, Negated_Inverse_Hessian=S)) #return them as a list
}

#2.5
#set up all the initial value
params <- c(1,1,1,1)
theta_start <- c(0,0,0,0)
x <- filament1$CAD_Weight # read data from filament1 with column name CAD_Weight
y <- filament1$Actual_Weight # read data from filament1 with column name Actual_Weight
#compute mode, Hessian, negated inverse of the Hessian
sol <- posterior_mode(theta_start, x, y, params)
#call out each component from the above list

#2.6
do_importance <- function(N, mu, S, x, y, params){
  sample <- rmvnorm(N, mean = mu, sigma = S)
  #compute the unnormalized weights
  log_weights <- log_posterior_density(sample, x, y, params)-dmvnorm(sample, mean = mu, sigma = S, log = TRUE)
  #compute the normalized weight
  log_weights_norm <- log_weights - log_sum_exp(log_weights) 
  log_weight_table <-data.frame(beta1 = sample[,1],
                                beta2 = sample[,2],
                                beta3 = exp(sample[,3]),
                                beta4 = exp(sample[,4]),
                                log_weights = log_weights_norm)
  log_weight_table
}

#2.7
mu <- sol$mode
S <- sol$Negated_Inverse_Hessian
data= do_importance(10000, mu, S, x, y, params)



  
#using wquantile to construct a 90% credible interval for all Beta
make_CI <- function(x, weights, prob){
  credible_interval <- wquantile(x, probs = c((1-prob)/2, 1-(1-prob)/2), weights = weights)
  credible_interval

}


#make the credible interval for Beta into a table
CI_table <- data.frame(beta=c("beta1","beta2","beta3","beta4"),
                       lower_bound = c(make_CI(data[,"beta1"], exp(data[,"log_weights"]), 0.9)[1],
                                       make_CI(data[,"beta2"], exp(data[,"log_weights"]), 0.9)[1],
                                       make_CI(data[,"beta3"], exp(data[,"log_weights"]), 0.9)[1],
                                       make_CI(data[,"beta4"], exp(data[,"log_weights"]), 0.9)[1]),
                       upper_bound =c(make_CI(data[,"beta1"], exp(data[,"log_weights"]), 0.9)[2],
                                      make_CI(data[,"beta2"], exp(data[,"log_weights"]), 0.9)[2],
                                      make_CI(data[,"beta3"], exp(data[,"log_weights"]), 0.9)[2],
                                      make_CI(data[,"beta4"], exp(data[,"log_weights"]), 0.9)[2]) )





