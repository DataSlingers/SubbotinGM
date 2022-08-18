#' Estimation algorithm for the Extreme Lasso problem.
#' 
#' @description Provides parameter estimates the Extreme Lasso problem, which is used for neighborhood selection.
#' Implements a proximal gradient descent with backtracking.
#' @param resp n by 1 response matrix.
#' @param pred_mat n by p predictor matrix.
#' @param sub_power even integer, power of Subbotin distribution.
#' @param abstol maximum absolute error.
#' @param reltol maximum relative error.
#' @param init_guess initialization of parameter estimates.
#' @param init_step initial step size scale parameter
#' @param shrink backtracking shrinkage parameter
#' @param lambda L1 regularization penalty.
#' @param maxit maximum iterations of gradient descent.
#' 
#' @return A list with the following elements:
#' \itemize{
#'   \item pars: p by 1 matrix, final parameter estimates
#'   \item sub_power: power of Subbotin distribution.
#'   \item iterations: number of iterations run.
#' }
extreme_lasso <- function(resp, pred_mat, sub_power, abstol = 0.001, reltol = 0.001,
                                  init_guess = NULL, init_step = 1, shrink = 0.5,
                                  lambda = 1, maxit = 100){
  
  # Checks
  if(lambda < 0) {stop("Invalid regularization parameter")}
  
  ## Initializations
  if(is.null(init_guess)){
    par_mat <- matrix(0, ncol(pred_mat), 1)
  } else {
    par_mat <- matrix(init_guess, ncol(pred_mat), 1)
  }
  
  resp <- matrix(resp, ncol = 1)
  
  # Keeps tack of objective value. 
  ### Change this.
  opt_path <- c(sum(abs((resp - pred_mat %*% par_mat)^sub_power))) + lambda * sum(abs(par_mat))
  delta <- Inf
  rel_delta <- Inf
  iter <- 0
  # Keeps track of all parameter estimates along gradient descent
  par_list <- list()
  par_list[[length(par_list) + 1]] <- par_mat
  
  ### Estimation
  while((delta > abstol && rel_delta > reltol && iter < maxit)){
    iter <- iter + 1 
    ## Reset
    step_size <- init_step / shrink
    lhs <- 1
    rhs <- 0
    
    ## Check step size with backtracking
    while(lhs > rhs){
      step_size <- step_size * shrink
      # Calculate gradient (dd), L1 prox of gradient (pdd)
      dd <- par_mat - step_size * -sub_power * (t(pred_mat) %*% ((resp - pred_mat %*% par_mat)^(sub_power - 1)))
      pdd <- matrix(NA, nrow(dd), 1)
      for(ii in 1:length(dd)){
        if(dd[ii] > lambda * step_size){
          pdd[ii] <- dd[ii] - lambda  * step_size
        } else if(dd[ii] < - lambda * step_size) {
          pdd[ii] <- dd[ii] + lambda  * step_size
        } else {
          pdd[ii] <- 0
        }
      }
      
      rr <- (par_mat - pdd) / step_size
      lhs <- sum(abs((resp - pred_mat %*% pdd)^sub_power))
      rhs <- sum(abs((resp - pred_mat %*% par_mat)^sub_power)) - 
        step_size * t(-sub_power * (t(pred_mat) %*% ((resp - pred_mat %*% par_mat)^(sub_power - 1)))) %*% rr +
        sum((rr)^2) * step_size / 2
    }
    
    ## Run gradient step
    # Calculate gradient (dd), L1 prox of gradient (pdd)
    dd <- par_mat - step_size * -sub_power * (t(pred_mat) %*% ((resp - pred_mat %*% par_mat)^(sub_power - 1)))
    pdd <- matrix(NA, nrow(dd), 1)
    for(ii in 1:length(dd)){
      if(dd[ii] > lambda * step_size){
        pdd[ii] <- dd[ii] - lambda * step_size 
      } else if(dd[ii] < -lambda * step_size) {
        pdd[ii] <- dd[ii] + lambda * step_size 
      } else {
        pdd[ii] <- 0
      }
    }
    
    ## Calculate likelihood value to check for convergence
    opt_path <- c(opt_path, sum((resp - pred_mat %*% pdd)^sub_power) + lambda * sum(abs(pdd)))
    delta <- opt_path[length(opt_path) - 1] - opt_path[length(opt_path)]
    rel_delta <- 1 - (opt_path[length(opt_path) - 1] / opt_path[length(opt_path)])
    par_list[[length(par_list) + 1]] <- pdd
    par_mat <- pdd
  }
  
  return(list(pars = par_list[[length(par_list)]], sub_power = sub_power, iterations = iter))
}
