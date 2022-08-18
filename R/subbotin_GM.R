#' Edge selection with the Subbotin graphical model using neighborhood selection.
#' @param dat n by p data matrix.
#' @param lambda value of regularization parameter.
#' @param nu power of Subbotin graphical model.
#' @param sym "and" vs. "or" rule for neighborhood selection.
#' @param ... variable arguments for neighborhood regression function.
#' 
#' @return A list with the following elements:
#' \itemize{
#'   \item edges: p by p adjacency matrix. 
#'   \item param: p by p asymetric parameter matrix from neighborhood regressions
#' }

subbotin_GM <- function(dat, lambda, nu, sym = "and", ...){
  
  if(any(class(dat) == "data.frame")){
    dat <- as.matrix(dat)
  }
  
  ## How to symmetrize edges of graph.
  if(sym == "or"){
    sym_num <- 1
  } else if (sym == "and") {
    sym_num <- 2
  } else {
    stop("Invalid sym argument.")
  }
  
  # Regularization parameter scaling
  lamb_nu <- (sqrt(log(ncol(dat)) * nrow(dat)) * lambda)^nu
  
  ## Fit individual neighborhood regressions.
  # Scale data a priori.
  dat <- apply(dat, 2, scale)
  for(ii in 1:ncol(dat)){
    sublasso_mat <- matrix(NA, ncol(dat), ncol(dat))
    y <- dat[, ii]
    x <- dat[, -ii]
    # Extreme lasso regression.
    reg_mod <- extreme_lasso(y, x, sub_power = nu, lambda = lamb_nu, ...)
    sublasso_mat[ii, -ii] <- as.numeric(c(reg_mod$pars != 0))
  }
  
  ## Edge selection
  diag(sublasso_mat) <- 0
  neighborhood_mat <- matrix(0, nrow(sublasso_mat), ncol(sublasso_mat))
  for(jj in 1:(nrow(sublasso_mat) - 1)){
    for(kk in (jj + 1):ncol(sublasso_mat)){
      if(sublasso_mat[jj, kk] + sublasso_mat[kk, jj] >= sym_num){
        neighborhood_mat[jj, kk] <- 1
        neighborhood_mat[kk, jj] <- 1
      }
    }
  }
  
  return(edges = neighborhood_mat, param = sublasso_mat)
}