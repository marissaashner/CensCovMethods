#' Inverse Probability Weighting for Censored Covariates
#'
#' Performs an inverse probability weighting estimation using least squares for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters
#' @param data a data frame containing columns for the censoring indicator and the variables in the formula
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param par_vec a character string indicating the parameter vector in the formula
#' @param starting_vals the starting values for the least squares algorithm. Must be a vector equal in length of the parameter vector
#' @param sandwich_se if \code{TRUE} (default), the empirical sandwich estimator for the standard error is calculated
#' @param weight_opt a character string indicating which method of weight calculation is to be done. One of "Cox", "AFT_lognormal", "MVN", "user" (if "user", then user provides weights)
#' @param weights_user if \code{weight_opt = "user"}, a vector of weights the same length as there are rows in \code{data}, otherwise \code{NULL} (default)
#' @param weights_cov if \code{weight_opt} one of \code{c("Cox", "AFT_lognormal", "MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model. Otherwise \code{NULL}
#'
#' @return A list with the following elements:
#' \item{beta_est}{a vector of the parameter estimates.}
#' \item{se_est}{if \code{sandwich_se = TRUE}, a vector of standard error estimates from the empirical sandwich estimator. Otherwise, \code{NULL}}
#'
#'
#' @import tidyverse
#' @import numDeriv
#'
#' @export
ipw_censored <- function(formula,
                        data,
                        cens_ind,
                        par_vec,
                        starting_vals,
                        sandwich_se = TRUE,
                        weight_opt,
                        weights_user = NULL,
                        weights_cov = NULL){

  # weights
  if(weight_opt == "user"){
    weights = weights_user
  }else if(weight_opt == "Cox"){
    weights = weights_cox()
  }else if(weight_opt == "AFT_lognormal"){
    weights = weights_aft()
  }else if(weight_opt == "MVN"){
    weights = weights_mvn()
  }

  # run nls
  start = list(temp = starting_vals)
  names(start) = par_vec
  beta_est <- summary(nls(formula,
                          data = data_cc,
                          start = start,
                          weights = get(cens_ind)*weights,
                          control = nls.control(minFactor = 1/5096,
                                                warnOnly = TRUE)))$coeff[,1] %>% t() %>% as.data.frame()

  # run sandwich estimator
  if(sandwich_se){
    se_est = ipw_sandwich(formula, data, cens_ind, par_vec, beta_est, weights)
  }else{
    se_est = NULL
  }

  # save beta estimates
  return(list(beta_est = beta_est,
              se_est = se_est))
}

ipw_sandwich <- function(formula,
                        data,
                        cens_ind,
                        par_vec,
                        beta_est,
                        weights){

  #convert beta_est to numeric
  beta_est = beta_est %>% as.numeric()

  # add weights to data frame
  data$weights = weights

  # extract variable names from formula
  varNames = all.vars(formula)
  form2 = formula
  form2[[2]] <- NULL
  varNamesRHS = all.vars(form2)
  Y = varNames[is.na(match(varNames, varNamesRHS))]

  do.call("<-", list(par_vec, beta_est))

  varNamesRHS = varNamesRHS[is.na(match(varNamesRHS, par_vec))]

  # turn formula into function

  cmd <- tail(as.character(formula),1)
  exp <- parse(text=cmd)
  exp_nobracket <- exp %>% as.character() %>% str_remove_all(., "\\[|\\]") %>%
    parse(text = .)
  m_func = function(p){
    with(as.list(p),
         eval(exp_nobracket))
  }

  # create "g" function for the sandwich estimator
  g = function(data, beta_est, m_func, par_vec, var_namesRHS, cens_ind){
    do.call("<-", list(varNamesRHS, data[varNamesRHS]))
    p = c(beta_est, lapply(varNamesRHS, get) %>% unlist()) %>% unlist()
    names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

    rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
      numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
      rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()
  }

  # jacobian g function
  g_jacobian = function(data, beta_est, m_func, par_vec, var_namesRHS, cens_ind){
    do.call("<-", list(varNamesRHS, data[varNamesRHS]))
    p = c(beta_est, lapply(varNamesRHS, get) %>% unlist()) %>% unlist()
    names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

    rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
      ((numDeriv::hessian(m_func, p)[1:length(beta_est), 1:length(beta_est)]*
          rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()) -
         outer(numDeriv::jacobian(m_func, p)[1:length(beta_est), 1:length(beta_est)])) %>% as.numeric()
  }

  # take the inverse first derivative of g
  first_der <- apply(data, 1, function(temp){
    g_jacobian(temp, beta_est, m_func, par_vec, var_namesRHS, cens_ind)
  })
  if(length(beta_est) > 1){
    first_der = first_der %>% rowMeans() %>% matrix(nrow = length(beta_est))
  }else{
    first_der = first_der %>% mean()
  }
  inv_first_der <- solve(first_der)

  # need to get the outer product of g at each observation and take the mean
  gs = apply(data, 1, function(temp)
    g(temp, beta_est, m_func, par_vec, var_namesRHS, cens_ind))
  if(length(beta_est) > 1){
    outer_prod = apply(gs, 2, function(g) g%*%t(g))
    outer_prod = outer_prod %>% rowMeans() %>% matrix(nrow = length(beta_est))
  }else{
    outer_prod = gs^2
    outer_prod = outer_prod %>% mean()
  }


  ## then need to put it all together
  se = sqrt((inv_first_der %*% outer_prod %*% t(inv_first_der) / nrow(data)) %>% diag())
  return(se)

}


weights_cox <- function(){

}

weights_aft <- function(){

}

weights_mvn <- function(){

}

