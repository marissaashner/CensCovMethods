#' Complete Case Analysis for Censored Covariates
#'
#' Performs a complete case analysis using least squares for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters.
#' @param data a data frame containing columns for the censoring indicator and the variables in the formula
#' @param cens_ind name of censoring indicator, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored.
#' @param par_vec name of parameter vector in the formula
#' @param starting_vals the starting values for the least squares algorithm. Must be a vector equal in length of the parameter vector
#' @param sandwich_se if \code{TRUE} (default), the empirical sandwich estimator for the standard error is calculated.
#'
#' @return A list with the following elements:
#' \item{beta_est}{a vector of the parameter estimates.}
#' \item{sigma_est}{the square root of the estimated variance of the random error.}
#' \item{se_est}{if \code{sandwich_se = TRUE}, a vector of standard error estimates from the empirical sandwich estimator. Otherwise, \code{NULL}}
#'
#'
#' @import tidyverse
#' @import rootSolve
#' @import numDeriv
#' @importFrom magrittr `%>%`
#'
#' @export
cc_censored <- function(formula,
                        data,
                        cens_ind,
                        par_vec,
                        starting_vals,
                        sandwich_se = TRUE){

  # Need to add error checks

  # subset to those that are not censored
  data_cc <- data %>% dplyr::filter(get(cens_ind) == 1)

  # run nls
  start = list(temp = starting_vals)
  names(start) = par_vec
  model_est <- nls(formula,
                          data = data_cc,
                          start = start,
                          control = nls.control(minFactor = 1/5096,
                                                warnOnly = TRUE))
  beta_est = summary(model_est)$coeff[,1] %>% t() %>% as.data.frame()
  sigma_est = summary(model_est)$sigma
  iteration_count = model_est$convInfo$finIter

  # run sandwich estimator
  if(sandwich_se){
    se_est = cc_sandwich(formula, data, cens_ind, par_vec, beta_est)
  }else{
    se_est = NULL
  }

  # save beta estimates
  return(list(beta_est = beta_est,
              sigma_est = sigma_est,
              se_est = se_est,
              iteration_count = iteration_count))
}

#' Sandwich Estimator for Complete Case for Censored Covariates
#'
#' Calculates the standard error sandwich estimate for a complete case estimator for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters
#' @param data a data frame containing columns for the censoring indicator and the variables in the formula
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param par_vec a character string indicating the parameter vector in the formula
#' @param beta_est the estimate from the complete case estimator
#'
#' @return A vector of the sandwich standard error estimates.
#'
#' @import tidyverse
#' @import rootSolve
#' @import survival
#' @import numDeriv
#'
#' @export
cc_sandwich <- function(formula,
                        data,
                        cens_ind,
                        par_vec,
                        beta_est){

  #convert beta_est to numeric
  beta_est = beta_est %>% as.numeric()

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
  exp_nobracket <- exp %>% as.character() %>% stringr::str_remove_all(., "\\[|\\]") %>%
    parse(text = .)
  m_func = function(p){
    with(as.list(p),
         eval(exp_nobracket))
  }

  # create "g" function for the sandwich estimator
  g = function(data, beta_est, m_func, par_vec, varNamesRHS, cens_ind){
    p = c(beta_est, data[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

    rep(data[cens_ind], length(beta_est)) %>% as.numeric()*
      numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
      rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
  }


  # first derivative function
  # calculates first derivative for subject i (data x)
  # f(x;beta) where beta is of length lb, x is a scalar
  # first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
  firstderivative <- function(beta, g, data, varNamesRHS, par_vec,
                              m_func, cens_ind){
    lb <- length(beta)
    derivs <- matrix(data = 0, nrow = lb, ncol = lb)
    delta <- beta * (10 ^ (- 4))
    betal <- betar <- beta
    for (i in 1:lb) {
      # Perturb the ith element of beta
      betal[i] <- beta[i] - delta[i]
      betar[i] <- beta[i] + delta[i]

      # Calculate function values
      yout1 <- g(data, betal, m_func, par_vec, varNamesRHS, cens_ind)
      yout2 <- g(data, betar, m_func, par_vec, varNamesRHS, cens_ind)

      # Calculate derivative and save in vector A
      derivs[i,] <- (yout2 - yout1) / (2 * delta[i])

      # Reset parameter vectors
      betal <- betar <- beta
    }
    #print("hi")
    return(derivs)
  }

  # take the inverse first derivative of g
  first_der <- apply(data, 1, function(temp){
    firstderivative(beta_est, g, temp, var_namesRHS, par_vec, m_func, cens_ind)
  })
  if(length(beta_est) > 1){
    first_der = first_der %>% rowMeans() %>% matrix(nrow = length(beta_est))
  }else{
    first_der = first_der %>% mean()
  }
  inv_first_der <- solve(first_der)

  # need to get the outer product of g at each observation and take the mean
  gs = apply(data, 1, function(temp)
    g(temp, beta_est, m_func, par_vec, varNamesRHS, cens_ind))
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







