#' Augmented Inverse Probability Weighting for Censored Covariates
#'
#' Performs an augmented inverse probability weighting estimation for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters
#' @param data a data frame containing columns for the censoring indicator and the variables in the formula
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
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
#' @import survival
#'
#' @export
aipw_censored <- function(formula,
                         data,
                         cens_ind,
                         cens_name,
                         par_vec,
                         starting_vals,
                         sandwich_se = TRUE,
                         weight_opt,
                         weights_user = NULL,
                         weights_cov = NULL){

  # Need to add error checks

  # add thresholding options

  # weights
  if(weight_opt == "user"){
    weights = weights_user
  }else if(weight_opt == "Cox"){
    weights = weights_cox(data, cens_ind, weights_cov, cens_name)
  }else if(weight_opt == "AFT_lognormal"){
    weights = weights_aft(data, cens_ind, weights_cov, cens_name)
  }else if(weight_opt == "MVN"){
    mvn_results = weights_mvn(data, cens_ind, weights_cov, cens_name)
    weights = mvn_results$weights
    mu_joint = mvn_results$mu_joint
    Sigma_joint = mvn_results$Sigma_joint
  }

  # add weights to data frame
  data$weights = weights

  # extract variable names from formula
  varNames = all.vars(formula)
  form2 = formula
  form2[[2]] <- NULL
  varNamesRHS = all.vars(form2)
  Y = varNames[is.na(match(varNames, varNamesRHS))]
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

  # Note: the A function is as follows
  # set p --- p = c(beta_temp, data[varNamesRHS])
  # set names for p ---- names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  # A = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]

  # CC to get estimates for starting values (this can be an option) and sigma2
  model_est_cc = cc_censored(formula, cens_ind, data, par_vec, starting_vals,
                             sandwich_se = FALSE)
  starting_vals = model_est_cc$beta_est %>% as.numeric()
  sigma2 = model_est_cc$sigma_est

  # set up multiroot function (the estimating equation we want to find the root of)
  multiroot_func = function(beta_temp, data,
                            Y, varNamesRHS, par_vec, cens_name, weights_cov, cens_ind,
                            m_func, integral_func_psi,
                            integral_func_denom, mu_joint, Sigma_joint, sigma2){
    pieces = apply(data, 1, function(temp){
      p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
      ipw_piece = rep(temp[cens_ind]*temp["weights"], length(beta_temp))*
        numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
        rep(temp[Y]-m_func(p), length(beta_temp))
      aipw_piece = rep(1 - temp[cens_ind]*temp["weights"], length(beta_temp))*
        psi_hat_i(temp, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                  beta_temp, m_func, integral_func_psi,
                  integral_func_denom, mu_joint, Sigma_joint, sigma2)
      ipw_piece + aipw_piece
    }) %>% unname()
    rowSums(pieces)
  }

  beta_est = rootSolve::multiroot(multiroot_func,
                                  data = data,
                                  Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                  cens_name = cens_name, weights_cov = weights_cov, cens_ind = cens_ind,
                                  m_func = m_func, integral_func_psi = integral_func_psi,
                                  integral_func_denom = integral_func_denom,
                                  mu_joint = mu_joint, Sigma_joint = Sigma_joint, sigma2 = sigma2,
                                  start = starting_vals)$root

  # run sandwich estimator (fix this)
  if(sandwich_se){
    se_est = aipw_sandwich(formula, data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                           beta_est, m_func, cens_ind, mu_joint, Sigma_joint, sigma2)
  }else{
    se_est = NULL
  }

  # save beta estimates
  return(list(beta_est = beta_est,
              se_est = se_est))
}

aipw_sandwich <- function(formula, data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                          beta_est, m_func, cens_ind, mu_joint, Sigma_joint, sigma2){

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
  exp_nobracket <- exp %>% as.character() %>% str_remove_all(., "\\[|\\]") %>%
    parse(text = .)
  m_func = function(p){
    with(as.list(p),
         eval(exp_nobracket))
  }

  # create "g" function for the sandwich estimator
  g = function(data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
               beta_est, m_func, cens_ind, mu_joint, Sigma_joint, sigma2){
    p = c(beta_est, data[varNamesRHS])
    names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

    ipw_piece = rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
      numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
      rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()
    aipw_piece = rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
      psi_hat_i(data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                beta_est, m_func, integral_func_psi,
                integral_func_denom, mu_joint, Sigma_joint, sigma2)

    ipw_piece + aipw_piece
  }

  # first derivative function
  # calculates first derivative for subject i (data x)
  # f(x;beta) where beta is of length lb, x is a scalar
  # first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
  firstderivative <- function(beta, g, data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                              m_func, cens_ind, mu_joint, Sigma_joint, sigma2){
    lb <- length(beta)
    derivs <- matrix(data = 0, nrow = lb, ncol = lb)
    delta <- beta * (10 ^ (- 4))
    betal <- betar <- beta
    for (i in 1:lb) {
      # Perturb the ith element of beta
      betal[i] <- beta[i] - delta[i]
      betar[i] <- beta[i] + delta[i]

      # Calculate function values
      yout1 <- g(data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                 betal, m_func, cens_ind, mu_joint, Sigma_joint, sigma2)
      yout2 <- g(data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                 betar, m_func, cens_ind, mu_joint, Sigma_joint, sigma2)

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
    firstderivative(beta_est, g, temp, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                    m_func, cens_ind, mu_joint, Sigma_joint, sigma2)
  })
  if(length(beta_est) > 1){
    first_der = first_der %>% rowMeans() %>% matrix(nrow = length(beta_est))
  }else{
    first_der = first_der %>% mean()
  }
  inv_first_der <- solve(first_der)

  # need to get the outer product of g at each observation and take the mean
  gs = apply(data, 1, function(temp)
    g(temp, Y, varNamesRHS, par_vec, cens_name, weights_cov,
      beta_est, m_func, cens_ind, mu_joint, Sigma_joint, sigma2))
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



