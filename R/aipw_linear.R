#' Augmented Inverse Probability Weighting for Censored Covariates with a Linear Model
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
#' @param weight_opt a character string indicating the method of weight calculation. One of "Cox", "AFT_lognormal", "MVN", "user" (if "user", then user provides weights).
#' @param weights_user if \code{weight_opt = "user"}, a vector of weights the same length as there are rows in \code{data}, otherwise \code{NULL} (default).
#' @param weights_cov if \code{weight_opt} one of \code{c("Cox", "AFT_lognormal", "MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model. Otherwise \code{NULL}.
#' @param weights_threshold the maximum weight for any one observation. If \code{NULL} (default), there is no thresholding.
#' @param cov_vars if \code{cov_dist_opt} one of \code{c("MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the covariate distribution Otherwise \code{NULL}.
#' @param ... additional arguments passed to function \code{multiroot}.
#'
#' @return A list with the following elements:
#' \item{beta_est}{a vector of the parameter estimates.}
#' \item{se_est}{if \code{sandwich_se = TRUE}, a vector of standard error estimates from the empirical sandwich estimator. Otherwise, \code{NULL}}
#' \item{iteration_count}{the number of iterations used in \code{multiroot}.}
#'
#' @import tidyverse
#' @import numDeriv
#' @import survival
#' @import rootSolve
#'
#' @export
aipw_censored_linear <- function(formula,
                          data,
                          cens_ind,
                          cens_name,
                          par_vec,
                          starting_vals,
                          sandwich_se = TRUE,
                          weight_opt,
                          weights_user = NULL,
                          weights_cov = NULL,
                          weights_threshold = NULL,
                          cov_vars,
                          ...){

  # Need to add error checks



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
  }

  # thresholding
  if(!is.null(weights_threshold)){
    weights = ifelse(weights > weights_threshold, weights_threshold, weights)
  }

  print("Weights complete!")

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

  # Find the moments for X | Y, Z
  ## want to estimate the parameters using AFT
  aft_formula <- as.formula(paste("survival::Surv(", cens_name, ", ", cens_ind, ") ~",
                                  paste(colnames(data %>% select(all_of(Y), all_of(cov_vars))),
                                        collapse = "+")))
  model_est_x_yz = survreg(aft_formula,
                          data = data,
                          dist = "lognormal")
  model_est_x_yz_coeff = model_est_x_yz$coefficients
  model_est_x_yz_sd = model_est_x_yz$scale
  x_yz_dist_params = list(model_est_x_yz_coeff = model_est_x_yz_coeff,
                         model_est_x_yz_sd = model_est_x_yz_sd)

  # Note: the A function is as follows
  # set p --- p = c(beta_temp, data[varNamesRHS])
  # set names for p ---- names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  # A = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]

  # CC to get estimates for starting values (this can be an option) and sigma2
  model_est_cc = cc_censored(formula, cens_ind, data, par_vec, starting_vals,
                             sandwich_se = FALSE)
  starting_vals = model_est_cc$beta_est %>% as.numeric()
  sigma2 = model_est_cc$sigma_est

  #### Find Psi for the CC Estimator
  # psi_all = apply(data, 1, function(temp){
  #   # print("hi")
  #   psi_hat_i_mvn(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
  #                 starting_vals, m_func, cov_dist_params$mu_joint,
  #                 cov_dist_params$Sigma_joint, sigma2)
  # }) %>% t() %>% as.data.frame()
  # colnames(psi_all) = paste0("psi", seq(1:length(starting_vals)))
  # data = cbind(data, psi_all)

  if(endsWith(cov_dist_opt, "MVN")){
    multiroot_results = rootSolve::multiroot(multiroot_func_mvn_linear,
                                             data = data,
                                             Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                             cens_name = cens_name, cov_vars = cov_vars, cens_ind = cens_ind,
                                             m_func = m_func,x_yz_dist_params = x_yz_dist_params,
                                             start = starting_vals, ...)
  }else if(cov_dist_opt == "AFT"){
    multiroot_results = rootSolve::multiroot(multiroot_func_aft,
                                             data = data,
                                             Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                             cens_name = cens_name, cov_vars = cov_vars, cens_ind = cens_ind,
                                             m_func = m_func, model_est_x_z_coeff = model_est_x_z_coeff,
                                             model_est_x_z_sd = model_est_x_z_sd, sigma2 = sigma2,
                                             start = starting_vals, ...)
  }

  beta_est = multiroot_results$root
  names(beta_est) = paste0(par_vec, seq(1:length(beta_est)))

  iteration_count = multiroot_results$iter

  print("Parameters estimated!")

  # run sandwich estimator
  if(sandwich_se){
    se_est = aipw_sandwich(formula, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                           beta_est, m_func, cens_ind, cov_dist_params, sigma2,
                           cov_dist_opt)
    print("Standard Errors estimated!")
  }else{
    se_est = NULL
  }



  # save beta estimates
  return(list(beta_est = beta_est,
              se_est = se_est,
              iteration_count = iteration_count))
}

aipw_sandwich <- function(formula, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                          beta_est, m_func, cens_ind, cov_dist_params, sigma2,
                          cov_dist_opt){

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
  if(endsWith(cov_dist_opt, "MVN")){
    g = function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 beta_est, m_func, cens_ind, cov_dist_params, sigma2){
      p = c(beta_est, data[varNamesRHS])
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      ipw_piece = rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
        rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()
      aipw_piece = rep(1 - data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
        psi_hat_i_mvn(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_est, m_func, cov_dist_params$mu_joint, cov_dist_params$Sigma_joint, sigma2)

      ipw_piece + aipw_piece
    }
  }else if(cov_dist_opt == "AFT"){
    g = function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 beta_est, m_func, cens_ind, cov_dist_params, sigma2){
      p = c(beta_est, data[varNamesRHS])
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      ipw_piece = rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
        rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()
      aipw_piece = rep(1 - data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
        psi_hat_i_aft(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_est, m_func, cov_dist_params$model_est_x_z_coeff,
                      cov_dist_params$model_est_x_z_sd, sigma2)

      ipw_piece + aipw_piece
    }

  }

  # first derivative function
  # calculates first derivative for subject i (data x)
  # f(x;beta) where beta is of length lb, x is a scalar
  # first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
  firstderivative <- function(beta, g, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              m_func, cens_ind, cov_dist_params, sigma2){
    lb <- length(beta)
    derivs <- matrix(data = 0, nrow = lb, ncol = lb)
    delta <- beta * (10 ^ (- 4))
    betal <- betar <- beta
    for (i in 1:lb) {
      # Perturb the ith element of beta
      betal[i] <- beta[i] - delta[i]
      betar[i] <- beta[i] + delta[i]

      # Calculate function values
      yout1 <- g(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 betal, m_func, cens_ind, cov_dist_params, sigma2)
      yout2 <- g(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 betar, m_func, cens_ind, cov_dist_params, sigma2)

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
    firstderivative(beta_est, g, temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                    m_func, cens_ind, cov_dist_params, sigma2)
  })
  if(length(beta_est) > 1){
    first_der = first_der %>% rowMeans() %>% matrix(nrow = length(beta_est))
  }else{
    first_der = first_der %>% mean()
  }
  inv_first_der <- solve(first_der)

  # need to get the outer product of g at each observation and take the mean
  gs = apply(data, 1, function(temp)
    g(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
      beta_est, m_func, cens_ind, cov_dist_params, sigma2))
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



