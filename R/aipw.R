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
#' @param se_opt a character string indicating if the weights are assumed to be known or if they are estimated. One of \code{c("known", "est_MVN")}.
#' @param weight_opt a character string indicating the method of weight calculation. One of "Cox", "AFT_lognormal", "MVN", "user" (if "user", then user provides weights).
#' @param weights_user if \code{weight_opt = "user"}, a vector of weights the same length as there are rows in \code{data}, otherwise \code{NULL} (default).
#' @param weights_cov if \code{weight_opt} one of \code{c("Cox", "AFT_lognormal", "MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model. Otherwise \code{NULL}.
#' @param weights_threshold the maximum weight for any one observation. If \code{NULL} (default), there is no thresholding.
#' @param weight_stabilize a character string indicating which method of weight stabilization is to be done (if any). One of \code{c("Mean", "KM", "None")}.
#' @param cov_dist_opt a character string indicating specification of the covariate distribution. One of "MVN", "user_MVN", "AFT_lognormal"
#' @param cov_vars if \code{cov_dist_opt} one of \code{c("MVN", "AFT_lognormal")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the covariate distribution Otherwise \code{NULL}.
#' @param cov_mean_user if \code{cov_dis_opt = "user_MVN"}, the mean of the multivariate normal distribution of \code{(log(X), log(C), Z)}.
#' @param cov_sigma_user if \code{cov_dis_opt = "user_MVN"}, the covariance matrix of the multivariate normal distribution of \code{(log(X), log(C), Z)}.
#' @param gh if \code{TRUE} (default), gauss-hermite quadrature will be run to estimate the integrals. Otherwise, the \code{integrate} function will be used.
#' @param gh_nodes number of nodes to use in gauss-hermite quadrature.
#' @param ... additional arguments passed to function \code{multiroot}.
#'
#' @return A list with the following elements:
#' \item{beta_est}{a vector of the parameter estimates.}
#' \item{se_est}{if \code{sandwich_se = TRUE}, a vector of standard error estimates from the empirical sandwich estimator. Otherwise, \code{NULL}}
#' \item{iteration_count}{the number of iterations used in \code{multiroot}.}
#'
#' @import tidyverse
#' @importFrom rootSolve multiroot
#' @import survival
#' @import numDeriv
#' @import statmod
#'
#' @export
aipw_censored <- function(formula,
                         data,
                         cens_ind,
                         cens_name,
                         par_vec,
                         starting_vals,
                         sandwich_se = TRUE,
                         se_opt = "known",
                         weight_opt,
                         weights_user = NULL,
                         weights_cov = NULL,
                         weights_threshold = NULL,
                         weight_stabilize = "None",
                         cov_dist_opt = "MVN",
                         cov_vars,
                         cov_mean_user = NULL,
                         cov_sigma_user = NULL,
                         gh = TRUE,
                         gh_nodes = 10,
                         ...){

  # Need to add error checks

# testing

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

  #stabilize weights
  if(weight_stabilize == "KM"){
    km_formula = as.formula(paste("survival::Surv(", cens_name, ", 1-", cens_ind, ") ~ 1"))
    km_fit = survival::survfit(km_formula, data = data)
    km_data <- data.frame(W = summary(km_fit, times = data[cens_name] %>% unlist(), extend = TRUE)$time,
                          surv_km = (summary(km_fit, times = data[cens_name] %>% unlist(), extend = TRUE)$surv))
    colnames(km_data)[1] = cens_name
    data <- data %>% dplyr::left_join(km_data %>% unique(), by = cens_name)
    weights = weights*data$surv_km
  }else if(weight_stabilize == "Mean"){
    weights = weights*mean(data[cens_ind] %>% unlist())
  }else if(weight_stabilize == "MVN"){
    weights = weights*pnorm(log(data[cens_name] %>% unlist()),
                            mean = mvn_results$params[2],
                            sd = sqrt(mvn_results$params[5]), lower.tail = FALSE)
  }


  # thresholding
  if(!is.null(weights_threshold)){
    weights = ifelse(weights > weights_threshold, weights_threshold, weights)
  }

  print("Weights complete!")

  # covariate distribution
  if(cov_dist_opt == "MVN" & weight_opt != "MVN"){
    mvn_results = weights_mvn(data, cens_ind, cov_vars, cens_name)
    mu_joint = mvn_results$mu_joint
    Sigma_joint = mvn_results$Sigma_joint
    cov_dist_params = list(mu_joint = mu_joint,
                           Sigma_joint = Sigma_joint)
  }else if(cov_dist_opt == "MVN" & weight_opt == "MVN"){
    mu_joint = mvn_results$mu_joint
    Sigma_joint = mvn_results$Sigma_joint
    cov_dist_params = list(mu_joint = mu_joint,
                           Sigma_joint = Sigma_joint)
  }else if(cov_dist_opt == "user_MVN"){
    mu_joint = cov_mean_user
    Sigma_joint = cov_sigma_user
    cov_dist_params = list(mu_joint = mu_joint,
                           Sigma_joint = Sigma_joint)
  }else if(cov_dist_opt == "AFT_lognormal"){
    ## want to estimate the parameters using AFT
    aft_formula <- as.formula(paste("survival::Surv(", cens_name, ", ", cens_ind, ") ~",
                                    paste(colnames(data %>% dplyr::ungroup() %>%
                                                     dplyr::select(all_of(cov_vars))),
                                          collapse = "+")))
    model_est_x_z = survreg(aft_formula,
                            data = data,
                            dist = "lognormal")
    model_est_x_z_coeff = model_est_x_z$coefficients
    model_est_x_z_sd = model_est_x_z$scale
    cov_dist_params = list(model_est_x_z_coeff = model_est_x_z_coeff,
                           model_est_x_z_sd = model_est_x_z_sd)
  }

  print("Covariate distribution parameters complete!")

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
  exp_nobracket <- exp %>% as.character() %>% stringr::str_remove_all(., "\\[|\\]") %>%
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
  model_est_cc = cc_censored(formula, data, cens_ind, par_vec, starting_vals,
                             sandwich_se = FALSE)
  starting_vals = model_est_cc$beta_est %>% as.numeric()
  sigma2 = model_est_cc$sigma_est^2

  if(endsWith(cov_dist_opt, "MVN") & !gh){
    multiroot_results = rootSolve::multiroot(multiroot_func_mvn,
                                             data = data,
                                             Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                             cens_name = cens_name, cov_vars = cov_vars, cens_ind = cens_ind,
                                             m_func = m_func, mu_joint = mu_joint,
                                             Sigma_joint = Sigma_joint, sigma2 = sigma2,
                                             start = starting_vals, ...)
  }else if(cov_dist_opt == "AFT_lognormal" & !gh){
    multiroot_results = rootSolve::multiroot(multiroot_func_aft,
                                             data = data,
                                             Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                             cens_name = cens_name, cov_vars = cov_vars, cens_ind = cens_ind,
                                             m_func = m_func, model_est_x_z_coeff = model_est_x_z_coeff,
                                             model_est_x_z_sd = model_est_x_z_sd, sigma2 = sigma2,
                                             start = starting_vals, ...)
  }else if(gh){
    multiroot_results = rootSolve::multiroot(multiroot_func_hermite_aipw,
                                             data = data,
                                             Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                             cens_name = cens_name, cov_vars = cov_vars, cens_ind = cens_ind,
                                             m_func = m_func, cov_dist_params = cov_dist_params,
                                             sigma2 = sigma2, gh_nodes = gh_nodes,
                                             start = starting_vals, ...)
  }

  beta_est = multiroot_results$root
  names(beta_est) = paste0(par_vec, seq(1:length(beta_est)))

  iteration_count = multiroot_results$iter

  print("Parameters estimated!")

  if(se_opt == "known"){
    params = NULL
  }else if(se_opt == "est_MVN"){
    params = mvn_results$params
  }

  # run sandwich estimator
  if(sandwich_se){
    if(!gh){
      se_est = aipw_sandwich(formula, data, par_vec, cens_name, cov_vars, weights_cov,
                             beta_est, cens_ind, cov_dist_params, sigma2,
                             cov_dist_opt)
    }else{
        se_est = aipw_sandwich_hermite(formula, data, par_vec, cens_name, cov_vars, weights_cov,
                                       beta_est, cens_ind, cov_dist_params, sigma2,
                                       cov_dist_opt, gh_nodes, params, se_opt)
    }

    print("Standard Errors estimated!")
  }else{
    se_est = NULL
  }



  # save beta estimates
  return(list(beta_est = beta_est,
              se_est = se_est,
              iteration_count = iteration_count))
}


#' Sandwich Estimator for Augmented Inverse Probability Weighting for Censored Covariates
#'
#' Calculates the standard error sandwich estimate for an AIPW estimator for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters
#' @param data a data frame containing columns for the censoring indicator and the variables in \code{formula}
#' @param par_vec a character string indicating the parameter vector in the formula
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#' @param cov_vars if \code{cov_dist_opt} one of \code{c("MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the covariate distribution Otherwise \code{NULL}.
#' @param weights_cov if \code{weight_opt} one of \code{c("Cox", "AFT_lognormal", "MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model. Otherwise \code{NULL}.
#' @param beta_est the estimate from the augmented complete case estimator
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param cov_dist_params a list of parameters for the conditional distribution of the censored covariate given the fully observed covariates
#' @param sigma2 the estimate of the error variance
#' @param cov_dist_opt a character string indicating specification of the covariate distribution. One of "MVN", "user_MVN", "AFT"
#'
#' @return A vector of the sandwich standard error estimates.
#'
#' @import tidyverse
#' @import rootSolve
#' @import survival
#' @import numDeriv
#'
#' @export
aipw_sandwich <- function(formula, data, par_vec, cens_name, cov_vars, weights_cov,
                          beta_est, cens_ind, cov_dist_params, sigma2,
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
  exp_nobracket <- exp %>% as.character() %>% stringr::str_remove_all(., "\\[|\\]") %>%
    parse(text = .)
  m_func = function(p){
    with(as.list(p),
         eval(exp_nobracket))
  }

  # create "g" function for the sandwich estimator
  if(endsWith(cov_dist_opt, "MVN")){
    g = function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 beta_est, m_func, cens_ind, cov_dist_params, sigma2){
      p = c(beta_est, data[varNamesRHS]) %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      ipw_piece = rep(as.numeric(data[[cens_ind]])*as.numeric(data[["weights"]]), length(beta_est)) %>% as.numeric()*
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
        rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
      aipw_piece = rep(1 - as.numeric(data[[cens_ind]])*as.numeric(data[["weights"]]), length(beta_est)) %>% as.numeric()*
        psi_hat_i_mvn(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                  beta_est, m_func, cov_dist_params$mu_joint, cov_dist_params$Sigma_joint, sigma2)

      ipw_piece + aipw_piece
    }
  }else if(cov_dist_opt == "AFT"){
    g = function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 beta_est, m_func, cens_ind, cov_dist_params, sigma2){
      p = c(beta_est, data[varNamesRHS])  %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      ipw_piece = rep(as.numeric(data[[cens_ind]])*as.numeric(data[["weights"]]), length(beta_est)) %>% as.numeric()*
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
        rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
      aipw_piece = rep(1 - as.numeric(data[[cens_ind]])*as.numeric(data[["weights"]]), length(beta_est)) %>% as.numeric()*
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

#' Sandwich Estimator for Augmented Inverse Probability Weighting for Censored Covariates with Gauss Hermite Quadrature
#'
#' Calculates the standard error sandwich estimate for an AIPW estimator for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters
#' @param data a data frame containing columns for the censoring indicator and the variables in \code{formula}
#' @param par_vec a character string indicating the parameter vector in the formula
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#' @param cov_vars if \code{cov_dist_opt} one of \code{c("MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the covariate distribution Otherwise \code{NULL}.
#' @param weights_cov if \code{weight_opt} one of \code{c("Cox", "AFT_lognormal", "MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model. Otherwise \code{NULL}.
#' @param beta_est the estimate from the augmented complete case estimator
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param cov_dist_params a list of parameters for the conditional distribution of the censored covariate given the fully observed covariates
#' @param sigma2 the estimate of the error variance
#' @param cov_dist_opt a character string indicating specification of the covariate distribution. One of "MVN", "user_MVN", "AFT"
#' @param gh_nodes number of nodes to use in gauss-hermite quadrature.
#' @param params if \code{se_opt = "est_MVN"}, a vector of nine parameter estimates representing those estimated from the maximum likelihood technique.
#' @param se_opt a character string indicating if the weights are assumed to be known or if they are estimated. One of \code{c("known", "est_MVN")}.
#'
#' @return A vector of the sandwich standard error estimates.
#'
#' @import tidyverse
#' @import rootSolve
#' @import survival
#' @import numDeriv
#' @import statmod
#'
#' @export
aipw_sandwich_hermite <- function(formula, data, par_vec, cens_name, cov_vars, weights_cov,
                          beta_est, cens_ind, cov_dist_params, sigma2,
                          cov_dist_opt, gh_nodes, params = NULL, se_opt = "known"){

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
    g = function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 beta_est, m_func, cens_ind, cov_dist_params, sigma2,
                 gh_nodes){
      p = c(beta_est, data[varNamesRHS]) %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      ipw_piece = rep(as.numeric(data[[cens_ind]])*as.numeric(data[["weights"]]), length(beta_est)) %>% as.numeric()*
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
        rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
      aipw_piece = rep(1 - as.numeric(data[[cens_ind]])*as.numeric(data[["weights"]]), length(beta_est)) %>% as.numeric()*
        psi_hat_i_hermite_aipw(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_est, m_func, cov_dist_params, sigma2, gh_nodes)

      ipw_piece + aipw_piece
  }

  # first derivative function
  # calculates first derivative for subject i (data x)
  # f(x;beta) where beta is of length lb, x is a scalar
  # first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
  firstderivative <- function(beta, g, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              m_func, cens_ind, cov_dist_params, sigma2, gh_nodes){
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
                 betal, m_func, cens_ind, cov_dist_params, sigma2, gh_nodes)
      yout2 <- g(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 betar, m_func, cens_ind, cov_dist_params, sigma2, gh_nodes)

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
                    m_func, cens_ind, cov_dist_params, sigma2, gh_nodes)
  })
  if(length(beta_est) > 1){
    first_der = first_der %>% rowMeans() %>% matrix(nrow = length(beta_est))
  }else{
    first_der = first_der %>% mean()
  }
  inv_first_der <- solve(first_der)

  if(se_opt == "known"){
    # need to get the outer product of g at each observation and take the mean
    gs = apply(data, 1, function(temp)
      g(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
        beta_est, m_func, cens_ind, cov_dist_params, sigma2, gh_nodes))
    if(length(beta_est) > 1){
      outer_prod = apply(gs, 2, function(g) g%*%t(g))
      outer_prod = outer_prod %>% rowMeans() %>% matrix(nrow = length(beta_est))
    }else{
      outer_prod = gs^2
      outer_prod = outer_prod %>% mean()
    }
  }else if(se_opt == "est_MVN"){
    # take the first derivative of g_gamma
    first_der_gamma <- apply(data, 1, function(temp){
      firstderivative_g_gamma_aipw(params, beta_est, g_gamma_aipw, temp, varNamesRHS, par_vec, m_func,
                              cens_ind, cens_name, weights_cov, Y,
                              cov_vars, cov_dist_params, sigma2, gh_nodes)
    })
    if(length(beta_est) > 1){
      first_der_gamma = first_der_gamma %>% rowMeans() %>% matrix(nrow = length(beta_est))
    }else{
      first_der_gamma = first_der_gamma %>% mean()
    }

    # take the first derivative of f
    first_der_f <- apply(data, 1, function(temp){
      firstderivative_f_gamma(params, f_gamma, temp, cens_ind, cens_name, weights_cov)
    })
    if(length(beta_est) > 1){
      first_der_f = first_der_f %>% rowMeans() %>% matrix(nrow = length(params))
    }else{
      first_der_f = first_der_f %>% mean()
    }

    # need to get the outer product of g at each observation and take the mean
    gs = apply(data, 1, function(temp)
      g(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
        beta_est, m_func, cens_ind, cov_dist_params, sigma2, gh_nodes) -
        first_der_gamma%*%solve(first_der_f)%*%f_gamma(params, cens_name, weights_cov, cens_ind, temp))
    if(length(beta_est) > 1){
      outer_prod = apply(gs, 2, function(g) g%*%t(g))
      outer_prod = outer_prod %>% rowMeans() %>% matrix(nrow = length(beta_est))
    }else{
      outer_prod = gs^2
      outer_prod = outer_prod %>% mean()
    }
  }



  ## then need to put it all together
  se = sqrt((inv_first_der %*% outer_prod %*% t(inv_first_der) / nrow(data)) %>% diag())
  return(se)

}











