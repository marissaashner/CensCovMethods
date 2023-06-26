#' MLE for Censored Covariates
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
#' @param cov_dist_opt a character string indicating which method of covariate distribution is to be done. One of "MVN", "user MVN", "AFT"
#' @param cov_vars if \code{cov_dist_opt} one of \code{c("MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the covariate distribution Otherwise \code{NULL}
#' @param cov_mean_user if \code{cov_dis_opt = "user MVN"}, the mean of the multivariate normal distribution of \code{(log(X), log(C), Z)}
#' @param cov_sigma_user if \code{cov_dis_opt = "user MVN"}, the covariance matrix of the multivariate normal distributionof \code{(log(X), log(C), Z)}
#' @param ... additional arguments passed to function \code{multiroot}.
#'
#' @return A list with the following elements:
#' \item{beta_est}{a vector of the parameter estimates.}
#' \item{se_est}{if \code{sandwich_se = TRUE}, a vector of standard error estimates from the empirical sandwich estimator. Otherwise, \code{NULL}}
#' \item{iteration_count}{the number of iterations used in \code{multiroot}.}
#'
#' @import tidyverse
#' @import rootSolve
#' @import survival
#' @import numDeriv
#'
#' @export
mle_censored <- function(formula,
                          data,
                          cens_ind,
                          cens_name,
                          par_vec,
                          starting_vals,
                          sandwich_se = TRUE,
                         cov_dist_opt = "MVN",
                         cov_vars,
                         cov_mean_user = NULL,
                         cov_sigma_user = NULL,
                         ...){

  # Need to add error checks

  # covariate distribution
  if(cov_dist_opt == "MVN"){
    mvn_results = weights_mvn(data, cens_ind, cov_vars, cens_name)
    mu_joint = mvn_results$mu_joint
    Sigma_joint = mvn_results$Sigma_joint
    x_cz_dist_params = list(mu_joint = mu_joint,
                            Sigma_joint = Sigma_joint)
  }else if(cov_dist_opt == "user_MVN"){
    mu_joint = cov_mean_user
    Sigma_joint = cov_sigma_user
    x_cz_dist_params = list(mu_joint = mu_joint,
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
    x_cz_dist_params = list(model_est_x_z_coeff = model_est_x_z_coeff,
                            model_est_x_z_sd = model_est_x_z_sd)
  }

  print("Covariate Distribution Estimated!")

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
  sigma2 = model_est_cc$sigma_est

  if(endsWith(cov_dist_opt,"MVN")){
    multiroot_results = rootSolve::multiroot(multiroot_func_mle_mvn,
                                             data = data,
                                             Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                             cens_name = cens_name, cov_vars = cov_vars, cens_ind = cens_ind,
                                             m_func = m_func, mu_joint = mu_joint,
                                             Sigma_joint = Sigma_joint, sigma2 = sigma2,
                                             start = starting_vals, ...)
  }else if(cov_dist_opt == "AFT_lognormal"){
    multiroot_results = rootSolve::multiroot(multiroot_func_mle_aft,
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

  print("Parameters Estimated!")

  # run sandwich estimator
  if(sandwich_se){
    se_est = mle_sandwich(formula, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                           beta_est, m_func, cens_ind, x_cz_dist_params, sigma2,
                          cov_dist_opt)
  }else{
    se_est = NULL
  }

  # save beta estimates
  return(list(beta_est = beta_est,
              se_est = se_est,
              iteration_count = iteration_count))
}


#' Sandwich Estimator for Maximum Likelihood for Censored Covariates with Gauss Hermite Quadrature
#'
#' Calculates the standard error sandwich estimate for an MLE estimator for a regression model with censored covariates.
#'
#' @param formula a linear or nonlinear model formula including variables and parameters
#' @param data a data frame containing columns for the censoring indicator and the variables in \code{formula}
#' @param Y a character string indicating the name of the outcome from \code{data}
#' @param varNamesRHS a vector of character strings indicating the names of the variables on the right side of \code{formula}
#' @param par_vec a character string indicating the parameter vector in the formula
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#' @param cov_vars if \code{cov_dist_opt} one of \code{c("MVN")}, a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the covariate distribution Otherwise \code{NULL}.
#' @param beta_est the estimate from the augmented complete case estimator
#' @param m_func the mean function of the regression model
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param x_cz_dist_params a list of parameters for the conditional distribution of the censored covariate given the censoring value and fully observed covariates
#' @param sigma2 the estimate of the error variance
#' @param cov_dist_opt a character string indicating specification of the covariate distribution. One of "MVN", "user MVN", "AFT"
#'
#' @return A vector of the sandwich standard error estimates.
#'
#' @import tidyverse
#' @importFrom rootSolve multiroot
#' @import survival
#' @import numDeriv
#'
#' @export
mle_sandwich <- function(formula, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                          beta_est, m_func, cens_ind, x_cz_dist_params, sigma2,
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
                 beta_est, m_func, cens_ind, x_cz_dist_params, sigma2){
      p = c(beta_est, data[varNamesRHS])  %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      if(data[cens_ind] == 1){
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
          rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
      }else{
        psi_hat_i_mle_mvn(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_est, m_func, x_cz_dist_params$mu_joint,
                      x_cz_dist_params$Sigma_joint, sigma2)
      }
    }
  }else if(cov_dist_opt == "AFT"){
    g = function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 beta_est, m_func, cens_ind, x_cz_dist_params, sigma2){
      p = c(beta_est, data[varNamesRHS])  %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

      if(data[cens_ind] == 1){
        numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
          rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
      }else{
        psi_hat_i_mle_aft(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_est, m_func, x_cz_dist_params$model_est_x_z_coeff,
                      x_cz_dist_params$model_est_x_z_sd, sigma2)
      }
    }
  }


  # first derivative function
  # calculates first derivative for subject i (data x)
  # f(x;beta) where beta is of length lb, x is a scalar
  # first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
  firstderivative <- function(beta, g, data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              m_func, cens_ind, x_cz_dist_params, sigma2){
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
                 betal, m_func, cens_ind, x_cz_dist_params, sigma2)
      yout2 <- g(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                 betar, m_func, cens_ind, x_cz_dist_params, sigma2)

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
                    m_func, cens_ind, x_cz_dist_params, sigma2)
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
      beta_est, m_func, cens_ind, x_cz_dist_params, sigma2))
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
