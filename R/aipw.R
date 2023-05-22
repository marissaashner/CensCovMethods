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

  print("TEST1")

  # add weights to data frame
  data$weights = weights

  print("TEST2")

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

  print("TEST3")

  # Note: the A function is as follows
  # set p --- p = c(beta_temp, data[varNamesRHS])
  # set names for p ---- names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  # A = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]

  # CC to get estimates for starting values (this can be an option) and sigma2
  model_est_cc = cc_censored(formula, cens_ind, data, par_vec, starting_vals,
                             sandwich_se = FALSE)
  starting_vals = model_est_cc$beta_est %>% as.numeric()
  sigma2 = model_est_cc$sigma_est

  print("TEST4")

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

  print("TEST5")

  beta_est = rootSolve::multiroot(multiroot_func,
                                  data = data,
                                  Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                  cens_name = cens_name, weights_cov = weights_cov, cens_ind = cens_ind,
                                  m_func = m_func, integral_func_psi = integral_func_psi,
                                  integral_func_denom = integral_func_denom,
                                  mu_joint = mu_joint, Sigma_joint = Sigma_joint, sigma2 = sigma2,
                                  start = starting_vals)$root

  print("TEST6")

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

# weighting options
weights_cox <- function(data, cens_ind, weights_cov, cens_name){
  cox_formula <- as.formula(paste("survival::Surv(", cens_name, ", 1-", cens_ind, ") ~",
                                  paste(colnames(data %>% select(all_of(weights_cov))),
                                        collapse = "+")))
  cox_fit <- survival::coxph(cox_formula, data = data)
  G <- survival::survfit(cox_fit, newdata = data)

  ## get Ghat estimates for each W in the observed data and join with the original df
  weights_data <- data.frame(W = summary(G, times = data[cens_name] %>% unlist(), extend = TRUE)$time,
                             weights = 1/(summary(G, times = data[cens_name] %>% unlist(), extend = TRUE)$surv %>% diag()))
  colnames(weights_data)[1] = cens_name
  data <- data %>% left_join(weights_data, by = cens_name)
  return(data$weights)
}

weights_aft <- function(data, cens_ind, weights_cov, cens_name){
  aft_formula <- as.formula(paste("survival::Surv(", cens_name, ", 1-", cens_ind, ") ~",
                                  paste(colnames(data %>% select(all_of(weights_cov))),
                                        collapse = "+")))
  model_est_c_z = survival::survreg(aft_formula,
                                    data = data,
                                    dist = "lognormal")
  model_est_c_z_coeff = model_est_c_z$coefficients
  model_est_c_z_sd = model_est_c_z$scale
  Z = data %>% select(all_of(weights_cov)) %>% as.matrix()
  Z = cbind(rep(1, ncol(Z)), Z)
  weights = 1/pnorm(log(data[cens_name] %>% unlist()), mean = Z %*% model_est_c_z_coeff,
                    sd = model_est_c_z_sd,
                    lower.tail = FALSE)
  return(weights)
}

weights_mvn <- function(data, cens_ind, weights_cov, cens_name){
  l_all = function(params, data){
    -(apply(data, 1, function(dat)
      l(params, log(dat[cens_name]), dat[weights_cov], dat[cens_ind])) %>% sum())
  }

  starting_vals_xcz = c(0,0,0,1,1,1,0,0,0)
  params_est <- optim(starting_vals_xcz, l_all, data = data,
                      method = "L-BFGS-B",
                      lower = c(-Inf, -Inf, -Inf, 1e-4,1e-4,1e-4,-Inf,-Inf,-Inf),
                      upper = rep(Inf, 9))
  params =params_est$par
  mu_joint = c(params[1], params[2], params[3])
  Sigma_joint = (matrix(c(params[4], params[7], params[8],
                          params[7], params[5], params[9],
                          params[8], params[9], params[6]),
                        nrow = 3))

  weights = 1/apply(data, 1, function(dat_row)
    condMVNorm::pcmvnorm(lower = log(dat_row[cens_name]), upper = Inf,
                         mean = mu_joint, sigma = Sigma_joint,
                         dependent.ind = 2,
                         given = c(1,3),
                         X.given = c(log(dat_row[cens_name]), dat_row[weights_cov])))
  return(list(weights = weights, mu_joint = mu_joint, Sigma_joint = Sigma_joint))
}

# want to make warnings go away on the MVN part

# helper functions for the MVN weights function
l = function(params, W, Z, D){
  mu = c(params[1], params[2], params[3])
  Sigma = (matrix(c(params[4], params[7], params[8],
                    params[7], params[5], params[9],
                    params[8], params[9], params[6]),
                  nrow = 3))

  if(D == 0){
    params_new = cond_normal_params(mu, Sigma, 1, c(W,Z))
    loglike = log(mvtnorm::dmvnorm(x=c(W, Z),
                                   mean = mu[2:3], sigma = Sigma[2:3,2:3])*
                    pnorm(W, mean = params_new$mu_new,
                          sd = sqrt(params_new$sig_new), lower.tail = FALSE))
    if(is.na(loglike) | loglike == -Inf){
      -10000
    }else{
      loglike
    }
  }else{
    params_new = cond_normal_params(mu, Sigma, 2, c(W,Z))
    loglike = log(mvtnorm::dmvnorm(x=c(W, Z),
                                   mean = mu[c(1,3)], sigma = Sigma[c(1,3),c(1,3)])*
                    pnorm(W, mean = params_new$mu_new,
                          sd = sqrt(params_new$sig_new), lower.tail = FALSE))
    if(is.na(loglike) | loglike == -Inf){
      -10000
    }else{
      loglike
    }
  }
}

cond_normal_params = function(mu, sigma, dependent.ind, X.given){
  sig_12 = sigma[-dependent.ind, dependent.ind] %>% as.matrix()
  sig_11 = sigma[dependent.ind, dependent.ind] %>% as.matrix()
  sig_22 = sigma[-dependent.ind, -dependent.ind] %>% as.matrix()
  mu_new = mu[dependent.ind] +
    t(sig_12)%*%solve(sig_22)%*%(X.given - mu[-dependent.ind])
  sig_new = sig_11 - t(sig_12)%*%solve(sig_22)%*%(sig_12)
  return(list(mu_new = mu_new, sig_new = sig_new))
}

# helper functions for AIPW
psi_hat_i <- function(data, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                      beta_temp, m_func, integral_func_psi,
                      integral_func_denom, mu_joint, Sigma_joint, sigma2){
  denominator =  integrate(integral_func_denom, 0, Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, weights_cov = weights_cov,
                           beta_temp = beta_temp, m_func = m_func,
                           mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-4){
    psi = 0
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi, 0, Inf, data_row = data, Y = Y,
                varNamesRHS = varNamesRHS, par_vec = par_vec,
                cens_name = cens_name, weights_cov = weights_cov,
                beta_temp = beta_temp, m_func = m_func,
                mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                sigma2 = sigma2, j = j,
                rel.tol = .Machine$double.eps^0.1,
                subdivisions = 1000)$value
    }) %>% unlist() %>% as.numeric()
  }
  psi = numerator/denominator
  #print(beta_temp)
  psi
}

integral_func_denom <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                                beta_temp, m_func,
                                mu_joint, Sigma_joint, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_z = condMVNorm::dcmvnorm(x = log(data_row[cens_name] %>% as.numeric())[1],
                                 mean = mu_joint,
                                 sigma = Sigma_joint,
                                 dependent.ind = 1,
                                 given.ind = c(3),
                                 X.given = c(data_row[weights_cov] %>% as.numeric()))
    value_ts[i] = f_y*f_x_z/t[i]
  }
  value_ts
}

integral_func_psi <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, weights_cov,
                              beta_temp, m_func,
                              mu_joint, Sigma_joint, j, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_z = condMVNorm::dcmvnorm(x = log(data_row[cens_name] %>% as.numeric())[1],
                                 mean = mu_joint,
                                 sigma = Sigma_joint,
                                 dependent.ind = 1,
                                 given.ind = c(3),
                                 X.given = c(data_row[weights_cov] %>% as.numeric()))
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y]-m_t)*f_y*f_x_z/t[i]
  }
  value_ts %>% unlist()
}


