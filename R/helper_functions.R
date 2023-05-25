
#############################
##### weighting options #####
#############################

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

  return_list = list(weights = weights, mu_joint = mu_joint, Sigma_joint = Sigma_joint)

  return(return_list)
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

#####################################
##### helper functions for AIPW #####
#####################################

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
    psi = numerator/denominator
  }
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

#####################################
##### helper functions for MLE ######
#####################################

psi_hat_i_mle <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_temp, m_func, mu_joint, Sigma_joint, sigma2){
  denominator =  integrate(integral_func_denom_mle, data[cens_name] %>% as.numeric(),
                           Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, cov_vars = cov_vars,
                           beta_temp = beta_temp, m_func = m_func,
                           C_val = data[cens_name]  %>% as.numeric(),
                           mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-4){
    psi = 0
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi_mle, data[cens_name] %>% as.numeric(),
                Inf, data_row = data, Y = Y,
                varNamesRHS = varNamesRHS, par_vec = par_vec,
                cens_name = cens_name, cov_vars = cov_vars,
                beta_temp = beta_temp, m_func = m_func,
                C_val = data[cens_name]  %>% as.numeric(),
                mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                sigma2 = sigma2, j = j,
                rel.tol = .Machine$double.eps^0.1,
                subdivisions = 1000)$value
    }) %>% unlist() %>% as.numeric()
    psi = numerator/denominator
  }
  #print(beta_temp)
  psi
}

integral_func_denom_mle <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                beta_temp, m_func, C_val,
                                mu_joint, Sigma_joint, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_cz = condMVNorm::dcmvnorm(x = log(data_row[cens_name] %>% as.numeric())[1],
                                 mean = mu_joint,
                                 sigma = Sigma_joint,
                                 dependent.ind = 1,
                                 given.ind = c(2, 3),
                                 X.given = c(log(C_val),
                                             data_row[cov_vars] %>% as.numeric()))
    value_ts[i] = f_y*f_x_cz/t[i]
  }
  value_ts
}

integral_func_psi_mle <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              beta_temp, m_func, C_val,
                              mu_joint, Sigma_joint, j, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_cz = condMVNorm::dcmvnorm(x = log(data_row[cens_name] %>% as.numeric())[1],
                                  mean = mu_joint,
                                  sigma = Sigma_joint,
                                  dependent.ind = 1,
                                  given.ind = c(2, 3),
                                  X.given = c(log(C_val),
                                              data_row[cov_vars] %>% as.numeric()))
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y]-m_t)*f_y*f_x_cz/t[i]
  }
  value_ts %>% unlist()
}











