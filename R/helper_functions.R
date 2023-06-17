
#############################
##### weighting options #####
#############################

weights_cox <- function(data, cens_ind, weights_cov, cens_name){
  cox_formula <- as.formula(paste("survival::Surv(", cens_name, ", 1-", cens_ind, ") ~",
                                  paste(colnames(data %>% dplyr::select(all_of(weights_cov))),
                                        collapse = "+")))
  cox_fit <- survival::coxph(cox_formula, data = data)
  G <- survival::survfit(cox_fit, newdata = data)

  ## get Ghat estimates for each W in the observed data and join with the original df
  weights_data <- data.frame(W = summary(G, times = data[cens_name] %>% unlist(), extend = TRUE)$time,
                             weights_cox = 1/(summary(G, times = data[cens_name] %>% unlist(), extend = TRUE)$surv %>% diag()))
  colnames(weights_data)[1] = cens_name
  data <- data %>% left_join(weights_data, by = cens_name)
  return(data$weights_cox)
}

weights_aft <- function(data, cens_ind, weights_cov, cens_name){
  aft_formula <- as.formula(paste("survival::Surv(", cens_name, ", 1-", cens_ind, ") ~",
                                  paste(colnames(data %>% dplyr::ungroup() %>%
                                                   dplyr::select(all_of(weights_cov))),
                                        collapse = "+")))
  model_est_c_z = survival::survreg(aft_formula,
                                    data = data,
                                    dist = "lognormal")
  model_est_c_z_coeff = model_est_c_z$coefficients
  model_est_c_z_sd = model_est_c_z$scale
  Z = data %>% dplyr::ungroup() %>%
    dplyr::select(all_of(weights_cov)) %>% as.matrix()
  Z = cbind(rep(1, nrow(Z)), Z)
  weights = 1/pnorm(log(data[cens_name] %>% unlist()), mean = Z %*% model_est_c_z_coeff,
                    sd = model_est_c_z_sd,
                    lower.tail = FALSE)
  return(weights)
}

weights_mvn <- function(data, cens_ind, weights_cov, cens_name){
  l_all = function(params, data){
    -(apply(data, 1, function(dat)
      l(params, log(dat[cens_name] %>% as.numeric()),
        dat[weights_cov] %>% as.numeric(),
        dat[cens_ind] %>% as.numeric())) %>% sum())
  }

  starting_vals_xcz = c(0,0,0,1,1,1,0,0,0)
  params_est <- optim(starting_vals_xcz, l_all, data = data,
                      method = "L-BFGS-B",
                      lower = c(-Inf, -Inf, -Inf, 1e-4,1e-4,1e-4,-Inf,-Inf,-Inf),
                      upper = rep(Inf, 9))
  params = params_est$par
  mu_joint = c(params[1], params[2], params[3])
  Sigma_joint = (matrix(c(params[4], params[7], params[8],
                          params[7], params[5], params[9],
                          params[8], params[9], params[6]),
                        nrow = 3))

  weights = 1/apply(data, 1, function(dat_row)
    condMVNorm::pcmvnorm(lower = log(dat_row[cens_name]%>% as.numeric()), upper = Inf,
                         mean = mu_joint, sigma = Sigma_joint,
                         dependent.ind = 2,
                         given = c(1,3),
                         X.given = c(log(dat_row[cens_name]%>% as.numeric()),
                                     dat_row[weights_cov]%>% as.numeric())))

  return_list = list(weights = weights, mu_joint = mu_joint, Sigma_joint = Sigma_joint)

  return(return_list)
}

weights_aft_acc <- function(data, cens_ind, weights_cov, Y, cens_name){
  ### For DELTA = 1
  aft_formula_c <- as.formula(paste("survival::Surv(", cens_name, ", 1-", cens_ind, ") ~",
                                    paste(colnames(data %>% dplyr::ungroup() %>%
                                                     dplyr::select(all_of(Y),
                                                                   all_of(weights_cov))),
                                          collapse = "+")))
  model_est_c_yz = survival::survreg(aft_formula_c,
                                     data = data,
                                     dist = "lognormal")
  model_est_c_yz_coeff = model_est_c_yz$coefficients
  model_est_c_yz_sd = model_est_c_yz$scale
  YZ = data %>% dplyr::ungroup() %>%
    dplyr::select(all_of(Y), all_of(weights_cov)) %>% as.matrix()
  YZ = cbind(rep(1, nrow(YZ)), YZ)
  weights_d1 = pnorm(log(data[cens_name] %>% unlist()), mean = YZ %*% model_est_c_yz_coeff,
                     sd = model_est_c_yz_sd,
                     lower.tail = FALSE)

  ### For DELTA = 0
  aft_formula_x <- as.formula(paste("survival::Surv(", cens_name, ", ", cens_ind, ") ~",
                                    paste(colnames(data %>% dplyr::select(all_of(Y), all_of(weights_cov))),
                                          collapse = "+")))
  model_est_x_yz = survival::survreg(aft_formula_x,
                                     data = data,
                                     dist = "lognormal")
  model_est_x_yz_coeff = model_est_x_yz$coefficients
  model_est_x_yz_sd = model_est_x_yz$scale
  weights_d0 = pnorm(log(data[cens_name] %>% unlist()), mean = YZ %*% model_est_x_yz_coeff,
                     sd = model_est_x_yz_sd,
                     lower.tail = TRUE)


  return(data$D*weights_d1 + (1-data$D)*weights_d0)
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


# first derivative function
# calculates first derivative for subject i (data x)
# f(x;beta) where beta is of length lb, x is a scalar
# first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
#firstderivative_mvn <- function(beta, l, data, cens_name, cov_vars, cens_ind){
#   lb <- length(beta)
#   derivs <- matrix(data = 0, nrow = lb, ncol = lb)
#   delta <- beta * (10 ^ (- 4))
#   betal <- betar <- beta
#   for (i in 1:lb) {
#     # Perturb the ith element of beta
#     betal[i] <- beta[i] - delta[i]
#     betar[i] <- beta[i] + delta[i]
#
#     # Calculate function values
#     yout1 <- l(betal, data[cens_name], data[cov_vars], data[cens_ind])
#     yout2 <- l(betar, data[cens_name], data[cov_vars], data[cens_ind])
#
#     # Calculate derivative and save in vector A
#     derivs[i,] <- (yout2 - yout1) / (2 * delta[i])
#
#     # Reset parameter vectors
#     betal <- betar <- beta
#   }
#   #print("hi")
#   return(derivs)
# }


# sandwich_ipw_estimate <- function(formula,
#                                  data,
#                                  cens_name,
#                                  weights_cov,
#                                  cens_ind,
#                                 par_vec,
#                                  beta_est,
#                                  weights,
#                                  weight_dist_params){
#
#   #convert beta_est to numeric
#   beta_est = beta_est %>% as.numeric()
#
#   # add weights to data frame
#   data$weights = weights
#
#   # extract variable names from formula
#   varNames = all.vars(formula)
#   form2 = formula
#   form2[[2]] <- NULL
#   varNamesRHS = all.vars(form2)
#   Y = varNames[is.na(match(varNames, varNamesRHS))]
#
#   do.call("<-", list(par_vec, beta_est))
#
#   varNamesRHS = varNamesRHS[is.na(match(varNamesRHS, par_vec))]
#
#   # turn formula into function
#   cmd <- tail(as.character(formula),1)
#   exp <- parse(text=cmd)
#   exp_nobracket <- exp %>% as.character() %>% str_remove_all(., "\\[|\\]") %>%
#     parse(text = .)
#   m_func = function(p){
#     with(as.list(p),
#          eval(exp_nobracket))
#   }
#
#
#   # create "g" function for the sandwich estimator
#   g = function(data, beta_est, m_func, par_vec, var_namesRHS, cens_ind){
#     p = c(beta_est, data[varNamesRHS])
#     names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)
#
#     rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
#       numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
#       rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()
#   }
#
#   # jacobian g function
#   g_jacobian = function(data, beta_est, m_func, par_vec, var_namesRHS, cens_ind){
#     p = c(beta_est, data[varNamesRHS])
#     names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)
#
#     j = numDeriv::jacobian(m_func, p)[1:length(beta_est)]
#
#     rep(data[cens_ind]*data["weights"], length(beta_est)) %>% as.numeric()*
#       ((rootSolve::hessian(m_func, p)[1:length(beta_est), 1:length(beta_est)]*
#           rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()) -
#          outer(j,j)) %>% as.numeric()
#   }
#
#
#   # take the inverse first derivative of g
#   first_der <- apply(data, 1, function(temp){
#     g_jacobian(temp, beta_est, m_func, par_vec, var_namesRHS, cens_ind)
#   })
#   if(length(beta_est) > 1){
#     first_der = first_der %>% rowMeans() %>% matrix(nrow = length(beta_est))
#   }else{
#     first_der = first_der %>% mean()
#   }
#   inv_first_der <- solve(first_der)
#
#   # need to get the outer product of g at each observation and take the mean
#   gs = apply(data, 1, function(temp)
#     g(temp, beta_est, m_func, par_vec, var_namesRHS, cens_ind))
#   if(length(beta_est) > 1){
#     outer_prod = apply(gs, 2, function(g) g%*%t(g))
#     outer_prod = outer_prod %>% rowMeans() %>% matrix(nrow = length(beta_est))
#   }else{
#     outer_prod = gs^2
#     outer_prod = outer_prod %>% mean()
#   }
#
#
#   ## then need to put it all together
#   se = sqrt((inv_first_der %*% outer_prod %*% t(inv_first_der) / nrow(data)) %>% diag())
#   return(se)
#
# }


#####################################
##### helper functions for AIPW #####
#####################################

#########################
##### Psi Functions #####
#########################

psi_hat_i_mvn <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_temp, m_func, mu_joint, Sigma_joint, sigma2){
  denominator =  integrate(integral_func_denom_mvn, 0, Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, cov_vars = cov_vars,
                           beta_temp = beta_temp, m_func = m_func,
                           mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-4){
    psi = rep(0, length(beta_temp))
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi_mvn, 0, Inf, data_row = data, Y = Y,
                varNamesRHS = varNamesRHS, par_vec = par_vec,
                cens_name = cens_name, cov_vars = cov_vars,
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

psi_hat_i_aft <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                          beta_temp, m_func, model_est_x_z_coeff,
                          model_est_x_z_sd, sigma2){
  denominator =  integrate(integral_func_denom_aft, 0, Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, cov_vars = cov_vars,
                           beta_temp = beta_temp, m_func = m_func,
                           model_est_x_z_coeff = model_est_x_z_coeff,
                           model_est_x_z_sd = model_est_x_z_sd,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-4){
    psi = rep(0, length(beta_temp))
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi_aft, 0, Inf, data_row = data, Y = Y,
                varNamesRHS = varNamesRHS, par_vec = par_vec,
                cens_name = cens_name, cov_vars = cov_vars,
                beta_temp = beta_temp, m_func = m_func,
                model_est_x_z_coeff = model_est_x_z_coeff,
                model_est_x_z_sd = model_est_x_z_sd,
                sigma2 = sigma2, j = j,
                rel.tol = .Machine$double.eps^0.1,
                subdivisions = 1000)$value
    }) %>% unlist() %>% as.numeric()
    psi = numerator/denominator
  }
  #print(beta_temp)
  psi
}

##########################################
##### Integral Denominator Functions #####
##########################################

integral_func_denom_mvn <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
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
                                 X.given = c(data_row[cov_vars] %>% as.numeric()))
    # f_c_xz = lapply(t, function(dummy_var)
    #   condMVNorm::pcmvnorm(lower = log(data_row[cens_name] %>% as.numeric()), upper = Inf,
    #                        mean = mu_joint, sigma = Sigma_joint,
    #                        dependent.ind = 2,
    #                        given = c(1,3),
    #                        X.given = c(log(data_row[cens_name] %>% as.numeric()),
    #                                    data_row[cov_vars] %>% as.numeric()))) %>% unlist()
    value_ts[i] = f_y*f_x_z/t[i]
  }
  value_ts
}

integral_func_denom_aft <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                beta_temp, m_func,
                                model_est_x_z_coeff, model_est_x_z_sd,
                                sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_z = dnorm(log(data_row[cens_name] %>% as.numeric()),
          mean = c(1, data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
          sd = model_est_x_z_sd)
    #print(f_x_z)
    value_ts[i] = f_y*f_x_z/t[i]
  }
  value_ts
}

########################################
##### Integral Numerator Functions #####
########################################

integral_func_psi_mvn <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
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
                                 X.given = c(data_row[cov_vars] %>% as.numeric()))
    # f_c_xz = lapply(t, function(dummy_var)
    #   condMVNorm::pcmvnorm(lower = log(data_row[cens_name] %>% as.numeric()), upper = Inf,
    #                        mean = mu_joint, sigma = Sigma_joint,
    #                        dependent.ind = 2,
    #                        given = c(1,3),
    #                        X.given = c(log(data_row[cens_name] %>% as.numeric()),
    #                                    data_row[cov_vars] %>% as.numeric()))) %>% unlist()
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_x_z/t[i]
  }
  value_ts %>% unlist()
}

integral_func_psi_aft <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                  beta_temp, m_func,
                                  model_est_x_z_coeff, model_est_x_z_sd, j, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_z = dnorm(log(data_row[cens_name] %>% as.numeric()),
                  mean = c(1, data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
                  sd = model_est_x_z_sd)
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_x_z/t[i]
  }
  value_ts %>% unlist()
}

###############################
##### Multiroot Functions #####
###############################

# set up multiroot function (the estimating equation we want to find the root of)

##### MVN
multiroot_func_mvn = function(beta_temp, data,
                          Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                          m_func, mu_joint, Sigma_joint, sigma2){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    ipw_piece = rep(as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    # psis = paste0("psi", seq(1:length(beta_temp)))
    aipw_piece = rep(1 - as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      psi_hat_i_mvn(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                beta_temp, m_func, mu_joint, Sigma_joint, sigma2)
    ipw_piece + aipw_piece
  }) %>% unname()
  rowSums(pieces)
}

##### AFT
multiroot_func_aft = function(beta_temp, data,
                          Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                          m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    ipw_piece = rep(as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    aipw_piece = rep(1 - as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      psi_hat_i_aft(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                beta_temp, m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2)
    ipw_piece + aipw_piece
  }) %>% unname()
  rowSums(pieces)
}

#####################################
##### helper functions for MLE ######
#####################################

#########################
##### Psi Functions #####
#########################

psi_hat_i_mle_mvn <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                      beta_temp, m_func, mu_joint, Sigma_joint, sigma2){
  denominator =  integrate(integral_func_denom_mle_mvn, data[cens_name] %>% as.numeric(),
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
    psi = rep(0, length(beta_temp))
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi_mle_mvn, data[cens_name] %>% as.numeric(),
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

psi_hat_i_mle_aft <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              beta_temp, m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2){
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
    psi = rep(0, length(beta_temp))
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

##########################################
##### Integral Denominator Functions #####
##########################################

integral_func_denom_mle_mvn <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
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

integral_func_denom_mle_aft <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                    beta_temp, m_func, C_val,
                                    model_est_x_z_coeff, model_est_x_z_sd, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_cz = dnorm(log(data_row[cens_name] %>% as.numeric()),
                   mean = (data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
                   sd = model_est_x_z_sd)
    value_ts[i] = f_y*f_x_cz/t[i]
  }
  value_ts
}

########################################
##### Integral Numerator Functions #####
########################################

integral_func_psi_mle_mvn <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
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
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_x_cz/t[i]
  }
  value_ts %>% unlist()
}

integral_func_psi_mle_aft <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                      beta_temp, m_func, C_val,
                                      model_est_x_z_coeff, model_est_x_z_sd, j, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_cz = dnorm(log(data_row[cens_name] %>% as.numeric()),
                   mean = (data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
                   sd = model_est_x_z_sd)
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_x_cz/t[i]
  }
  value_ts %>% unlist()
}

###############################
##### Multiroot Functions #####
###############################

# set up multiroot function (the estimating equation we want to find the root of)

##### MVN
multiroot_func_mle_mvn = function(beta_temp, data,
                          Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                          m_func, mu_joint, Sigma_joint, sigma2){
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    if(temp[cens_ind] == 1){
      piece = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
        rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    }else{
      piece = psi_hat_i_mle_mvn(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                            beta_temp, m_func, mu_joint, Sigma_joint, sigma2)
    }
    # print(piece)
    piece
  }) %>% unname()
  rowSums(pieces)
}

##### AFT
multiroot_func_mle_aft = function(beta_temp, data,
                              Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                              m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2){
    pieces = apply(data, 1, function(temp){
      p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
      if(temp[cens_ind] == 1){
        piece = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
          rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
      }else{
        piece = psi_hat_i_mle_aft(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              beta_temp, model_est_x_z_coeff, model_est_x_z_sd, sigma2)
      }
      # print(piece)
      piece
    }) %>% unname()
    rowSums(pieces)
  }




#####################################
##### helper functions for ACC  #####
#####################################

#########################
##### Psi Functions #####
#########################

psi_hat_i_mvn_acc <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                          beta_temp, m_func, mu_joint, Sigma_joint, sigma2){
  denominator =  integrate(integral_func_denom_mvn_acc, 0, Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, cov_vars = cov_vars,
                           beta_temp = beta_temp, m_func = m_func,
                           mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-4){
    psi = rep(0, length(beta_temp))
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi_mvn_acc, 0, Inf, data_row = data, Y = Y,
                varNamesRHS = varNamesRHS, par_vec = par_vec,
                cens_name = cens_name, cov_vars = cov_vars,
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

psi_hat_i_aft_acc <- function(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                          beta_temp, m_func, model_est_x_z_coeff,
                          model_est_x_z_sd, sigma2){
  denominator =  integrate(integral_func_denom_aft, 0, Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, cov_vars = cov_vars,
                           beta_temp = beta_temp, m_func = m_func,
                           model_est_x_z_coeff = model_est_x_z_coeff,
                           model_est_x_z_sd = model_est_x_z_sd,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-4){
    psi = rep(0, length(beta_temp))
  }else{
    numerator = lapply(1:length(beta_temp), function(j) {
      integrate(integral_func_psi_aft, 0, Inf, data_row = data, Y = Y,
                varNamesRHS = varNamesRHS, par_vec = par_vec,
                cens_name = cens_name, cov_vars = cov_vars,
                beta_temp = beta_temp, m_func = m_func,
                model_est_x_z_coeff = model_est_x_z_coeff,
                model_est_x_z_sd = model_est_x_z_sd,
                sigma2 = sigma2, j = j,
                rel.tol = .Machine$double.eps^0.1,
                subdivisions = 1000)$value
    }) %>% unlist() %>% as.numeric()
    psi = numerator/denominator
  }
  #print(beta_temp)
  psi
}

##########################################
##### Integral Denominator Functions #####
##########################################

integral_func_denom_mvn_acc <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
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
                                 X.given = c(data_row[cov_vars] %>% as.numeric()))
    f_c_xz = condMVNorm::pcmvnorm(lower = log(data_row[cens_name] %>% as.numeric()), upper = Inf,
                           mean = mu_joint, sigma = Sigma_joint,
                           dependent.ind = 2,
                           given = c(1,3),
                           X.given = c(log(data_row[cens_name] %>% as.numeric()),
                                       data_row[cov_vars] %>% as.numeric()))
    value_ts[i] = f_y*f_x_z*f_c_xz/t[i]
  }
  value_ts
}

integral_func_denom_aft_acc <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                    beta_temp, m_func,
                                    model_est_x_z_coeff, model_est_x_z_sd,
                                    sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_z = dnorm(log(data_row[cens_name] %>% as.numeric()),
                  mean = (data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
                  sd = model_est_x_z_sd)
    value_ts[i] = f_y*f_x_z/t[i]
  }
  value_ts
}

########################################
##### Integral Numerator Functions #####
########################################

integral_func_psi_mvn_acc <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
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
                                 X.given = c(data_row[cov_vars] %>% as.numeric()))
    f_c_xz = condMVNorm::pcmvnorm(lower = log(data_row[cens_name] %>% as.numeric()), upper = Inf,
                                  mean = mu_joint, sigma = Sigma_joint,
                                  dependent.ind = 2,
                                  given = c(1,3),
                                  X.given = c(log(data_row[cens_name] %>% as.numeric()),
                                              data_row[cov_vars] %>% as.numeric()))
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_x_z*f_c_xz/t[i]
  }
  value_ts %>% unlist()
}

integral_func_psi_aft_acc <- function(t, data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                  beta_temp, m_func,
                                  model_est_x_z_coeff, model_est_x_z_sd, j, sigma2){
  value_ts = vector("numeric", length(t))
  for(i in 1:length(t)){
    data_row[cens_name] = t[i]
    p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    m_t = m_func(p)
    f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
    f_x_z = dnorm(log(data_row[cens_name] %>% as.numeric()),
                  mean = (data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
                  sd = model_est_x_z_sd)
    value_ts[i] =
      numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_x_z/t[i]
  }
  value_ts %>% unlist()
}

###############################
##### Multiroot Functions #####
###############################

# set up multiroot function (the estimating equation we want to find the root of)

##### MVN
multiroot_func_mvn_acc = function(beta_temp, data,
                              Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                              m_func, mu_joint, Sigma_joint, sigma2){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    cc_piece = rep(temp[cens_ind], length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    acc_piece = rep(temp[cens_ind]-temp["weights"], length(beta_temp))*
      psi_hat_i_mvn_acc(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                    beta_temp, m_func, mu_joint, Sigma_joint, sigma2)
    cc_piece + acc_piece
  }) %>% unname()
  rowSums(pieces)
}

##### AFT
#multiroot_func_aft_acc = function(beta_temp, data,
#                              Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
#                              m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2){
#   print(beta_temp)
#   pieces = apply(data, 1, function(temp){
#     p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
#     names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
#     ipw_piece = rep(temp[cens_ind]*temp["weights"], length(beta_temp))*
#       numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
#       rep(temp[Y]-m_func(p), length(beta_temp))
#     aipw_piece = rep(1 - temp[cens_ind]*temp["weights"], length(beta_temp))*
#       psi_hat_i_aft(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
#                     beta_temp, m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2)
#     ipw_piece + aipw_piece
#   }) %>% unname()
#   rowSums(pieces)
# }

