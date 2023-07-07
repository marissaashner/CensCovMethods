
############################
##### Helper Functions #####
############################

###########################
##### For MVN Weights #####
###########################

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
  if(abs(denominator) < 10e-16){
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
  if(abs(denominator) < 10e-16){
    psi = rep(0, length(beta_temp))
  }else{
    tryCatch(
      expr = {
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
      },
      error = function(e){
        print("ERROR")
      }
    )
    if(!exists("numerator")){
      numerator = lapply(1:length(beta_temp), function(j) {
        integrate(integral_func_psi_aft, 0, 1000, data_row = data, Y = Y,
                  varNamesRHS = varNamesRHS, par_vec = par_vec,
                  cens_name = cens_name, cov_vars = cov_vars,
                  beta_temp = beta_temp, m_func = m_func,
                  model_est_x_z_coeff = model_est_x_z_coeff,
                  model_est_x_z_sd = model_est_x_z_sd,
                  sigma2 = sigma2, j = j,
                  rel.tol = .Machine$double.eps^0.1,
                  subdivisions = 1000)$value
      }) %>% unlist() %>% as.numeric()
    }
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
  if(abs(denominator) < 10e-16){
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
  denominator =  integrate(integral_func_denom_mle_aft, data[cens_name] %>% as.numeric(),
                           Inf, data_row = data, Y = Y,
                           varNamesRHS = varNamesRHS, par_vec = par_vec,
                           cens_name = cens_name, cov_vars = cov_vars,
                           beta_temp = beta_temp, m_func = m_func,
                           C_val = data[cens_name]  %>% as.numeric(),
                           model_est_x_z_coeff = model_est_x_z_coeff,
                           model_est_x_z_sd =model_est_x_z_sd,
                           sigma2 = sigma2,
                           rel.tol = .Machine$double.eps^0.1,
                           subdivisions = 1000)$value
  if(abs(denominator) < 10e-16){
    psi = rep(0, length(beta_temp))
  }else{
    tryCatch(
      expr = {
        numerator = lapply(1:length(beta_temp), function(j) {
          integrate(integral_func_psi_mle_aft, data[cens_name] %>% as.numeric(),
                    Inf, data_row = data, Y = Y,
                    varNamesRHS = varNamesRHS, par_vec = par_vec,
                    cens_name = cens_name, cov_vars = cov_vars,
                    beta_temp = beta_temp, m_func = m_func,
                    C_val = data[cens_name]  %>% as.numeric(),
                    model_est_x_z_coeff = model_est_x_z_coeff,
                    model_est_x_z_sd =model_est_x_z_sd,
                    sigma2 = sigma2, j = j,
                    rel.tol = .Machine$double.eps^0.1,
                    subdivisions = 1000)$value
        }) %>% unlist() %>% as.numeric()
      },
      error = function(e){
        print("ERROR")
        numerator = lapply(1:length(beta_temp), function(j) {
          integrate(integral_func_psi_mle_aft, data[cens_name] %>% as.numeric(),
                    1000, data_row = data, Y = Y,
                    varNamesRHS = varNamesRHS, par_vec = par_vec,
                    cens_name = cens_name, cov_vars = cov_vars,
                    beta_temp = beta_temp, m_func = m_func,
                    C_val = data[cens_name]  %>% as.numeric(),
                    model_est_x_z_coeff = model_est_x_z_coeff,
                    model_est_x_z_sd =model_est_x_z_sd,
                    sigma2 = sigma2, j = j,
                    rel.tol = .Machine$double.eps^0.1,
                    subdivisions = 1000)$value
        }) %>% unlist() %>% as.numeric()
      }
    )
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
                   mean =  c(1, data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
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
                   mean =  c(1, data_row[cov_vars] %>% as.numeric()) %*% model_est_x_z_coeff,
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
  print(beta_temp)
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
  print(beta_temp)
    pieces = apply(data, 1, function(temp){
      p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
      if(temp[cens_ind] == 1){
        piece = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
          rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
      }else{
        piece = psi_hat_i_mle_aft(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                              beta_temp, m_func, model_est_x_z_coeff, model_est_x_z_sd, sigma2)
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
  if(abs(denominator) < 10e-16){
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
  if(abs(denominator) < 10e-16){
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


multiroot_func_mvn_linear_x_yz = function(beta_temp, data,
                                     Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                                     m_func, x_yz_dist_params){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    ipw_piece = rep(temp[cens_ind]*temp["weights"], length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y]-m_func(p), length(beta_temp))
    aipw_piece = rep(1 - temp[cens_ind]*temp["weights"], length(beta_temp))*
      psi_hat_i_mvn_linear(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                           beta_temp, m_func, x_yz_dist_params)
    ipw_piece + aipw_piece
  }) %>% unname()
  rowSums(pieces)
}





