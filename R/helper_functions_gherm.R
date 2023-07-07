
################
##### AIPW #####
################

hermite_numerator_aipw <- function(x, data_row, cens_name, beta_temp, par_vec,
                                   varNamesRHS, m_func, Y, sigma2, j){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y
}


hermite_denominator_aipw <- function(x, data_row, cens_name, beta_temp, par_vec,
                                     varNamesRHS, m_func, Y, sigma2){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  f_y
}

psi_hat_i_hermite_aipw = function(data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                             beta_temp, m_func, cov_dist_params, sigma2, gh_nodes){
  if(!is.null(cov_dist_params$model_est_x_z_coeff)){
    mu <- c(1, data_row[cov_vars] %>% as.numeric()) %*% cov_dist_params$model_est_x_z_coeff
    sigma <- cov_dist_params$model_est_x_z_sd
  }else if(!is.null(cov_dist_params$mu_joint)){
    mu <- cov_dist_params$mu_joint[1] +
      cov_dist_params$Sigma_joint[1,3]*cov_dist_params$Sigma_joint[3,3]^(-1)*
      (data_row[cov_vars] %>% as.numeric() - cov_dist_params$mu_joint[3])
    sigma <- sqrt(cov_dist_params$Sigma_joint[1,1] -
      cov_dist_params$Sigma_joint[1,3]^2*cov_dist_params$Sigma_joint[3,3]^(-1))
  }

  gherm <- statmod::gauss.quad.prob(gh_nodes, dist="normal", mu = mu %>% as.numeric(),
                                    sigma = sigma)

  numerator = lapply(1:length(beta_temp), function(j) {
    sum(gherm$weights * lapply(gherm$nodes, function(node){
      hermite_numerator_aipw(node, data_row, cens_name, beta_temp,
                        par_vec, varNamesRHS, m_func, Y, sigma2, j)
    }) %>% unlist())
  }) %>% unlist()

  denominator = sum(gherm$weights * lapply(gherm$nodes, function(node){
    hermite_denominator_aipw(node, data_row, cens_name, beta_temp,
                        par_vec, varNamesRHS, m_func, Y, sigma2)
  }) %>% unlist())

  numerator/denominator
}

multiroot_func_hermite_aipw = function(beta_temp, data,
                                  Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                                  m_func, cov_dist_params, sigma2, gh_nodes){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    ipw_piece = rep(as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    aipw_piece = rep(1 - as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      psi_hat_i_hermite_aipw(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                        beta_temp, m_func, cov_dist_params, sigma2, gh_nodes)
    ipw_piece + aipw_piece
  }) %>% unname()
  rowSums(pieces)
}


###############
##### ACC #####
###############

hermite_numerator_acc <- function(x, data_row, cens_name, beta_temp, par_vec,
                                   varNamesRHS, m_func, Y, sigma2, j, mu_c, sigma_c,
                                  cov_dist_params){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  if(is.function(mu_c)){
    mu_c = mu_c(x, cov_dist_params)
  }
  f_c_xz = pnorm(x, mean = mu_c, sd = sigma_c, lower.tail = FALSE)
  numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y*f_c_xz
}


hermite_denominator_acc <- function(x, data_row, cens_name, beta_temp, par_vec,
                                     varNamesRHS, m_func, Y, sigma2, mu_c, sigma_c,
                                    cov_dist_params){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  if(is.function(mu_c)){
    mu_c = mu_c(x, cov_dist_params)
  }
  f_c_xz = pnorm(x, mean = mu_c, sd = sigma_c, lower.tail = FALSE)
  f_y*f_c_xz
}

psi_hat_i_hermite_acc = function(data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                  beta_temp, m_func, cov_dist_params, sigma2, gh_nodes){
  if(!is.null(cov_dist_params$model_est_x_z_coeff)){
    mu <- c(1, data_row[cov_vars] %>% as.numeric()) %*% cov_dist_params$model_est_x_z_coeff
    sigma <- cov_dist_params$model_est_x_z_sd
    mu_c <- c(1, data_row[cov_vars] %>% as.numeric()) %*% cov_dist_params$model_est_c_z_coeff
    sigma_c <- cov_dist_params$model_est_c_z_sd
  }else if(!is.null(cov_dist_params$mu_joint)){
    mu <- cov_dist_params$mu_joint[1] +
      cov_dist_params$Sigma_joint[1,3]*cov_dist_params$Sigma_joint[3,3]^(-1)*
      (data_row[cov_vars] %>% as.numeric() - cov_dist_params$mu_joint[3])
    sigma <- sqrt(cov_dist_params$Sigma_joint[1,1] -
      cov_dist_params$Sigma_joint[1,3]^2*cov_dist_params$Sigma_joint[3,3]^(-1))
    mu_c <- function(x, cov_dist_params) {cov_dist_params$mu_joint[2] +
      cov_dist_params$Sigma_joint[2, c(1,3)]%*%solve(cov_dist_params$Sigma_joint[c(1, 3), c(1, 3)])%*%
      (c(x, data_row[cov_vars] %>% as.numeric())  - cov_dist_params$mu_joint[c(1,3)])}
    sigma_c <- sqrt(cov_dist_params$Sigma_joint[2,2] -
      cov_dist_params$Sigma_joint[2, c(1,3)]%*%solve(cov_dist_params$Sigma_joint[c(1, 3), c(1, 3)])%*%
      cov_dist_params$Sigma_joint[c(1,3), 2])
  }

  gherm <- statmod::gauss.quad.prob(gh_nodes, dist="normal", mu = mu %>% as.numeric(),
                                    sigma = sigma)

  numerator = lapply(1:length(beta_temp), function(j) {
    sum(gherm$weights * lapply(gherm$nodes, function(node){
      hermite_numerator_acc(node, data_row, cens_name, beta_temp,
                             par_vec, varNamesRHS, m_func, Y, sigma2, j, mu_c, sigma_c,
                            cov_dist_params)
    }) %>% unlist())
  }) %>% unlist()

  denominator = sum(gherm$weights * lapply(gherm$nodes, function(node){
    hermite_denominator_acc(node, data_row, cens_name, beta_temp,
                             par_vec, varNamesRHS, m_func, Y, sigma2, mu_c, sigma_c,
                            cov_dist_params)
  }) %>% unlist())

  numerator/denominator
}

multiroot_func_hermite_acc = function(beta_temp, data,
                                       Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                                       m_func, cov_dist_params, sigma2, gh_nodes){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    cc_piece = rep(as.numeric(temp[[cens_ind]]), length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    acc_piece = rep(as.numeric(temp[[cens_ind]]) - as.numeric(temp[["weights"]]), length(beta_temp))*
      psi_hat_i_hermite_acc(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                             beta_temp, m_func, cov_dist_params, sigma2, gh_nodes)
    cc_piece - acc_piece
  }) %>% unname()
  rowSums(pieces)
}

###############
##### MLE #####
###############

hermite_numerator_mle <- function(x, data_row, cens_name, beta_temp, par_vec,
                                   varNamesRHS, m_func, Y, sigma2, j){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y
}


hermite_denominator_mle <- function(x, data_row, cens_name, beta_temp, par_vec,
                                     varNamesRHS, m_func, Y, sigma2){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  f_y
}

psi_hat_i_hermite_mle = function(data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                  beta_temp, m_func, cov_dist_params, sigma2, gh_nodes){
  if(!is.null(cov_dist_params$model_est_x_z_coeff)){
    mu <- c(1, data_row[cov_vars] %>% as.numeric()) %*% cov_dist_params$model_est_x_z_coeff
    sigma <- cov_dist_params$model_est_x_z_sd
  }else if(!is.null(cov_dist_params$mu_joint)){
    mu <- cov_dist_params$mu_joint[1] +
      cov_dist_params$Sigma_joint[1, c(2,3)]%*%solve(cov_dist_params$Sigma_joint[c(2,3),c(2,3)])%*%
      (c(log(data_row[cens_name]) %>% as.numeric(),
         data_row[cov_vars] %>% as.numeric()) - cov_dist_params$mu_joint[c(2,3)])
    sigma <- sqrt(cov_dist_params$Sigma_joint[1,1] -
                    cov_dist_params$Sigma_joint[1,c(2,3)]%*%solve(cov_dist_params$Sigma_joint[c(2,3),c(2,3)])%*%
                    cov_dist_params$Sigma_joint[c(2,3), 1])
  }

  gherm <- statmod::gauss.quad.prob(gh_nodes, dist="normal", mu = mu %>% as.numeric(),
                                    sigma = sigma %>% as.numeric())

  numerator_gh = lapply(1:length(beta_temp), function(j) {
    sum(gherm$weights * lapply(gherm$nodes, function(node){
      hermite_numerator_mle(node, data_row, cens_name, beta_temp,
                             par_vec, varNamesRHS, m_func, Y, sigma2, j)
    }) %>% unlist())
  }) %>% unlist()

  denominator_gh = sum(gherm$weights * lapply(gherm$nodes, function(node){
    hermite_denominator_mle(node, data_row, cens_name, beta_temp,
                             par_vec, varNamesRHS, m_func, Y, sigma2)
  }) %>% unlist())

  x_sintegral = seq(0, data_row[cens_name] %>% as.numeric(), length = 100)

  numerator_s = lapply(1:length(beta_temp), function(j){
    sintegral(x = x_sintegral,
              fx = integral_func_psi_mle_mvn(x_sintegral, data_row = data_row, Y = Y,
                                               varNamesRHS = varNamesRHS, par_vec = par_vec,
                                               cens_name = cens_name, cov_vars = cov_vars,
                                               beta_temp = beta_temp, m_func = m_func,
                                               C_val = data_row[cens_name]  %>% as.numeric(),
                                               mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                                               j = j, sigma2 = sigma2))$int
  }) %>% unlist()


  denominator_s = sintegral(x = x_sintegral,
                            fx = integral_func_denom_mle_mvn(x_sintegral, data_row = data_row, Y = Y,
                                                             varNamesRHS = varNamesRHS, par_vec = par_vec,
                                                             cens_name = cens_name, cov_vars = cov_vars,
                                                             beta_temp = beta_temp, m_func = m_func,
                                                             C_val = data_row[cens_name]  %>% as.numeric(),
                                                             mu_joint = mu_joint, Sigma_joint = Sigma_joint,
                                                             sigma2 = sigma2))$int

  (numerator_gh-numerator_s)/(denominator_gh-denominator_s)
}

multiroot_func_hermite_mle = function(beta_temp, data,
                                       Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                                       m_func, cov_dist_params, sigma2, gh_nodes){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    if(temp[cens_ind] == 1){
      piece = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
        rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    }else{
      piece = psi_hat_i_hermite_mle(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                                beta_temp, m_func, cov_dist_params, sigma2, gh_nodes)
    }
    # print(piece)
    piece
  }) %>% unname()
  rowSums(pieces)
}





