hermite_numerator_aipw <- function(x, data_row, cens_name, beta_temp, par_vec,
                                   varNamesRHS, m_func, Y, sigma2, j){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  numDeriv::jacobian(m_func, p)[j]*(data_row[Y] %>% as.numeric()-m_t)*f_y/exp(x)
}


hermite_denominator_aipw <- function(x, data_row, cens_name, beta_temp, par_vec,
                                     varNamesRHS, m_func, Y, sigma2){
  data_row[cens_name] = exp(x)
  p = c(beta_temp, data_row[varNamesRHS]) %>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  m_t = m_func(p)
  f_y = dnorm(data_row[Y] %>% as.numeric(), mean = m_t, sd = sqrt(sigma2))
  f_y/exp(x)
}

psi_hat_i_hermite_aipw = function(data_row, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                             beta_temp, m_func, cov_dist_params, sigma2, gherm){
  if(!is.null(cov_dist_params$model_est_x_z_coeff)){
    mu <- c(1, data_row[cov_vars] %>% as.numeric()) %*% cov_dist_params$model_est_x_z_coeff
    sigma <- cov_dist_params$model_est_x_z_sd
  }else if(!is.null(cov_dist_params$mu_joint)){
    mu <- cov_dist_params$mu_joint[1] +
      cov_dist_params$Sigma_joint[1,3]*cov_dist_params$Sigma_joint[3,3]^(-1)*
      (data_row[cov_vars] %>% as.numeric() - cov_dist_params$mu_joint[3])
    sigma <-cov_dist_params$Sigma_joint[1,1] +
      cov_dist_params$Sigma_joint[1,3]^2*cov_dist_params$Sigma_joint[3,3]^(-1)
  }


  numerator = lapply(1:length(beta_temp), function(j) {
    sum(gherm$weights * sigma * lapply(gherm$nodes, function(node){
      hermite_numerator_aipw(sqrt(2)*sigma*node + mu, data_row, cens_name, beta_temp,
                        par_vec, varNamesRHS, m_func, Y, sigma2, j)
    }) %>% unlist())
  }) %>% unlist()

  denominator = sum(gherm$weights * sigma * lapply(gherm$nodes, function(node){
    hermite_denominator_aipw(sqrt(2)*sigma*node + mu, data_row, cens_name, beta_temp,
                        par_vec, varNamesRHS, m_func, Y, sigma2)
  }) %>% unlist())

  numerator/denominator
}

multiroot_func_hermite_aipw = function(beta_temp, data,
                                  Y, varNamesRHS, par_vec, cens_name, cov_vars, cens_ind,
                                  m_func, cov_dist_params, sigma2, gherm){
  print(beta_temp)
  pieces = apply(data, 1, function(temp){
    p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
    names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
    ipw_piece = rep(as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
      rep(temp[Y] %>% as.numeric()-m_func(p), length(beta_temp))
    aipw_piece = rep(1 - as.numeric(temp[[cens_ind]])*as.numeric(temp[["weights"]]), length(beta_temp))*
      psi_hat_i_hermite_aipw(temp, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                        beta_temp, m_func, cov_dist_params, sigma2, gherm)
    ipw_piece + aipw_piece
  }) %>% unname()
  rowSums(pieces)
}




