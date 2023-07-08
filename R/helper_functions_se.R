##### helper functions for MVN est outer product #####

############################
##### For IPW and AIPW #####
############################

# estimating function for gamma
f_gamma = function(params, cens_name, weights_cov, cens_ind, data){
  lb <- length(params)
  derivs <- matrix(data = 0, nrow = lb, ncol = 1)
  delta <- params * (10 ^ (- 4))
  paramsl <- paramsr <- params
  for (i in 1:lb) {
    # Perturb the ith element of beta
    paramsl[i] <- params[i] - delta[i]
    paramsr[i] <- params[i] + delta[i]

    # Calculate function values
    yout1 <- l(paramsl, data[[cens_name]], data[[weights_cov]], data[[cens_ind]])
    yout2 <- l(paramsr, data[[cens_name]], data[[weights_cov]], data[[cens_ind]])

    # Calculate derivative and save in vector A
    derivs[i,] <- (yout2 - yout1) / (2 * delta[i])

    # Reset parameter vectors
    paramsl <- paramsr <- params
  }
  #print("hi")
  return(derivs)
}

firstderivative_f_gamma <- function(params, f_gamma, data, cens_ind, cens_name, weights_cov){
  lb <- length(params)
  derivs <- matrix(data = 0, nrow = lb, ncol = lb)
  delta <- params * (10 ^ (- 4))
  paramsl <- paramsr <- params
  for (i in 1:lb) {
    # Perturb the ith element of beta
    paramsl[i] <- params[i] - delta[i]
    paramsr[i] <- params[i] + delta[i]

    # Calculate function values
    yout1 <- f_gamma(paramsl, cens_name, weights_cov, cens_ind, data)
    yout2 <- f_gamma(paramsr, cens_name, weights_cov, cens_ind, data)

    # Calculate derivative and save in vector A
    derivs[i,] <- (yout2 - yout1) / (2 * delta[i])

    # Reset parameter vectors
    paramsl <- paramsr <- params
  }
  #print("hi")
  return(derivs)
}


###############
##### IPW #####
###############

# g as a function of gamma
g_gamma = function(params, data, beta_est, m_func, par_vec, varNamesRHS, cens_ind,
                   cens_name, weights_cov, Y){
  mu_joint = c(params[1], params[2], params[3])
  Sigma_joint = (matrix(c(params[4], params[7], params[8],
                          params[7], params[5], params[9],
                          params[8], params[9], params[6]),
                        nrow = 3))

  weights = pnorm(log(data[cens_name]%>% as.numeric()),
                  mean = mu_joint[2],
                  sd = sqrt(Sigma_joint[2,2]), lower.tail = FALSE)/
    condMVNorm::pcmvnorm(lower = log(data[cens_name]%>% as.numeric()), upper = Inf,
                         mean = mu_joint, sigma = Sigma_joint,
                         dependent.ind = 2,
                         given = c(1,3),
                         X.given = c(log(data[cens_name]%>% as.numeric()),
                                     data[weights_cov]%>% as.numeric()))

  p = c(beta_est, data[varNamesRHS])%>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

  rep(as.numeric(data[[cens_ind]])*as.numeric(weights),
      length(beta_est)) %>% as.numeric()*
    numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
    rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
}

firstderivative_g_gamma <- function(params, beta, g_gamma, data, varNamesRHS, par_vec,
                                    m_func, cens_ind, cens_name, weights_cov, Y){
  lb <- length(params)
  derivs <- matrix(data = 0, nrow = length(beta), ncol = lb)
  delta <- params * (10 ^ (- 4))
  paramsl <- paramsr <- params
  for (i in 1:lb) {
    # Perturb the ith element of beta
    paramsl[i] <- params[i] - delta[i]
    paramsr[i] <- params[i] + delta[i]

    # Calculate function values
    yout1 <- g_gamma(paramsl, data, beta, m_func, par_vec, varNamesRHS,
                     cens_ind, cens_name, weights_cov, Y)
    yout2 <- g_gamma(paramsr, data, beta, m_func, par_vec, varNamesRHS,
                     cens_ind, cens_name, weights_cov, Y)

    # Calculate derivative and save in vector A
    derivs[,i] <- (yout2 - yout1) / (2 * delta[i])

    # Reset parameter vectors
    paramsl <- paramsr <- params
  }
  #print("hi")
  return(derivs)
}


################
##### AIPW #####
################

# g as a function of gamma
g_gamma_aipw = function(params, data, beta_est, m_func, par_vec, varNamesRHS, cens_ind,
                        cens_name, weights_cov, Y, cov_vars, cov_dist_params, sigma2, gh_nodes){
  mu_joint = c(params[1], params[2], params[3])
  Sigma_joint = (matrix(c(params[4], params[7], params[8],
                          params[7], params[5], params[9],
                          params[8], params[9], params[6]),
                        nrow = 3))

  weights = pnorm(log(data[cens_name]%>% as.numeric()),
                  mean = mu_joint[2],
                  sd = sqrt(Sigma_joint[2,2]), lower.tail = FALSE)/
    condMVNorm::pcmvnorm(lower = log(data[cens_name]%>% as.numeric()), upper = Inf,
                         mean = mu_joint, sigma = Sigma_joint,
                         dependent.ind = 2,
                         given = c(1,3),
                         X.given = c(log(data[cens_name]%>% as.numeric()),
                                     data[weights_cov]%>% as.numeric()))

  p = c(beta_est, data[varNamesRHS])%>% as.numeric()
  names(p) = c(paste0(par_vec, seq(1:length(beta_est))), varNamesRHS)

  ##### CHANGE THIS TO AIPW G
  ipw_piece = rep(as.numeric(data[[cens_ind]])*as.numeric(weights),
                  length(beta_est)) %>% as.numeric()*
    numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
    rep(data[Y]  %>% as.numeric()-m_func(p), length(beta_est)) %>% as.numeric()
  aipw_piece = rep(1 - as.numeric(data[[cens_ind]])*as.numeric(weights), length(beta_est)) %>% as.numeric()*
    psi_hat_i_hermite_aipw(data, Y, varNamesRHS, par_vec, cens_name, cov_vars,
                           beta_est, m_func, cov_dist_params, sigma2, gh_nodes)

  ipw_piece + aipw_piece
}


firstderivative_g_gamma_aipw <- function(params, beta, g_gamma, data, varNamesRHS, par_vec,
                                         m_func, cens_ind, cens_name, weights_cov, Y,
                                         cov_vars, cov_dist_params, sigma2, gh_nodes){
  lb <- length(params)
  derivs <- matrix(data = 0, nrow = length(beta), ncol = lb)
  delta <- params * (10 ^ (- 4))
  paramsl <- paramsr <- params
  for (i in 1:lb) {
    # Perturb the ith element of beta
    paramsl[i] <- params[i] - delta[i]
    paramsr[i] <- params[i] + delta[i]

    # Calculate function values
    yout1 <- g_gamma(paramsl, data, beta, m_func, par_vec, varNamesRHS,
                     cens_ind, cens_name, weights_cov, Y,
                     cov_vars, cov_dist_params, sigma2, gh_nodes)
    yout2 <- g_gamma(paramsr, data, beta, m_func, par_vec, varNamesRHS,
                     cens_ind, cens_name, weights_cov, Y,
                     cov_vars, cov_dist_params, sigma2, gh_nodes)

    # Calculate derivative and save in vector A
    derivs[,i] <- (yout2 - yout1) / (2 * delta[i])

    # Reset parameter vectors
    paramsl <- paramsr <- params
  }
  #print("hi")
  return(derivs)
}







