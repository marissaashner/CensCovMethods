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
mle_censored <- function(formula,
                          data,
                          cens_ind,
                          cens_name,
                          par_vec,
                          starting_vals,
                          sandwich_se = TRUE){

  # Need to add error checks

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

  # Note: the A function is as follows
  # set p --- p = c(beta_temp, data[varNamesRHS])
  # set names for p ---- names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
  # A = numDeriv::jacobian(m_func, p)[1:length(beta_temp)]

  # CC to get estimates for starting values (this can be an option) and sigma2
  model_est_cc = cc_censored(formula, cens_ind, data, par_vec, starting_vals,
                             sandwich_se = FALSE)
  starting_vals = model_est_cc$beta_est %>% as.numeric()
  sigma2 = model_est_cc$sigma_est

  # set up multiroot function (the estimating equation we want to find the root of)
  multiroot_func = function(beta_temp, data,
                            Y, varNamesRHS, par_vec, cens_name, cens_ind,
                            m_func, mu_joint, Sigma_joint, sigma2){
    pieces = apply(data, 1, function(temp){
      p = c(beta_temp, temp[varNamesRHS]) %>% as.numeric()
      names(p) = c(paste0(par_vec, seq(1:length(beta_temp))), varNamesRHS)
      if(temp[cens_ind] == 1){
        print("1")
        numDeriv::jacobian(m_func, p)[1:length(beta_temp)]*
          rep(temp[Y]-m_func(p), length(beta_temp))
      }else{
        print("0")
        psi_hat_i_mle(temp, Y, varNamesRHS, par_vec, cens_name,
                  beta_temp, m_func, mu_joint, Sigma_joint, sigma2)
      }
    }) %>% unname()
    rowSums(pieces)
  }

  beta_est = rootSolve::multiroot(multiroot_func,
                                  data = data,
                                  Y = Y, varNamesRHS = varNamesRHS, par_vec = par_vec,
                                  cens_name = cens_name, cens_ind = cens_ind,
                                  m_func = m_func, mu_joint = mu_joint,
                                  Sigma_joint = Sigma_joint, sigma2 = sigma2,
                                  start = starting_vals)$root
  names(beta_est) = paste0(par_vec, seq(1:length(beta_est)))

  # run sandwich estimator
  if(sandwich_se){
    se_est = mle_sandwich(formula, data, Y, varNamesRHS, par_vec, cens_name,
                           beta_est, m_func, cens_ind, mu_joint, Sigma_joint, sigma2)
  }else{
    se_est = NULL
  }

  # save beta estimates
  return(list(beta_est = beta_est,
              se_est = se_est))
}

mle_sandwich <- function(formula, data, Y, varNamesRHS, par_vec, cens_name,
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

    if(data[cens_ind] == 1){
      numDeriv::jacobian(m_func, p)[1:length(beta_est)]*
        rep(data[Y]-m_func(p), length(beta_est)) %>% as.numeric()
    }else{
      psi_hat_i_mle(data, Y, varNamesRHS, par_vec, cens_name,
                beta_est, m_func, mu_joint, Sigma_joint, sigma2)
    }
  }

  # first derivative function
  # calculates first derivative for subject i (data x)
  # f(x;beta) where beta is of length lb, x is a scalar
  # first derivative = ( f(beta+ delta) - f(beta-delta) )/ (2 * delta)
  firstderivative <- function(beta, g, data, Y, varNamesRHS, par_vec, cens_name,
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
      yout1 <- g(data, Y, varNamesRHS, par_vec, cens_name,
                 betal, m_func, cens_ind, mu_joint, Sigma_joint, sigma2)
      yout2 <- g(data, Y, varNamesRHS, par_vec, cens_name,
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
    firstderivative(beta_est, g, temp, Y, varNamesRHS, par_vec, cens_name,
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
    g(temp, Y, varNamesRHS, par_vec, cens_name,
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
