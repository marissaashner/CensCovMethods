
#############################
##### weighting options #####
#############################

#' Cox Weighting for IPW/AIPW
#'
#' Performs a cox proportional hazards regression model to find the inverse probability of being censored given fully observed covariates.
#'
#' @param data a data frame containing columns for the censoring indicator and the all covariates involved in the weights
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param weights_cov a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#'
#' @return A vector of weights
#'
#' @import tidyverse
#' @import survival
#'
#' @export
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



#' AFT Weighting for IPW/AIPW
#'
#' Performs a log-normal accelerated failure time regression model to find the inverse probability of being censored given fully observed covariates.
#'
#' @param data a data frame containing columns for the censoring indicator and the all covariates involved in the weights
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param weights_cov a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#'
#' @return A vector of weights
#'
#' @import tidyverse
#' @import survival
#'
#' @export
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

#' Multivariate Normal Weighting for IPW/AIPW
#'
#' Calculates the inverse probability of being censored given one fully observed covariate and the censoring value by maximizing the likelihood of a multivariate normal. The multivariate nornmal is the joint distribution of the joint multivariate normal distribution of the log of the censored covariate, the log of the censored value, and the fully observed covariate.
#'
#' @param data a data frame containing columns for the censoring indicator and the all covariates involved in the weights
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param weights_cov a list of character string indicating the names of the variables from \code{data} to be used as a predictor in the weights model
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#'
#' @return A list with the following elements:
#' \item{weights}{a vector of weights}
#' \item{mu_joint}{the mean of the joint multivariate normal distribution of the log of the censored covariate, the log of the censored value, and the fully observed covariate.}
#' \item{Sigma_joint}{the covariance matrix of the joint multivariate normal distribution of the log of the censored covariate, the log of the censored value, and the fully observed covariate.}
#' \item{params}{a vector of nine parameter estimates representing those estimated from the maximum likelihood technique.}
#'
#' @import tidyverse
#' @import condMVNorm
#' @import mvtnorm
#'
#' @export
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

  return_list = list(weights = weights, mu_joint = mu_joint,
                     Sigma_joint = Sigma_joint, params = params)

  return(return_list)
}


#' AFT Weighting for ACC
#'
#' Performs log-normal accelerated failure time regression models to find the inverse probability of being censored given fully observed covariates.
#'
#' @param data a data frame containing columns for the censoring indicator, the outcome, and the all covariates involved in the weights
#' @param cens_ind a character string indicating the name of censoring indicator from \code{data}, defined to be \code{=1} if observation is uncensored and \code{=0} if observation is censored
#' @param weights_cov a list of character strings indicating the names of the variables from \code{data} to be used as predictors in the weights model
#' @param Y a character string indicating the name of the outcome from \code{data}
#' @param cens_name a character string indicating the name of censored covariate from \code{data}
#'
#' @return A vector of weights
#'
#' @import tidyverse
#' @import survival
#'
#' @export
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
                                    paste(colnames(data %>% dplyr::ungroup() %>%
                                                     dplyr::select(all_of(Y), all_of(weights_cov))),
                                          collapse = "+")))
  model_est_x_yz = survival::survreg(aft_formula_x,
                                     data = data,
                                     dist = "lognormal")
  model_est_x_yz_coeff = model_est_x_yz$coefficients
  model_est_x_yz_sd = model_est_x_yz$scale
  weights_d0 = pnorm(log(data[cens_name] %>% unlist()), mean = YZ %*% model_est_x_yz_coeff,
                     sd = model_est_x_yz_sd,
                     lower.tail = TRUE)


  return(data[cens_ind]*weights_d1 + (1-data[cens_ind])*weights_d0)
}
