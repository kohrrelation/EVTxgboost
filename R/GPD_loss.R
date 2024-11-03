
# log_xi <- log_xi_out <- log(0.1)
# log_sigma <- log_sigma_out <- log(1)

d_prob <- function(x, preds, log_xi, exponential=FALSE, orthogonal=FALSE) {
  if (exponential==TRUE){
    -( log(1/exp(preds)) - x/exp(preds) )
  } else {
    if (orthogonal==FALSE){
      -log( 1/exp(preds) * (1+exp(log_xi)*x/exp(preds))^(-(exp(log_xi)+1)/exp(log_xi)) )
    } else {
      (1+1/exp(log_xi))*log(1+ exp(log_xi)*(exp(log_xi)+1)*x/exp(preds)) +((preds)) -log(exp(log_xi)+1)
    }
  }
}

# d_prob(0.6,log(1.5),log(0.1))
# fExtremes::dgpd(x=0.6,mu=0,beta=(1.5),xi=0.1, log=TRUE)
# d_prob(0.6,log((1.5)*((0.1)+1)),log(0.1), orthogonal=TRUE)


grad_gpd <- function(x, preds, log_xi, exponential=FALSE, orthogonal=FALSE) {
  if (exponential==TRUE){
    exp(preds)/exp(preds)^2/(1/exp(preds)) - x * exp(preds)/exp(preds)^2
  } else {
    if (orthogonal==FALSE){
      (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                              1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                        (exp(log_xi) * x * exp(preds)/exp(preds)^2))) + exp(preds)/exp(preds)^2 *
         (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)))/(1/exp(preds) *
                                                                               (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)))
    } else {
      1 - (1 + 1/exp(log_xi)) * (exp(log_xi) * (exp(log_xi) + 1) *
                                   x * exp(preds)/exp(preds)^2/(1 + exp(log_xi) * (exp(log_xi) +
                                                                                     1) * x/exp(preds)))
    }
  }
}

grad_gpd_xi <- function(x, preds, log_xi, exponential=FALSE, orthogonal=FALSE) {
  if (orthogonal==FALSE){
  -(1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                           1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                     (exp(log_xi) * x/exp(preds))) - (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) +
                                                                                                                                                           1)/exp(log_xi)) * (log((1 + exp(log_xi) * x/exp(preds))) *
                                                                                                                                                                                (exp(log_xi)/exp(log_xi) - (exp(log_xi) + 1) * exp(log_xi)/exp(log_xi)^2)))/(1/exp(preds) *
                                                                                                                                                                                                                                                               (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi))))
  } else {
    (1 + 1/exp(log_xi)) * ((exp(log_xi) * (exp(log_xi) + 1) + exp(log_xi) *
                              exp(log_xi)) * x/exp(preds)/(1 + exp(log_xi) * (exp(log_xi) +
                                                                                1) * x/exp(preds))) - exp(log_xi)/exp(log_xi)^2 * log(1 +
                                                                                                                                        exp(log_xi) * (exp(log_xi) + 1) * x/exp(preds)) - exp(log_xi)/(exp(log_xi) +
                                                                                                                                                                                                         1)
  }
}


hess_gpd <- function(x, preds, log_xi,exponential=FALSE, orthogonal=FALSE) {
  if (exponential==TRUE){
    (exp(preds)/exp(preds)^2 - exp(preds) * (2 * (exp(preds) * exp(preds)))/(exp(preds)^2)^2)/(1/exp(preds)) +
      exp(preds)/exp(preds)^2 * (exp(preds)/exp(preds)^2)/(1/exp(preds))^2 -
      (x * exp(preds)/exp(preds)^2 - x * exp(preds) * (2 * (exp(preds) *
                                                              exp(preds)))/(exp(preds)^2)^2)
  } else {
    if (orthogonal==FALSE){
      (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                              1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                        (exp(log_xi) * x * exp(preds)/exp(preds)^2 - exp(log_xi) *
                                                                                           x * exp(preds) * (2 * (exp(preds) * exp(preds)))/(exp(preds)^2)^2)) -
                         (1 + exp(log_xi) * x/exp(preds))^(((-(exp(log_xi) + 1)/exp(log_xi)) -
                                                              1) - 1) * (((-(exp(log_xi) + 1)/exp(log_xi)) - 1) * (exp(log_xi) *
                                                                                                                     x * exp(preds)/exp(preds)^2)) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                                                                                        (exp(log_xi) * x * exp(preds)/exp(preds)^2))) - exp(preds)/exp(preds)^2 *
         ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) + 1)/exp(log_xi)) -
                                              1) * ((-(exp(log_xi) + 1)/exp(log_xi)) * (exp(log_xi) *
                                                                                          x * exp(preds)/exp(preds)^2))) + ((exp(preds)/exp(preds)^2 -
                                                                                                                               exp(preds) * (2 * (exp(preds) * exp(preds)))/(exp(preds)^2)^2) *
                                                                                                                              (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)) -
                                                                                                                              exp(preds)/exp(preds)^2 * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                                                                                                                                                                1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                                                                                                                                                          (exp(log_xi) * x * exp(preds)/exp(preds)^2)))))/(1/exp(preds) *
                                                                                                                                                                                                                                                                             (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi))) +
        (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                                1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                          (exp(log_xi) * x * exp(preds)/exp(preds)^2))) + exp(preds)/exp(preds)^2 *
           (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi))) *
        (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                                1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                          (exp(log_xi) * x * exp(preds)/exp(preds)^2))) + exp(preds)/exp(preds)^2 *
           (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) +
                                                 1)/exp(log_xi)))/(1/exp(preds) * (1 + exp(log_xi) *
                                                                                     x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)))^2
    } else {
      -((1 + 1/exp(log_xi)) * ((exp(log_xi) * (exp(log_xi) + 1) * x *
                                  exp(preds)/exp(preds)^2 - exp(log_xi) * (exp(log_xi) + 1) *
                                  x * exp(preds) * (2 * (exp(preds) * exp(preds)))/(exp(preds)^2)^2)/(1 +
                                                                                                        exp(log_xi) * (exp(log_xi) + 1) * x/exp(preds)) + exp(log_xi) *
                                 (exp(log_xi) + 1) * x * exp(preds)/exp(preds)^2 * (exp(log_xi) *
                                                                                      (exp(log_xi) + 1) * x * exp(preds)/exp(preds)^2)/(1 + exp(log_xi) *
                                                                                                                                          (exp(log_xi) + 1) * x/exp(preds))^2))
    }
  }
}


hess_gpd_xi <- function(x, preds, log_xi, orthogonal=FALSE) {
  if (orthogonal==FALSE){
    (1/exp(preds) * (((1 + exp(log_xi) * x/exp(preds))^(((-(exp(log_xi) +
                                                            1)/exp(log_xi)) - 1) - 1) * (((-(exp(log_xi) + 1)/exp(log_xi)) -
                                                                                            1) * (exp(log_xi) * x/exp(preds))) - (1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                                                                                                                                        1)/exp(log_xi)) - 1) * (log((1 + exp(log_xi) * x/exp(preds))) *
                                                                                                                                                                                                  (exp(log_xi)/exp(log_xi) - (exp(log_xi) + 1) * exp(log_xi)/exp(log_xi)^2))) *
                     ((-(exp(log_xi) + 1)/exp(log_xi)) * (exp(log_xi) * x * exp(preds)/exp(preds)^2)) +
                     (1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) + 1)/exp(log_xi)) -
                                                         1) * ((-(exp(log_xi) + 1)/exp(log_xi)) * (exp(log_xi) *
                                                                                                     x * exp(preds)/exp(preds)^2) - (exp(log_xi)/exp(log_xi) -
                                                                                                                                       (exp(log_xi) + 1) * exp(log_xi)/exp(log_xi)^2) * (exp(log_xi) *
                                                                                                                                                                                           x * exp(preds)/exp(preds)^2))) + exp(preds)/exp(preds)^2 *
     ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) + 1)/exp(log_xi)) -
                                          1) * ((-(exp(log_xi) + 1)/exp(log_xi)) * (exp(log_xi) *
                                                                                      x/exp(preds))) - (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) +
                                                                                                                                             1)/exp(log_xi)) * (log((1 + exp(log_xi) * x/exp(preds))) *
                                                                                                                                                                  (exp(log_xi)/exp(log_xi) - (exp(log_xi) + 1) * exp(log_xi)/exp(log_xi)^2))))/(1/exp(preds) *
                                                                                                                                                                                                                                                  (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi))) -
    (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                            1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                      (exp(log_xi) * x * exp(preds)/exp(preds)^2))) + exp(preds)/exp(preds)^2 *
       (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi))) *
    (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                            1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                      (exp(log_xi) * x/exp(preds))) - (1 + exp(log_xi) *
                                                                                                                         x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)) *
                       (log((1 + exp(log_xi) * x/exp(preds))) * (exp(log_xi)/exp(log_xi) -
                                                                   (exp(log_xi) + 1) * exp(log_xi)/exp(log_xi)^2))))/(1/exp(preds) *
                                                                                                                        (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)))^2
  } else {
    (1 + 1/exp(log_xi)) * ((exp(log_xi) * (exp(log_xi) + 1) + exp(log_xi) *
                              exp(log_xi) + (exp(log_xi) * exp(log_xi) + exp(log_xi) *
                                               exp(log_xi))) * x/exp(preds)/(1 + exp(log_xi) * (exp(log_xi) +
                                                                                                  1) * x/exp(preds)) - (exp(log_xi) * (exp(log_xi) + 1) + exp(log_xi) *
                                                                                                                          exp(log_xi)) * x/exp(preds) * ((exp(log_xi) * (exp(log_xi) +
                                                                                                                                                                           1) + exp(log_xi) * exp(log_xi)) * x/exp(preds))/(1 + exp(log_xi) *
                                                                                                                                                                                                                              (exp(log_xi) + 1) * x/exp(preds))^2) - exp(log_xi)/exp(log_xi)^2 *
      ((exp(log_xi) * (exp(log_xi) + 1) + exp(log_xi) * exp(log_xi)) *
         x/exp(preds)/(1 + exp(log_xi) * (exp(log_xi) + 1) * x/exp(preds))) -
      ((exp(log_xi)/exp(log_xi)^2 - exp(log_xi) * (2 * (exp(log_xi) *
                                                          exp(log_xi)))/(exp(log_xi)^2)^2) * log(1 + exp(log_xi) *
                                                                                                   (exp(log_xi) + 1) * x/exp(preds)) + exp(log_xi)/exp(log_xi)^2 *
         ((exp(log_xi) * (exp(log_xi) + 1) + exp(log_xi) * exp(log_xi)) *
            x/exp(preds)/(1 + exp(log_xi) * (exp(log_xi) + 1) *
                            x/exp(preds)))) - (exp(log_xi)/(exp(log_xi) + 1) -
                                                 exp(log_xi) * exp(log_xi)/(exp(log_xi) + 1)^2)
  }
}


myobjective_gpd <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi <- EVTxgboost_env$log_xi
  grad <- grad_gpd(labels, preds, log_xi, exponential=FALSE, orthogonal=FALSE)
  hess <- hess_gpd(labels, preds, log_xi, exponential=FALSE, orthogonal=FALSE)
  return(list(grad = grad, hess = hess))
}


myobjective_gpd_ortho <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi <- EVTxgboost_env$log_xi
  grad <- grad_gpd(labels, preds, log_xi, exponential=FALSE, orthogonal=TRUE)
  hess <- hess_gpd(labels, preds, log_xi, exponential=FALSE, orthogonal=TRUE)
  return(list(grad = grad, hess = hess))
}


myobjective_gpd_xi <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_sigma <- EVTxgboost_env$log_sigma
  grad <- grad_gpd_xi(labels, log_sigma, preds, orthogonal=FALSE)
  hess <- hess_gpd_xi(labels, log_sigma, preds, orthogonal=FALSE)
  return(list(grad = grad, hess = hess))
}


myobjective_gpd_xi_ortho <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_sigma <- EVTxgboost_env$log_sigma
  grad <- grad_gpd_xi(labels, log_sigma, preds, orthogonal=TRUE)
  hess <- hess_gpd_xi(labels, log_sigma, preds, orthogonal=TRUE)
  return(list(grad = grad, hess = hess))
}


myobjective_exponential <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi <- EVTxgboost_env$log_xi
  grad <- grad_gpd(labels, preds, log_xi, exponential=TRUE)
  hess <- hess_gpd(labels, preds, log_xi, exponential=TRUE)
  return(list(grad = grad, hess = hess))
}



evalerror_gpd <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi <- EVTxgboost_env$log_xi
  p <-  d_prob(labels, preds, log_xi=log_xi , exponential=FALSE, orthogonal=FALSE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}



evalerror_gpd_ortho <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi <- EVTxgboost_env$log_xi
  p <-  d_prob(labels, preds, log_xi=log_xi , exponential=FALSE, orthogonal=TRUE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}


evalerror_exponential <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  p <-  d_prob(labels, preds, exponential=TRUE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}

evalerror_gpd_xi <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_sigma <- EVTxgboost_env$log_sigma
  p <-  d_prob(labels, log_sigma, log_xi=preds, orthogonal=FALSE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}


evalerror_gpd_xi_ortho <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_sigma <- EVTxgboost_env$log_sigma
  p <-  d_prob(labels, log_sigma, log_xi=preds, orthogonal=TRUE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}

evalerror_gpd_cv <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi_out <- EVTxgboost_env$log_xi_out
  p <-  d_prob(labels, preds, log_xi=log_xi_out , exponential=FALSE, orthogonal=FALSE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}


evalerror_gpd_cv_ortho <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_xi_out <- EVTxgboost_env$log_xi_out
  p <-  d_prob(labels, preds, log_xi=log_xi_out , exponential=FALSE, orthogonal=TRUE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}


evalerror_gpd_xi_cv <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_sigma_out <- EVTxgboost_env$log_sigma_out
  p <-  d_prob(labels, log_sigma_out, log_xi=preds, orthogonal=FALSE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}


evalerror_gpd_xi_cv_ortho <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  log_sigma_out <- EVTxgboost_env$log_sigma_out
  p <-  d_prob(labels, log_sigma_out, log_xi=preds, orthogonal=TRUE)
  err <- mean( (p))
  return(list(metric = "GPD_loss", value = err))
}



EVTxgboost_env <- new.env(parent = emptyenv())

update_log_sigma <- function(xgb_model, dtrain_all) {
  # Calculate the new value for log_sigma
  log_sigma_value <- predict(xgb_model, dtrain_all)

  # Assign the calculated value to the environment
  EVTxgboost_env$log_sigma <- log_sigma_value
}


update_log_sigma_out <- function(xgb_model, dtrain_all) {
  # Calculate the new value for log_sigma
  log_sigma_value <- predict(xgb_model, dtrain_all)

  # Assign the calculated value to the environment
  EVTxgboost_env$log_sigma_out <- log_sigma_value
}


update_log_xi <- function(xgb_model, dtrain_all) {
  # Calculate the new value for log_sigma
  log_xi_value <- predict(xgb_model, dtrain_all)

  # Assign the calculated value to the environment
  EVTxgboost_env$log_xi <- log_xi_value
}


update_log_xi_out <- function(xgb_model, dtrain_all) {
  # Calculate the new value for log_sigma
  log_xi_value <- predict(xgb_model, dtrain_all)

  # Assign the calculated value to the environment
  EVTxgboost_env$log_xi_out <- log_xi_value
}

get_log_sigma <- function() {
  return(EVTxgboost_env$log_sigma)
}

get_log_xi <- function() {
  return(EVTxgboost_env$log_xi)
}


get_log_sigma_out <- function() {
  return(EVTxgboost_env$log_sigma_out)
}

get_log_xi_out <- function() {
  return(EVTxgboost_env$log_xi_out)
}



