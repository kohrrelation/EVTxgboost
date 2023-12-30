


d_prob <- function(x, preds, log_xi, exponential=FALSE) {
  if (exponential==TRUE){
    -( log(1/exp(preds)) - x/exp(preds) )
  } else {
    -log( 1/exp(preds) * (1+exp(log_xi)*x/exp(preds))^(-(exp(log_xi)+1)/exp(log_xi)) )
  }
}

grad_gpd <- function(x, preds, log_xi, exponential=FALSE) {
  if (exponential==TRUE){
    exp(preds)/exp(preds)^2/(1/exp(preds)) - x * exp(preds)/exp(preds)^2
  } else {
    (1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                            1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                      (exp(log_xi) * x * exp(preds)/exp(preds)^2))) + exp(preds)/exp(preds)^2 *
       (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)))/(1/exp(preds) *
                                                                             (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi)))
  }
}

grad_gpd_xi <- function(x, preds, log_xi) {
  -(1/exp(preds) * ((1 + exp(log_xi) * x/exp(preds))^((-(exp(log_xi) +
                                                           1)/exp(log_xi)) - 1) * ((-(exp(log_xi) + 1)/exp(log_xi)) *
                                                                                     (exp(log_xi) * x/exp(preds))) - (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) +
                                                                                                                                                           1)/exp(log_xi)) * (log((1 + exp(log_xi) * x/exp(preds))) *
                                                                                                                                                                                (exp(log_xi)/exp(log_xi) - (exp(log_xi) + 1) * exp(log_xi)/exp(log_xi)^2)))/(1/exp(preds) *
                                                                                                                                                                                                                                                               (1 + exp(log_xi) * x/exp(preds))^(-(exp(log_xi) + 1)/exp(log_xi))))
}

hess_gpd <- function(x, preds, log_xi,exponential=FALSE) {
  if (exponential==TRUE){
    (exp(preds)/exp(preds)^2 - exp(preds) * (2 * (exp(preds) * exp(preds)))/(exp(preds)^2)^2)/(1/exp(preds)) +
      exp(preds)/exp(preds)^2 * (exp(preds)/exp(preds)^2)/(1/exp(preds))^2 -
      (x * exp(preds)/exp(preds)^2 - x * exp(preds) * (2 * (exp(preds) *
                                                              exp(preds)))/(exp(preds)^2)^2)
  } else {
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
  }
}



hess_gpd_xi <- function(x, preds, log_xi) {
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
}


myobjective_gpd <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  grad <- grad_gpd(labels, preds, log_xi, exponential)
  hess <- hess_gpd(labels, preds, log_xi, exponential)
  return(list(grad = grad, hess = hess))
}


myobjective_gpd_xi <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  grad <- grad_gpd_xi(labels, log_sigma, preds)
  hess <- hess_gpd_xi(labels, log_sigma, preds)
  return(list(grad = grad, hess = hess))
}



evalerror_gpd <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  p <-  d_prob(labels, preds, log_xi=log_xi )
  err <- mean( (p))
  return(list(metric = "MyError", value = err))
}

evalerror_gpd_xi <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  p <-  d_prob(labels, log_sigma, log_xi=preds)
  err <- mean( (p))
  return(list(metric = "MyError", value = err))
}