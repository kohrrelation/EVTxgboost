

grad_gpd <- function(x, log_sigma, xi) {
  (1/exp(log_sigma) * ((1 + xi * x/exp(log_sigma))^((-(xi + 1)/xi) -
                                                      1) * ((-(xi + 1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2))) +
     exp(log_sigma)/exp(log_sigma)^2 * (1 + xi * x/exp(log_sigma))^(-(xi +
                                                                        1)/xi))/(1/exp(log_sigma) * (1 + xi * x/exp(log_sigma))^(-(xi +
                                                                                                                                     1)/xi))
}




hess_gpd <- function(x, log_sigma, xi){
  (1/exp(log_sigma) * ((1 + xi * x/exp(log_sigma))^((-(xi + 1)/xi) -
                                                      1) * ((-(xi + 1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2 -
                                                                                xi * x * exp(log_sigma) * (2 * (exp(log_sigma) * exp(log_sigma)))/(exp(log_sigma)^2)^2)) -
                         (1 + xi * x/exp(log_sigma))^(((-(xi + 1)/xi) - 1) - 1) *
                         (((-(xi + 1)/xi) - 1) * (xi * x * exp(log_sigma)/exp(log_sigma)^2)) *
                         ((-(xi + 1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2))) -
     exp(log_sigma)/exp(log_sigma)^2 * ((1 + xi * x/exp(log_sigma))^((-(xi +
                                                                          1)/xi) - 1) * ((-(xi + 1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2))) +
     ((exp(log_sigma)/exp(log_sigma)^2 - exp(log_sigma) * (2 *
                                                             (exp(log_sigma) * exp(log_sigma)))/(exp(log_sigma)^2)^2) *
        (1 + xi * x/exp(log_sigma))^(-(xi + 1)/xi) - exp(log_sigma)/exp(log_sigma)^2 *
        ((1 + xi * x/exp(log_sigma))^((-(xi + 1)/xi) - 1) * ((-(xi +
                                                                  1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2)))))/(1/exp(log_sigma) *
                                                                                                                             (1 + xi * x/exp(log_sigma))^(-(xi + 1)/xi)) + (1/exp(log_sigma) *
                                                                                                                                                                              ((1 + xi * x/exp(log_sigma))^((-(xi + 1)/xi) - 1) * ((-(xi +
                                                                                                                                                                                                                                        1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2))) +
                                                                                                                                                                              exp(log_sigma)/exp(log_sigma)^2 * (1 + xi * x/exp(log_sigma))^(-(xi +
                                                                                                                                                                                                                                                 1)/xi)) * (1/exp(log_sigma) * ((1 + xi * x/exp(log_sigma))^((-(xi +
                                                                                                                                                                                                                                                                                                                  1)/xi) - 1) * ((-(xi + 1)/xi) * (xi * x * exp(log_sigma)/exp(log_sigma)^2))) +
                                                                                                                                                                                                                                                              exp(log_sigma)/exp(log_sigma)^2 * (1 + xi * x/exp(log_sigma))^(-(xi +
                                                                                                                                                                                                                                                                                                                                 1)/xi))/(1/exp(log_sigma) * (1 + xi * x/exp(log_sigma))^(-(xi +
                                                                                                                                                                                                                                                                                                                                                                                              1)/xi))^2

}




myobjective_gpd <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")

  grad <- grad_gpd(labels, preds, xi)

  hess <- hess_gpd(labels, preds, xi)

  return(list(grad = grad, hess = hess))
}



evalerror_gpd <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")

  p <-  d_prob(labels, log_sigma=(preds), xi=xi )

  err <- mean( (p))
  return(list(metric = "MyError", value = err))
}


