
n <- 10000
X <- cbind(rnorm(n,0,2), rnorm(n,0,2))
y <- fExtremes::rgpd(n, xi=0.1, beta = exp(1+X[,1]) )

hist(exp(-1+(X[,1])/15), breaks=1000)
y <- numeric()
for (i in 1:n){
  y[i] <- fExtremes::rgpd(1, xi=exp(-1+(X[i,1])/15), beta = exp(1+X[i,1]) )
}

init.sigma <- 1
stepsize=0.1
tree_depth=3
rounds=100
xi=0.1
exponential=FALSE

GPDxgb.train <- function(y, X, xi=0.1, simultaneous=FALSE, init.sigma, stepsize=0.1, tree_depth=5, rounds){

  # if (!inherits(X, "matrix"))
  #   stop("second argument X must be xgb.DMatrix")

  dtrain_all <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
  xgboost::setinfo(dtrain_all, 'base_margin', rep(init.sigma, nrow(X)))

  log_xi <- log(xi)

  watchlist <- list(eval = dtrain_all)

  param_list <- list(booster = 'gbtree',
                     objective= myobjective_gpd ,
                     tree_method = 'hist',
                     learning_rate=stepsize
                     , max_depth=tree_depth,
                     eval_metric= evalerror_gpd
  )


 if (simultaneous==FALSE){
  xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                            nrounds = rounds, verbose=TRUE, maximize=FALSE)
 } else {

   xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                   nrounds = 1, verbose=TRUE, maximize=FALSE)

   log_sigma <- predict(xgb_model, dtrain_all)

   dtrain_all_xi <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
   xgboost::setinfo(dtrain_all_xi, 'base_margin', rep(log_xi, nrow(X)))
   watchlist_xi <- list(eval = dtrain_all_xi)

   param_list_xi <- list(booster = 'gbtree',
                      objective= myobjective_gpd_xi ,
                      tree_method = 'hist',
                      learning_rate=stepsize
                      , max_depth=tree_depth,
                      eval_metric= evalerror_gpd_xi
   )

   xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                   nrounds = 1, verbose=TRUE, maximize=FALSE)
   log_xi <- predict(xgb_model_xi, dtrain_all_xi)

   if (rounds>1){
    for (iter in 2:rounds){
      xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                      nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model)
      log_sigma <- predict(xgb_model, dtrain_all)
      xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                         nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model_xi)
      log_xi <- predict(xgb_model_xi, dtrain_all_xi)
    }
   }
 }

}

boxplot(exp(predict(xgb_model_xi, dtrain_all_xi)))
abline(h=0.1)
plot((predict(xgb_model, dtrain_all)), x= (1+X[,1]), cex=0.2, pch=19)
abline(a=0,b=1)
plot((predict(xgb_model_xi, dtrain_all_xi)), x= (-1+(X[,1])/15), cex=0.2, pch=19)
abline(a=0,b=1)

plot((predict(xgb_model, dtrain_all)), x= X[,2], cex=0.2, pch=19)
plot(exp(predict(xgb_model_xi, dtrain_all_xi)), x= X[,2], cex=0.2, pch=19)

