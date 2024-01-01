

#' Title
#'
#' @param y response
#' @param X
#' @param xi
#' @param simultaneous
#' @param init.sigma
#' @param stepsize
#' @param tree_depth
#' @param tree_depth_xi
#' @param rounds
#' @param exponential
#'
#' @return
#' @export
#'
#' @examples
GPDxgb.train <- function(y, X, xi=0.1, simultaneous=FALSE, init.sigma, stepsize=0.1, tree_depth=5,
                         tree_depth_xi=3, rounds, exponential=FALSE){

  # if (!inherits(X, "matrix"))
  #   stop("second argument X must be xgb.DMatrix")

  if (xi<=0 & exponential==FALSE)
    stop("This function only runs for heavy-tailed responses with GPD shape xi>0.
         For the exponential-tailed case with xi=0, set exponential=TRUE")
  if (xi>0 & exponential==TRUE)
    warning("Eventhough xi>0 was supplied, The function will run with xi=0, as exponential=TRUE was supplied")
  if (xi==0 & exponential==FALSE){
    warning("Eventhough xi=0 was supplied, The function will run with xi>0 with initial xi equal to 0.1,
            as exponential=FALSE was supplied. To run with xi=0, set exponential=TRUE")
    xi=0.1
  }
  if (exponential==TRUE & simultaneous==TRUE){
    warning("The function will run with xi=0, simultaneous=TRUE supplied will be ignored")
    simultaneous=FALSE
  }
  if (init.sigma<0)
    stop("Sigma cannot be negative")
  if (!inherits(rounds, "numeric") | !exists('rounds'))
    stop("Please provide the number of boosting iterations (numeric)")


  dtrain_all <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
  xgboost::setinfo(dtrain_all, 'base_margin', rep(log(init.sigma), nrow(X)))

  #log_xi <- log(xi)

  assignInNamespace(x="log_xi", value=log(xi), ns="EVTxgboost")


  watchlist <- list(eval = dtrain_all)

  if (exponential==FALSE){
    param_list <- list(booster = 'gbtree',
                       objective= myobjective_gpd ,
                       tree_method = 'hist',
                       learning_rate=stepsize
                       , max_depth=tree_depth,
                       eval_metric= evalerror_gpd
    )
  } else {
    param_list <- list(booster = 'gbtree',
                       objective= myobjective_exponential ,
                       tree_method = 'hist',
                       learning_rate=stepsize
                       , max_depth=tree_depth,
                       eval_metric= evalerror_exponential
    )
  }

 if (simultaneous==FALSE){
   xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                            nrounds = rounds, verbose=TRUE, maximize=FALSE)
   xgb_model_xi <- NULL

 } else {

   xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                   nrounds = 1, verbose=TRUE, maximize=FALSE)

   assignInNamespace(x="log_sigma", value=predict(xgb_model, dtrain_all), ns="EVTxgboost")

   dtrain_all_xi <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
   xgboost::setinfo(dtrain_all_xi, 'base_margin', rep(log_xi, nrow(X)))
   watchlist_xi <- list(eval = dtrain_all_xi)

   param_list_xi <- list(booster = 'gbtree',
                      objective= myobjective_gpd_xi ,
                      tree_method = 'hist',
                      learning_rate=stepsize
                      , max_depth=tree_depth_xi,
                      eval_metric= evalerror_gpd_xi
   )

   xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                   nrounds = 1, verbose=TRUE, maximize=FALSE)

   assignInNamespace(x="log_xi", value=predict(xgb_model_xi, dtrain_all_xi), ns="EVTxgboost")


   if (rounds>1){
    for (iter in 2:rounds){
      xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                      nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model)
      assignInNamespace(x="log_sigma", value=predict(xgb_model, dtrain_all), ns="EVTxgboost")
      xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                         nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model_xi)
      assignInNamespace(x="log_xi", value=predict(xgb_model_xi, dtrain_all_xi), ns="EVTxgboost")
    }
   }
 }

  result_list <- list()

  result_list$xgb_model <- xgb_model
  result_list$xgb_model_xi <- xgb_model_xi


  return(result_list)
}






