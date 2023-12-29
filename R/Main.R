
GPDxgb.train <- function(y, X, xi){

  if (!inherits(X, "xgb.DMatrix"))
    stop("second argument X must be xgb.DMatrix")

  dtrain_all <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
  xgboost::setinfo(dtrain_all, 'label_lower_bound', max_data[exceed_indx])

  watchlist <- list(eval = dtrain_all)

  best_param <- list(booster = 'gbtree',
                     objective= myobjective_gpd ,
                     tree_method = 'hist',
                     learning_rate=0.2
                     , max_depth=5,
                     eval_metric= evalerror_gpd
  )


  model_intensity <- xgb_model <- xgb.train(params = best_param, data = dtrain_all, watchlist=watchlist,
                                            nrounds = best.nround, verbose=TRUE, maximize=FALSE)
}
