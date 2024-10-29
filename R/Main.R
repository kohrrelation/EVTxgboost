

#' Fit a Generalized Pareto (GPD) boosting model with shape parameter xi>=0
#'
#' @param y a vector of observations of length \code{n}, or a matrix with \code{n} rows.
#' @param X a matrix or data.frame object. Design matrix of dimension \code{n * p}.
#' @param xi a non-negative numeric. The shape parameter of the GPD. Also indicates the starting value of xi in the boosting algorithm if \code{simultaneous=TRUE}. Note, if \code{xi=0}, the distribution is of type exponential
#' @param simultaneous a logical. Should simultaneous boosting of the the sigma and xi parameters be performed? If \code{FALSE}, the xi parameter will be fixed. By default, \code{FALSE}.
#' @param init.sigma a positive numeric. Indicates the starting value of sigma in the boosting algorithm.
#' @param stepsize a positive scalar \code{learning_rate} parameter for the boosting model for sigma.
#' See \link[xgboost]{xgb.train} for more details.
#' @param stepsize_xi a positive numeric. \code{learning_rate} parameter for the boosting model for xi.
#' See \link[xgboost]{xgb.train} for more details.
#' @param tree_depth a positive numeric. \code{tree_depth} parameter for the boosting model for sigma.
#' See \link[xgboost]{xgb.train} for more details.
#' @param tree_depth_xi a positive numeric. \code{tree_depth} parameter for the boosting model for xi.
#' See \link[xgboost]{xgb.train} for more details.
#' @param nrounds a positive numeric. \code{nrounds} parameter for xgboost calibration.
#' See \link[xgboost]{xgb.train} for more details.
#' @param X_test a matrix or data.frame object. Design matrix with the same row names as \code{X} to perform predictions on. By default set to \code{X}.
#' @param orthogonal a logical. Should the two parameters be orthogonalized. By default, \code{FALSE}.
#'
#' @return a list
#' @export
#'
#' @examples
GPDxgb.train <- function(y, X, xi, simultaneous=FALSE, init.sigma, stepsize=0.1, stepsize_xi=0.1, tree_depth=5,
                         tree_depth_xi=3, nrounds, X_test=NULL, orthogonal=FALSE){

  if (!inherits(X, "matrix") & !inherits(X, "data.frame"))
    stop("argument X must be a matrix or a data.frame")
  if (xi<0)
    stop("This function only runs for heavy-tailed and exponential-tailed responses with GPD shape xi>=0")
  if (xi>0){
    exponential=FALSE
  }
  if (xi==0){
    exponential=TRUE
  }
  if (xi==0 & simultaneous==TRUE){
    warning("The function will run with xi=0. The input simultaneous=TRUE supplied will be ignored, as simultaneous boosting will only run with positive xi.
            If you would like to boost the two parameters xi and sigma simultaneously, set simultaneous=TRUE and xi>0")
    simultaneous=FALSE
  }

  if (init.sigma<0)
    stop("Sigma cannot be negative")

  if (is.null(X_test)){
    X_test <- X
  }

  if (orthogonal==TRUE){
    init.sigma=(init.sigma*(xi+1))
  }


  dtrain_all <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
  xgboost::setinfo(dtrain_all, 'base_margin', rep(log(init.sigma), nrow(X)))

  #log_xi <- log(xi)

  assignInNamespace(x="log_xi", value=log(xi), ns="EVTxgboost")

  watchlist <- list(eval = dtrain_all)




  if (exponential==FALSE){
    if (orthogonal==FALSE){
      param_list <- list(booster = 'gbtree',
                         objective= myobjective_gpd ,
                         tree_method = 'hist',
                         learning_rate=stepsize
                         , max_depth=tree_depth,
                         eval_metric= evalerror_gpd)
    } else {
      param_list <- list(booster = 'gbtree',
                         objective= myobjective_gpd_ortho ,
                         tree_method = 'hist',
                         learning_rate=stepsize
                         , max_depth=tree_depth,
                         eval_metric= evalerror_gpd_ortho)
    }
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
                                            nrounds = nrounds, verbose=TRUE, maximize=FALSE)
   xgb_model_xi <- NULL

   dtrain_test <- xgboost::xgb.DMatrix(data=as.matrix(X_test))
   xgboost::setinfo(dtrain_test, 'base_margin', rep(log(init.sigma), nrow(X_test)))
   sigma.pred <- exp(predict(xgb_model, dtrain_test))
   xi.pred <- NULL

   plot(xgb_model$evaluation_log$eval_GPD_loss, x=xgb_model$evaluation_log$iter,
        xlab='Number of trees/boosting iterations', ylab='GPD loss likelihood (training)', type='b',
        pch=19)


 } else {

   xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                   nrounds = 1, verbose=TRUE, maximize=FALSE)

   assignInNamespace(x="log_sigma", value=predict(xgb_model, dtrain_all), ns="EVTxgboost")

   dtrain_all_xi <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)
   xgboost::setinfo(dtrain_all_xi, 'base_margin', rep(log_xi, nrow(X)))
   watchlist_xi <- list(eval = dtrain_all_xi)

   if (orthogonal==FALSE){
     param_list_xi <- list(booster = 'gbtree',
                        objective= myobjective_gpd_xi ,
                        tree_method = 'hist',
                        learning_rate=stepsize_xi
                        , max_depth=tree_depth_xi,
                        eval_metric= evalerror_gpd_xi)
   } else {
     param_list_xi <- list(booster = 'gbtree',
                           objective= myobjective_gpd_xi_ortho ,
                           tree_method = 'hist',
                           learning_rate=stepsize_xi
                           , max_depth=tree_depth_xi,
                           eval_metric= evalerror_gpd_xi_ortho)
   }

   xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                   nrounds = 1, verbose=TRUE, maximize=FALSE)

   assignInNamespace(x="log_xi", value=predict(xgb_model_xi, dtrain_all_xi), ns="EVTxgboost")


   if (nrounds>1){
    for (iter in 2:nrounds){
      xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                      nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model)
      assignInNamespace(x="log_sigma", value=predict(xgb_model, dtrain_all), ns="EVTxgboost")
      xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                         nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model_xi)
      assignInNamespace(x="log_xi", value=predict(xgb_model_xi, dtrain_all_xi), ns="EVTxgboost")
    }
   }

   dtrain_test <- xgboost::xgb.DMatrix(data=as.matrix(X_test))
   xgboost::setinfo(dtrain_test, 'base_margin', rep(log(init.sigma), nrow(X_test)))
   sigma.pred <- exp(predict(xgb_model, dtrain_test))
   dtrain_test <- xgboost::xgb.DMatrix(data=as.matrix(X_test))
   xgboost::setinfo(dtrain_test, 'base_margin', rep(log(xi), nrow(X_test)))
   xi.pred <- exp(predict(xgb_model_xi, dtrain_test))

   plot(matrix(rbind(xgb_model$evaluation_log$eval_GPD_loss ,
                     xgb_model_xi$evaluation_log$eval_GPD_loss)),
        x=1:(2*length(xgb_model$evaluation_log$eval_GPD_loss)) ,
        xlab='Number of trees/boosting iterations', ylab='GPD loss likelihood (training)', type='b',
        pch=19)


 }

  if (orthogonal==TRUE){
    sigma.pred <- sigma.pred/(xi.pred+1)
  }

  result_list <- list()

  result_list$xgb_model <- xgb_model
  result_list$xgb_model_xi <- xgb_model_xi
  result_list$sigma.pred <- sigma.pred
  result_list$xi.pred <- xi.pred
  result_list$xi <- xi
  result_list$init.sigma <- init.sigma

  return(result_list)
}









#' Cross-validation for a Generalized Pareto (GPD) boosting model with shape parameter xi>=0. Grid search for xi.
#'
#' @param y a vector of observations of length \code{n}, or a matrix with \code{n} rows.
#' @param X a matrix or data.frame object. Design matrix of dimension \code{n * p}.
#' @param xi a non-negative numeric. The shape parameter of the GPD. Also indicates the starting value of xi in the boosting algorithm if \code{simultaneous=TRUE}. Note, if \code{xi=0}, the distribution is of type exponential
#' @param simultaneous a logical. Should simultaneous boosting of the the sigma and xi parameters be performed? If \code{FALSE}, the xi parameter will be fixed. By default, \code{FALSE}.
#' @param init.sigma a positive numeric. Indicates the starting value of sigma in the boosting algorithm.
#' @param stepsize a positive numeric. \code{learning_rate} parameter for the boosting model for sigma.
#' See \link[xgboost]{xgb.train} for more details.
#' @param stepsize_xi a positive numeric. \code{learning_rate} parameter for the boosting model for xi.
#' See \link[xgboost]{xgb.train} for more details.
#' @param tree_depth a positive numeric. \code{tree_depth} parameter for the boosting model for sigma.
#' See \link[xgboost]{xgb.train} for more details.
#' @param tree_depth_xi a positive numeric. \code{tree_depth} parameter for the boosting model for xi.
#' See \link[xgboost]{xgb.train} for more details.
#' @param nrounds a positive scalar \code{nrounds} parameter for xgboost calibration.
#' See \link[xgboost]{xgb.train} for more details.
#' @param xis a non-negative vector. Will be ignored if \code{simultaneous=TRUE}. If not supplied when \code{simultaneous=FALSE}, the argument \code{xi} will be used instead.
#' @param cv.nfold a positive scalar. Determines the number of folds for cross-validation.
#' @param nthread a positive scalar. Number of thread used for cross-validation. If not supplied, all available cores are used
#' @param orthogonal a logical. Should the two parameters be orthogonalized. By default, \code{FALSE}.
#'
#' @return a list
#' @export
#' @importFrom foreach %dopar% %do%
#'
#' @examples
GPDxgb.cv <- function(y, X, xi, simultaneous=FALSE, init.sigma, stepsize=0.1, stepsize_xi=0.1, tree_depth=5,
                         tree_depth_xi=3, nrounds, xis=NULL, cv.nfold = 5, nthread=NULL, orthogonal=FALSE){

  if (!inherits(X, "matrix") & !inherits(X, "data.frame"))
    stop("argument X must be a matrix or a data.frame")
  if (xi<0)
    stop("This function only runs for heavy-tailed and exponential-tailed responses with GPD shape xi>=0")
  if (xi>0){
    exponential=FALSE
  }
  if (xi==0){
    exponential=TRUE
  }
  if (xi==0 & simultaneous==TRUE){
    warning("The function will run with xi=0. The input simultaneous=TRUE supplied will be ignored.
            If you would like to boost the two parameters xi and sigma simultaneously, set simultaneous=TRUE and xi>0")
    simultaneous=FALSE
  }


  if (init.sigma<=0)
    stop("Sigma must be positive")


  if (orthogonal==TRUE){
    init.sigma=(init.sigma*(xi+1))
  }

  # if (nthread<=0)
  #   stop("nthread must be positive")

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

    if (is.null(xis)){
      xis <- xi
      message('No vector of xis supplied for cross-validation. Only xi will be used.')
    } else {
      if (any(xis<0)){
        stop('Supply a vector of xi for cross-validation that has non-negative entries')
      }
    }

    model_list_cv <- list()

    for (j in 1:length(xis)){
      assignInNamespace(x="log_xi", value=log(xis[j]), ns="EVTxgboost")

      if (xis[j]!=0){
        if (orthogonal==FALSE){
          param_list <- list(booster = 'gbtree',
                             objective= myobjective_gpd ,
                             tree_method = 'hist',
                             learning_rate=stepsize
                             , max_depth=tree_depth,
                             eval_metric= evalerror_gpd
        )
        } else {

          param_list <- list(booster = 'gbtree',
                             objective= myobjective_gpd_ortho ,
                             tree_method = 'hist',
                             learning_rate=stepsize
                             , max_depth=tree_depth,
                             eval_metric= evalerror_gpd_ortho
          )
        }
      } else {
        param_list <- list(booster = 'gbtree',
                           objective= myobjective_exponential ,
                           tree_method = 'hist',
                           learning_rate=stepsize
                           , max_depth=tree_depth,
                           eval_metric= evalerror_exponential
        )
      }

      if (is.null(nthread)){
        mdcv <- xgboost::xgb.cv(data=dtrain_all, params = param_list,
                       nfold=cv.nfold, nrounds=nrounds,
                       verbose = T, maximize=FALSE)
      } else {
        mdcv <- xgboost::xgb.cv(data=dtrain_all, params = param_list, nthread=nthread,
                                nfold=cv.nfold, nrounds=nrounds,
                                verbose = T, maximize=FALSE)
      }

      model_list_cv[[j]] <- mdcv

    }


    table_values <- data.frame(x = rep(1:nrounds , length(xis)),
                            y = as.numeric(sapply(1:length(xis), function(x) model_list_cv[[x]]$evaluation_log$test_GPD_loss_mean )),
                            ll95 = as.numeric(sapply(1:length(xis), function(x) model_list_cv[[x]]$evaluation_log$test_GPD_loss_mean + 1*model_list_cv[[x]]$evaluation_log$test_GPD_loss_std/sqrt(cv.nfold) )),
                            ul95 = as.numeric(sapply(1:length(xis), function(x) model_list_cv[[x]]$evaluation_log$test_GPD_loss_mean - 1*model_list_cv[[x]]$evaluation_log$test_GPD_loss_std/sqrt(cv.nfold) )),
                            xi=rep(round(xis,2),each=nrounds))



    tab.cv <- cbind(table_values$x,table_values$y, table_values$xi, table_values$ll95)
    chosen_xi <- tab.cv[which.min(tab.cv[,2]),3]
    which.xi <- which(tab.cv[,3]==chosen_xi)
    indx.min <- tab.cv[which.min(tab.cv[,2]),1]
    indx.1se <- which(tab.cv[which.xi,2]<
                        tab.cv[which.xi,4][indx.min] )[1]
    if (is.na(indx.1se)){
      indx.1se <- indx.min
    }

    table_values_2 <- table_values[which(table_values[,5]==chosen_xi),]


    gplt_xi <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y, colour = factor(xi)), linetype = 1, data = table_values) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = factor(xi)), size=0.2, data = table_values)+
      ggplot2::geom_ribbon(ggplot2::aes(x = x, y = y, ymin = ll95, ymax = ul95), colour = 'grey',
                           alpha = 0.2, data = table_values_2)  +
      ggplot2::xlab('Number of trees/boosting iterations') +ggplot2::ylab('Validation score') +
      ggplot2::theme_bw()+
      ggplot2::geom_hline(yintercept=tab.cv[which.xi,2][indx.1se], linetype='dashed', color='black', size=0.2) +
      ggplot2::geom_vline(xintercept=indx.1se, linetype='dashed', color='black', size=0.2) +
      ggplot2::labs(colour = 'xi')
    print(gplt_xi)

  } else {

    folds <- caret::createFolds(1:length(y), k = cv.nfold, list = TRUE, returnTrain = FALSE)

    runXG <- function(i){

      out_obs <- folds[[i]]


      assignInNamespace(x="log_xi", value=log(xi), ns="EVTxgboost")
      assignInNamespace(x="log_xi_out", value=log(xi), ns="EVTxgboost")
      assignInNamespace(x="log_sigma", value=log(init.sigma), ns="EVTxgboost")
      assignInNamespace(x="log_sigma_out", value=log(init.sigma), ns="EVTxgboost")

      if (orthogonal==FALSE){

        param_list <- list(booster = 'gbtree',
                           objective= myobjective_gpd ,
                           tree_method = 'hist',
                           learning_rate=stepsize
                           , max_depth=tree_depth,
                           eval_metric= evalerror_gpd_cv)

        param_list_xi <- list(booster = 'gbtree',
                              objective= myobjective_gpd_xi ,
                              tree_method = 'hist',
                              learning_rate=stepsize_xi
                              , max_depth=tree_depth_xi,
                              eval_metric= evalerror_gpd_xi_cv)
      } else {

        param_list <- list(booster = 'gbtree',
                           objective= myobjective_gpd_ortho ,
                           tree_method = 'hist',
                           learning_rate=stepsize
                           , max_depth=tree_depth,
                           eval_metric= evalerror_gpd_cv)

        param_list_xi <- list(booster = 'gbtree',
                              objective= myobjective_gpd_xi ,
                              tree_method = 'hist',
                              learning_rate=stepsize_xi
                              , max_depth=tree_depth_xi,
                              eval_metric= evalerror_gpd_xi_cv)
      }


      dtrain_all <- xgboost::xgb.DMatrix(data=as.matrix(X[-out_obs,]),label=y[-out_obs])
      xgboost::setinfo(dtrain_all, 'base_margin', rep(log(init.sigma), nrow(as.matrix(X[-out_obs,]))))

      dtrain_out <- xgboost::xgb.DMatrix(data=as.matrix(X[out_obs,]),label=y[out_obs])
      xgboost::setinfo(dtrain_out, 'base_margin', rep(log(init.sigma), nrow(as.matrix(X[out_obs,]))))

      watchlist <- list(eval = dtrain_out)

      xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                      nrounds = 1, verbose=TRUE, maximize=FALSE)

      assignInNamespace(x="log_sigma", value=predict(xgb_model, dtrain_all), ns="EVTxgboost")
      assignInNamespace(x="log_sigma_out", value=predict(xgb_model, dtrain_out), ns="EVTxgboost")

      dtrain_all_xi <- xgboost::xgb.DMatrix(data=as.matrix(X[-out_obs,]),label=y[-out_obs])
      xgboost::setinfo(dtrain_all_xi, 'base_margin', rep(log(xi), nrow(X[-out_obs,])))

      dtrain_out_xi <- xgboost::xgb.DMatrix(data=as.matrix(X[out_obs,]),label=y[out_obs])
      xgboost::setinfo(dtrain_out_xi, 'base_margin', rep(log(xi), nrow(X[out_obs,])))
      watchlist_xi <- list(eval = dtrain_out_xi)


      xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                         nrounds = 1, verbose=TRUE, maximize=FALSE)

      assignInNamespace(x="log_xi", value=predict(xgb_model_xi, dtrain_all_xi), ns="EVTxgboost")
      assignInNamespace(x="log_xi_out", value=predict(xgb_model_xi, dtrain_out_xi), ns="EVTxgboost")

      #log_xi <- predict(xgb_model_xi, dtrain_all_xi)


      if (nrounds>1){
        for (iter in 2:nrounds){
          xgb_model <- xgboost::xgb.train(params = param_list, data = dtrain_all, watchlist=watchlist,
                                          nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model)
          assignInNamespace(x="log_sigma", value=predict(xgb_model, dtrain_all), ns="EVTxgboost")
          assignInNamespace(x="log_sigma_out", value=predict(xgb_model, dtrain_out), ns="EVTxgboost")
          xgb_model_xi <- xgboost::xgb.train(params = param_list_xi, data = dtrain_all_xi, watchlist=watchlist_xi,
                                             nrounds = 1, verbose=TRUE, maximize=FALSE, xgb_model = xgb_model_xi)
          assignInNamespace(x="log_xi", value=predict(xgb_model_xi, dtrain_all_xi), ns="EVTxgboost")
          assignInNamespace(x="log_xi_out", value=predict(xgb_model_xi, dtrain_out_xi), ns="EVTxgboost")
        }
      }

      return(matrix(rbind(xgb_model$evaluation_log$eval_GPD_loss ,
                          xgb_model_xi$evaluation_log$eval_GPD_loss)))
    }

    if (is.null(nthread)){
      cl <- parallel::makeCluster(parallel::detectCores())
      doParallel::registerDoParallel(cl)
      table_values_par <- foreach::foreach(i=1:cv.nfold, .combine=c,
                                       .packages = 'EVTxgboost') %dopar% {
        runXG(i)
      }
      parallel::stopCluster(cl)
    } else if (nthread==1){
    # RUN FUNCTION
      table_values_par <- foreach::foreach(i=1:cv.nfold, .combine=c,
                                           .packages = 'EVTxgboost') %do% {
                                             runXG(i)
                                           }
    } else {
      cl <- parallel::makeCluster(nthread)
      doParallel::registerDoParallel(cl)
      #parallel::clusterCall(cl, function() {library(EVTxgboost)})
      table_values_par <- foreach::foreach(i=1:cv.nfold, .combine=c,
                                       .packages = 'EVTxgboost') %dopar% {
        runXG(i)
      }
      parallel::stopCluster(cl)
    }

    mean_error <- apply(matrix(table_values_par, ncol=cv.nfold), 1, mean)
    std_error <- apply(matrix(table_values_par, ncol=cv.nfold), 1, function(x) sqrt(var(x)) )

    table_values <- data.frame(x = 1:(nrounds*2),
                               y = mean_error,
                               ll95 = mean_error + 1*std_error/sqrt(cv.nfold) ,
                               ul95 = mean_error- 1*std_error/sqrt(cv.nfold) )

    tab.cv <- cbind(table_values$x,table_values$y, table_values$ll95)
    indx.min <- tab.cv[which.min(tab.cv[,2]),1]
    indx.1se <- which(tab.cv[,2]<
                        tab.cv[indx.min,3] )[1]
    if (is.na(indx.1se)){
      indx.1se <- indx.min
    }

    gplt_xi <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y), linetype = 1, data = table_values) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y), size=0.2, data = table_values)+
      ggplot2::geom_ribbon(ggplot2::aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = table_values)  +
      ggplot2::xlab('Number of trees/boosting iterations') +ggplot2::ylab('Validation score') +
      ggplot2::theme_bw()+
      ggplot2::geom_hline(yintercept=tab.cv[,2][indx.1se], linetype='dashed', color='black', size=0.2) +
      ggplot2::geom_vline(xintercept=indx.1se, linetype='dashed', color='black', size=0.2) +
      ggplot2::labs(colour = 'xi')
    print(gplt_xi)


    chosen_xi=NULL

  }

  result_list <- list()
  result_list$table_values <- table_values
  result_list$chosen_xi <- chosen_xi
  result_list$indx.min <- indx.min
  result_list$indx.1se <- indx.1se


  return(result_list)
}

