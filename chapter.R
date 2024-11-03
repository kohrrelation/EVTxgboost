#library(devtools); install_github("kohrrelation/EVTxgboost")
library(EVTxgboost); library(mev)
library(scoringRules); library(plotrix)
data(GB_DF)


head(colnames(GB_DF)); DF <- GB_DF

DF <- data_train_DF

indx.train <- which(DF$year<2011)
indx.test <- which(DF$year>=2011)
seqs <- seq(quantile(DF$BA[indx.train],0.85),
            quantile(DF$BA[indx.train],0.995), length.out=40)
par_est <- st_est <- numeric()
for (i in 1:length(seqs)){
  par_est[i] <- gp.fit(DF$BA,seqs[i])$estimate[2]
  st_est[i] <- gp.fit(DF$BA,seqs[i])$std.err[2]
}
par(mfrow=c(1,3))
hist(log10(DF$BA[indx.train]+1), breaks=30, xlab='log(BA+1)')
s <- NC.diag(x=DF$BA[indx.train], u=seqs, my.xlab='BA (ac)')
u <- seqs[which(s$e.p.values>0.1)[1]]
abline(v=u, lty='dashed')
plotCI(x=seqs, y=par_est, li=par_est-1.96*st_est,
       ui=par_est+ 1.96*st_est, pch=19, cex=0.5, ylab='xi')
abline(v=u, lty='dashed')

fit_gpd <- mev::gp.fit(DF$BA[indx.train], u)$estimate
init.sigma <- fit_gpd[1]
y <- DF$BA[indx.train] - u; which.excess <- which(y - u>0)
y <- y[which.excess]; X <- DF[indx.train,3:37][which.excess,]

model_cv <- GPDxgb.cv(y, X, xi=0.4, init.sigma=init.sigma,
                      stepsize=0.02, nrounds=200, simultaneous=FALSE,
                      xis=seq(0.3,0.6,by=0.05), cv.nfold = 5)
c(model_cv$indx.1se, model_cv$chosen_xi, model_cv$indx.min)

y_test <- DF$BA[indx.test] - u
indx.excess.test <- which(y_test - u>0)
X_test <- DF[indx.test,3:37][indx.excess.test,]
model_fit <- GPDxgb.train(y, X, xi=model_cv$chosen_xi,
                          init.sigma=init.sigma, stepsize=0.02,
                          nrounds=model_cv$indx.1se, simultaneous=FALSE,
                          X_test=X_test)


model_cv_2 <- GPDxgb.cv(y, X, xi=0.5, init.sigma=init.sigma,
                        stepsize=0.01, stepsize_xi=0.01, nrounds=400,
                        simultaneous=TRUE, cv.nfold = 5)


model_fit_2 <- GPDxgb.train(y, X, xi=0.5, init.sigma=init.sigma,
                            stepsize=0.01, stepsize_xi=0.01,
                            nrounds=round(model_cv_2$indx.1se/2),
                            simultaneous = TRUE, X_test=X_test)


dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
model_cv_3 <- xgboost::xgb.cv(data=dtrain, eta = 0.01, nrounds = 200, nfold=5)
plot(model_cv_3$evaluation_log$test_rmse_mean)
model_fit_3 <- xgboost::xgboost(data = dtrain, eta = 0.01, nrounds = 100)
dtest <- xgboost::xgb.DMatrix(data = as.matrix(X_test), label = y_test[indx.excess.test])
mean_basic <- rep(fit_gpd[1], n.exc)/(1-rep(fit_gpd[2], n.exc))
mean_model <- model_fit$sigma.pred/rep(model_fit$xi, n.exc)
mean_model_2 <- model_fit_2$sigma.pred/model_fit_2$xi.pred
mean_sqloss <- predict(model_fit_3, dtest)

sqrt(mean((mean_basic-test_excess)^2))
sqrt(mean((mean_model-test_excess)^2))
sqrt(mean((mean_model_2-test_excess)^2))
sqrt(mean((mean_sqloss-test_excess)^2))



test_excess <- y_test[indx.excess.test]; n.exc <- length(test_excess)
crps_basic <- crps_gpd(test_excess, scale=rep(fit_gpd[1], n.exc),
                       shape=rep(fit_gpd[2], n.exc))
crps_model <- crps_gpd(test_excess, scale=model_fit$sigma.pred,
                       shape=rep(model_fit$xi, n.exc))
crps_model_2 <- crps_gpd(test_excess,
                         scale=model_fit_2$sigma.pred,
                         shape=(model_fit_2$xi.pred))
crps_model_3 <- scoringRules::crps_norm(test_excess,
                                        mean=predict(model_fit_3, dtest),sd=1)
c(mean(crps_basic), mean(crps_model), mean(crps_model_2), mean(crps_model_3))

xgboost::xgb.importance(model=model_fit$xgb_model)
xgboost::xgb.importance(model=model_fit_2$xgb_model_xi)



