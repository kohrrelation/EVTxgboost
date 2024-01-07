
n <- 10000
X <- cbind(rnorm(n,0,2), rnorm(n,0,2), rnorm(n,0,2), rnorm(n,0,2))
y <- fExtremes::rgpd(n, xi=0.2, beta = exp(1+X[,1]) )

hist(exp(-1+(X[,1]*X[,2])/200), breaks=1000)

y <- numeric()
for (i in 1:n){
  y[i] <- fExtremes::rgpd(1, xi=exp(-1+(X[i,1]*X[i,2])/200), beta = exp(1+X[i,1]) )
}

test_cv <- GPDxgb.cv(y,X, xi=0.25, init.sigma=10, stepsize=0.1, tree_depth = 5 , nrounds=100,
                     simultaneous = TRUE, cv.nfold = 5, nthread = 5)

test_cv$indx.1se
test_cv$table_values
test_cv$indx.min

test_result <- GPDxgb.train(y,X, xi=0.25, init.sigma=10, stepsize=0.1, tree_depth = 5,
                            nrounds=test_cv$indx.1se/2,
                            simultaneous = TRUE)

meds <- meds_truth <- numeric()

for (j in 1:n){
  meds[j] <- fExtremes::qgpd(p=0.5, xi = test_result$xi.pred[j],
                          beta = test_result$sigma.pred[j])
  meds_truth[j] <- fExtremes::qgpd(p=0.5, beta = exp((1+X[j,1])),
                                xi = exp(-1+(X[j,1]*X[j,2])/200))
}

plot(x=meds, y=meds_truth, pch=19, cex=0.2)
abline(a=0,b=1, lty='dashed')

plot( y=(test_result$sigma.pred), x=exp(1+X[,1]) )
abline(a=0,b=1, lty='dashed')








# second example


test_cv_2 <- GPDxgb.cv(y,X, xi=0.25, init.sigma=10, stepsize=0.1, tree_depth = 5 , nrounds=40,
                       simultaneous = FALSE, cv.nfold = 5)
test_cv_2$table_values

test_cv_2 <- GPDxgb.cv(y,X, xi=0.25, init.sigma=10, stepsize=0.1, tree_depth = 5 , nrounds=40,
                            simultaneous = FALSE, xis=seq(0.15,0.25,by=0.05), cv.nfold = 5)

test_cv$table_values[[1]]

test_cv$indx.1se
test_cv$table_values
chosen_xi <- test_cv$chosen_xi
test_cv$indx.min


test_result <- GPDxgb.train(y,X, xi=chosen_xi, init.sigma=10, stepsize=0.1, tree_depth = 5,
                            nrounds=test_cv$indx.min,
                            simultaneous = FALSE)

meds <- fExtremes::qgpd(p=0.5, xi = chosen_xi,
                           beta = test_result$sigma.pred)
meds_truth <- fExtremes::qgpd(p=0.5, beta = exp((1+X[,1])),
                                 xi = 0.2)
plot(x=meds, y=meds_truth, pch=19, cex=0.2)
abline(a=0,b=1, lty='dashed')

plot( y=(test_result$sigma.pred), x=exp(1+X[,1]) )
abline(a=0,b=1, lty='dashed')




init.sigma <- 1
stepsize=0.1
tree_depth=3
rounds=100
xi=0.1
exponential=FALSE

#boxplot(exp(predict(xgb_model_xi, dtrain_all_xi)))
#abline(h=0.1)
plot((predict(xgb_model, dtrain_all)), x= (1+X[,1]), cex=0.2, pch=19)
abline(a=0,b=1)
plot(exp(predict(xgb_model_xi, dtrain_all_xi)), x= exp(-1+(X[,1]*X[,2])/100), cex=0.2, pch=19)
abline(a=0,b=1)

meds <- numeric()
meds_truth <- numeric()

xis <- exp(predict(xgb_model_xi, dtrain_all_xi))
betas <- exp(predict(xgb_model, dtrain_all))

for (i in 1:n){
  meds[i] <- fExtremes::qgpd(p=0.5, xi = xis[i],
                             beta = betas[i])
  meds_truth[i] <- fExtremes::qgpd(p=0.5, beta = exp((1+X[i,1])),
                                   xi = exp(-1+(X[i,1]*X[i,2])/100))
}

plot((meds), x=(meds_truth), pch=19, cex=0.1 )
abline(a=0,b=1)

plot((predict(xgb_model, dtrain_all)), x= X[,2], cex=0.2, pch=19)
plot(exp(predict(xgb_model_xi, dtrain_all_xi)), x= X[,2], cex=0.2, pch=19)
