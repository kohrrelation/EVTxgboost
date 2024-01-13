

n <- 30000
X <- cbind(rnorm(n,0,2), rnorm(n,1,2), rnorm(n,4,2), rnorm(n,8,2))

beta_truth <- exp(1+X[,1]/10 + (X[,4]>8)/10 )
xi_truth <- exp(-1+(X[,1]*X[,2])/100 +(X[,3]>4)/3)

hist(xi_truth, breaks=1000)
hist(beta_truth, breaks=1000)

y <- numeric()
for (i in 1:n){
  y[i] <- fExtremes::rgpd(1, xi=xi_truth[i],
                          beta = beta_truth[i]  )
}

mean(beta_truth)
mean(xi_truth)

test_cv <- GPDxgb.cv(y,X, xi=0.4, init.sigma=3, stepsize=0.01, stepsize_xi=0.01,
                     tree_depth = 3 ,
                     tree_depth_xi = 3,
                     nrounds=150,
                     simultaneous = TRUE, cv.nfold = 5)

test_cv$indx.1se
test_cv$table_values
test_cv$indx.min

test_result <- GPDxgb.train(y,X, xi=0.4, init.sigma=3, stepsize=0.01, stepsize_xi=0.01,
                            tree_depth = 3,
                            tree_depth_xi = 3,
                            nrounds=test_cv$indx.1se/2,
                            simultaneous = TRUE)

beta_preds <- test_result$sigma.pred
xi_preds <- test_result$xi.pred

meds <- meds_truth <- numeric()

for (j in 1:n){
  meds[j] <- fExtremes::qgpd(p=0.5, beta = beta_preds[j],
                             xi = xi_preds[j])
  meds_truth[j] <- fExtremes::qgpd(p=0.5, beta = beta_truth[j],
                                xi = xi_truth[j])
}

plot(x=meds, y=meds_truth, pch=19, cex=0.2)
abline(a=0,b=1, lty='dashed')





y <- fExtremes::rgpd(n, xi=0.2, beta = exp(1+X[,1] + (X[,4]>0)/10) )



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






## Fit to data
.
colnames(Cali_DF)

quantile(Cali_DF$BA, 0.95)

which.excess <- which(y-quantile(Cali_DF$BA, 0.95)>0)

y <- Cali_DF$BA-quantile(Cali_DF$BA, 0.95)
y <- y[which.excess]
X <- Cali_DF[which.excess,3:38]

test_cv <- GPDxgb.cv(y,X, xi=0.25, init.sigma=10, stepsize=0.1, tree_depth = 5 , nrounds=40,
                       simultaneous = FALSE, xis=seq(0.95,1.25,by=0.05), cv.nfold = 5)

test_cv$table_values[[1]]

test_cv$indx.1se
test_cv$table_values
chosen_xi <- test_cv$chosen_xi
test_cv$indx.min

chosen_xi

test_result <- GPDxgb.train(y,X, xi=chosen_xi, init.sigma=10, stepsize=0.1, tree_depth = 5,
                            nrounds=test_cv$indx.1se,
                            simultaneous = FALSE)

test_result$xgb_model







