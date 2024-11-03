

n <- 1000
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


test_result <- GPDxgb.train(y,X, xi=0.4, init.sigma=3, stepsize=0.02, stepsize_xi=0.01,
                            tree_depth = 3,
                            tree_depth_xi = 3,
                            nrounds=300,
                            simultaneous = TRUE, orthogonal = FALSE)

boxplot((beta_truth-test_result$sigma.pred))
abline(h=0)
boxplot(xi_truth-test_result$xi.pred)
abline(h=0)

test_cv <- GPDxgb.cv(y,X, xi=0.2, xis=seq(0.5,0.6,by=0.05),
                     init.sigma=1, stepsize=0.01, stepsize_xi=0.01,
                     tree_depth = 3 ,
                     tree_depth_xi = 3,
                     nrounds=100,
                     simultaneous = FALSE, cv.nfold = 5)

test_cv <- GPDxgb.cv(y,X, xi=0.2, xis=seq(0.5,0.6,by=0.05),
                     init.sigma=2, stepsize=0.01, stepsize_xi=0.05,
                     tree_depth = 3 ,
                     tree_depth_xi = 3,
                     nrounds=10,
                     simultaneous = FALSE, cv.nfold = 5, orthogonal=TRUE)

test_cv$chosen_xi
test_cv$indx.1se

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





###########################
## Fit to data . For wildfire chapter
###########################

data(Cali_DF)
colnames(Cali_DF)

data(USA)


data("GB_DF")
DF <- rbind(GB_DF)
# DF <- rbind(Cali_DF)
#DF <- data_train_DF

indx.train <- which(DF$year<2011)
indx.test <- which(DF$year>=2011)

DF$BA <- DF$BA/247.1
hist(log10(DF$BA[indx.train]+1), breaks=150, xlab='log(BA+1)', main='')


library(mev)

seqs <- seq(quantile(DF$BA,0.8),quantile(DF$BA,0.99), length.out=35)
par_est <- numeric()
st_est <- numeric()

for (i in 1:length(seqs)){
  par_est[i] <- mev::gp.fit(DF$BA,seqs[i])$estimate[2]
  st_est[i] <- mev::gp.fit(DF$BA,seqs[i],show=FALSE)$std.err[2]
}

pdf(file='ba_dist.pdf', width=5.5, height=2.5)
par(mfrow=c(1,3))
hist(log10(DF$BA+1), xlab='log(BA+1)', main='')
s <- NC.diag( x=DF$BA, u=seqs, my.xlab='BA (ac)')
u <- seqs[which(s$e.p.values>0.05)[1]]
abline(v=u, lty='dashed')
plotrix::plotCI(x=seqs,y=par_est,li=par_est- 1.96*st_est,
                ui=par_est+ 1.96*st_est, pch=19, cex=0.5, ylab='xi',
                xlab='BA (ac)')
abline(v=u, lty='dashed')
dev.off()

ecdf(DF$BA)(u)

# u <- quantile(DF$BA[which(DF$BA>0)], 0.85)

y <- DF$BA[indx.train]-u
indx.excess <- which(y-u>0)
y <- y[indx.excess]
X <- DF[indx.train,3:37][indx.excess,][,-c(3,4)]

fit_gpd <- mev::gp.fit(DF$BA[indx.train],u)$estimate
init.sigma <- fit_gpd[1]
xi.start <- fit_gpd[2]

pdf(file='cv_xi.pdf', width=6, height=4)
model_cv <- GPDxgb.cv(y, X, xi=xi.start, init.sigma=init.sigma,
                      stepsize=0.01, nrounds=200, simultaneous=FALSE,
                      xis=seq(0.2,0.8,by=0.05), cv.nfold = 3)
dev.off()

model_cv$indx.1se
model_cv$chosen_xi
model_cv$indx.min




y_test <- DF$BA[indx.test]-u
indx.excess.test <- which(y_test-u>0)
X_test <- DF[indx.test,3:37][indx.excess.test,][,-c(3,4)]


model_fit <- GPDxgb.train(y,X, xi=model_cv$chosen_xi, init.sigma=init.sigma, stepsize=0.01,
                            nrounds=model_cv$indx.1se,
                            simultaneous = FALSE, X_test=X_test )

model_fit$xgb_model
model_fit$sigma.pred
model_fit$xi


xgboost::xgb.importance(model=model_fit$xgb_model)

hist(X_test$clim6, breaks=10)
hist(X_test$clim2, breaks=10)




model_cv_2 <- GPDxgb.cv(y,X, xi=model_cv$chosen_xi,init.sigma=fit_gpd[1],
                        stepsize=0.01,stepsize_xi=0.01, nrounds=200,
                      simultaneous = TRUE, cv.nfold = 3, orthogonal = FALSE)

model_cv_2$indx.1se

model_fit_2 <- GPDxgb.train(y,X, xi=model_cv$chosen_xi, init.sigma=fit_gpd[1], stepsize=0.01,stepsize_xi=0.01,
                          nrounds=round(model_cv_2$indx.1se/2),
                          simultaneous = TRUE, X_test=X_test , orthogonal = FALSE)

max(model_fit_2$xi.pred)
min(model_fit_2$xi.pred)
plot(model_fit_2$xi.pred)


dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
model_cv_3 <- xgboost::xgb.cv(data=dtrain, eta = 0.01, nrounds = 200, nfold=5)
plot(model_cv_3$evaluation_log$test_rmse_mean)
model_fit_3 <- xgboost::xgboost(data = dtrain, eta = 0.01, nrounds = 100)
dtest <- xgboost::xgb.DMatrix(data = as.matrix(X_test), label = y_test[indx.excess.test])

fit_gpd <- mev::gp.fit(DF$BA[indx.train],u)$estimate
#init.sigma <- fit_gpd[1]

test_excess <- y_test[indx.excess.test]; n.exc <- length(test_excess)
crps_basic <- scoringRules::crps_gpd(test_excess, scale=rep(fit_gpd[1], n.exc ),
                                     shape=rep(fit_gpd[2], n.exc ) )
crps_model <- scoringRules::crps_gpd(test_excess, scale=model_fit$sigma.pred,
                                     shape=rep(model_fit$xi, n.exc ) )
crps_model_2 <- scoringRules::crps_gpd(test_excess, scale=model_fit_2$sigma.pred,
                                       shape=(model_fit_2$xi.pred) )
crps_model_3 <- scoringRules::crps_norm(test_excess, mean = predict(model_fit_3, dtest),
                                        sd=1)

hist(evd::pgpd(q=test_excess, scale = model_fit$sigma.pred, shape = model_fit$xi), breaks=25)
hist(evd::pgpd(q=test_excess, scale = fit_gpd[1], shape = fit_gpd[2]), breaks=25)

mean(crps_basic)
mean(crps_model)
mean(crps_model_2)
mean(crps_model_3)

level <- 1000

plot(test_excess, x=1:length(test_excess), pch=19, cex=0.5)
abline(h=level, lty='dashed')
probs <- evd::pgpd(q=level, scale = model_fit$sigma.pred, shape = model_fit$xi,lower.tail=FALSE )
par(new=TRUE)
plot(y=probs, x=1:length(test_excess), pch=4, xlim=c(0,length(test_excess)),
     axes=FALSE, col='red', cex=0.5, xlab='', ylab='', type='l')
par(new=TRUE)
plot(y=probs[which(test_excess>=level)], x=c(which(test_excess>=level)), pch=4, xlim=c(0,length(test_excess)),
     axes=FALSE, col='red', cex=0.5, xlab='', ylab='')

abline(h=evd::pgpd(q=level, scale = fit_gpd[1], shape = fit_gpd[2],lower.tail=FALSE ), pch=19, col='red', cex=0.5,
       lty='dashed')
mtext("Cell Density",side=4,col="red",line=4)
axis(4, col="red",col.axis="red",las=1)

test_excess[39]
X_test[39,]
probs[39]
probs_2

plot(evd::dgpd(seq(0.1,5000,length.out=100), scale = fit_gpd[1], shape = fit_gpd[2], log=TRUE), x=seq(0.1,5000,length.out=100),
     type='l')
lines(evd::dgpd(seq(0.1,5000,length.out=100), scale =  model_fit$sigma.pred[39], shape =  model_fit$xi, log=TRUE), x=seq(0.1,5000,length.out=100))
abline(v=test_excess[39], col='red')


plot(evd::dgpd(seq(0.1,5000,length.out=100), scale = fit_gpd[1], shape = fit_gpd[2], log=TRUE), x=seq(0.1,5000,length.out=100),
     type='l')
lines(evd::dgpd(seq(0.1,5000,length.out=100), scale =  model_fit$sigma.pred[which(test_excess>=level)[1]], shape =  model_fit$xi, log=TRUE), x=seq(0.1,5000,length.out=100))
abline(v=test_excess[which(test_excess>=level)[1]], col='red')
hist(X_test$clim6, breaks=10)
abline(v=X_test$clim6[which(test_excess>=level)[1]])
hist(X_test$clim2, breaks=10)
abline(v=X_test$clim2[which(test_excess>=level)[1]])


probs_2 <- evd::pgpd(q=500, scale = fit_gpd[1], shape = fit_gpd[2],lower.tail=FALSE )






save(model_cv,crps_basic, crps_model, crps_model_2, file='model_test.rda')



crps_basic <- scoringRules::logs_gpd(y_test[indx.excess.test], scale=rep(fit_gpd[1], length(indx.excess.test) ),
                                     shape=rep(fit_gpd[2], length(indx.excess.test) ) )
crps_model <- scoringRules::logs_gpd(y_test[indx.excess.test], scale=model_fit$sigma.pred,
                                     shape=rep(model_fit$xi, length(indx.excess.test) ) )
crps_model_2 <- scoringRules::logs_gpd(y_test[indx.excess.test], scale=model_fit_2$sigma.pred,
                                       shape=(model_fit_2$xi.pred) )

mean(crps_basic)
mean(crps_model)
mean(crps_model_2)

u_vec <- numeric()

for (j in 1:length(indx.excess.test)){
  u_vec[j] <- evd::pgpd(y_test[indx.excess.test][j] , scale=model_fit$sigma.pred[j],
                        shape=model_fit$xi)
}
hist(u_vec)

fit_gpd <- mev::gp.fit(DF$BA[indx.train],u)$estimate

for (j in 1:length(indx.excess.test)){
  u_vec[j] <- evd::pgpd(y_test[indx.excess.test][j] , scale=fit_gpd[1],
                        shape=fit_gpd[2])
}

hist(u_vec)


for (j in 1:length(indx.excess.test)){
  u_vec[j] <- evd::pgpd(y_test[indx.excess.test][j] , scale=model_fit_2$sigma.pred[j],
                        shape=model_fit_2$xi.pred[j])
}

hist(u_vec)


mean(crps_model)
mean(crps_basic)






devtools::install_github("ModelOriented/shapviz", dependencies = TRUE)
library(shapviz)  # Needs development version 0.9.0 from github

dmat_train <- xgboost::xgb.DMatrix(data=as.matrix(X),label=y)

sv <- shapviz(test_result$xgb_model_xi,
              X_pred =  dmat_train, X=X )
sv_dependence(
  sv,
  v = c("clim1", "clim4"),
  alpha = 0.2
)

sv_dependence2D(sv, x = "lon", y = "lat") +
  coord_equal()

plot(DF$lon, DF$lat)


test_result <- GPDxgb.train(y,X, xi=chosen_xi, init.sigma=10, stepsize=0.01, tree_depth = 5,
                            nrounds=test_cv$indx.1se,stepsize_xi = 0.02,
                            simultaneous = TRUE)
test_result$xgb_model_xi

plot(test_result$xi.pred)


x.dat <- data.frame(lon=X$lon, lat=X$lat, xi=(test_result$xi.pred) )

plot.dat <- x.dat %>%
  group_by(lon, lat) %>%
  summarise(avg = mean(xi))


gplt <- ggplot(data = ca_county) +
  theme_bw()+
  coord_quickmap() +
  geom_tile(data=plot.dat, aes(x=lon, y=lat, fill= avg), colour='white',
            show.legend = TRUE, inherit.aes = FALSE) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "transparent", color = "black")




library(ggmap)
states <- map_data("state")

west_coast <- states %>%
  filter(region %in% c("california", "oregon", "washington"))

test_plt <- sv_dependence2D(sv, x = "lon", y = "lat")
test_plt$data

counties <- map_data("county")
ca_county <- counties %>%
  filter(region == "california")


gplt <- ggplot(data = ca_county) +
  theme_bw()+
  coord_quickmap() +
  geom_tile(data=test_plt$data, aes(x=lon, y=lat, fill= SHAP), colour='white',
            show.legend = TRUE, inherit.aes = FALSE) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "transparent", color = "black")


print(gplt)


ggplot(data = west_coast) + geom_tile(data=test_plt$data, aes(x=lon, y=lat, fill= SHAP), colour='white', show.legend = TRUE, inherit.aes = FALSE)+
  theme_bw() +  coord_fixed()



