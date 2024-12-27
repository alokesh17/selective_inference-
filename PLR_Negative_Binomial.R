library(glmnet)
library(MASS)
n = 2000
p = 40
Cov = diag(p)
 mu <- rep(1,p)##Made it otherwise mean will be negative
 beta <- rep(0,(p+1))
 beta[1:11] <- c(.5,c(1:10)/15)#c(1:10)/15
 
m <- length(beta)
X <- MASS::mvrnorm(n,mu,Cov)
Xb <-  cbind(1,X) %*% beta
pr <- 1/(1+exp(-Xb))
theta=4.5

Y <- rnegbin(n=n, mu = exp(Xb), theta = theta)#rbinom(n = n, prob = pr, size = 1)
# loading <- MASS::mvrnorm(1,mu,Cov/2)
# loading[11:p] <- loading[11:p]/25

 
# cvfit <- cv.glmregNB( Y~.,data=data.frame(Y,X),alpha=0.1)
# lambda_chosen=cvfit$lambda.optim #based on cv which gives minimum deviance
# coefficients_estimated=as.numeric(c(coef(cvfit)))
# beta.lasso.approx1=coefficients_estimated
experiment= 50
result_grand_nb=matrix(rep(0,experiment*(p+1)),nrow=p+1)
diff_grand_nb=matrix(rep(0,experiment*(p+1)),nrow=p+1)

for( l in 1:experiment){
simulations=50

betamatrix=matrix(rep(0,simulations*(p+1)),ncol=simulations)

for(j in 1:simulations){
samp=sample(1:n,n,rep=T)
Xsamp <- X[samp,]
Ysamp <- Y[samp]

theta_est=glm.nb(Ysamp~.,data=data.frame(Ysamp,Xsamp))$theta

cvfit_samp <-suppressWarnings(mpath:: cv.glmreg( Ysamp~.,data=data.frame(Ysamp,Xsamp),alpha=1,family="negbin", theta = theta_est,nfolds=5,n.cores=3,plot.it=FALSE,parallel=TRUE))#cv.glmnet(Xsamp, Ysamp, negative.binomial(theta =theta ))
#lambda_chosen=cvfit_samp$lambda.min #based on cv which gives minimum deviance
pblasso=as.numeric(c(coef(cvfit_samp)))##Bootstrap lasso estimate pblasso
s_pblasso=which(pblasso!=0)               ##which are nonzero in pblasso
s_pblasso_complement=which(pblasso==0) ##which are zero in pblasso
penalty=rep(0,(p+1))#0 means always included, no shrinkage
penalty[s_pblasso_complement]=1# only overshrunk 0 variables are done with ridge
if(sum(penalty)==1) penalty[s_pblasso[2]]=1 #in case if just one is zero. penalty.factor did not work in that case.So just select randomly another one.
cvfit_pblpr=suppressWarnings(mpath::cv.glmreg( Ysamp~.,data=data.frame(Ysamp,Xsamp),alpha=0.001,penalty.factor=penalty[2:(p+1)],family="negbin", theta =theta_est,nfolds=5,n.cores=3,plot.it=FALSE,parallel=TRUE))
#cv.glmnet(Xsamp, Ysamp, family = "binomial",alpha = 0,penalty.factor=penalty[2:(p+1)])##Fit a ridge regression ##Alpha cannot be zero here it gives error

beta_pblpr=as.numeric(coef(cvfit_pblpr))


beta.lasso.approx <- beta_pblpr
betamatrix[,j]=c(beta.lasso.approx)
#print(j)
}

#adjustments=theta_hat%*%sigma%*%t(theta_hat)*(1/n)
lower_CI=sapply(1:(p+1), function(k){
 min(quantile(c(betamatrix[k,]),prob=c(.025)), quantile(c(betamatrix[k,]),prob=c(.925))) 
  
})
#b_hat-1.96*sqrt(diag(adjustments))
upper_CI=sapply(1:(p+1), function(k){
    max(quantile(c(betamatrix[k,]),prob=c(.025)), quantile(c(betamatrix[k,]),prob=c(.925))) 

})#b_hat-1.96*sqrt(diag(adjustments))

diff=upper_CI-lower_CI

result=c(sapply(1:(p+1), function(i){ifelse((c(beta[i])>c(lower_CI[i])) && (c(beta[i])<c(upper_CI[i])),1,0)}))
#print(cbind(beta,lower_CI,upper_CI,diff,result))

  
result_grand_nb[,l]=result
diff_grand_nb[,l]=diff
print(c("experiment number:",l))
}
