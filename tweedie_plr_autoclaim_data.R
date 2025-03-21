library(glmnet)
## LOAD AUTO CLAIM DATA
data(AutoClaim, package = "cplm")

# ## REMOVE FEATURES
# x <- AutoClaim[, -c(1:5, 10, 16, 29)]
# 
# ## INTEGER CODING TABLE FOR CATEGORICAL FEATURES
# coding <- sapply(x, function(x) {
#   if(is.factor(x))
#     data.frame(
#       code = seq_along(levels(x)),
#       level = levels(x)
#     )
# })
# coding[sapply(coding, is.null)] <- NULL
# print(coding)
# 
# ## ENCODE CATEGORICAL FEATURES INTO INTEGERS
# index <- sapply(x, is.factor)
# x[index] <- lapply(x[index], as.integer)
# 
# dt <- list(data = as.matrix(x), label = AutoClaim$CLM_AMT5 / 1000)
# 
# X=as.matrix(x)
# Y=AutoClaim$CLM_AMT5 / 1000
# n=length(Y)
# p=dim(X)[2]
# 
# impute_median <- function(x) {
#   x[is.na(x)] <- median(x, na.rm = TRUE)
#   return(x)
#}

## REMOVE FEATURES
x <- AutoClaim[, -c(1:5, 10, 16, 29)]

str(x)
# Apply the function to each column of the matrix
X_imputed <- apply(X, 2, impute_median)

X_imputed <- apply(x, 2, impute_median)

library(dplyr)
library(tidyr)
x <- x %>% 
  mutate(across(where(is.numeric), ~replace_na(., median(., na.rm = TRUE))))


X_imputed <- x

#X_imputed <- as.matrix(X_imputed)

str(X_imputed); str(x)

install.packages("fastDummies")
library(fastDummies)



X_imputed <- fastDummies::dummy_cols(x, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

sum(is.na(X_imputed))  # Total count of NA values

X_imputed <- data.frame(lapply(X_imputed, function(x) {
  if (is.factor(x)) as.numeric(as.factor(x)) else x
}))

X_imputed <- as.matrix(X_imputed)
mode(X_imputed) <- "numeric"  # Ensure all values are numeric


fit_1 <- cv.glmnet(X_imputed, Y, family = statmod::tweedie, alpha = 1)

fit_2 <- cv.glmnet(X_imputed, Y, family = statmod::tweedie, alpha = 0)


coef(fit_1)

coef(fit_2)


simulations=100
 Y=AutoClaim$CLM_AMT5 / 1000
 n=length(Y)
p=dim(X_imputed)[2]

betamatrix=matrix(rep(0,simulations*(p+1)),ncol=simulations)
for(j in 1:simulations){
  samp=sample(1:n,n,rep=T)
  Xsamp <- X_imputed[samp,]
  Ysamp <- Y[samp]
  cvfit_samp <-suppressWarnings(cv.glmnet( Xsamp,Ysamp,data=data.frame(Ysamp,Xsamp),alpha=1,family=statmod::tweedie))#cv.glmnet(Xsamp, Ysamp, negative.binomial(theta =theta ))
  #lambda_chosen=cvfit_samp$lambda.min #based on cv which gives minimum deviance
  pblasso=as.numeric(as.vector(coef(cvfit_samp)))##Bootstrap lasso estimate pblasso
  s_pblasso=which(pblasso!=0)               ##which are nonzero in pblasso
  s_pblasso_complement=which(pblasso==0) ##which are zero in pblasso
  penalty=rep(0,(p+1))#0 means always included, no shrinkage
  penalty[s_pblasso_complement]=1##only zero variables to be included in 
  if(sum(penalty)==1) penalty[s_pblasso[2]]=1 #in case if just one is zero. penalty.factor did not work in that case.So just select randomly another one.
  #cvfit_pblpr=suppressWarnings(cv.glmnet(Xsamp, Ysamp,data=data.frame(Ysamp,Xsamp),alpha=0,
  #                                       family=statmod::tweedie,penalty.factor=penalty)) ##Change
  cvfit_pblpr=suppressWarnings(cv.glmnet(Xsamp, Ysamp,data=data.frame(Ysamp,Xsamp),alpha=0,
                                        family=statmod::tweedie,penalty.factor=penalty[2:(p+1)]))
  #cv.glmnet(Xsamp, Ysamp, family = "binomial",alpha = 0,penalty.factor=penalty[2:(p+1)])##Fit a ridge regression ##Alpha cannot be zero here it gives error
  
  beta_pblpr=as.numeric(coef(cvfit_pblpr))
  
  
  beta.lasso.approx <- beta_pblpr
  betamatrix[,j]=c(beta.lasso.approx)
  print(j)
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


df=data.frame(variable=c("intercept",colnames(X_imputed)))|>
   bind_cols(lower_CI=lower_CI)|>
   bind_cols(upper_CI=upper_CI)|>
   bind_cols(length_of_CI=diff)
library(xtable)
xtable(df,caption="Confidence interval of coefficients",digits=10)



#####````
library(lightgbm)
library(tidyr)
library(dplyr)
library(tweedie)
library(HDtweedie)
rm(list = ls())
library(lightgbm)
library(tweedie)
library(ggplot2)


## TWEEDIE NEGATIVE LOG-LIKELIHOOD
dtweedie_nlogl <- function(phi, y, mu, power) {
  ## `y` and `mu` MUST BE OF THE SAME LENGTH
  ans <- -2 * mean(log(
    dtweedie(y = y, mu = mu, phi = phi,  power = power)
  ))
  if (is.infinite(ans)) {
    ans <- mean(
      tweedie.dev(y = y, mu = mu, power = power)
    )
  }
  attr(ans, "gradient") <- dtweedie.dldphi(
    y = y, mu = mu, phi = phi, power = power
  )
  return(ans)
}


## LOAD AUTO CLAIM DATA
data(AutoClaim, package = "cplm")

## REMOVE FEATURES
x <- AutoClaim[, -c(1:5, 10, 16, 29)]

## INTEGER CODING TABLE FOR CATEGORICAL FEATURES
coding <- sapply(x, function(x) {
  if(is.factor(x))
    data.frame(
      code = seq_along(levels(x)),
      level = levels(x)
    )
})
coding[sapply(coding, is.null)] <- NULL
print(coding)

## ENCODE CATEGORICAL FEATURES INTO INTEGERS
index <- sapply(x, is.factor)
x[index] <- lapply(x[index], as.integer)

dt <- list(data = as.matrix(x), label = AutoClaim$CLM_AMT5 / 1000)

X=as.matrix(x)
Y=AutoClaim$CLM_AMT5 / 1000
n=length(Y)

alpha=.05
eps=.001
# percentage of 0 is 0.6323
sum(Y==0)/length(Y)
hist(Y)
simulations=100

pow <- seq(1.2, 1.5, by = 0.025)
nrounds <- 2000L
num_leaves <- 10L
learning_rate <- 0.005
nfold <- 5L

asymetric_pearson_residuals_rates=rep(0,simulations)
pearson_residuals_rates=rep(0,simulations)
pearson_residuals_rates_n3=rep(0,simulations)

deviance_residuals_rates=rep(0,simulations)
anscombe_residuals_rates=rep(0,simulations)


asymetric_pearson_residuals_gapOFci=rep(0,simulations)

pearson_residuals_gapOFci_n3=rep(0,simulations)

pearson_residuals_gapOFci=rep(0,simulations)##for test
deviance_residuals_gapOFci=rep(0,simulations)
anscombe_residuals_gapOFci=rep(0,simulations)

##we evaluate the best model gradually.

best_model_evaluation=function(
    dtrain=train_data ,
    dt = list(data = as.matrix(Xsamp), label = Ysamp
    ),
    nrounds=100){
  
  loglik <- rep(NA, length(pow))
  model <- vector("list", length(pow))
  set.seed(2024)
  for (k in seq_along(pow)) {
    pp <- pow[k]
    params <- list(
      num_leaves = num_leaves,
      learning_rate = learning_rate,
      objective = "tweedie",
      tweedie_variance_power = pp,
      verbose = 0L,
      use_missing = TRUE
    )
    
    fit <- lgb.cv(
      params = params,
      data = dtrain,
      nfold = nfold,
      nrounds = nrounds
    )
    best_iter <- fit$best_iter
    # print(best_iter <- fit$best_iter)
    # print(fit$best_score)
    # 
    # ## CROSS-VALIDATION ERROR
    # metric <- unlist(fit$record_evals$valid$tweedie$eval)
    # metric_se <- unlist(fit$record_evals$valid$tweedie$eval_err)
    # 
    # ## U-CURVE FOR CROSS-VALIDATION ERROR
    # data.frame(
    #   metric = metric,
    #   lower = metric - 2 * metric_se,
    #   upper = metric + 2 * metric_se,
    #   boost = rep(1:nrounds, 3)
    # ) |> ggplot() + 
    #   geom_line(aes(x = boost, y = metric), linewidth = 2) +
    #   geom_ribbon(aes(x = boost, ymin = lower, ymax = upper),
    #               alpha = 0.5, fill = "lightblue") +
    #   labs(x = "Number of boosting iterations",
    #        y = "Cross-validation error")
    
    ## REFIT MODEL WITH FULL DATA USING BEST NROUNDS
    fit <- lgb.train(
      params = params,
      data = dtrain,
      nrounds = best_iter
    )
    model[[k]] <- fit
    preds <- predict(fit, dt$data)
    
    ##===== ESTIMATION OF DISPERSION =====##
    
    # ## PEARSON ESTIMATION
    # phi_moment <- mean((dt$label - preds)^2 / (preds^pp))
    # phi <- phi_moment
    
    ## DEVIANCE ESTIMATION
    phi_saddle <- mean(
      tweedie.dev(y = dt$label, mu = preds, power = pp)
    )
    # phi <- phi_saddle
    
    ## MAXIMUM LIKELIHOOD ESTIMATION
    lower_limit <- min(0.001, 0.5 * phi_saddle)
    upper_limit <- 10 * phi_saddle
    ans <- optimize(
      f = dtweedie_nlogl, maximum = FALSE,
      interval = c(lower_limit, upper_limit),
      power = pp, mu = preds, y = dt$label
    )
    phi <- ans$minimum
    
    print(sprintf(
      "Power parameter = %.3f, best iteration = %4d, dispersion parameter = %.4f",
      pp, best_iter, phi
    ))
    loglik[k] <- mean(log(
      dtweedie(y = dt$label, mu = preds, phi = phi, power = pp)
    ))
  }
  
  data.frame(loglik = loglik, power = pow) |>
    ggplot() + geom_line(aes(x = power, y = loglik)) +
    labs(x = "Power parameter", y = "Log-likelihood")
  
  kmax <- which.max(loglik)
  print(pp <- pow[kmax])
  fit <- model[[kmax]]
  
  return(list(pp,phi,fit))
}



for(j in 1:simulations){
  
  n1= 4000
  n2= 4000
  n3= n-n1-n2
  #samp=sample(1:n1+n2,n1,rep=F)
  ##
  samp=sample(1:(n1+n2),n1,rep=F)
  Xsamp <- X[samp,]
  Ysamp <- Y[samp]
  
  #####
  
  
  train_data <- lgb.Dataset(data = as.matrix(Xsamp),label = Ysamp,
                            categorical_feature = colnames(x)[which(index)])
  
  test_data <- lgb.Dataset(data = as.matrix(X[8001:n,]),label = Y[8001:n],
                           categorical_feature = colnames(x)[which(index)])
  
  train_data_bigger <- lgb.Dataset(data = as.matrix(X[1:8000,]),label = Y[1:8000],
                                   categorical_feature = colnames(x)[which(index)])
  
  
  # Train LightGBM model
  modd <- best_model_evaluation(
    dtrain=train_data ,
    dt = list(data = as.matrix(Xsamp), label = Ysamp
    ),
    nrounds=50*j)
  
  model= modd[[3]]
  power_parameter_calculated=modd[[1]]
  phi_calculated=modd[[2]]
  
  
  predictions <- predict(model, as.matrix(Xsamp))%>%
    data.frame()
  
  
  #cvfit_pblpr=cv.glmnet(Xsamp[,-1],Ysamp,family=tweedie(link.power=0,var.power=1.7),alpha=0)##Fit a ridge regression 
  ##?tweedie #a value of q=link power=0 is log(\mu)=x^{'}beta
  Ypred_n1=predict(model,Xsamp
                   #, s=cvfit_pblpr$lambda.min,type = "response" 
  )  
  n2_set=setdiff(1:(n1+n2),samp) ## the n2 set
  Ypred_n2=predict(model,X[n2_set,]
                   #, s=cvfit_pblpr$lambda.min,type = "response" 
  )  
  
  
  
  
  n3_set=c(8001:n)
  Ypred_n3=predict(model,X[n3_set,]
                   #, s=cvfit_pblpr$lambda.min,type = "response" 
  )  
  
  #Obtain Pearson, Deviance, and Anscombe residuals
  pearson_residuals <- (Y[n2_set]-Ypred_n2)/sqrt(phi_calculated*(Ypred_n2)^power_parameter_calculated)
  deviance_residuals <- 2*(((Y[n2_set])^(2-power_parameter_calculated)/((1-power_parameter_calculated)*(2-power_parameter_calculated)))-
                             (((Y[n2_set])*(Ypred_n2)^(1-power_parameter_calculated))/(1-power_parameter_calculated))+
                             (((Ypred_n2)^(2-power_parameter_calculated))/(2-power_parameter_calculated))  )
  
  #2*(Y[n2_set]*log((Y[n2_set]+eps)/(Ypred_n2+eps))-(Y[n2_set]-Ypred_n2)) #just to avoid NAN
  anscombe_residuals <- sqrt(abs(Y[n2_set]-Ypred_n2))*sign(Y[n2_set]-Ypred_n2)
  
  q=ceiling((1-alpha)*(n2+1))/(n2)
  
  asymmetric_upper_pearson_residuals_quantile=quantile(pearson_residuals,ceiling((1-(.05/2))*(n2+1))/(n2))
  asymmetric_lower_pearson_residuals_quantile=quantile(pearson_residuals,floor(((.05/2))*(n2+1))/(n2))
  
  
  # asymmetric_upper_pearson_residuals_quantile=quantile(pearson_residuals,.975)
  # asymmetric_lower_pearson_residuals_quantile=quantile(pearson_residuals,.025)
  # 
  pearson_residuals_quantile=c(quantile(abs(pearson_residuals),q))
  deviance_residuals_quantile=quantile(deviance_residuals,q)
  anscombe_residuals_quantile=quantile(anscombe_residuals,q)
  
  asymetric_pearson_residuals_CI=cbind((Ypred_n1/sqrt(phi_calculated*Ypred_n1^power_parameter_calculated)) +asymmetric_lower_pearson_residuals_quantile,
                                       (Ypred_n1/sqrt(phi_calculated*Ypred_n1^power_parameter_calculated)) +asymmetric_upper_pearson_residuals_quantile)
  
  pearson_residuals_CI=cbind((Ypred_n1/sqrt(phi_calculated*Ypred_n1^power_parameter_calculated)) -pearson_residuals_quantile,
                             (Ypred_n1/sqrt(phi_calculated*Ypred_n1^power_parameter_calculated)) +pearson_residuals_quantile)
  
  ##a completely new test set
  pearson_residuals_CI_n3=cbind((Ypred_n3/sqrt(phi_calculated*Ypred_n3^power_parameter_calculated)) -pearson_residuals_quantile,
                                (Ypred_n3/sqrt(phi_calculated*Ypred_n3^power_parameter_calculated)) +pearson_residuals_quantile)
  
  
  
  #pearson_residuals_CI=cbind(Ypred_n1 -pearson_residuals_quantile,Ypred_n1 +pearson_residuals_quantile)
  deviance_residuals_CI=cbind(Ypred_n1 -deviance_residuals_quantile,Ypred_n1 +deviance_residuals_quantile)
  anscombe_residuals_CI=cbind(Ypred_n1 -anscombe_residuals,Ypred_n1 + anscombe_residuals)
  
  colnames(asymetric_pearson_residuals_CI)=c("L","U")
  
  colnames(pearson_residuals_CI)=c("L","U")
  colnames(pearson_residuals_CI_n3)=c("L","U")
  
  colnames(deviance_residuals_CI)=c("L","U")
  colnames(anscombe_residuals_CI)=c("L","U")
  
  asymetric_pearson_residuals_rates[j]=asymetric_pearson_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%
    rowwise()%>%
    mutate(result=if(((y/sqrt(phi_calculated*y^power_parameter_calculated +eps))>L)&&((y/sqrt(phi_calculated*y^power_parameter_calculated +eps))<U)){1}else{0})%>%dplyr::select(result)%>%sum()/n1
  
  
  pearson_residuals_rates[j]= pearson_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%
    rowwise()%>%
    mutate(result=if(((y/sqrt(phi_calculated*y^power_parameter_calculated +eps))>L)&&((y/sqrt(phi_calculated*y^power_parameter_calculated +eps))<U)){1}else{0})%>%dplyr::select(result)%>%sum()/n1
  
  
  pearson_residuals_rates_n3[j]= pearson_residuals_CI_n3%>%as.data.frame()%>%mutate(y=Y[8001:n])%>%
    rowwise()%>%
    mutate(result=if(((y/sqrt(phi_calculated*y^power_parameter_calculated +eps))>L)&&((y/sqrt(phi_calculated*y^power_parameter_calculated +eps))<U)){1}else{0})%>%dplyr::select(result)%>%sum()/n3
  
  # pearson_residuals_rates[j]= pearson_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
  #   mutate(result=if((y>L)&&(y<U)){1}else{0})%>%dplyr::select(result)%>%sum()/n1
  # 
  
  asymetric_pearson_residuals_gapOFci[j]=asymetric_pearson_residuals_CI%>%
    as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
    mutate(diff=U-L)%>%dplyr::select(diff)%>%sum()/n1
  
  pearson_residuals_gapOFci[j]=pearson_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
    mutate(diff=U-L)%>%dplyr::select(diff)%>%sum()/n1
  
  pearson_residuals_gapOFci_n3[j]=pearson_residuals_CI_n3%>%as.data.frame()%>%mutate(y=Y[8001:n])%>%rowwise()%>%
    mutate(diff=U-L)%>%dplyr::select(diff)%>%sum()/n3
  
  deviance_residuals_rates[j]= deviance_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
    mutate(result=if((y>L)&&(y<U)){1}else{0})%>%dplyr::select(result)%>%sum()/n1
  
  deviance_residuals_gapOFci[j]=deviance_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
    mutate(diff=U-L)%>%dplyr::select(diff)%>%sum()/n1
  
  anscombe_residuals_rates[j]= anscombe_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
    mutate(result=if((y>L)&&(y<U)){1}else{0})%>%dplyr::select(result)%>%sum()/n1
  
  anscombe_residuals_gapOFci[j]= anscombe_residuals_CI%>%as.data.frame()%>%mutate(y=Y[samp])%>%rowwise()%>%
    mutate(diff=U-L)%>%dplyr::select(diff)%>%sum()/n1
  
  ##
  
  print(j)
}

asymetric_pearson_residuals_rates

pearson_residuals_rates
# deviance_residuals_rates
# anscombe_residuals_rates
asymetric_pearson_residuals_gapOFci
pearson_residuals_gapOFci
pearson_residuals_gapOFci_n3

# deviance_residuals_gapOFci
# anscombe_residuals_gapOFci


plot(1:simulations,asymetric_pearson_residuals_gapOFci,ty="l",xlab="boosting rounds*50",
     ylab="average_gap_of_CI_residual_scale",
     col = "red",ylim=c(3.5,5),
     main="red:asymetric,blue:symetric pearson residual")
lines(1:simulations,pearson_residuals_gapOFci,col="blue")


plot(1:simulations,pearson_residuals_gapOFci,ty="l",xlab="boosting rounds")




plot(1:simulations,asymetric_pearson_residuals_rates,ty="l",xlab="boosting rounds*50",ylab="rates",
     col = "red",ylim=c(.9,1.1),
     main="red:asymetric,blue:symetric pearson residual")

lines(1:simulations,pearson_residuals_rates,col="blue")

ggplot(data=data.frame(claim_amount=Y), aes(claim_amount)) + 
  geom_histogram(binwidth=1,color="darkblue", fill="lightblue")+
  labs(title = "",
       x = "claim amount",
       y = "Frequency") +
  theme_minimal() 

library(reshape)
meltData <- melt(data.frame(Y,x))
boxplot(data=meltData, value~variable)


ggplot(data=data.frame(Y,x), aes(x=category, y=Y, fill=category)) +
  geom_bar(stat="identity") +
  scale_fill_viridis_d() +
  facet_grid(. ~ variable)



x

cat_cols <- sapply(x, is.factor)

for (col in names(cat_cols)[cat_cols]) {
  barplot(table(x[[col]]), main = paste("Barplot of", col), xlab = col, ylab = "Frequency")
}


# Plot each categorical column

library(gridExtra)
library(ggplot2)
p=list()
# Plot each categorical column
for (colm in names(cat_cols)[cat_cols]) {
  p[[colm]] <- ggplot(x, aes_string(x = colm,fill=colm)) +
    geom_bar() +
    ggtitle(paste("Barplot of", colm)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.03, hjust=1))+
    xlab(colm) +
    ylab("Frequency")+
    theme_minimal()+
    coord_flip()
  
  #print(p)
}

do.call(grid.arrange,p)

