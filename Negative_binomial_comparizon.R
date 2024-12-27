#Data
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

y <- rnegbin(n=n, mu = exp(Xb), theta = theta)#rbinom(n = n, prob = pr, size = 1)
# loading <- MASS::mvrnorm(1,mu,Cov/2)
# loading[11:p] <- loading[11:p]/25



library(MASS)

library(MASS)

lasso_glm_bootstrap_nb <- function(data, formula, alpha = 1, n_boot = 1000) {
  # Step 1: Fit the initial Negative Binomial Model
  nb_fit <- glm.nb(formula, data = data)
  theta <- nb_fit$theta  # Dispersion parameter
  y <- model.response(model.frame(formula, data))  # Observed response
  x <- model.matrix(formula, data)  # Design matrix
  
  # Step 2: Determine fitted values and Pearson residuals
  y_hat <- predict(nb_fit, type = "response")
  var_hat <- y_hat + (y_hat^2 / theta)  # Variance for Negative Binomial
  pearson_residuals <- (y - y_hat) / sqrt(var_hat)
  
  # Step 3: Bootstrap procedure
  bootstrap_coefs <- replicate(n_boot, {
    # Generate new response values by bootstrapping residuals
    resampled_residuals <- sample(pearson_residuals, size = length(y), replace = TRUE)
    y_star <- sqrt(var_hat) * resampled_residuals + y_hat
    y_star <- pmax(y_star, 0)  # Ensure non-negative response
    y_star <- ceiling(y_star)  # Ensure integer response

    # Refit Negative Binomial Model
    data_star <- data
    data_star[[as.character(formula[[2]])]] <- y_star  # Replace response in the data
    nb_fit_star <- glm.nb(formula, data = data_star)  # Refit model
    coef(nb_fit_star)  # Extract coefficients
  })
  
  # Step 4: Calculate confidence intervals
  ci_lower <- apply(bootstrap_coefs, 1, quantile, probs = 0.025)
  ci_upper <- apply(bootstrap_coefs, 1, quantile, probs = 0.975)
  
  # Combine results
  list(
    beta_hat = coef(nb_fit),
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    theta = theta
  )
}

example_data<- data.frame(X,y)

# Fit bootstrap Lasso GLM
results <- lasso_glm_bootstrap_nb(
  data = example_data,
  #formula = y ~ x1 + x2,
  formula = y ~ .,
  #family = "negative_binomial",
  n_boot = 500
)

# Display results
results

library(MASS)

# Fit the initial Negative Binomial Model
nb_fit <- glm.nb(y ~ ., data = example_data)

# Calculate the 95% confidence intervals for the model coefficients
ci_nb <- confint(nb_fit, level = 0.95)

# Display the results
#ci

confint_diff_nb=ci_nb[,2]-ci_nb[,1]

# Define bootstrap function
bootstrap_nb <- function(data, indices) {
  boot_data <- data[indices, ]
  fit <- glm.nb(y ~ ., data = boot_data)
  coef(fit)  # Extract coefficients
}

# Perform bootstrapping
set.seed(123)
boot_results <- boot(data.frame(y,X), statistic = bootstrap_nb, R = 1000)

# Calculate 95% confidence intervals for coefficients
ci_boot <- t(apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975)))

df2=data.frame(diff_grand_nb[2:41,-1]|>rowMeans(),results_nb[-1],(ci_boot[,2]-ci_boot[,1])[-1])


colnames(df1)=c("PLR","Algorithm 2", "Bootstrap Lasso")
library(xtable)
xtable(df1[1:10,],caption="Confidence interval comparizon for nonzero coefficients",digits =3)


