library(glmnet)
library(MASS)
library(parallel)
library(foreach)
library(mpath)

n = 2000
p = 40
Cov = diag(p)
#mu <- rep(1,p)##Made it otherwise mean will be negative
 beta <- rep(0,(p+1))
 beta[1:11] <- c(.5,c(1:10)/15)#c(1:10)/15
 
m <- length(beta)
X <- MASS::mvrnorm(n,mu=runif(p,-2,2),Cov)
Xb <-  cbind(1,X) %*% beta
mu=exp(Xb)
pr <- 1/(1+exp(-Xb))
theta=4.5

y <- sapply(1:n, function(i) rpois(n=1, lambda = exp(Xb)[i]))

# Load necessary libraries
library(glmnet) # For Lasso penalization
library(boot)   # For bootstrapping

# Define the bootstrap function
lasso_glm_bootstrap <- function(data, formula, alpha = 1, n_boot = 1000, family = "gaussian", lambda = NULL) {
  
  # Step 1: Fit Lasso GLM and determine optimal lambda
  x <- model.matrix(formula, data)#[, -1]  # Design matrix (excluding intercept)
  y <- model.response(model.frame(formula, data))
  
  # Cross-validation to find optimal lambda if not provided
  if (is.null(lambda)) {
    cv_fit <- cv.glmnet(x, y, alpha = alpha, family = family)
    lambda_opt <- cv_fit$lambda.min
  } else {
    lambda_opt <- lambda
  }

  lasso_fit <- glmnet(x, y, alpha = alpha, lambda = lambda_opt, family = family)
  beta_hat <- coef(lasso_fit)#coef(lasso_fit)[-1]  # Exclude intercept

  # Step 2: Calculate Pearson residuals
  y_hat <- predict(lasso_fit, x, type = "response")
  var_hat <- family_variance(y_hat, family)
  phi_hat <- mean((y - y_hat)^2 / var_hat)  # Estimate dispersion
  pearson_residuals <- (y - y_hat) / sqrt(var_hat)

  # Bootstrap procedure
  bootstrap_coefs <- replicate(n_boot, {
    # Step 3: Bootstrap residuals and generate new response values
    resampled_residuals <- sample(pearson_residuals, size = length(y), replace = TRUE)
    y_star <- sqrt(var_hat) * resampled_residuals + y_hat
    
    y_star <- ifelse(y_star>=0,y_star,0)
    
    # Step 4: Refit Lasso GLM
    lasso_fit_star <- glmnet(x, y_star, alpha = alpha, lambda = lambda_opt, family = family)
    coef(lasso_fit_star)[-1]  # Exclude intercept
  })

  # Step 5: Calculate confidence intervals
  ci_lower <- apply(bootstrap_coefs, 1, quantile, probs = 0.025)
  ci_upper <- apply(bootstrap_coefs, 1, quantile, probs = 0.975)

  # Combine results
  list(
    bootstrap_coefs,
    beta_hat = beta_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}

# Helper function to calculate variance for specific GLM families
family_variance <- function(mu, family, theta = NULL) {
  if (family == "gaussian") {
    return(rep(1, length(mu)))
  } else if (family == "poisson") {
    return(mu)
  } else if (family == "binomial") {
    return(mu * (1 - mu))
  } else if (family == "negative_binomial") {
    if (is.null(theta)) {
      stop("For negative binomial family, the dispersion parameter 'theta' must be provided.")
    }
    return(mu + (mu^2 / theta))
  } else {
    stop("Unsupported family")
  }
 }
# Example usage
set.seed(13)
# example_data <- data.frame(
#   x1 = rnorm(100),
#   x2 = rnorm(100),
#   y = rbinom(100, 1, 0.5)
# )

example_data<- data.frame(X,y)

# Fit bootstrap Lasso GLM
results <- lasso_glm_bootstrap(
  data = example_data,
  #formula = y ~ x1 + x2,
  formula = y ~ .,
  family = "poisson",
  n_boot = 500
)

# Display results
#results

results_poisson=results$ci_upper-results$ci_lower
#results_poisson

# Fit the Poisson regression model
poisson_fit <- glm(y~X, family = "poisson", data = data)

# Get the confidence intervals
ci_poisson <- confint.default(poisson_fit, level = 0.95)

# Print the confidence intervals
confint_diff_poisson=ci_poisson[,2]-ci_poisson[,1]

library(glmnet)

# Example data
set.seed(42)
# x <- matrix(rnorm(100 * 5), 100, 5)  # 100 observations, 5 predictors
# y <- rpois(100, lambda = exp(1 + x[, 1]))  # Poisson response variable

# Fit LASSO with glmnet
fit <- glmnet(X, y, family = "poisson")

# Bootstrapping function
bootstrap_lasso <- function(data, indices) {
  x_boot <- data[indices, -1]
  y_boot <- data[indices, 1]
  fit <- glmnet(x_boot, y_boot, family = "poisson", lambda = cv_fit$lambda.min)
  as.numeric(coef(fit))  # Extract coefficients
}

# Perform cross-validation to determine best lambda
cv_fit <- cv.glmnet(x, y, family = "poisson")

# Combine response and predictors for bootstrapping
data <- data.frame(y, X)
boot_results <- boot::boot(data, statistic = bootstrap_lasso, R = 1000)

# Calculate 95% confidence intervals
ci_boot <- t(apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975)))

# Print confidence intervals
print(ci_boot)

df1=data.frame(diff_grand_pois[2:41,-1]|>rowMeans(),results_poisson[-1],(ci_boot[,2]-ci_boot[,1])[-1])

colnames(df1)=c("PLR","Algorithm 2", "Bootstrap Lasso")
library(xtable)
xtable(df1[1:10,],caption="Confidence interval comparizon for nonzero coefficients",digits =3)
