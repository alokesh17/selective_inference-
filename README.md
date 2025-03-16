# selective_inference-
Interval Estimation of Coefficients in Penalized Regression Models of Insurance Data

The Tweedie exponential dispersion family is a popular choice among many to model insurance
losses that consist of zero-inflated semicontinuous data. In such data, it is often important to
obtain credibility (inference) of the most important features that describe the endogenous variables.
Post-selection inference is the standard procedure in statistics to obtain confidence intervals of
model parameters after performing a feature extraction procedure. For a linear model, the lasso
estimate often has non-negligible estimation bias for large coefficients corresponding to exogenous
variables. To have valid inference on those coefficients, it is necessary to correct the bias of the
lasso estimate. Traditional statistical methods, such as hypothesis testing or standard confidence
interval construction might lead to incorrect conclusions during post-selection, as they are generally
too optimistic. Here we discuss a few methodologies for constructing confidence intervals of the
coefficients after feature selection in the Generalized Linear Model (GLM) family with application
to insurance data.
