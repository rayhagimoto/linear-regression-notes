#import "@local/templates:1.0.0" : * 

== Autocorrelation

In this section we will investigate the effect of autocorrelated errors on regression estimates. 
First, we consider the case that the residuals $epsilon_t$ follow an $"AR"(1)$ process (discussed in more detail in #sref(<sec:tsa-ar1>)).
We will see that the OLS estimator $hat(beta)_"OLS"$ remains _unbiased_, but its variance is not appropriately estimated by #eref(<eq:ols-sigmahat>). 
Therefore, significance testing via the $t$-statistic defined in #eref(<eq:t-statistic>) is unreliable. 
We will then demonstrate that this can be addressed using _robust standard errors_, e.g. _Newey-West standard errors_. 
Following our discussion of AR(1) autocorrelated errors we will see how autocorrelated errors can arise via ommitted variables and misspecified models. 
Through that exercise we will also demonstrate how autocorrelation can be diagnosed using techniques from time series analysis like the Durbin-Watson test (or Ljung-Box).

=== AR(1) errors

Suppose the data generating process is 

$
  y_t = alpha + beta x_t + epsilon_t \
  epsilon_t = rho epsilon_(t-1) + u_t, quad u_t ~ N(0, (1 - rho^2) thin sigma^2) \
$

where $abs(rho) < 1$ and $u_t$ is a white noise process ($u_t$ are iid samples).
Then one can show that $Var(epsilon) = sigma^2$ (which is the reason for the factor of $(1-rho^2)$ I wrote in the distribution for $u_t$).