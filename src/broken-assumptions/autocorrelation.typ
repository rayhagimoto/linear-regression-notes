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
$ <eq:ar1-error-model>

where $abs(rho) < 1$ and $u_t$ is a white noise process ($u_t$ are iid samples).
Then one can show that $Var(epsilon) = sigma^2$ (which is the reason for the factor of $1-rho^2$ that appears in the distribution for $u_t$).

@fig:ar1-ols-errors-data-comparison[Figure] depicts 100 samples of $epsilon$ generated via white noise (WN) (left), AR(1) with $rho = 0.95$ (middle), and AR(1) with $rho = -0.95$ (right).
In all cases $sigma = 0.64$. 
Notice how the _overall_ spread appears to be the same in all cases.
However, the WN has no discernable pattern, the +AR(1) process seems to follow a trend, and the -AR(1) process seems to alternate.

$Var(epsilon^("AR1")) = Var(epsilon^("WN")) = sigma^2 = 0.01$.
#figure(
  image("figs/autocorrelation/example_data.png"),
) <fig:ar1-ols-errors-data-comparison>

The examples shown here are kind of extreme but the point I want to get across is that when there is positive autocorrelation you are more likely, by chance, to obtain data where there is apparently a 'statistically significant' trend. 
Positive autocorrelation leads to more false positives. 
Negative autocorrelation leads to more false negatives.

#figure(
  image("figs/autocorrelation/beta_ols_histogram.svg"),
) <fig:ar1-ols-beta-comparison>

