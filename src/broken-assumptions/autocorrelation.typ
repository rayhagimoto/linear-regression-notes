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
Then one can show that $sigma^2_epsilon equiv Var(epsilon) = sigma^2$ (which is the reason for the factor of $1-rho^2$ that appears in the distribution for $u_t$). Let's look at some examples of this and compare it to data that satisyf the OLS assumptions. 
@fig:ar1-ols-errors-data-comparison[Figure] depicts 50 samples of $epsilon$ generated via white noise (WN) (left), AR(1) with $rho = 0.95$ (middle), and AR(1) with $rho = -0.95$ (right).
In all cases $sigma^2_epsilon = 0.176$
Notice how the _overall_ spread appears to be the same in all cases.
However, the WN has no discernable pattern, the +AR(1) process seems to follow a trend, and the $-$AR(1) process seems to alternate.

The examples shown in @fig:ar1-ols-beta-comparison[figure] are somewhat extreme but the point I want to get across is that when there is positive autocorrelation you are more likely, by chance, to obtain data where there is apparently a 'statistically significant' trend. 
Positive autocorrelation leads to more false positives. 
Negative autocorrelation leads to more false negatives.

#figure(
  image("figs/autocorrelation/example_data.png"),
  caption:[]
) <fig:ar1-ols-errors-data-comparison>

The tendency for false positives when $rho > 0$ is illustrated by the fat tail of the subplots in the middle column of @fig:ar1-ols-beta-comparison[figure], whereas the tendency for false negatives when $rho < 0$ is illustrated by the fat tail of the subplots in the rightmost column. 
The red line overlaying the subplots is the standard normal pdf, which is what we assume to obtain $p$-values during significance testing.

#figure(
  image("figs/autocorrelation/t_stats.png"),
  caption:[]
) <fig:ar1-ols-beta-comparison>

If we suspect that our residuals exhibit autocorrelation we can do better by using _Newey-West standard errors_.
This basic idea here is to go back to #eref(<eq:multiple-lr-var-deriv-1>) and notice that the equation
$
  Var(hat(beta)) = (X^T X)^(-1) X^T Cov(epsilon) X (X^T X)^(-1)
$
still holds even when $Cov(epsilon)$ isn't equal to $sigma^2 bb(1)$.
This seems to imply we can have reliable significance testing if we can get a good estimate of $Cov(epsilon)$, but that's hard. 
As a next step, we might think to approximate the combination $X^T Cov(epsilon) X$ using $X^T hat(epsilon) hat(epsilon)^T X$ where $hat(epsilon)$ are the residual errors.
However in practice this doesn't work well because it weights higher order lags with the same importance as lower order lags (which have more effective samples, and are more reliable). 
Standard practice is to use Newey-West standard errors which are calculated as
$
  hat(V)_"NW" = (X^T X)^(-1) (hat(S)_0 + 2 sum_(ell = 1)^L w_ell hat(S)_ell) (X^T X)^(-1) \
  hat(S)_ell = sum_(t=ell+1)^T hat(epsilon)_t hat(epsilon)_(t - ell) x_t x_(t-ell)^T \
  w_ell = "weight kernel, often chosen to be" 1 - ell (L + 1) "for some max lag" L quad
$
Notice that the principle difference here is the introduction of a weight kernel $w_ell$ which discounts the importance of the higher order lags. 
In `statsmodels` one can use NW standard errors by doing `sm.OLS(y, X).fit(cov_type='HAC', cov_kwds=dict(maxlags=L))`.
HAC stands for "heteroskedasticity and autocorrelation consistent".
In the limit of infinite samples the $t$-stat computed with NW standard errors is $N(0,1)$ so the two-sided tail is $2(1 - Phi(abs(t)))$.
$t$-statistics computed with NW standard errors for various choices of $N_"samples"$ are depicted in @fig:ar1-nw-t-stat[figure].
We observe that while the NW errors are less reliable for the white noise samples than than the usual OLS standard errors for a small-ish sample size of $N_"samples" = 100$ (cf. top left of @fig:ar1-nw-t-stat[figure] vs. top left of @fig:ar1-ols-beta-comparison[figure]) the results are more reliable in the presence of autocorrelated errors. 
Moreover, as $N_"samples" -> infinity$ the distributions asymptote to a standard normal. The takeaway is that NW standard errors improve things when autocorrelation is present, but not when OLS assumptions are satisfied. Hence, it's important to diagnose the presence of autocorrelation first before using this technique.

#figure(
  image("figs/autocorrelation/t_stats_NW.png"),
  caption:[
    Distributions of $t$-statistics computed using Newey-West standard errors. 
    Critical values are obtained by assuming the distributions asymptote to a standard normal distribution whose pdf is overlaid in red. 
    Distributions of $t$-stats are derived using a suite of $N_"sims" =$ 10,000 simulations of length $N_"samples"$ (shown on the left hand side of the plot). 
    The plots confirm that the $t$-stat approaches a standard normal random variable. 
    For small samples, significance tests are more reliable than under OLS assumptions, but still not perfect.
  ]
) <fig:ar1-nw-t-stat>

=== Diagnosing autocorrelation

To do this we can use Durbin-Watson or Ljung-Box. 
The Durbin-Watson statistic is 
$
  "DW" = (sum_(t=2)^T (hat(epsilon)_t - hat(epsilon)_(t-1))^2) / (sum_(t=1)^T hat(epsilon)_t^2) med .
$ <eq:durbin-watson-stat>
This statistic compares the summed squared change series of $hat(epsilon)_t$ against the summed square of the series itself. 
It's basically $Var(Delta epsilon) slash Var(epsilon)$. 
Notice that if $hat(epsilon)_t$ and $hat(epsilon)_(t-1)$ are uncorrelated then this reduces to $(Var(epsilon_t) + Var(epsilon_(t-1))) slash Var(epsilon_t) = 2$ (under OLS assumptions).
In the presence of positive autocorrelation $"DW" < 2$ since residuals will tend to cluster more. 
Otherwise, if there is negative autocorrelation $"DW" > 2$ since residuals will swing more wildly.
It's worth going back to @fig:ar1-ols-errors-data-comparison[figure] to confirm this makes sense.

The Ljung-Box statistic tests for the significance of a set of autocorrelations up to a maximum lag $m$. 
Under the null hypothesis that none of the lags are nonzero, it is asymptotically a chi-square distribution with $m$ degrees of freedom.
$
  "Ljung-Box" = Q(m) = T(T+2) sum_(ell=1)^m hat(rho)^2_ell / (T - ell) med .
$ <eq:ljung-box-stat>

