#import "@local/templates:1.0.0" : *
#import "/utils.typ" : *

= Consequences of violating Gauss Markov assumptions

== Weak exogeneity

Some terms:
#defenv[
  _Exogeneity_ is the assumption that measurement errors are uncorrelated with the covariate $x$. In other words, $Cov(x,epsilon) = 0.$ We often write $EE[epsilon|x] = 0$
]
#defenv[
  _Endogeneity_ refers to the errors in measurement of $Y$ being correlated with measurements of $x$. 
]

In this section we consider the single-variable model in #eref(<eq:reg-model>).
We have assumed that there are no errors in our observations of the covariate $x$, but it's possible there actually are errors. If we naiively use the OLS estimator for $hat(beta)$ how does the estimate relate to the true value?
Violation of weak exogeneity is sometimes referred to as errors-in-variables.
In OLS regression it leads to _attenuation bias_, where $hat(beta)$ becomes biased towards 0. 

First, let's arrive at the effect using intuition. The OLS estimator for  $beta$ is  $S_(x y) slash S_x$.
If there are no errors in the measurement of $x$ then the only thing obscuring our ability to see the true covariance between $x$ and $y$ are the errors in $y$ that we assume in OLS regression.
Adding errors to $x$ has the effect of reducing the observed covariance between $x$ and $y$, so we should expect that if we use the OLS estimator in this case, our estimate would be biased towards zero than the same estimator when used in the case when there are no errors in $x$.

Now some maths. Denote the true value of $x$ by $x^*$ and let the error in measurements of $x_*$ be $eta$.
The model is given by taking #eref(<eq:reg-model>) and replacing $x -> x_*$,
$
  y = beta_0 + x_* beta_1 + epsilon med ,
$
but since we can only measure $x = x_* + eta$ we have, in practice,
$
  y &= beta_0 + (x - eta) beta_1 + epsilon med #no-num \
  &= beta_0 + x beta_1 + (epsilon - beta_1 eta) \ 
  &equiv beta_0 + x beta_1 + tilde(epsilon) med,
$
where $tilde(epsilon) = epsilon - beta_1 eta$ is identified as the "new" residual, which is now correlated with $x$. 
The OLS estimator for $beta_1$ then converges to
$
  hat(beta)_1 = (sum_i (x_i - xbar)(y_i - ybar)) / (sum_i (x_i - xbar)^2) arrow Cov(x,y) / Var(x) = (sigma_(x_*)^2) / (sigma_(x_*)^2 + sigma_eta^2) thin beta_1,
$
which is less than or equal to $beta_1$. This effect is called _attenuation damping_. 
In deriving this expression I used the fact that we condition on the observed $x$ but are uncertain about the true value $x_*$ and the noise $eta$.
We have,
$
  Cov(x,y) 
  &= Cov(x_* + eta, beta_0 + beta_1 x_* + epsilon) #no-num \
  &= Cov(x_*, beta_1, x_*) + Cov(eta, beta_1 x_*) + Cov(eta, epsilon) #no-num \
  &= beta_1 Var(x_*) + 0 + 0 equiv beta_1 sigma^2_(x*) #no-num
$


_Note:_ The first time I encountered this I was very confused about the meaning of $Cov(x,y)$ because I had the perspective that $x$ is not a random variable and $y$ is. 