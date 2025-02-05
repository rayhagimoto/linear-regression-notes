#import "@local/templates:1.0.0" : * 

== Homoscedasticity

- Significance tests become unreliable
- OLS estimator is no longer the BLUE. The intuitive explanation for this is that it weights more noisy terms the same as less noisy terms. Therefore the strategy should be to downweight the importance of the more noisy samples compared to other samples. This line of reasoning leads us to weighted least squares regression. 

Suppose the variance of the residuals is not constant. Assuming there is still no autocorrelation of errors we can write the general case as 
$Var(epsilon_i) = sigma^2_i$. 
Then if we take the linear model
$
  y_i = beta X_i + epsilon_i
$
(where $X_i in RR^(p+1)$) and normalise both sides by $1 / sigma_i$ we can define
$
  y_i / sigma_i = beta (X_i / sigma_i) + epsilon_i / sigma_i med .
$
Since $Var(epsilon_i slash sigma_i) = 1$ the residuals have constant variance. 
Moreover, the model is still linear, and $beta$ is unchanged. 
So if we do OLS estimation on the augmented data $(X_i slash sigma_i, y_i slash sigma_i)$ we no longer have problems with heteroscedasticity and can use the usual OLS estimator methods.
In matrix form, let
$
  W equiv "diag"(1 slash sigma_i^2) in RR^(n times n) med .
$
Then the weighted least squares estimator can be written

#bluebox[
  _Weighted least squares (WLS) estimator:_
  $
    hat(beta) = (X^T W X)^(-1) X^T W y med ,
  $
  where 
  $
    W = "diag"(sigma_i^(-2)) med .
  $
]

Of course, we rarely know $sigma^2_i$ precisely, so this also needs to be estimated. 
If we have reason to believe that the errors are dependent on $X$ one can fit another model to estimate $sigma^2(X)$, e.g. using another linear model, a decision tree, or a neural network.
This might seem to violate the assumption of weak exogeneity, but this is not necessarily the case. 
You could have $sigma^2 = sigma^2(X)$ without violating $EE[epsilon|X] = 0$.

=== WLS Example

Consider the following setup:
$
  x ~ U(1, 5) quad 
  epsilon|x ~ "N"(0, 1 / 4 x^2) \
  y|x,epsilon = beta x + epsilon med .
$
That is, the error term has a variance that is explicitly dependent on $x$. 

As an experiment, I generated $N_"samples" = 100$ data using this setup and estimated $beta$ using OLS and WLS regression. 
For WLS regression I used weights $W = "diag"(1 slash x_i^2)$.
I repeated this process $N_"sims" = 10,000$ times, storing the estimated $hat(beta)_"OLS"$ and $hat(beta)_"WLS"$ in each run. 
@fig:ols-wls-beta-hist[Figure] shows the distribution of the estimates from OLS and WLS regression as blue and orange histograms, respectively. 
The true $beta$ is plotted as a vertical dashed red line. 
The OLS histogram is slightly broader than the WLS histogram, indicating that the WLS estimator is more efficient.
Furthermore, both histograms have their mode close to the true value.

#figure(
  image("figs/homoscedasticity/beta-hist.svg", width:68%),
  caption: [
    Distribution of $hat(beta)$ obtained by performing OLS (blue) and WLS (orange) regression. 
    The WLS histogram is slightly narrower, indicating that it is more efficient than the OLS estimator on the same size sample data.
  ]
) <fig:ols-wls-beta-hist>

@fig:ols-wls-pred[Figure] shows an example of the predicted trend lines from OLS and WLS regression on a set of 100 points. 
In this particular example I cherry-picked it so that the OLS estimator actually does a bit better because although the OLS estimator has higher variance it can sometimes beat the WLS estimator _by pure chance_. 
The idea is that the WLS estimator is generally more reliable. 
However, it is only more reliable if the weights we have chosen are good. 
In this example we chose ideal weights that use the known $sigma^2(x) = 0.5 thin x^2$ relationship.

#figure(
  image("figs/homoscedasticity/trendlines.svg", width:68%),
  caption: [
    Example of the 
  ]
) <fig:ols-wls-pred>





#red[
  *Need to add:*
  + Show that the $hat(t)$ statistic is not $t$-distributed when conditional homoscedasticity is violated.
]
