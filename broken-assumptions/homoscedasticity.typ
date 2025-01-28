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
  W equiv "diag"(1 slash sigma_i) in RR^(n times n) med .
$
Then the weighted least squares estimator can be written

#bluebox[
  _Weighted least squares (WLS) estimator:_
  $
    hat(beta) = (X^T W^T W X)^(-1) X^T W^T (W y) med ,
  $
  where 
  $
    W = "diag"(sigma_i^(-1)) med .
  $
]

Of course, we rarely know $sigma^2_i$ precisely, so this also needs to be estimated. 
If we have reason to believe that the errors are dependent on $X$ one can fit another model to estimate $sigma^2(X)$, e.g. using another linear model, a decision tree, or a neural network.
This might seem to violate the assumption of weak exogeneity, but this is not necessarily the case. 
You could have $sigma^2 = sigma^2(X)$ without violating $EE[epsilon|X] = 0$.

#red[
  *Need to add:*
  + Diagnosing heteroscedasticity (look at residuals vs predictions/vs each feature)
  + Analytic example where $sigma^2(x) = beta x + alpha$
  + Example where you fit $sigma^2$ with a model
]
