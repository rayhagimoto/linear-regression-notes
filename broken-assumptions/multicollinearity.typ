#import "@local/templates:1.0.0" : *
#import "/utils.typ" : *

== Multicollinearity

#red[
  _Note:_ My understanding of the maths of this section are fuzzy.
  I was following #link("https://www.statlect.com/glossary/variance-inflation-factor")[#underline[these notes]] for much of this section, but it seems like they do not include the constant $beta_0$ term in their regression model. 
  Of course, this can be achieved by standardising the target and the covariates e.g. $y -> y - ybar$, $x_i -> x_i - xbar_i$.
  Wikipedia suggests that the formulas here still stand up when you include $beta_0$. 
  In particular, in the expression below make the replacement $(X^T X)^(-1) -> [sum (x - xbar)^2]^(-1)$ to get Wikipedia's expression.
]

This topic is a favourite in quant finance interviews.  
Multicollinearity means that two or more variables are linearly dependent (in practice, approximately linearly dependent) so that the covariance matrix $X^T X$ becomes (approximately) singular and some of the regression parameter estimates are undefined (blow up).

The most significant consequence of multicollinearity is _variance inflation_. 
The basic idea is that since the regression model is effectively trying to find how the target changes with each covariate while holding all but one constant, it isn't able to pick up on degeneracies. 
We can illustrate this with a simple example. 
Suppose we have just two covariates $x_1$ and $x_2$, and $x_1 = x_2$. 
Then our regression model is 
$
  y = beta_0 + beta_1 x_1 + beta_2 x_2 + epsilon med .
$
which can be rewritten as 
$
  y = beta_0 + (beta_1 + beta_2) x_1 + epsilon med . #<eq:multicollinearity-1>
$
Now notice that an increase in $beta_1$ can be compensated by a decrease in $beta_2$ and the equation remains unchanged. 
Our regression model estimates $beta_1$ and $beta_2$ separately, but they are not uniquely determined, even in the limit of infinite data. 
So the coefficients $beta_1$ and $beta_2$ will have _high variance_.
Yet another way to think about #eref(<eq:multicollinearity-1>) is that the loss function $L = sum_(i=1)^N (y_i - sum_(j=1)^p x_(i j) beta_j)$ will have a flat ridge minima.

_Variance inflation factor_

Start with #eref(<eq:multiple-lr-var>) (repeated below for convenience)
$
  Var(hat(beta)) = sigma^2 (X^T X)^(-1) #no-num 
$
Then by #eref(<eq:multiple-lr-var>) we have 
$
  Var(hat(beta)_k) = sigma^2 (X_(dot,k)^T X_(dot,k))^(-1) 1 / (1 - R^2_k) med,
$
where 
$X_(dot,k)$ is the $k$-th column of $X$ and $R^2_k$ is the R-squared obtained by regressing the $k$-th regressor on all the other regressors.
The rightmost factor is known as the _variance inflation factor (VIF)_ and it's important enough that I'll enshrine it in a blue box. 
#bluebox[
  $
    "VIF"_k = 1 / (1 - R^2_k)
  $
]
Clearly, if $R^2_k approx 1$, then the variance will be large. 
The important thing is that the VIF characterises how much of $x_k$ can be explained by the other variables. 
If most of it can, then $x_k$ may be worth removing.



