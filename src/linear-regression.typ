#import "@preview/physica:0.9.4" : *
#import "@local/templates:1.0.0" : *
#import "utils.typ" : *

= Introduction

The most general relationship between variables $x$ and $y$ is a statistical one. 
Every data point $(x,y)$ is generated by sampling from the joint distribution between $x$ and $y$, denoted by $p(x,y)$.
It is useful to write this relationship in terms of the distribution of $y$ conditioned on $x$, since often we care about predicting $y$ given observations of $x$. 
We therefore write

// #bluebox[
  $ (x,y) &~ p(y|x) thin p(x) med , #<eq:conditional-sample> $
// ]

where $p(x)$ is the marginal distribution of $x$.
In general we want to learn $p(y|x)$ from observed data $cal(D) equiv {(x_i, y_i) : i = 1, dots, N}$.
However, we are often limited to learning the conditional mean $EE[y|x]$ (as in the case of minimising an $L_2$ loss), or median (as in the case of minimising an $L_1$ loss).


== Linear model

The simplest model is a linear one that assumes $y$ depends linearly on the model parameters $beta$. One example, for the univariate case is
$ y = beta_0 + beta_1 x + epsilon med , $ <eq:reg-model>
where $epsilon$ is the residual, a random variable which captures the uncertainty in measurements of $y$. 
Another example is 

$ y = beta_0 + beta_1 x^2 + epsilon med $
The only difference is that $x$ has been replaced with $x^2$, which makes the model non-linear in $x$. However, since the model is still linear in $y$ and the model parameters $beta_0, beta_1$, this is still considered a linear model. *Linearity, in this context, means linear w.r.t $y$ and $beta$*. 

Without loss of generality we take #eref(<eq:reg-model>) to be our model.
Having chosen a model the next obvious question is how we fit the model parameters (in this case $beta_0$ and $beta_1$) given some data? 
A common approach is to do _ordinary least squares (OLS)_ regression, where one quantifies the performance of a set of parameters by the sum of squared differences between predictions and observed values, $L = sum_i [y_i - beta_0 - beta_1 x_i]^2$.

It is certainly reasonable to consider this loss function, but why not the sum of absolute values or sum of 4-th power residuals, or something else entirely? Does it even matter? It turns out it does matter. Since it matters, it's important to motivate this loss function to see what implicit assumptions are being made.
We will do this in the next section.

== Deriving the least-squares loss

We start by specifying the conditional distribution $f(y|x)$. Given $x$, the randomness in $y$ is sourced by the residual $epsilon$. If we assume $epsilon ~ N(0, sigma^2)$ then for a single observation we get the log-likelihood
$
  - 2 log cal(L)(y|x,beta) = 1 / (sigma^2) (y - beta_0 - beta_1 x)^2 med + "constants" .
$
When we have $N$ data this becomes
$
  - 2 log cal(L)(y|x,beta) = 1 / (sigma^2) sum_(i=1)^N (y_i - beta_0 - beta_1 x_i)^2 med .
$ <eq:least-square>
In #eref(<eq:least-square>) the $y$ appearing in the LHS of the equation represents the collection $(y_1, y_2, dots, y_N)$, and similarly for $x$. 

We can derive statistical estimators for $beta_0$ and $beta_1$ by finding the values where they maximise the likelihood, or equivalently, minimise the negative log-likelihood. 
From #eref(<eq:least-square>) we identify that the loss given by the negative loss-likelihood is precisely the least-squares loss function.

In other words, *assuming Gaussian residuals $epsilon_i ~ N(0, sigma^2)$ leads to the least squares loss*.


== Aside: what do we learn by minimising the least-squares loss function?

Suppose we have a flexible model $f(x;theta)$ with parameters $theta$ that we wish to train to predict $y$ given measurements of $x$.
If we identify the best-fit parameters $theta^*$ as those that minimise the squared difference between our model and true values on our dataset $cal(D)$.
That is, by minimising
$
  L[f] =  sum_i [y_i - f(x_i;theta)]^2 med ,
$
where I've written $L[f]$ to emphasize that the loss function can be interpreted as a functional in terms of the model $f$.
In the limit of infinite data the sum over $i$ becomes an average weighted by the joint distribution $f(x,y)$.
$
  L[f]
  &-> 
  integral.double [y - f(x;theta)]^2 p(x,y) dif x dif y #no-num\ 
  &=  
  integral.double [y - f(x;theta)]^2 p(y|x) p(x) dif x dif y #<eq:ls-general>
$
Now we minimise $L$ by varing $f$ (we could equivalently varying $theta$, but doing it this way is more clean, and more fun). 
Setting $delta L = 0$ yields
$
 (delta L) / (delta f) &= integral (y - f(x;theta)) thin p(y|x) thin p(x) dif y = (EE[y|x] - f(x;theta)) thin p(x) = 0 #no-num\
 &=> f(x;theta^*) = EE[y|x] med.
$
This is an important result. 
It tells us that even if we have _infinite data_ and an _arbitrarily flexible model_, *the best we can do by minimising a least-squares loss is to learn the conditional expectation of $y$ given $x$*.

#align(center)[
  _Note: A really nice reference for the content in this section is the  introduction of ref._ @bishop1994mixture.
]

== OLS estimators

The maximum likelihood estimator is obtained by taking the derivative of #eref(<eq:least-square>) w.r.t. to $beta_0$ and $beta_1$, setting their results equal to zero, and rearranging. 
The results are simply,

#bluebox[
$
  hat(beta)_1 &= (sum_i (x_i - xbar)(y_i - ybar)) / (sum_i (x_i - xbar^&)^2) = S_(x y)^2 / S_x^2 = rho_(x y) S_y / S_x, #<eq:beta-0-est> \
  
  hat(beta)_0 &= ybar - hat(beta)_1 xbar med, #<eq:beta-1-est> 
$
]
where I have introduced the estimators for standard error $S$ and correlation $rho$,
$
  S_(x y)^2 &equiv 1 / (N-1) sum_i (x_i - xbar)(y_i - ybar) \ 
  S_(x)^2 &equiv 1 / (N-1) sum_i (x_i - xbar)^2 \
  rho_(x y) &equiv S_(x y)^2  / (S_x S_y) med , \

$
and overlines denote sample means, e.g. $xbar = 1 / N sum_(i=1)^N x_i$.

== Properties of estimators 

_*Bias*_

The first property of the OLS estimators is that they are _unbiased_, when we condition on $x$. 
This can be shown with a straightforward calculation that I will carry out below.
Note that in the following all expectation values are _conditional on $x$_. Hence, when I write $EE(y)$ I really mean $EE[y|x]$ (the expectation of $y$ conditioned on $x$). 
This implies that, the expectation of any arbitrary function of $x = (x_1, x_2, dots, x_N)$ is itself when conditioned on $x$, e.g. $EE[norm(x)^2] = norm(x)^2$.

First, let's evaluate the expectation of $hat(beta)_1$. We have,
$
  EE hat(beta)_1 
  &= (sum_i (x_i - xbar) EE (y_i - ybar)) / (sum_i (x_i - xbar)^2) med ,
$
but,
$
  EE(y_i - ybar) = (beta_0 + beta_1 x_i - beta_0 - beta_1 xbar) = beta_1(x_i - xbar) med ,
$
and so,
$
  EE hat(beta)_1 = beta_1 med .
$
Thus, for $hat(beta)_0 = ybar - hat(beta)_1 xbar$ we have 
$
  EE hat(beta)_0 = EE ybar - xbar thin EE hat(beta)_1 = beta_0 + beta_1 xbar - xbar hat(beta)_1 = beta_0 med .
$

*_Variance_*\

I'll present the formulas for quick reference then derive the formula for $Var(hat(beta)_1)$.

#bluebox[
$
  Var(hat(beta)_1) = sigma^2 / (sum_i (x_i - xbar)^2) #<eq:slr-var-beta1> \
  Var(hat(beta)_0) = (sigma^2 sum_i x_i^2) / (n sum_i (x_i - xbar)^2) #<eq:slr-var-beta0>
$
]

// To derive the formula for $Var(hat(beta)_1)$, I actually find it easier to solve the multivariable case with $p$ regressors $x_1, x_2, ..., x_p$ (see @sec:multiple-regressors[sec.] for details) but I think it would be useful to include a brute force derivation. 
// I won't present one for $Var(hat(beta)_0)$ since it's more annoying and I genuinely believe the best way is the linear algebra method explored in @sec:multiple-regressors[sec.].

The variance of $hat(beta)_1$ may be written as 
$
  Var(hat(beta)_1) = [sum_i (x_i - xbar)^2]^(-2) Var[sum_i (x_i - xbar)(y_i - ybar) ] med ,
$
where I have used the fact that we are conditioning on $x$ to factor the denonimator out a-la the identity $Var(k Y) = k^2 Var(Y)$ for constant $k$ and random variable $Y$.
Now let's focus on the variance factor on the right. 
For convenience, introduce the notation $x^*_i equiv x_i - xbar$ and notice that $y_i - ybar = beta_1 x^*_i + epsilon_i - overline(epsilon)$, so that when we take the variance only the $epsilon$ terms will be relevant:
$
  Var[sum_i (x_i - xbar)(y_i - ybar) ] 
  &= Var(sum_i x^*_i [epsilon_i - overline(epsilon)]) \
  &= EE[ sum_(i,j) x^*_i x^*_j (epsilon_i - overline(epsilon))(epsilon_j - overline(epsilon)) ] \
  &= sum_(i,j) x^*_i x^*_j { EE[epsilon_i epsilon_j- epsilon_i overline(epsilon) - epsilon_j overline(epsilon) + overline(epsilon)^2]} med , #<eq:simple-var-deriv-1>
$
where in the second equality we used the fact that $EE[epsilon_i] = EE[overline(epsilon)] = 0$. 
Using the linearity of expectation to expand #eref(<eq:simple-var-deriv-1>)  yields
$
  sum_(i,j) x^*_i x^*_j { EE[epsilon_i epsilon_j] - EE[epsilon_i overline(epsilon)] - EE[epsilon_j overline(epsilon)] + EE[overline(epsilon)^2]} med .
$
However, since $EE[epsilon_i epsilon_j] = delta_(i,j) sigma^2$ and $overline(epsilon) = 1/ n sum_k epsilon_k $ all of these expectation values can be simplified. 
$
  &quad sum_(i,j) x^*_i x^*_j { sigma^2 delta_(i j) - sigma^2 / n - sigma^2 / n + sigma^2 / n} \
  &= sigma^2 sum_i (x^*_i)^2 - sigma^2 / n sum_(i,j) x^*_i x^*_j med .
$
In the rightmost term we recognise that $sum_(i,j) x^*_i x^*_j = (sum_i x^*_i)^2$, and furthermore, $sum_i x^*_i = 0$ by definition, so the term vanishes and we're left with 
$
  Var(hat(beta)_1) = (sigma^2 sum_i (x^*_i)^2) / ((sum_i (x^*_i)^2)^2) = sigma^2 / (sum_i (x_i - xbar)^2) med ,
$
which exactly matches #eref(<eq:slr-var-beta1>).


== Signifcance testing

The linear correlation between $x$ and $y$ is typically assessed via the $t$-statistic,\
#bluebox[
  $
    hat(t) = hat(beta)_1 / Var(hat(beta)_1) =  hat(beta)_1 / (hat(sigma) slash S_x ) med ,
  $
]
where $hat(sigma)$ is the estimator for the standard deviation of the residuals and is given by
$
  hat(sigma)^2 = 1 / (n-1) sum_i (y_i - hat(beta)_0 - hat(beta)_1 x_i )^2 med .
$
If the $epsilon_i$ are assumed to (1) be *Gaussian* with mean zero (2) have *no autocorrelation* (3) exhibit *weak exogeneity*, then the $t$-statistic follows a #box(fill:cyan)[$t$ distribution with $n - p$ degrees of freedom]. 
This can be used to calculate $p$-values for significance testing.
However, if any of these assumptions are violated you can't use the standard $p$-values. 
This happens basically all the time in financial time series analysis where, for example, you may model the next time step $y_t$ as a linear combination of lagged values.
This introduces autocorrelation in the residuals.
The Dickey-Fuller test takes this into account when calculating $p$-values for the presence of a unit root.

== $t$-test

Let's explore the $t$-statistics properties in more detail.
Under the model given in #eref(<eq:reg-model>) we are assuming that the observed values of $y$ fluctuate around the 'true trend' $beta x$ due to Gaussian noise
#footnote[#eref(<eq:reg-model>) doesn't require the noise to be Gaussian, but this is the most common assumption. ]
.
If there is no relationship between $x$ and $y$, then under #eref(<eq:reg-model>) this means $beta = 0$.
However, even if $beta = 0$ our OLS estimate of $beta$, $hat(beta)$ will generally be nonzero in a given sample.
The question is how do we determine if an obtained nonzero $hat(beta)$ is statistically significant?
Assume the null hypothesis $beta = 0$ and $epsilon ~ N(0, sigma^2)$.
Then,
$
  hat(beta) = (X^T X)^(-1)X^T y ~ 
$


== Multiple regressors 
<sec:multiple-regressors>

Suppose we want to use $p$ covariates to predict the variate $y$. We can write down this model as
$
  y_i = beta_0 + sum_(k=1)^p x_(k,i) beta_k  + epsilon_k "for" i = 1, 2, dots, N. #<eq:multiple-lr-model>
$
It's conventional to define the so-called _design matrix_ $X in RR^(N times (p + 1))$ as#footnote[In the literature I mostly see people saying the design matrix is $N times p$]
$
  X equiv mat(
    1, x_(1,1), x_(2,1), dots, x_(N,1);
    1, x_(1,2), x_(2,2), dots, x_(N,2);
    dots.v, dots.v, dots.v, dots.v, dots.v ;
    1, x_(1,N), x_(2,N), dots, x_(N,2);
  ) = mat(
    bar, bar, bar, bar, bar;
    1, x_1, x_2, dots, x_N;
    bar, bar, bar, bar, bar;
  )
$

In this case the estimator for the regression parameters is\
#bluebox[$
  hat(beta) = (X^T X)^(-1) X^T y med .
$ <eq:multiple-lr-est>]
Its variance-covariance matrix is #footnote[this can be easily derived by using the identity $Var(A x) = A Var(x) A^T$, and using the assumption of homoscedasticity to write $Var(epsilon) = sigma^2 bb(1)_(p+1)$.]\
#bluebox[
$
  Var(hat(beta)) = sigma^2 (X^T X)^(-1) med . #<eq:multiple-lr-var>
$ 
] 
#Eref(<eq:multiple-lr-var>) can be derived easily by using the identity $Var(A arrow(x)) = A Var(arrow(x)) thin A^T$ where $Var(arrow(x))$ denotes the variance-covariance matrix of the elements of $arrow(x)$: $[Var(arrow(x))]_(i j) = Cov(x_i, x_j)$.
Here's a quick derivation of the identity, and then how it can be applied to #eref(<eq:multiple-lr-est>).
I'm going to use the Einstein summation convention and denote the $i$-th element of $x$ by $x^i$ just for this derivation.
$
  [Var(A arrow(x))]_(i j) 
  &= EE[ (A_(i k) x^k)(A_(j l) x^l) ] - EE[A_(i k) x^k] EE[A_(j l) x^l] #no-num \
  &= A_(i k) EE [x^k x^l] A_(j l) - A_(i k) EE [x^k] EE[x^l] A_(j l) #no-num \
  &= A_(i k) EE [x^k x^l] (A^T)_(l j) - A_(i k) EE [x^k] EE[x^l] (A^T)_(l j) #no-num \
  &= A_(i k) (EE [x^k x^l] - EE [x^k] EE[x^l]) (A^T)_(l j) #no-num \
  &= A_(i k) Cov(x^k, x^l) (A^T)_(l j) #no-num \
  &= [A Var(x) A^T]_(i j) med . #no-num
$
Applying this identity to #eref(<eq:multiple-lr-est>) we get 
$
  Var(hat(beta)) 
  &= Var((X^T X)^(-1) X^T y) = (X^T X)^(-1) X^T Var(y) X (X^T X)^(-1) med . #<eq:multiple-lr-var-deriv-1>
$
The variance-covariance matrix $Var(y)$ can be written as $Var(X beta + epsilon) = Var(epsilon) = sigma^2 bb(1)_(p+1)$. Substituting this into #eref(<eq:multiple-lr-var-deriv-1>) immediately yields #eref(<eq:multiple-lr-var>).

== Assumptions

Up until this point I haven't gone into much detail about the assumptions we have made. 
I've just blitzed through the derivation of the estimators. Here we enumerate the assumptions and give them fancy names which I think were popularised by econometrics.
Memorising the assumptions is important because they are almost always violated. 
If they're violated a little then you're probably fine proceeding as usual, but when they're violated a lot we need to understand their implications.
This will help us recognize the fingerprints of each assumption violation. 

#bluebox[
  #underline[*Linear Regression Assumptions*]
  + *Linearity*: the model is linear in the variate and parameters.
  + *Random sampling*: the data $(x_i, y_i)$ are i.i.d., ensuring that the sample is representative of the population.
  + *No perfect multicollinearity*: the $p$ covariates are linearly independent. This imposes $"Rank"(X) = p$. 
  + *Weak exogeneity*: no information loss in $Y$ when conditioned on $X$, $EE[epsilon|X] = 0$.
  + *Homoscedasticity* the variance of errors is constant across all values of $X$.
  + *No autocorrelation*: errors are uncorrelated, $EE[epsilon_i, epsilon_j] = 0 space forall space i != j.$.
  + *Errors follow a distribution (optional)*: here we assumed Gaussian, but they could've been $t$-distributed
  + *Model specification*: basically "the model is correct". This assumption is often violated if for example there are additional features which have not been included in the model.
]

== Gauss Markov theorem

One of the most famous results is that the estimator #eref(<eq:multiple-lr-est>) is the _best linear unbiased estimator (BLUE)_, where best means lowest variance. The derivation is pretty straightforward so I will present it here.
First we define an arbitrary linear estimator of $beta$ as an estimator of the form 
$
  tilde(beta) = A y med ,
$
where $A in RR^((p + 1) times N)$. 
If it's unbiased then,
$
  EE(tilde(beta)) = beta med .\
$ <eq:gauss-markov-1>
On the other hand substituting #eref(<eq:multiple-lr-model>) for $y$ and using the fact that $EE(epsilon) = 0$ yields
$
  EE(tilde(beta)) = A EE(X beta + epsilon) = A X beta med .
$ <eq:gauss-markov-2>
Combining eqs. #ref(<eq:gauss-markov-1>) and @eq:gauss-markov-2 gives
$
  A X beta = beta => A X = bb(1)_(p + 1) med .
$ <eq:gauss-markov-3>
#Eref(<eq:gauss-markov-3>) motivates us to decompose $A$ as 
$
  A = (X^T X)^(-1) X^T + C med ,
$
where $C in RR^((p+1) times N)$ is in the null space of $X$, i.e., $C X = 0$.
The first term can't simply be $X^(-1)$ since we need a matrix with the shape $(p + 1) times N$, and $X^(-1)$ would be $N times (p + 1)$. 

The variance of $tilde(beta)$ can be written as 
$
  Var(tilde(beta)) &= Var(A(X beta + epsilon)) = Var(A epsilon) = A Var(epsilon) A^T \
  &= [(X^T X)^(-1) X^T + C] sigma^2 [(X^T X)^(-1) X^T + C]^T "(using Var"(epsilon) = sigma^2")" \
  &= sigma^2 { (X^T X)^(-1) + X^T X^(-1) X^T C^T + C X (X^T X)^(-1) + C C^T} #<eq:gauss-markov-4c> \
  &= Var(hat(beta)) + sigma^2 C C^T med . #<eq:gauss-markov-4d> 
$
To go from #eref(<eq:gauss-markov-4c>) to #eref(<eq:gauss-markov-4d>) I eliminated the cross terms in the middle via the fact that $C X = 0$ and rewrote the first term using #eref(<eq:multiple-lr-var>).
Since $C$ is a positive semi-definite matrix we have shown that $Var(tilde(beta))$ exceeds $Var(hat(beta))$ by a positive semi-definite matrix
#footnote[To verify that $C C^T$ is positive semi-definite simply write $v = C^T x$, then for any $x$ $abs(v)^2 = x^T C C^T x >= 0.$], $sigma^2 C C^T$.

== Nested model comparison using F-test

An alternative way of determining whether a particular covariate is significant is using an $F$-test. 