#import "@local/templates:1.0.0" : *

= Time series analysis 

#red[
  #emph[
    This section is basically just scratch notes at this point.
    I haven't really figured out structure or presentation. 
    I'm just throwing down useful definitions for quick reference, and derivations that I found useful but weren't in my reference textbook.
  ]
]

The foundation of time series analysis is stationarity.
#definition[
  A time series ${r_t}$ is said to be weakly stationary if $EE[r_t] = mu$ is independent of $t$ and the autocovariance $EE[(r_t - mu)(r_(t-ell) - mu)] = gamma_ell$ is only a function of the lag. 
]

Stationarity can be checked using the Dickey-Fuller (or augmented Dickey-Fuller) test.
To explain this test let's look at a simple example of a linear time series.

== AR(2) model

The AR(2) model for $r_t$ is defined as
$
  r_t = phi_0 + phi_1 r_(t-1) + phi_2 r_(t-2) + epsilon_t med .
$
Assuming weak stationarity we can derive the mean $EE[r_t]$ as
$
  EE[r_t] 
  &= phi_0 + phi_1 EE[r_(t-1)] + phi_2 EE[r_(t-2)] + EE[epsilon_t] \
  => mu &=  phi_0 + phi_1 mu + phi_2 mu \
  => mu &= phi_0 / (1 - phi_1 - phi_2) med .
$
Next we can compute its autocovariance and autocorrelation functions. 
Rewriting the time series (in terms of deviations from the mean) we get
$
  r_t - mu = phi_1 (r_(t-1) - mu) + phi_2 (r_(t-2) - mu) + epsilon_t med .
$
Now we multiply this by the lagged values on both sides:
$
  (r_(t - ell) - mu)(r_t - mu) = phi_1 (r_(t - ell)- mu)(r_(t-1) - mu) + phi_2 (r_(t - ell) - mu)(r_(t-2) - mu) + epsilon_t
$
Taking the expectation,
$
  EE[(r_(t - ell) - mu)(r_t - mu)] 
  &= phi_1 EE[(r_(t - ell)- mu)(r_(t-1) - mu)] #no-num \
  &quad + phi_2 EE[(r_(t - ell) - mu)(r_(t-2) - mu)] + EE[epsilon_t] #no-num \
$
Apply stationarity to the expectation values,
$
  => gamma_ell &= phi_1 gamma_(ell - 1) + phi_2 gamma_(ell - 2) med . #<eq:autocovf-recursion>
$
Divide #eref(<eq:autocovf-recursion>) by $sqrt(Var(r_(t-ell))Var(r_(t))) = Var(r_t) = gamma_0$ to convert to autocorrelation function.
$
  rho_ell &= phi_1 rho_(ell - 1) + phi_2 rho_(ell - 2) med .
$
This gives us a second order recursive relation for the autocorrelation function. 
It's a second-order difference equation.
Introducing the lag operator $L rho_ell equiv rho_(ell - 1)$ we can write this as,
$
  (1 - phi_1 L - phi_2 L^2) rho_ell = 0 med . #<eq:ar2-diff>
$
Introduce an ansatz of the form $rho_ell = z^ell$, then
$
  (1 - phi_1 L - phi_2 L^2) z^ell = 0 #no-num \
  => (z^ell - phi_1 z^(ell-1) - phi_2 z^(ell-2)) = 0 #no-num \ 
  => z^(ell - 2) (z^2 - phi_1 z - phi_2) = 0 med .
$
Assume $z != 0$ to find non-trivial solutions. 
This yields the characteristic equation 
$
  z^2 - phi_1 z - phi_2 = 0 med .
$
Roots of this polynomial determine the asymptotic properties of the autocovariance.
The presence of a unit root implies $rho_ell$ grows exponentially with $ell$.


== AR(p) model

#bluebox[
  An $A R(p)$ time series is stationary _if and only if_ its characteristic equation 
  $
    z^p - phi_1 z^(p-1) - dots.c - phi_(p-1) z - phi_p = 0
  $
  has no unit roots $|z_*|^2 < 1$.
]

_Derivation of AIC for AR(p) models_

The likelihood of an AR(p) model, to generate $T$ samples ${r_t : t = 1, 2, ..., T}$, given $p$ previous values ${r_t : t = 0, -1, ..., p - 1}$, and assuming the noise term is Gaussian with zero mean and variance $sigma^2$, is
$
  cal(L(phi)) = product_(t=1)^(T) 1 / sqrt(2 pi sigma^2)
  exp[ 1 / (2 sigma^2) ((1 - phi[L]) r_t)^2 ] med .
$
So the log-likelihood is 
$
  ln cal(L(phi)) = - T / 2 ln (2 pi) - T / 2 ln sigma^2 - 1 / (2 sigma^2) sum_(t=1)^T [(1 - phi[L])r_t]^2 med .
$
Substituting $sigma^2$ and $phi$ with their MLEs $hat(sigma)^2 = "SSE" slash T$, and $hat(phi)$ yields
$
  ln cal(L) 
  &= - T / 2 ln hat(sigma)^2 - 1 / (2 hat(sigma)^2) sum_(t=1)^T [(1 - hat(phi)[L])r_t]^2 med \
  &= - T / 2 ln hat(sigma)^2 - 1 / (2 hat(sigma)^2) "SSE" med \
  &= - T / 2 ln hat(sigma)^2 - 1 / (2 hat(sigma)^2) [(T - ell) hat(sigma)^2] med \
  &= - T / 2 ln hat(sigma)^2 - 1 / (2) (T - ell) med .
$
The last term is a constant and can be dropped, since when we use AIC to compare models the constant terms will be the same across models.
Hence,
#definition[
  _Akaike Information Criterion (AIC)_ (smaller is better)
  $
    "AIC" equiv 2 times ("number of parameters") - 2 ln("likelihood") med .
  $
]
#bluebox[
  For AR$(p)$ models the AIC is given by
  $
    "AIC"(p) = 2 / T p - 2 / T cal(L)
    = (2p) / T + ln hat(sigma)^2 med .
  $
]

#red[Note: here, I'm using SSE to mean "sum of squared errors" $= sum (y - "model")^2$]