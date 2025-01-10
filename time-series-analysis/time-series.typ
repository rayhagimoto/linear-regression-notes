#import "@local/templates:1.0.0" : *

= Time series analysis 

The foundation of time series analysis is stationarity.
#definition[
  A time series ${r_t}$ is said to be weakly stationary if $EE[r_t] = mu$ is independent of $t$ and the autocovariance $EE[(r_t - mu)(r_(t-ell) - mu)] = gamma_ell$ is only a function of the lag. 
]

Stationarity can be checked using the Dickey-Fuller (or augmented Dickey-Fuller) test.
To explain this test let's look at a simple example of a linear time series.

== AR(2) model

$
  r_t = phi_0 + phi_1 r_(t-1) + phi_2 r_(t-2) + epsilon_t
$

$
  EE[r_t] 
  &= phi_0 + phi_1 EE[r_(t-1)] + phi_2 EE[r_(t-2)] + EE[epsilon_t] \
  => mu &=  phi_0 + phi_1 mu + phi_2 mu \
  => mu &= phi_0 / (1 - phi_1 - phi_2) med .
$
Rewriting the time series (in terms of deviations from the mean):
$
  r_t - mu = phi_1 (r_(t-1) - mu) + phi_2 (r_(t-2) - mu) + epsilon_t med .
$
$
  (r_(t - ell) - mu)(r_t - mu) = phi_1 (r_(t - ell)- mu)(r_(t-1) - mu) + phi_2 (r_(t - ell) - mu)(r_(t-2) - mu) + epsilon_t
$
Taking the expectation,
$
  EE[(r_(t - ell) - mu)(r_t - mu)] 
  &= phi_1 EE[(r_(t - ell)- mu)(r_(t-1) - mu)] #no-num \
  &quad + phi_2 EE[(r_(t - ell) - mu)(r_(t-2) - mu)] + EE[epsilon_t] #no-num \
  => gamma_ell &= phi_1 gamma_(ell - 1) + phi_2 gamma_(ell - 2) med . #<eq:autocovf-recursion>
$
Divide #eref(<eq:autocovf-recursion>) by $gamma_0$ to convert to autocorrelation function.
$
  rho_ell &= phi_1 rho_(ell - 1) + phi_2 rho_(ell - 2) med .
$
Can be written as a second-order difference equation
$
  (1 - phi_1 L - phi_2 L^2) rho_ell = 0 med .
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

#bluebox[
  An $A R(p)$ time series is stationary _if and only if_ its characteristic equation 
  $
    z^p - phi_1 z^(p-1) - dots.c - phi_(p-1) z - phi_p = 0
  $
  has no unit roots $|z_*|^2 < 1$.
]