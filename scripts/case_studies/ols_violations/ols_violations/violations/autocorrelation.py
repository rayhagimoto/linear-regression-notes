import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t, norm
from tqdm import tqdm
from dataclasses import dataclass
import statsmodels.api as sm
from scipy.stats import t
if 'axion' in plt.style.available:
    plt.style.use('axion')

from pathlib import Path
MODULE = Path(__file__).name.strip(".py")
OUT_PATH = list(Path(__file__).parents)[5] / f"src/broken-assumptions/figs/{MODULE}"

# from .base_violation import OLSViolationStudy

@dataclass
class SimulationParams:

    sigma : float
    rho : float
    n_simulations : int
    n_samples : int
    beta : float
    alpha : float

    def __iter__(self):
        for (k, v) in self.__dict__.items():
            yield k, v

def generate_ar1_errors(params):
    u_t = np.random.normal(scale=params.sigma, size=(params.n_simulations, params.n_samples))
    et = np.zeros_like(u_t)
    et[:, 0] = u_t[:, 0]
    for i in range(1, params.n_samples):
        et[:, i] = params.rho * et[:, i - 1] + u_t[:, i]
    return et

def generate_iid_errors(params):
    scale_adjusted = params.sigma / np.sqrt(1 - params.rho ** 2)
    return np.random.normal(
        scale=scale_adjusted, size=(params.n_simulations, params.n_samples)
    )

def get_beta_ols_estimates(params, error_generator):

    xt = np.linspace(-1.0, 1.0, params.n_samples)
    et = error_generator(params)
    yt = params.alpha + params.beta * xt + et

    X = sm.add_constant(xt)

    beta_ols_vals = []
    t_stats = []
    false_positives = 0
    for i in range(params.n_simulations):
        mdl = sm.OLS(yt[0], X).fit()
        beta = mdl.params[1]
        beta_ols_vals.append(beta)
        df = params.n_samples - len(mdl.params)
        sigma_hat_sq = sm.ssr / (df)
        ssx = np.var(xt) * xt.shape[0]
        sdev_beta = np.sqrt(sigma_hat_sq / ssx)
        t_stat = beta / sdev_beta
        p_value = 1 - (t.cdf(t_stat, df=df) - t.cdf(-t_stat, df=df))
        t_stats.append(t_stat)

        if p_value < 0.05:
            false_positives += 1
        
    print(f"False positive rate:\t{false_positives / params.n_simulations}")

    results = {
        "yt" : yt,
        "xt" : xt,
        "beta_ols_vals" : beta_ols_vals,
        "false_positive_rate" : false_positives / params.n_simulations,
        "t_stats" : t_stats
    }

    return results

def plot_example_data(params, results=None):

    if results:
        results
        xt = results['iid']['xt']
        yt_iid = results["iid"]["yt"]
        yt_ar1_pos = results["ar1_pos"]["yt"]
        yt_ar1_neg = results["ar1_neg"]["yt"]
    if results is None:
        params.n_simulations = 1
        params.rho = abs(params.rho)

        # Shared
        xt = np.linspace(-1, 1, params.n_samples)

        # WN
        et = generate_iid_errors(params)
        yt_iid = params.alpha + params.beta * xt + et

        # AR(1) pos
        et = generate_ar1_errors(params)
        yt_ar1_pos = params.alpha + params.beta * xt + et

        # AR(1) neg
        params.rho = - abs(params.rho)
        et = generate_ar1_errors(params)
        yt_ar1_neg = params.alpha + params.beta * xt + et

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(6,2), sharey=True)
    ax1.scatter(xt, yt_iid, s=5)
    ax2.scatter(xt, yt_ar1_pos, s=5)
    ax3.scatter(xt, yt_ar1_neg, s=5)

    ax1.set_ylim(-2, 2)
    for ax in (ax1, ax2, ax3):
        ax.set_xlim(-1.3, 1.3)
    ax1.set_ylabel("$\\varepsilon$")
    ax1.set_xlabel("$x$")
    ax2.set_xlabel("$x$")
    ax3.set_xlabel("$x$")
    ax1.set_title(f'white noise')
    ax2.set_title(f"$\\rho = {abs(params.rho)}$")
    ax3.set_title(f"$\\rho = {-abs(params.rho)}$")
    plt.tight_layout()
    filepath = OUT_PATH / 'example_data.png'
    plt.savefig(str(filepath.absolute()), format='png')
    plt.close()


def autocorrelation_study():

    rho = 0.95

    params = SimulationParams(
        sigma=0.2,
        rho=rho,
        n_simulations=1000,
        n_samples = 100,
        beta=0.0,
        alpha=0.0
    )


    # White noise
    results = {}
    results['iid'] = get_beta_ols_estimates(params, generate_iid_errors)
    results['ar1_pos'] = get_beta_ols_estimates(params, generate_iid_errors)
    params.rho = - rho
    results['ar1_neg'] = get_beta_ols_estimates(params, generate_iid_errors)

    plot_example_data(params)

class AutocorrelationStudy(OLSViolationStudy):
    """
    Investigates the effect of autocorrelation on the significance testing under OLS assumptions.
    
    For each simulation the predictor is fixed (np.linspace(-1,1)) and the response is generated in two ways:
      - AR(1) errors generated recursively using a specified ρ.
      - IID errors with an inflated variance (i.e. sigma / sqrt(1-ρ^2)).
      
    In the rendered figure, two subplots show the distribution of OLS estimates (β̂) for:
      ρ = 0.7   and   ρ = -0.7.
    Each subplot overlays the histogram for the AR(1) simulation (blue) and the IID simulation (red),
    plus a dashed black line showing the reference normal PDF (based on the AR(1) average variance).
    """

    def simulate(self, seed=None):
        """
        Run simulations for two values of ρ (0.7 and -0.7). In each case set up:
          - AR(1) error simulations.
          - IID error simulations with inflated variance.
        We then compute the OLS beta estimates, their estimated variance, and the false positive rate.
        """
        if seed is not None:
            np.random.seed(seed)

        # Store simulation results per autoregressive coefficient.
        self.results = {}

        # Define sigma and the predictor (X)
        sigma = 0.25
        xt = np.linspace(-1, 1, self.n_samples)
        self.xt = xt  # in case it is needed later
        self.true_beta = 0.0  # Under the null, no relationship
        self.rho_vals = [0.95, -0.95]

        # Run simulation for each value of ρ
        for rho in self.rho_vals:
            # AR(1) simulation for errors:
            yt_ar1 = self._generate_ar1_errors(rho, sigma)
            beta_ar1, avg_beta_var_ar1, false_rate_ar1 = self._run_ols_simulation(xt, yt_ar1)

            # IID simulation (variance inflated to match AR(1)):
            yt_iid = self._generate_iid_errors(rho, sigma)
            beta_iid, avg_beta_var_iid, false_rate_iid = self._run_ols_simulation(xt, yt_iid)

            # Save the results by ρ value.
            self.results[rho] = {
                "yt_ar1" : yt_ar1,
                "yt_iid" : yt_iid,
                "beta_ar1": beta_ar1,
                "avg_beta_var_ar1": avg_beta_var_ar1,
                "false_rate_ar1": false_rate_ar1,
                "beta_iid": beta_iid,
                "avg_beta_var_iid": avg_beta_var_iid,
                "false_rate_iid": false_rate_iid,
            }

    

    def _run_ols_simulation(self, xt, yt):
        """
        For a given predictor xt and a 2D array of responses yt (one row per simulation),
        run OLS regressions to obtain:
          - The array of beta (slope) estimates.
          - The average estimated variance of beta.
          - The false positive rate based on a two-sided t-test at 5% significance.
        """
        beta_ols_list = []
        false_positives = 0
        sum_beta_var_hat = 0.0

        # Denom used in the variance estimator for beta
        xt_mean = np.mean(xt)
        denom = np.sum((xt - xt_mean) ** 2)

        # Loop over simulations – you can wrap with tqdm if desired.
        for i in range(self.n_simulations):
            cov_val = np.cov(xt, yt[i, :], ddof=self.n_samples - 1)[0, 1]
            var_x = np.var(xt, ddof=self.n_samples - 1)
            beta_ols = cov_val / var_x
            beta_ols_list.append(beta_ols)

            # Estimate sigma squared from residuals and compute beta variance estimate.
            residuals = yt[i, :] - beta_ols * xt
            sigma_hat_sq = np.var(residuals, ddof=1)
            beta_var_hat = sigma_hat_sq / denom
            sum_beta_var_hat += beta_var_hat

            # t-statistic for testing beta = 0:
            t_stat = beta_ols / np.sqrt(beta_var_hat)
            p_value = 1 - (t.cdf(t_stat, df=self.n_samples - 1) -
                           t.cdf(-t_stat, df=self.n_samples - 1))
            if p_value < 0.05:
                false_positives += 1

        avg_beta_var_hat = sum_beta_var_hat / self.n_simulations
        false_positive_rate = false_positives / self.n_simulations

        return np.array(beta_ols_list), avg_beta_var_hat, false_positive_rate

    def render_plots(self):
        """
        Render a single figure with two subplots. Each subplot overlays histogram
        of OLS beta estimates from AR(1) and IID errors for a fixed value of ρ
        (either 0.7 or -0.7). The subplot title shows the corresponding ρ value.
        """

        # 1) Example data

        # Fetch the iid data:
        results = self.results[self.rho_vals[0]]
        yt_iid = results['yt_iid'][0] # get the first one
        yt_ar1_pos = results['yt_ar1'][0] 
        results = self.results[self.rho_vals[1]]
        yt_ar1_neg = results['yt_ar1'][0] 

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(6,2.0), sharey=True)
        ax1.scatter(self.xt, yt_iid, label=r'WN', s=5)
        ax2.scatter(self.xt, yt_ar1_pos, label=f"$\\rho = {self.rho_vals[0]}$", s=5)
        ax3.scatter(self.xt, yt_ar1_neg, label=f"$\\rho = {self.rho_vals[1]}$", s=5)
        ax1.set_ylabel(r'$\epsilon$')
        ax1.set_xlabel('$x$')
        ax2.set_xlabel('$x$')
        ax3.set_xlabel('$x$')
        ax1.text(0, 2.5, r"WN", horizontalalignment='center')
        ax2.text(0, 2.5, f"$\\rho = {self.rho_vals[0]}$", horizontalalignment="center")
        ax3.text(0, 2.5, f"$\\rho = {self.rho_vals[1]}$", horizontalalignment="center")
        plt.ylim(-4, 4)
        plt.tight_layout()
        self.save_figure("example_data")
        plt.close()

        # 2) Histogram of beta_ols estimates compared to gaussian
        fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharey=True)
        bins = 50

        i = 0
        for ax, rho in zip(axes, self.rho_vals):
            results = self.results[rho]

            # Histogram for AR(1) errors (blue)
            ax.hist(
                results["beta_ar1"],
                bins,
                histtype="step",
                density=True,
                color="blue",
                label="AR(1) Errors",
            )
            # Histogram for IID errors (red)
            ax.hist(
                results["beta_iid"],
                bins,
                histtype="step",
                density=True,
                color="red",
                label="IID Errors",
            )
            # Overlay a normal PDF based on the AR(1) simulation's average beta variance
            x_vals = np.linspace(-1.5, 1.5, 200)
            pdf_vals = norm(
                loc=0, scale=np.sqrt(results["avg_beta_var_ar1"])
            ).pdf(x_vals)
            ax.plot(x_vals, pdf_vals, "k--", lw=2, label="Normal PDF")

            ax.set_title(r"$\rho = {}$".format(rho))
            ax.set_xlabel(r"OLS $\hat{\beta}$")
            if rho == self.rho_vals[0]:
                ax.set_ylabel("Density")
            if i == 0:
                ax.legend()
                i += 1

            ax.set_xlim(-1.5, 1.5)

        plt.tight_layout()
        self.save_figure("beta_ols_histogram")
        plt.close()


if __name__ == "__main__":
    # Adjust n_simulations and n_samples as needed (here we use a reduced number for speed)
    study = AutocorrelationStudy(n_simulations=1000, n_samples=100)
    study.simulate(seed=123)
    study.render_plots()
