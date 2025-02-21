import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t, norm
from tqdm import tqdm

from .base_violation import OLSViolationStudy


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
        sigma = 1.0
        xt = np.linspace(-1, 1, self.n_samples)
        self.xt = xt  # in case it is needed later
        self.true_beta = 0.0  # Under the null, no relationship

        # Run simulation for each value of ρ
        for rho in [0.7, -0.7]:
            # AR(1) simulation for errors:
            yt_ar1 = self._generate_ar1_errors(rho, sigma)
            beta_ar1, avg_beta_var_ar1, false_rate_ar1 = self._run_ols_simulation(xt, yt_ar1)

            # IID simulation (variance inflated to match AR(1)):
            yt_iid = self._generate_iid_errors(rho, sigma)
            beta_iid, avg_beta_var_iid, false_rate_iid = self._run_ols_simulation(xt, yt_iid)

            # Save the results by ρ value.
            self.results[rho] = {
                "beta_ar1": beta_ar1,
                "avg_beta_var_ar1": avg_beta_var_ar1,
                "false_rate_ar1": false_rate_ar1,
                "beta_iid": beta_iid,
                "avg_beta_var_iid": avg_beta_var_iid,
                "false_rate_iid": false_rate_iid,
            }

    def _generate_ar1_errors(self, rho, sigma):
        """
        Generate AR(1) errors of shape (n_simulations, n_samples) using
        
            e[0] = u[0]
            e[t] = ρ * e[t-1] + u[t]
        
        where u ~ N(0, sigma).
        """
        u_t = np.random.normal(scale=sigma, size=(self.n_simulations, self.n_samples))
        et = np.zeros_like(u_t)
        et[:, 0] = u_t[:, 0]
        for i in range(1, self.n_samples):
            et[:, i] = rho * et[:, i - 1] + u_t[:, i]
        return et

    def _generate_iid_errors(self, rho, sigma):
        """
        Generate IID errors of shape (n_simulations, n_samples) with standard deviation
        sigma / sqrt(1 - ρ^2) so that their overall variance matches the AR(1) process.
        """
        scale_adjusted = sigma / np.sqrt(1 - rho ** 2)
        return np.random.normal(
            scale=scale_adjusted, size=(self.n_simulations, self.n_samples)
        )

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
        fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharey=True)
        bins = 50

        i = 0
        for ax, rho in zip(axes, [0.7, -0.7]):
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
            ax.set_ylabel("Density")
            if i == 0:
                ax.legend()
                i += 1

        plt.tight_layout()
        self.save_figure("autocorrelation-study")
        plt.show()


if __name__ == "__main__":
    # Adjust n_simulations and n_samples as needed (here we use a reduced number for speed)
    study = AutocorrelationStudy(n_simulations=1000, n_samples=100)
    study.simulate(seed=123)
    study.render_plots()