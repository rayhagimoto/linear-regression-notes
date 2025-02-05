from .base_violation import OLSViolationStudy
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from utils.plotting import plot_distribution, plot_residuals, plot_loghist


class HomoscedasticityStudy(OLSViolationStudy):
    """
    Investigates heteroscedasticity in a model where sigma = 0.5 * X.
    """

    def simulate(self, seed=None):
        rng = np.random.default_rng(seed=seed)

        self.true_beta = 2.0
        self.ols_betas = []
        self.wls_betas = []
        self.ols_alphas = []
        self.wls_alphas = []
        for _ in range(self.n_simulations):
            X = rng.uniform(1, 5, self.n_samples)
            X = sm.add_constant(X)

            # Heteroscedastic errors
            errors = rng.normal(loc=0, scale=0.5 * X[:, 1], size=self.n_samples)
            Y = self.true_beta * X[:, 1] + errors

            # Fit OLS
            ols_model = sm.OLS(Y, X).fit()
            self.ols_alphas.append(ols_model.params[0])
            self.ols_betas.append(ols_model.params[1])

            # Compute residuals for first dataset
            if _ == 0:
                residuals = ols_model.resid
                self.residuals = residuals

            # Fit WLS
            weights = 1 / (0.5 * X[:, 1])**2
            wls_model = sm.WLS(Y, X, weights=weights).fit()
            self.wls_alphas.append(wls_model.params[0])
            self.wls_betas.append(wls_model.params[1])
        
        # Store the last simulated data
        self.ols_model = ols_model.model
        self.wls_model = wls_model.model
        self.X = X
        self.Y = Y

    def render_plots(self):
        # First plot: Distribution of beta estimates
        plt.figure(figsize=(6, 3))
        plot_loghist(self.ols_betas, 50, label='OLS fit', color='C0', density=True)
        plot_loghist(self.wls_betas, 50, label='WLS fit', color='C1', density=True)
        plt.axvline(self.true_beta, label='true', linestyle='dashed', color='red')
        plt.legend()
        plt.xlabel("$\\hat{\\beta}$")
        plt.ylabel("Density")
        self.save_figure("beta-hist")
        plt.close()

        # Second plot: Residuals vs X
        plot_residuals(self.X[:, 1], self.residuals)
        self.save_figure("resid-vs-x")
        plt.close()

        # Third plot: WLS Preds, OLS Preds, True trend, Data
        x = np.linspace(1, 5, 25)
        X = self.X
        Y = self.Y
        ols_pred = self.ols_alphas[-1] + self.ols_betas[-1] * x
        wls_pred = self.wls_alphas[-1] + self.wls_betas[-1] * x
        true_trend = self.true_beta * x

        plt.figure(figsize=(6,3))
        plt.scatter(X[:, 1], Y, alpha=0.3, c='C0', label='Data')
        plt.plot(x, ols_pred, label=f'OLS $\\hat{{\\beta}} = {self.ols_betas[-1]:.2f}$', linestyle='dashed')
        plt.plot(x, wls_pred, label=f'WLS $\\hat{{\\beta}} = {self.wls_betas[-1]:.2f}$', linestyle='dashed')
        plt.plot(x, true_trend, label=f'True $\\hat{{\\beta}} = {self.true_beta:.2f}$', c='black')
        plt.ylim(1e-1, 13 )
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.legend()
        self.save_figure("trendlines")
