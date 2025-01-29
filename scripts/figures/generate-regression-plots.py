import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from pathlib import Path

# Global constants and settings
OUT_PATH = Path('../../src/')
plt.style.use('axion')  
plt.rcParams['figure.figsize'] = (6, 3)

def _check_1d(arr):
    """Check if array is 1-dimensional."""
    return len(arr.shape) == 1

def filepath(subdir, filename):
    """Generates the full file path for saving plots, creating directories if necessary."""
    path = OUT_PATH / subdir / filename
    path.parent.mkdir(parents=True, exist_ok=True)  # Create subdir if needed
    return path

def set_xy_labels(ax, x_label="X", y_label="Y"):
    """
    Helper function to set consistent X, Y labels on a given Axes object.
    Modify these defaults if you prefer different labels in specific plots.
    """
    ax.set_xlabel(f"${x_label}$")
    ax.set_ylabel(f"${y_label}$")

class LinearModel:
    """Wrapper around sklearn's LinearRegression model."""
    def __init__(self, **model_kwargs):
        """
        Initialize the LinearRegression model with any additional keyword arguments,
        e.g. LinearModel(fit_intercept=False, etc.).
        """
        self.model = LinearRegression(**model_kwargs)

    def fit(self, x, y, **kwargs):
        """
        Fit the linear model to the data, passing any keyword arguments
        (like `sample_weight`) to .fit() directly.
        """
        x = x.reshape(-1, 1) if _check_1d(x) else x
        self.model.fit(x, y, **kwargs)

    def predict(self, x):
        """Predict using the fitted model."""
        x = x.reshape(-1, 1) if _check_1d(x) else x
        return self.model.predict(x)

def linear_fit(x, y, **kwargs):
    """
    Convenient wrapper to:
    1. create a LinearModel,
    2. fit it on (x, y) with all keyword arguments passed to .fit().

    Examples of kwargs: sample_weight, etc.
    """
    model = LinearModel()
    model.fit(x, y, **kwargs)
    return model

# Function to simulate and plot weak exogeneity violation (Errors in Variables)
def render_weak_exogeneity(subdir, filename):
    """Generates a plot to demonstrate weak exogeneity violation."""
    
    # Simulating data (true x, y) with a linear relationship
    n = 100
    true_x = np.random.normal(size=n)
    true_y = 2.5 * true_x + np.random.normal(scale=0.5, size=n)  # y = 2.5 * x + error
    
    # Introduce measurement error in x
    eta = np.random.normal(scale=0.5, size=n)  # Measurement error
    observed_x = true_x + eta  # Observed x with error
    
    # Fit linear model on true data (no errors in x)
    model_true = linear_fit(true_x, true_y)
    pred_true = model_true.predict(true_x)
    
    # Fit linear model on observed data (with errors in x)
    x = np.linspace(true_x.min(), true_x.max(), 25)
    model_observed = linear_fit(observed_x, true_y)
    pred_observed = model_observed.predict(x)

    # Plotting the results
    fig, ax = plt.subplots()
    
    ax.scatter(true_x, true_y, label='true data', alpha=0.2, color='C1', s=3)
    ax.plot(x, 2.5 * x, label='true model', linestyle='--', color='C1')
    
    ax.scatter(observed_x, true_y, label='observed data (with error)', alpha=0.5, color='C0', s=3)
    ax.plot(x, pred_observed, label='fitted model', linestyle='--', color='C0')
    
    # Use helper function to set labels consistently
    set_xy_labels(ax, x_label="$x$", y_label="$y$")
    ax.legend()
    
    # Save and show
    plt.savefig(filepath(subdir=subdir, filename=filename), bbox_inches='tight')
    plt.show()

# Function to simulate and plot heteroscedasticity violation
def render_heteroscedasticity(subdir, filename):
    """
    Demonstrate heteroscedasticity where errors' scale is drawn from a scaled chi-square 
    distribution with mean=1 and variance=0.25, ensuring positivity.
    We'll compare standard OLS vs. weighted OLS (WLS).
    """
    n = 100
    x = np.random.normal(0, 1, size=n)

    # Chi-square parameters (k=8, alpha=0.125) => E[sigma]=1, Var(sigma)=0.25, sigma>0
    k = 8
    alpha = 2.0
    sigma = alpha * np.random.chisquare(k, size=n)

    # True relationship: y = 2*x + noise, where noise ~ Normal(0, sigma_i)
    y = 2.0 * x + np.random.normal(scale=sigma, size=n)

    # Fit OLS (no weights)
    model_ols = linear_fit(x, y)

    # Fit Weighted OLS, weights = 1 / sigma^2
    w = 1.0 / (sigma**2)
    model_wls = linear_fit(x, y, sample_weight=w)

    # Create a sparser domain of x-values for prediction
    domain_x = np.linspace(x.min(), x.max(), 25)
    pred_ols = model_ols.predict(domain_x)
    pred_wls = model_wls.predict(domain_x)

    # Plot
    fig, ax = plt.subplots()

    ax.scatter(x, y, label='data', alpha=0.5, color='C0', s=3)
    ax.plot(domain_x, pred_ols, label='OLS fit', linestyle='--', color='C0')
    ax.plot(domain_x, pred_wls, label='WLS fit', linestyle='--', color='C1')

    # Use helper to set lower-case axis labels
    set_xy_labels(ax, x_label="$x$", y_label="$y$")

    ax.legend()
    plt.savefig(filepath(subdir, filename), bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    render_weak_exogeneity(subdir="broken-assumptions/figs", filename="weak_exogeneity_example.svg")
    render_heteroscedasticity(subdir="broken-assumptions/figs", filename="heteroscedasticity_example.svg")

