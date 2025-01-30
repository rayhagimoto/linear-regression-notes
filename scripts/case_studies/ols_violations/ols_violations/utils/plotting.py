import matplotlib.pyplot as plt
import numpy as np

def _default_kwargs(dict):
    if dict:
        return dict
    else:
        return {}

def plot_distribution(data, true_value, xlabel, title, label, color):
    plt.hist(data, bins=50, alpha=0.6, density=True, label=label, color=color)
    if true_value:
        plt.axvline(true_value, color='red', linestyle='dashed', linewidth=2, label=f"True {xlabel}")
    plt.xlabel(xlabel)
    plt.ylabel("Density")
    plt.legend()
    plt.title(title)

def plot(x, y, ax=None, figsize=(10,6), xlabel=None, ylabel=None, title=None, plot_kwargs=None):

    # Ensure scatter_kwargs is at least an empty dictionary
    plot_kwargs = _default_kwargs(plot_kwargs)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x, y, **plot_kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

def scatter(x, y, ax=None, figsize=(10,6), xlabel=None, ylabel=None, title=None, scatter_kwargs=None):

    # Ensure scatter_kwargs is at least an empty dictionary
    scatter_kwargs = _default_kwargs(scatter_kwargs)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x, y, **scatter_kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

def plot_residuals(x, residuals):
    kwargs = dict(
        xlabel="X",
        ylabel="Residuals",
        title="Residuals vs X (Heteroskedasticity Pattern)",
        scatter_kwargs = dict(alpha=0.6, label="Residuals"),
    )
    
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter(x, residuals, ax=ax, **kwargs)
    plt.axhline(0, color='red', linestyle='dashed', linewidth=2, label="Zero residual line")
    plt.legend()

def plot_loghist(x, n : int, **kwargs):
    "Plot a histogram of x where x is binned in a logspace."
    x = np.array(x).flatten()
    bins = np.geomspace(max(x.min(), 1e-13), max(x.max(), 1e-13), n)
    _ = plt.hist(x, bins, **kwargs)
