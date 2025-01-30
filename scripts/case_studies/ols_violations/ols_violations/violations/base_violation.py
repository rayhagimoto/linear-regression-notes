import os
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

class OLSViolationStudy(ABC):
    def __init__(self, n_simulations=10000, n_samples=100):
        """
        Base class for OLS violation studies.

        Args:
            filename (str): Base filename for saving figures.
            n_simulations (int): Number of Monte Carlo simulations.
            n_samples (int): Number of samples per simulation.
            true_beta (float): True beta coefficient.
        """
        # Automatically derive subdir from class name
        self.subdir = self.__class__.__name__.replace("Study", "").lower()
        self.n_simulations = n_simulations
        self.n_samples = n_samples

    @abstractmethod
    def simulate(self, seed=None):
        """Method to run simulations (implemented in subclasses)."""
        pass

    @abstractmethod
    def render_plots(self):
        """Method to generate and display plots (implemented in subclasses)."""
        pass

    def save_figure(self, filename):
        """
        Saves the current figure to the inferred subdirectory inside `/src/broken-assumptions/figs/`.

        Args:
            index (int): The figure number (e.g., 1, 2, 3) to append to the filename.
        """
        # Define the full path to the subdirectory
        base_dir = os.path.join("src", "broken-assumptions", "figs", self.subdir)

        # Create the directory if it doesn't exist
        os.makedirs(base_dir, exist_ok=True)

        # Construct full file path with index
        file_path = os.path.join(base_dir, f"{filename}.svg")

        # Save the figure
        plt.savefig(file_path, format="svg", bbox_inches="tight")
        print(f"Figure saved: {file_path}")
