import sys
from pathlib import Path
abs_path = Path(__file__).parent / "case_studies/ols_violations/ols_violations"
print("exists", (abs_path / "violations").exists())
print(f"Module Path {abs_path}")
sys.path.insert(0, str(abs_path))

import matplotlib.pyplot as plt
if 'axion' in plt.style.available:
    plt.style.use('axion')

from violations import HomoscedasticityStudy, AutocorrelationStudy

def render(s):
    "Wrapper that runs simulation and renders a plot"
    study = s()
    study.simulate(seed=12345)
    study.render_plots()

if __name__ == '__main__':
    
    studies = [
        AutocorrelationStudy,
    ]

    for study in studies:
        render(study)
