from pathlib import Path

for parent in Path(".").absolute().parents:
    if Path(parent / "src").exists():
        workspaceDirectory = parent
        pysrc_path = str(workspaceDirectory / "scripts/case_studies/ols_violations")


# Create the YAML string with the absolute PYTHONPATH
env_yaml = f"""
name: ols_env
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - numpy
  - matplotlib
  - statsmodels
  - scipy
  - jupyter
  - ipykernel
  - nbdime
variables:
  PYTHONPATH: {pysrc_path}
"""

# Write the YAML content to 'environment.yml'
SCRIPT_PATH = str(Path(__file__).parent.absolute())
ENV_PATH = f"{SCRIPT_PATH}/case_studies/ols_violations/environment.yml"

print(f"Writing to", ENV_PATH)

with open(ENV_PATH, "w") as file:
    file.write(env_yaml)
