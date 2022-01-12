# marc-tools
Collection of tools for post-processing of marc.

Tools are based on the newest hdf-file format

# Installation

1. Clone the repo from github

## Install in developer mode

2. Create a virtual environment for the project. In a terminal, move to
repository root (git\flex-fatigue-tools) and run deploy:
```
python deploy.py
```
This creates a .venv folder with a virtual environment, installing all packages lited in requirements.txt

3. Activate environment
```
.venv\Scripts\activate
```

4. Install necessary notebook extensions:
```
python -m ipykernel install --user --name=flex-fatigue-tools
```

5. Enable developer mode:
```
pip install -e .
```

## Using the package

You can now import the package in Jupyter and Python (you may have to restart 
the Jupyter session)

In your notebook, run
```
import flex_fatigue_tools as fft
```
to import the module


## Testing

Put tests in the *tests* folder and run `tox`. Tox is install in the base python
environment so you don't have to activate your project environment.

You can look at test coverage by looking at *.tox/py38/tmp/cov/index.html*

Testing.
