# marc-tools
Collection of tools for post-processing of marc.

Tools are based on the newest hdf-file format

## Install in developer mode

1. Clone the repo from github

2. Move to repository root in terminal:
```
cd git\marc-tools
```

3. Run deploy script:
```
python deploy.py
```
This creates a .venv folder with a virtual environment, installing all packages lited in requirements.txt

4. Activate environment
```
.venv\Scripts\activate
```

5. Install necessary notebook extensions:
```
python -m ipykernel install --user --name=marc-tools
```

6. Enable developer mode:
```
pip install -e .
```

## Using the package

You can now import the package in Jupyter and Python (you may have to restart 
the Jupyter session)

In your notebook, run
```
import marc_tools as mt
```
to import the module


## Testing

Put tests in the *tests* folder and run `tox`. Tox is install in the base python
environment so you don't have to activate your project environment.

You can look at test coverage by looking at *.tox/py38/tmp/cov/index.html*

Testing.
