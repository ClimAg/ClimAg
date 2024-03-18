# ClimAg: Multifactorial causes of fodder crises in Ireland and risks due to climate change

![ClimAg project logos](https://raw.githubusercontent.com/ClimAg/.github/main/images/logos.png)

[![Documentation Status](https://readthedocs.org/projects/climag/badge/?version=latest)](https://climag.readthedocs.io/?badge=latest)

ClimAg is examining past fodder crises such as the 2018 dry summer and placing them in the context of long-term climate change.

ClimAg seeks to identify the multifactorial drivers of fodder crises by:

- developing a detailed understanding of the multiple interlinked drivers of previous fodder crises affecting the Irish agricultural sector
- combining datasets from 21st century climate simulations with grass growth models to predict the frequency and severity of fodder crisis events under future climate change scenarios

## Repositories

Documentation is available at: <https://climag.readthedocs.io>.
All repositories can be found in the [ClimAg GitHub organisation](https://github.com/ClimAg).
[This repository](https://github.com/ClimAg/ClimAg) hosts Python code for the grass growth model and scripts to perform data preparation, model simulations, and analysis.

## Acknowledgements

ClimAg is a three-year research project funded by the [Environmental Protection Agency (EPA)](https://www.epa.ie/) under the Climate Change Research Programme grant number 2018-CCRP-MS.50, with additional funding provided under the COVID-19 research support scheme of the [Higher Education Authority](https://hea.ie/).

The Python implementation of the [ModVege](https://code.europa.eu/agri4cast/modvege) pasture model adapted for use in this project was translated from Java to Python by Y. Chemin.
This Python implementation was originally published as public domain software on GitHub under the [Unlicence license](https://github.com/ClimAg/modvege).
The Java model was provided by R. Martin of INRAE UREP Clermont-Ferrand for the original Python implementation.
The original ModVege pasture model was developed by [Jouven et al.](https://doi.org/10.1111/j.1365-2494.2006.00515.x).

## Installation

This project uses [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) with [Python](https://www.python.org/) 3.10.

> [!NOTE]
> Windows users should use Conda within [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install), as some packages (e.g. [CDO](https://code.mpimet.mpg.de/projects/cdo)) are unavailable for Windows.

Clone the ClimAg repository including submodules:

```sh
git clone --recurse-submodules https://github.com/ClimAg/ClimAg.git
```

Navigate to the directory of the cloned repository:

```sh
cd ClimAg
```

Create a virtual environment and install all requirements:

```sh
conda env create
```

Activate the virtual environment:

```sh
conda activate ClimAg
```

To run tests:

```sh
python -m pytest --cov
```

To update the virtual environment:

```sh
conda env update
```

To build the documentation locally:

```sh
cd doc && make html && cd ../docs && cp ../doc/_build/html/objects.inv . && make html
```

To clean build the documentation locally:

```sh
cd doc && make clean html && cd ../docs && cp ../doc/_build/html/objects.inv . && make clean html
```

## Licence

Copyright 2022-2024 N. Streethran

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  <https://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
