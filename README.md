[![pipeline status](https://gitlab.pyadf.org/vib/LocVib/badges/install_requires_sphinx_doc/pipeline.svg)](https://gitlab.pyadf.org/vib/LocVib/-/commits/install_requires_sphinx_doc)
[![Custom Badge](https://gitlab.pyadf.org/vib/LocVib/badges/install_requires_sphinx_doc/coverage.svg?min_medium=30&min_acceptable=50&min_good=75&key_text=Unittest+Coverage&key_width=130)](https://gitlab.pyadf.org/vib/LocVib/-/commits/install_requires_sphinx_doc)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![Numpy 1.23.4](https://img.shields.io/badge/numpy-1.23.4-blue.svg)](https://numpy.org/devdocs/release/1.23.4-notes.html)
[![Matplotlib 3.6.1](https://img.shields.io/badge/matplotlib-3.6.1-blue.svg)](https://matplotlib.org/stable/users/installing/index.html)
[![pip - Openbabel 3.1.1](https://img.shields.io/badge/openbabel-3.1.1-cornflowerblue.svg)](https://bioweb.pasteur.fr/packages/pack@openbabel@3.1.1)
[![conda - Openbabel 3.1.1.16](https://img.shields.io/badge/openbabel--wheel-3.1.1.16-blueviolet.svg)](https://pypi.org/project/openbabel-wheel/3.1.1.16/)




# LocVib/VibTools

Python tools for localizing normal modes.

Copyright (C) 2009-2023 by Christoph R. Jacob and others.

In scientific publications using the LocVib tools, please cite:
  Ch. R. Jacob, J. Chem. Phys 130 (2009), 084106 (https://doi.org/10.1063/1.3077690).

We can find the documentation for LocVib/VibTools on the website

https://vib.gitlab-pages.pyadf.org/LocVib/ 

Download-Link:

https://github.com/chjacob-tubs/LocVib

## Requirements

LocVib is an independent code that for running needs only Python standard
packages, extended with NumPy, Matplotlib and Openbabel.
More details are in the documentation.

## Installation

Just clone this repository and update `$PYTHONPATH` environment variable accordingly.
Or install it using pip and run 

`LocVib/ % pip install . ` 

in the folder where setup.py is located.
More details are in the documentation.

Verify the installation with running pytest (must be installed before):

`/LocVib/tests % pytest - v`

Other installation options can be found in the further documentation.

## Documentation

We can find the documentation, as mentioned above,
on the Homepage (https://vib.gitlab-pages.pyadf.org/LocVib/) 
or we can generate the HTML documentation  ourselves with Sphinx.

Build documentation via Sphinx and extension packages (must be installed before)

`Locvib/doc/ % sphinx-build . Build`.

Opening the index.html file takes you to the home directory of the code documentation:

`Locvib/doc/Build/ % open index.html`

## Usage

See `example/` directory for some examples of typical runs.

### Any suggestions and improvements are welcome.
