************
Installation
************

The *LocVib* *Python* package releases with its *VibTools* modules 
are available as a Github repository
https://github.com/chjacob-tubs/LocVib.

LocVib relies on additional *Python* packages, that have to be installed on our system:

* **Python 3** (https://www.python.org/).

* **Numpy** (https://numpy.org/).

* **Matplotlib** (https://matplotlib.org/).

* **Openbabel 3** (https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python.html).

Additional (optinal) extensions:

* **Pytest** (https://pytest.org) (recommended).

* **Sphinx** (https://www.sphinx-doc.org)

For the installation we just need to download the file and follow the instructions below.

.. note::
   We highly recommend using a **Conda environment** (see `Conda Environment`_).

Download the Code
=================

Download the Github repository as a zip file from the link below:

https://github.com/chjacob-tubs/LocVib

Unzip the LocVib-Package-Zip-File:

   >>> unzip file.zip

We should now find the following folder structure:

.. code-block:: bash

   LocVib/
   ├── COPYING
   ├── doc/
   ├── Dockerfile
   ├── example/
   ├── README.md
   ├── requirements/
   ├── setup.py
   ├── src/VibTools/
   └── tests/

Install with pip (easiest method)
=================================

Execute the following command in the main code directory
(Python-3 and a current pip version must already be installed):


    >>> LocVib/ % pip install .

.. _automatic pip installation:


Install with Conda
==================

We highly recommend using the *Conda* package manager(https://conda.io/) for the environment 
and the use of *Anaconda* (https://www.anaconda.com/) for using *Python*.
The above *Python* dependencies (*Numpy*, *Matplotlib*, *Openbabel*) must be installed for *LocVib* to work. 

With *Pytest* we can determine the correct executability of the program.

If we are interested in the further development of *LocVib* ourselfs, 
it makes sense to install *Sphinx* for the documentation.

.. _Conda Environment:

Installation Conda Environment
------------------------------

First we install *Conda* on our system.
The best way to do this is to follow the instructions 
on the *Conda* homepage (https://conda.io/).

Here we show in short form which *Conda* commands are necessary 
to install the necessary *Python* packages, 
provided that our conda installation worked.

**Conda initialization:**

   >>> % conda init zsh

Or usage without initialization only with `source activate`:

   >>> % source activate
   >>> (base)%

Pip Installation in Conda Environment (recommended)
---------------------------------------------------

Create and activate the LocVib Conda environment:

   >>> (base)% conda create --name VibToolsCondaENV python=3.11.4
   >>> (base)% conda activate VibToolsCondaENV
   >>> (VibToolsCondaENV)% 

As a prerequisite we still need the pip package:

   >>> (VibToolsCondaENV)% conda install pip

Select the *LocVib* folder and run the *pip* installation:

   >>> (VibToolsCondaENV)% cd LocVib
   >>> (VibToolsCondaENV)/LocVib% pip install .

.. code-block:: console

   Processing ~/LocVib
     Preparing metadata (setup.py) ... done
   Requirement already satisfied: numpy in /home/name/.conda/envs/LVenv/lib/python3.10/site-packages (from VibTools==0.0.1) (1.22.3)
   Requirement already satisfied: matplotlib in /home/name/.conda/envs/LVenv/lib/python3.10/site-packages (from VibTools==0.0.1) (3.5.2)
   Requirement already satisfied: pytest in /home/name/.conda/envs/LVenv/lib/python3.10/site-packages (from VibTools==0.0.1) (7.1.2)
   Requirement already satisfied: python-dateutil>=2.7 in /home/name/.conda/envs/LVenv/lib/python3.10/site-packages (from matplotlib->VibTools==0.0.1) (2.8.2)
   Requirement already satisfied: cycler>=0.10 in /home/name/.conda/envs/LVenv/lib/python3.10/site-packages (from matplotlib->VibTools==0.0.1) (0.11.0)
   .
   .
   .
   Building wheels for collected packages: VibTools
     Building wheel for VibTools (setup.py) ... done
     Created wheel for VibTools: filename=VibTools-0.0.1-py3-none-any.whl size=96726 sha256=0e2110eaffbb70ba64c2e3cf5bf1dd724387a642a73ead3b244b454dea79ff9b
     Stored in directory: /tmp/pip-ephem-wheel-cache-jijb06l6/wheels/70/2a/ae/c4a6afe46f78a2dd633299e079f4d909310bc94ec529e1388d
   Successfully built VibTools
   Installing collected packages: VibTools
   Successfully installed VibTools-0.0.1 

.. note:: 
   We can also install with pip in developer mode (editable). 
      >>> (VibToolsCondaENV)/LocVib% pip install -e .


Here everything is done regarding the installation.

**Optional:** For verification (see `Verify Installation with Pytest`_) of successful installation install **pytest** additionally:

   >>> (VibToolsCondaENV)/LocVib% conda install -c conda-forge pytest=7.2.0

.. note::
   The automatic Pip installation also works without Conda environment but Pip must be installed anyway.

Semi-Automatic Installation with Conda
--------------------------------------

Go to the LocVib folder and perform the creation of the appropriate environment:

   >>> (base)% cd LocVib/
   >>> (base)LocVib/ % conda env create -f requirements/environment.yml
   Collecting package metadata (repodata.json): done
   Downloading and Extracting Packages
   kiwisolver-1.4.4     | 70 KB     | ################################################## | 100% 
   libllvm14-14.0.6     | 33.4 MB   | ################################################## | 100% 
   openssl-1.1.1v       | 3.7 MB    | ################################################## | 100% 
   .
   .
   .
   certifi-2023.7.22    | 154 KB    | ################################################## | 100% 
   Solving environment: done
   Preparing transaction: done
   Verifying transaction: done
   Executing transaction: done
   #
   # To activate this environment, use
   #
   #     $ conda activate VibToolsCondaENV
   #
   # To deactivate an active environment, use
   #
   #     $ conda deactivate

Now you only need to activate the environment and add LocVib:

   >>> (base)LocVib/ % conda activate VibToolsCondaENV
   >>> (VibToolsCondaENV)LocVib% conda develop src/
   added /home/User/LocVib/src
   completed operation for: /home/User/LocVib/src

Manual Installation
-------------------

**Creating suitable LocVib environment:**

   >>> conda create -n VibToolsCondaENV

   >>> conda activate VibToolsCondaENV

**Installation of the necessary packages:**

   >>> (VibToolsCondaENV)% conda install anaconda
   >>> (VibToolsCondaENV)% conda install -c conda-forge python~=3.11.4
   >>> (VibToolsCondaENV)% conda install -c conda-forge numpy~=1.23.4
   >>> (VibToolsCondaENV)% conda install -c conda-forge matplotlib`=3.6.1
   >>> (VibToolsCondaENV)% conda install -c conda-forge openbabel~=3.1.1
   >>> (VibToolsCondaENV)% conda install -c conda-forge pytest~=7.2.0

optional for developing the Documentation:

   >>> (VibToolsCondaENV)% conda install -c conda-forge sphinx~=5.3.0
   >>> (VibToolsCondaENV)% pip install sphinx_rtd_theme~=0.4.3
   >>> (VibToolsCondaENV)% pip install sphinx_mdinclude~=0.5.3

Installation of LocVib itself:

   >>> (VibToolsCondaENV)LocVib% conda develop src/
   added /home/User/LocVib/src
   completed operation for: /home/User/LocVib/src

.. warning:: 
   If we use an up-to-date conda version, unfortunately the develop command is no longer supported
   and the pip editable mode is recommended(`pip install -e .`).

.. note::
   In principle, other versions of the respective packages are also usable, 
   but with the specified versions, the runnability is guaranteed in any case.

.. _download LocVib:




Manual Installation (PYTHONPATH)
================================

We have to include the subdirectory 'LocVib/src/VibTools' in our
PYTHONPATH environment variable.

We can modify your *.zprofile* file with adding:

   >>> export PYTHONPATH="${PYTHONPATH}:/home/yourname/LocVib/src/"

Or we use the following lines of code in our scripts for importing LocVib:

   >>> import sys
   >>> sys.path.append('/home/yourname/LocVib/src/')



Verify Installation with Pytest
===============================

.. _Verify Installation with Pytest:

The prerequisite for the check is that we have pytest installed.

Go to the appropiate test folder:

   >>> % cd LocVib/tests/

Run the test:

   >>> LocVib/tests/% pytest -v

If everything runs correctly, we will get the following output:

.. code-block:: console

    ============================================== test session starts ==============================================
    platform linux -- Python 3.10.4, pytest-7.1.2, pluggy-1.0.0 -- /home/yourname/.conda/envs/LVenv/bin/python
    cachedir: .pytest_cache
    rootdir: /home/yourname/LocVib
    collected 104 items                                                                                             
    
    test_VibTools.py::test_read_from_coord PASSED                                                             [  0%]
    test_VibTools.py::test_get_fragment PASSED                                                                [  1%]
    test_VibTools.py::test_reset_molcache PASSED                                                              [  2%]
    test_VibTools.py::test_write_and_read PASSED                                                              [  3%]
    test_VibTools.py::test_get_natoms PASSED                                                                  [  4%]
    test_VibTools.py::test_get_atmasses PASSED                                                                [  5%]
    test_VibTools.py::test_get_atnums PASSED                                                                  [  6%]
    test_VibTools.py::test_get_coordinates PASSED                                                             [  7%]
    test_VibTools.py::test_add_atoms PASSED                                                                   [  8%]
    test_VibTools.py::test_residue_groups PASSED                                                              [  9%]
    .
    .
    .
    test_VibTools.py::test_get_gaussian_spectrum PASSED                                                       [ 94%]
    test_VibTools.py::test_get_rect_spectrum PASSED                                                           [ 95%]
    test_VibTools.py::test_scale_range PASSED                                                                 [ 96%]
    test_VibTools.py::test_get_band_maxima PASSED                                                             [ 97%]
    test_VibTools.py::test_get_band_minima PASSED                                                             [ 98%]
    test_VibTools.py::test_get_plot PASSED                                                                    [ 99%]
    test_VibTools.py::test_get_rect_plot PASSED                                                               [100%]
    
    ============================================= 104 passed in 14.38s ==============================================

.. note::
   Another possibility to check the executability of the program is to calculate the code examples. 
   See

   .. toctree::

      examples
