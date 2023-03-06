************
Installation
************

The *LocVib* *Python* package releases with its *VibTools* modules 
are available as a zip file on website 
https://www.tu-braunschweig.de/pci/agjacob/software.

LocVib relies on additional *Python* packages, that have to be installed on your system:

* **Python 3** (https://www.python.org/).

* **Numpy** (https://numpy.org/).

* **Matplotlib** (https://matplotlib.org/).

* **Openbabel 3** (https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python.html).

Additional (optinal) extensions:

* **Pytest** (https://pytest.org) (recommended).

* **Sphinx** (https://www.sphinx-doc.org)

For the installation you just need to download the file and follow the instructions below.

.. note::
   After downloading *LocVib* (`download LocVib`_) 
   you can start an **automatic installation using pip** (`automatic pip installation`_).
   But openbabel must be properly integrated/installed in your system.

.. note::
   We highly recommend using a **Conda environment** (see `Conda Environment`_).

.. warning::
    A pip installation for **Openbabel** on linux/ios is not available at the moment, 
    use the conda install or sudo-apt get options instead.

Suitable Python Environment (Linux/IOs)
=======================================

We highly recommend using the *Conda* package (https://conda.io/) for the environment 
and the use of *Anaconda* (https://www.anaconda.com/) for using *Python*.
The above *Python* dependencies (*Numpy*, *Matplotlib*, *Openbabel*) must be installed for *LocVib* to work. 

With *Pytest* you can determine the correct executability of the program.

If you are interested in the further development of *LocVib* yourself, 
it makes sense to install *Sphinx* for the documentation.

.. _Conda Environment:

Conda Environment (recommended)
-------------------------------

First you install *Conda* on your system.
The best way to do this is to follow the instructions 
on the *Conda* homepage (https://conda.io/).

Here we show in short form which *Conda* commands are necessary 
to install the necessary *Python* packages, 
provided that your conda installation worked.

Manual Installation
^^^^^^^^^^^^^^^^^^^

**Conda initialization:**

   >>> conda init zsh

**Or usage without initialization:**

   >>> source activate LVenv

or

   >>> conda create -n LVenv
   >>> conda activate LVenv

**Creating suitable LocVib environment:**

   >>> conda install anaconda
   >>> conda install -c conda-forge numpy
   >>> conda install -c conda-forge matplotlib
   >>> conda install -c conda-forge openbabel
   >>> conda install -c conda-forge pytest

optional for developing the Documentation:

   >>> conda install -c conda-forge sphinx
   >>> pip install sphinx_theme


.. _download LocVib:

Download the Code
=================

Download the zip file from the link below:

https://www.tu-braunschweig.de/pci/agjacob/software

Unzip the LocVib-Package-Zip-File:

   >>> unzip file.zip

You should now find the following folder structure:

.. code-block:: bash

   LocVib
   ├── COPYING
   ├── doc
   ├── environment.yml
   ├── example
   ├── README.md
   ├── requirements.txt
   ├── setup.py
   ├── src
   └── tests

.. _automatic pip installation:

(Automatic) Pip Installation
----------------------------

As a prerequisite we still need the pip package:

   >>> conda install -c anaconda pip

Select the *LocVib* folder and run the *pip* installation:

   >>> cd LocVib
   >>> /LocVib/pip install .

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



Manual Installation
-------------------

Using Pythonpath
^^^^^^^^^^^^^^^^

you have to include the subdirectory 'LocVib/src/VibTools' in your
PYTHONPATH environment variable.

You can modify your *.zprofile* file with adding:

   >>> export PYTHONPATH="${PYTHONPATH}:/home/yourname/LocVib/src/"

Or you use the following lines of code in your scripts for importing LocVib:

   >>> import sys
   >>> sys.path.append('/home/yourname/LocVib/src/')


Add LocVib to Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installing *LocVib* in your conda environment:

    >>> conda develop /home/name/LocVib/src/


Verify LocVib Installation with Pytest
======================================

The prerequisite for the check is that you have pytest installed.

Go to the appropiate test folder:

    >>> cd LocVib/test/

Run the test:

    >>> LocVib/test/pytet -v test_VibTools

If everything runs correctly, you will get the following output:

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

If an error occurred, check your installation or contact support (email: jacob_software_support@tu-bs.de) 
or ask Stackoverflow.

.. note::
   Another possibility to check the executability of the program is to calculate the code examples. 
   See

   .. toctree::

      examples
