import setuptools

setuptools.setup(
      include_package_data = True,
      name         = 'VibTools',
      version      = '0.0.1',
      description  = 'LocVib/VibTools: Python tools for localizing normal modes',
      author       = 'Christoph Jacob, Julia Brueggemann, Mario Wolter, Michael Welzel and others',
      url          = 'https://www.tu-braunschweig.de/pci/agjacob/software',
      license      = 'GPLv3',
      package_dir  = {'': 'src/'},
      python_requires = '>=3.10.4',
      install_requires = ['numpy==1.23.4','matplotlib==3.6.1','openbabel-wheel==3.1.1.16','pytest==7.2.0'],
      classifiers  = ["Programming Language :: Python :: 3",
                      "Operating System :: OS Independent"],
     )
