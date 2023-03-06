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
      install_requires = ['numpy','matplotlib','pytest'],
      classifiers  = ["Programming Language :: Python :: 3",
                      "Operating System :: OS Independent"],
     )
