.. image:: http://mkoeppe.github.io/cutgeneratingfunctionology/graphics-nonfree/Z11_058_github_template.jpg
   :width:  100%
   :target: https://github.com/mkoeppe/cutgeneratingfunctionology
   :alt:    mkoeppe/cutgeneratingfunctionology: Python code for computation and experimentation with cut-generating functions

Most of the code is for the 1-dimensional Gomory-Johnson infinite
group problem, including an electronic compendium of extreme
functions.

See the survey "Light on the Infinite Group Relaxation" 
(http://www.optimization-online.org/DB_HTML/2014/10/4620.html)
for the mathematical background and a table of functions in the 
electronic compendium.  See also the paper "An electronic compendium 
of extreme functions for the Gomory--Johnson infinite group problem"
(http://www.optimization-online.org/DB_HTML/2014/11/4646.html) for 
a discussion of several functions in the compendium.

See http://www.sagemath.org/doc/tutorial/ for information on how to
use Sage.

.. image:: https://img.shields.io/travis/mkoeppe/cutgeneratingfunctionology
   :alt: Travis CI
   :target: https://travis-ci.org/mkoeppe/cutgeneratingfunctionology/
.. image:: https://img.shields.io/pypi/l/cutgeneratingfunctionology
   :alt: License: GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
.. image:: https://img.shields.io/pypi/v/cutgeneratingfunctionology
   :alt: PyPI package
   :target: https://pypi.org/project/cutgeneratingfunctionology/
.. https://img.shields.io/pypi/pyversions/cutgeneratingfunctionology
.. image:: https://mybinder.org/badge_logo.svg
   :alt: Run it on mybinder.org
   :target: https://mybinder.org/v2/gh/mkoeppe/cutgeneratingfunctionology/master?filepath=demo.ipynb
.. image:: https://img.shields.io/github/last-commit/mkoeppe/cutgeneratingfunctionology/gh-pages?label=sphinx%20doc%20built
   :alt: Sphinx documentation built
   :target: http://mkoeppe.github.io/cutgeneratingfunctionology/doc/html/
.. image:: https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fmkoeppe%2Fcutgeneratingfunctionology
   :alt: Twitter
   :target: https://twitter.com/intent/tweet?text=%23cutgeneratingfunctionology:&url=https%3A%2F%2Fgithub.com%2Fmkoeppe%2Fcutgeneratingfunctionology


Authors
-------

See file `<AUTHORS.rst>`_ and also `<THANKS.rst>`_

License
-------

The code is released under the GNU General Public License, version 2,
or any later version as published by the Free Software Foundation. 

Documentation
-------------

http://mkoeppe.github.io/cutgeneratingfunctionology/doc/html/

Using the cutgeneratingfunctionology package
--------------------------------------------

There are many ways to run this package.

A. Run it online on mybinder.org
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://mybinder.org/badge_logo.svg
           :target: https://mybinder.org/v2/gh/mkoeppe/cutgeneratingfunctionology/master?filepath=demo.ipynb

B. Install released version from PyPI and run it within conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://img.shields.io/pypi/v/cutgeneratingfunctionology
   :alt: PyPI package
   :target: https://pypi.org/project/cutgeneratingfunctionology/

- Install Miniconda from https://docs.conda.io/en/latest/miniconda.html

- Set up the conda environment described in https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/environment.yml::

    curl https://raw.githubusercontent.com/mkoeppe/cutgeneratingfunctionology/master/environment.yml
    conda env create -n sage-cgf -f environment.yml
    conda activate sage-cgf

  This takes a while; it installs SageMath, which has many dependencies.

- Install PyPI package::

    pip install cutgeneratingfunctionology

- Start Sage.  You can either use the terminal (IPython)::

    sage

  or a Jupyter notebook::

    sage -n jupyter

- At the Sage prompt, type::

    import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

- Follow the instructions and examples in `<demo.rst>`_.


C.  Clone from GitHub and run it within conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Install Miniconda from https://docs.conda.io/en/latest/miniconda.html

- Clone the GitHub repository https://github.com/mkoeppe/cutgeneratingfunctionology.git::

    git clone https://github.com/mkoeppe/cutgeneratingfunctionology.git
    cd cutgeneratingfunctionology

- Set up the conda environment described in https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/environment.yml::

    conda env create -n sage-cgf -f environment.yml
    conda activate sage-cgf

  This takes a while; it installs SageMath which has many dependencies.

- (Optional:) Install the cutgeneratingfunctionology package using pip::

    pip install .

- Start Sage.  You can either use the terminal (IPython)::

    sage

  or a Jupyter notebook::

    sage -n jupyter

- At the Sage prompt, type::

    import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

- Follow the instructions and examples in `<demo.rst>`_.


D.  Run in a standalone installation of the SageMath distribution (no conda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Install the SageMath distribution:

   a) Either from source from http://www.sagemath.org/

   b) or with a binary from http://www.sagemath.org/

  The SageMath distribution brings its own installation of Python and many packages.

- Clone the GitHub repository https://github.com/mkoeppe/cutgeneratingfunctionology.git::

    git clone https://github.com/mkoeppe/cutgeneratingfunctionology.git
    cd cutgeneratingfunctionology

- (Optional:) Install optional SageMath distribution packages::

    sage -i lrslib pynormaliz

- Install the cutgeneratingfunctionology package using pip::

    sage -pip install -r requirements.txt

- Start SageMath.  You can either use the terminal (IPython)::

    sage

  or a Jupyter notebook::

    sage -n jupyter

- At the Sage prompt, type::

    import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

- Follow the instructions and examples in `<demo.rst>`_.
