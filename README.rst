.. image:: http://mkoeppe.github.io/cutgeneratingfunctionology/graphics-nonfree/Z11_058_github_template.jpg
   :width:  100%
   :target: https://github.com/mkoeppe/cutgeneratingfunctionology
   :alt:    mkoeppe/cutgeneratingfunctionology: Python code for computation and experimentation with cut-generating functions

Most of the code is for the 1-dimensional Gomory-Johnson infinite
group problem, including an electronic compendium of extreme
functions.

See http://www.sagemath.org/doc/tutorial/ for information on how to
use Sage.

.. badges

.. image:: https://img.shields.io/travis/mkoeppe/cutgeneratingfunctionology
   :alt: Travis CI
   :target: https://travis-ci.org/mkoeppe/cutgeneratingfunctionology/

.. image:: https://img.shields.io/pypi/l/cutgeneratingfunctionology
   :alt: License: GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
   :target: https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/COPYING

.. image:: https://img.shields.io/pypi/v/cutgeneratingfunctionology
   :alt: PyPI package
   :target: https://pypi.org/project/cutgeneratingfunctionology/

.. image:: https://mybinder.org/badge_logo.svg
   :alt: Run it on mybinder.org
   :target: https://mybinder.org/v2/gh/mkoeppe/cutgeneratingfunctionology/master?filepath=demo.ipynb

.. image:: https://img.shields.io/github/last-commit/mkoeppe/cutgeneratingfunctionology/gh-pages?label=sphinx%20doc%20built
   :alt: Sphinx documentation built
   :target: http://mkoeppe.github.io/cutgeneratingfunctionology/doc/html/

.. image:: https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fmkoeppe%2Fcutgeneratingfunctionology
   :alt: Twitter
   :target: https://twitter.com/intent/tweet?text=%23cutgeneratingfunctionology:&url=https%3A%2F%2Fgithub.com%2Fmkoeppe%2Fcutgeneratingfunctionology

.. add later: https://img.shields.io/pypi/pyversions/cutgeneratingfunctionology

Authors
-------

See https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/AUTHORS.rst and also https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/THANKS.rst

License
-------

The code is released under the GNU General Public License, version 2,
or any later version as published by the Free Software Foundation. 

Documentation
-------------

http://mkoeppe.github.io/cutgeneratingfunctionology/doc/html/

Using the cutgeneratingfunctionology package
--------------------------------------------
.. how_to_run

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

    e=environment.yml; curl -o $e https://raw.githubusercontent.com/mkoeppe/cutgeneratingfunctionology/master/$e
    conda env create -n sage-cgf -f $e
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

- Follow the instructions and examples in https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/demo.rst or https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/demo.ipynb .


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

- Follow the instructions and examples in https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/demo.rst or https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/demo.ipynb .


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

    sage -pip install .

- Start SageMath.  You can either use the terminal (IPython)::

    sage

  or a Jupyter notebook::

    sage -n jupyter

- At the Sage prompt, type::

    import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

- Follow the instructions and examples in https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/demo.rst or https://github.com/mkoeppe/cutgeneratingfunctionology/blob/master/demo.ipynb .
