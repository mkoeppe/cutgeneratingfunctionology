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

Run it on mybinder.org
----------------------

.. image:: https://mybinder.org/badge_logo.svg
           :target: https://mybinder.org/v2/gh/mkoeppe/cutgeneratingfunctionology/master?filepath=demo.ipynb


How to run the code in a local copy of Sage
-------------------------------------------

1. Download the code from
   https://github.com/mkoeppe/cutgeneratingfunctionology.git

2. Install SageMath:

   a) Either from source from http://www.sagemath.org/

   b) or with a binary from http://www.sagemath.org/

   c) or via conda (recommended):

      - Install Miniconda from https://docs.conda.io/en/latest/miniconda.html

      - From the directory "cutgeneratingfunctionology", set up the environment::

          conda env create -n sage-cgf -f environment.yml
          conda activate sage-cgf

3. From the directory "cutgeneratingfunctionology", start
   Sage.  You can either use the terminal (IPython)::

        sage

   or a Jupyter notebook::

        sage -n jupyter

4. At the Sage prompt, type::

    import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

5. Follow the instructions and examples in `<demo.rst>`_.


How to run the code online via cloud.sagemath.com
-------------------------------------------------

1. Create a user account at https://cloud.sagemath.com

2. Log in at https://cloud.sagemath.com

3. Create a new project "Group relaxation" (or any name)

4. Open the project

5. Create a directory: 
   Paste in the weblink: https://github.com/mkoeppe/cutgeneratingfunctionology.git
   and hit enter

6. Enter that directory

7. Click "+ New", select "Sage worksheet"

8. Type::

    import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

   and hit shift+enter

9. Follow the instructions and examples in `<demo.sage>`_.


To update the code to the latest version:

1. In the project "Group relaxation", open the directory "cutgeneratingfunctionology".
   
2. In the line "Terminal command...", enter::
     
    git pull 


