#! /bin/bash -x
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e

if [ -n "${SAGE_AGE}" ]; then
  SAGE_IMAGE=`python2 -c "import sage_version; print sage_version.get_all_version_names('${SAGE_SERVER}',${SAGE_AGE})"`
  cd $HOME
  echo "Obtaining Sage image:" ${SAGE_IMAGE}
  if [ ! -x SageMath/sage ] ; then
      rm -f SageMath.tar.bz2
      wget --progress=dot:giga ${SAGE_SERVER}${SAGE_IMAGE} -O SageMath.tar.bz2
      tar xf SageMath.tar.bz2
  fi
  MAKE="make -j4"
  export MAKE
  # Install packages
  SageMath/sage -i lrslib pynormaliz
  # To initialize matplotlib font manager
  $HOME/SageMath/sage -python -c 'import matplotlib.pyplot'
  $HOME/SageMath/sage -pip install --user sphinxcontrib-websupport
else
    #https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#using-pip-in-an-environment
    pip -i pynormaliz
fi
