#! /bin/bash -x
#
# Best if sourced! It sets the PATH where it installs sage.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e

if [ -n "${SAGE_AGE}" ]; then
    ##
    ## Install a Sage binary.
    ##
  SAGE_IMAGE=`python2 -c "import sage_version; print sage_version.get_all_version_names('${SAGE_SERVER}',${SAGE_AGE})"`
  save_dir=`pwd`
  cd $HOME
  echo "Obtaining Sage image:" ${SAGE_IMAGE}
  if [ ! -x SageMath/sage ] ; then
      rm -f SageMath.tar.bz2
      wget --progress=dot:giga ${SAGE_SERVER}${SAGE_IMAGE} -O SageMath.tar.bz2
      tar xf SageMath.tar.bz2
  fi
  # Disable recompiles of sagelib after installing packages, which times out on Travis CI
  sed -i.bak $'s/^sage:/sage:\\\nrebuild-sage-lib:/' "$HOME/SageMath/src/Makefile"
  sed -i.bak $'s/^sage:/sage:\\\nrebuild-sage-lib:/' "$HOME/SageMath/src/Makefile.in"
  # Back to the correct directory.
  cd "$save_dir"
  export PATH="$HOME/SageMath/:$PATH"
  MAKE="make -j4"
  export MAKE
  # Install packages
  sage -i lrslib pynormaliz
  # To initialize matplotlib font manager
  sage -python -c 'import matplotlib.pyplot'
  sage -pip install --user sphinxcontrib-websupport
fi
# Display the banner, so we can be sure what version we ended up with!
sage < /dev/null
