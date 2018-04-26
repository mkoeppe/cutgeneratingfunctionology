#!/bin/sh
set -e
if [ "${SAGE_AGE}" = "-1" ]; then
    sudo sage -pip install --user --upgrade -v -i https://pypi.python.org/pypi sagemath
    sudo sage -pip install --upgrade --no-index -v .
    sudo sage ./setup.py test
    (cd docs && sudo sage -sh -c "make html")
    sudo sage -pip uninstall .
else
    ${HOME}/SageMath/sage -pip install --upgrade -v -i https://pypi.python.org/pypi sagemath
    ${HOME}/SageMath/sage -pip install --upgrade --no-index -v .
    ${HOME}/SageMath/sage setup.py test
    (cd docs && ${HOME}/SageMath/sage -sh -c "make html")
    ${HOME}/SageMath/sage -pip uninstall .
fi
