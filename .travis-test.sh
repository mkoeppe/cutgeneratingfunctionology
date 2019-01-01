#!/bin/sh
set -e
$HOME/SageMath/sage -pip install --upgrade --no-index -v .
make check-long SAGE="$HOME/SageMath/sage" #SAGE_CHECK_FLAGS="--verbose"
make doc SAGE="$HOME/SageMath/sage"
