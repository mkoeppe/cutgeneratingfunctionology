#! /bin/sh
set -e
make check SAGE="$HOME/SageMath/sage"
make doc SAGE="$HOME/SageMath/sage"
