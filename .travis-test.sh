#! /bin/sh
set -e
make check-long SAGE="$HOME/SageMath/sage"
make doc SAGE="$HOME/SageMath/sage"
