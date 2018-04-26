#!/bin/sh
set -e
make check-long SAGE="$HOME/SageMath/sage" SAGE_CHECK_FLAGS="--verbose"
make doc SAGE="$HOME/SageMath/sage"
