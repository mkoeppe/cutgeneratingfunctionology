#!/bin/sh
set -e
export PATH="$HOME/SageMath/:$PATH"
sage -pip -vvv install --upgrade . || pip -vvv install --upgrade .
make check-long #SAGE_CHECK_FLAGS="--verbose"
make doc
