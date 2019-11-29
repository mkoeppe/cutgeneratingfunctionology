#!/bin/sh
set -e
sage -pip -vvv install --upgrade . || pip -vvv install --upgrade .
make check
