#!/bin/sh
set -e
sage -pip -vvv install --upgrade -r requirements.txt || pip -vvv install --upgrade -r requirements.txt
make check
