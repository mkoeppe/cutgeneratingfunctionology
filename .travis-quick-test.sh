#!/bin/sh
set -e
sage -pip -vvv install --upgrade-strategy only-if-needed -r requirements.txt || pip -vvv install --upgrade-strategy only-if-needed -r requirements.txt
make check
