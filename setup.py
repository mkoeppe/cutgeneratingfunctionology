## -*- encoding: utf-8 -*-
import os
import sys
from setuptools import setup
from codecs import open # To open the README file with proper encoding
from setuptools.command.test import test as TestCommand # for tests


# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

# For the tests
class SageTest(TestCommand):
    def run_tests(self):
        errno = os.system("make check")
        if errno != 0:
            sys.exit(1)

setup(
    packages = ['cutgeneratingfunctionology', 'cutgeneratingfunctionology.igp', 'cutgeneratingfunctionology.multirow', 'cutgeneratingfunctionology.dff', 'cutgeneratingfunctionology.spam', 'cutgeneratingfunctionology.igp.subadditivity_slack_diagrams', 'cutgeneratingfunctionology.igp.procedures'],
    include_package_data=True,     # to install the .sage files too
    cmdclass = {'test': SageTest}, # adding a special setup command for tests
)
