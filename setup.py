from setuptools import setup

setup(
    packages = ['cutgeneratingfunctionology', 'cutgeneratingfunctionology.igp', 'cutgeneratingfunctionology.multirow', 'cutgeneratingfunctionology.dff', 'cutgeneratingfunctionology.spam', 'cutgeneratingfunctionology.igp.subadditivity_slack_diagrams', 'cutgeneratingfunctionology.igp.procedures'],
    include_package_data=True,     # to install the .sage files too
)
