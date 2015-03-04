import os
import git
from setuptools import setup, find_packages

import seqpy

#repo = git.Repo(os.path.dirname(os.path.realpath(__file__)))

README = open('README.rst').read()

setup(
    #version = repo.commit('master').id,
    version = 0.1,
    name='seqpy',
    packages = find_packages(),
    description = "Tools for analyzing high-throughput sequencing data",
    author='Charles Murphy',
    license='Non-commercial',
    long_description=README
)
