import os
from setuptools import setup, find_packages

README = open('README.rst').read()

setup(
    #version = repo.commit('master').id,
    version = 0.1,
    name='seqpy',
    packages = find_packages(),
    description = "Python and R tools for analyzing high-throughput sequencing data",
    author='Charles Murphy',
    license='Non-commercial',
    long_description=README
)
