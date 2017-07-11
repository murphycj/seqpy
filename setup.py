import os
from setuptools import setup, find_packages

README = open('README.md').read()

setup(
    version = 0.1,
    name='seqpy',
    packages = find_packages(),
    description = "Python and R tools for analyzing high-throughput sequencing data",
    author='Charles Murphy',
    include_package_data=True,
    license='Non-commercial',
    long_description=README
)
