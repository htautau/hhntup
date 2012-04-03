#!/usr/bin/env python

from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from glob import glob

from distutils.core import Extension

jet_cleaning = Extension('_libcleaning',
                    sources = ['_libcleaning.cpp',
                               '_cleaning.cpp'])

setup(name='higgstautau',
      author='Noel Dawe',
      author_email='noel.dawe@cern.ch',
      packages=find_packages(),
      ext_modules = [jet_cleaning]
     )
