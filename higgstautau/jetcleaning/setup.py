#!/usr/bin/env python

from setuptools import setup, find_packages, Extension

jet_cleaning = Extension('_libcleaning',
    sources = [
        '_libcleaning.cpp',
        '_cleaning.cpp'])

setup(name='higgstautau',
      author='Noel Dawe',
      author_email='noel.dawe@cern.ch',
      packages=find_packages(),
      ext_modules = [jet_cleaning])
