#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:09:32 2017

@author: mattjohnson
"""

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

#setup(
#    ext_modules=[
#        Extension("lines.line", ["lines/line.pyx"],
#                  include_dirs=[numpy.get_include()]),
#    ],
#)

# Or, if you use cythonize() to make the ext_modules list,
# include_dirs can be passed to setup()

setup(
    ext_modules=cythonize("*.pyx"),
    include_dirs=[numpy.get_include()]
) 