#!/usr/bin/env python3

from distutils.core import Extension
from Cython.Build import cythonize
import numpy

NAME = "bootstrap_resample"

extensions = [Extension(NAME, [f'{NAME}/{NAME}.pyx'], extra_compile_args=['-O2'] )]
cythonize(extensions, compiler_directives={'language_level' : "3"})
