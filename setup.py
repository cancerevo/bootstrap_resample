#!/usr/bin/env python3
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
import os

setup(
    name="bootstrap_resampling",
    description="Bootstrap resampling class for Pandas.DataFrames",
    url="https://github.com/cancerevo/bootstrap_resample",
    author="Christopher McFarland",
    author_email="christopher.mcfarland@case.edu",
    license="MIT",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=["bootstrap_resampling"],
    ext_modules=cythonize(Extension("*", ["bootstrap_resample/*.pyx"], include_dirs=[numpy.get_include()])),
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "progressbar2",
        "cython",
        "tzlocal",
        "joblib",
    ],
    python_requires="~=3.2",
    version="0.2",
)
