#!/usr/bin/env python3
from distutils.core import setup, Extension
import numpy
import os

NAME = "bootstrap_resample"
VERSION = 0.1
ext = '.c'

extensions = [Extension(NAME, [NAME+ext], extra_compile_args=['-O2'] )]


setup(
    name=NAME,
    description="Bootstrap resampling class for Pandas.DataFrames",
    version=VERSION,
    url="https://github.com/cancerevo/bootstrap_resample",
    author="Christopher McFarland",
    author_email="christopher.mcfarland@case.edu",
    license="MIT",
    classifiers=[
		'Programming Language :: Python :: 3',
		'Programming Language :: Cython',
        "Development Status :: 2 - Pre-Alpha",
		'License :: OSI Approved :: MIT License',
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=[NAME],
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "progressbar2",
#        "cython",
#        "tzlocal",
        "joblib",
    ],
    python_requires="~=3.2",
)
