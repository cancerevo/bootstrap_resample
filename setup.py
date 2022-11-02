#!/usr/bin/env python3
from distutils.core import setup, Extension
import numpy
import os

NAME = "bootstrap_resample"

if NAME+'.c'

if RE_CYTHONIZE:
    from Cython.Build import cythonize

# Delete fast_prng.c to re-cythonize code, this script will automatically regenerated fast_prng.c using your cython package.
USE_CYTHON = False if Name+'.c' in os.listdir(os.getcwd()) else True

ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension(Name, [Name+ext], extra_compile_args=['-O2'] )]

if USE_CYTHON:
	from Cython.Build import cythonize
	extensions = cythonize(extensions)


with open(Name+'.pyx', 'r') as f:
    for line in f:
        if '__version__' in line:
            exec(line)
setup(
    name=name,
    description="Bootstrap resampling class for Pandas.DataFrames",
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
    packages=["bootstrap_resampling"],
    ext_modules=cythonize(Extension("*", [f"{NAME}/*.pyx"], include_dirs=[numpy.get_include()])),
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
#        "matplotlib",
        "progressbar2",
#        "cython",
#        "tzlocal",
        "joblib",
    ],
    python_requires="~=3.2",
    version="0.2",
)
