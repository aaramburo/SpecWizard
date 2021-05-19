from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
ext_modules = cythonize("projectdata.pyx"),
include_dirs=[numpy.get_include()]
)




# setup(
# ext_modules = cythonize("computeib.pyx"),
# include_dirs=[numpy.get_include()]
# )
