from distutils.core import setup, Extension
import numpy as np

# define the extension module
generalied_chi_squared = Extension(
    'generalized_chi_squared',
    include_dirs=[np.get_include()],
    sources=['generalized_chi_squared.cpp'],
    language='c++',
    extra_compile_args=['-std=c++11', '-Wunused-function']
    )

# run the setup
setup(name='generalized_chi_squared', ext_modules=[generalied_chi_squared])

