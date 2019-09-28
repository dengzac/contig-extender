#!/usr/bin/python
from setuptools import setup, find_packages, Extension
from Cython.Compiler import Options
import subprocess
Options.embed = 'main'
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np


# flags = subprocess.run(['python3-config', '--cflags', '--ldflags'], stdout=subprocess.PIPE).stdout.decode()
# flags = flags.strip().replace('\n', ' -I' + np.get_include() + ' extender/extender.c ') + ' -fPIC -o extender/extend'
# print(flags)
# p  =subprocess.run(['gcc ' + flags ], shell=True)

setup(
    name = "contig-extender",
    version="0.0.1",
    ext_modules = [Extension("extender", ["extender/extender.pyx"])], 
    include_dirs = [np.get_include()],
    cmdclass={'build_ext': build_ext},
    scripts=['extender/extender_wrapper.py'],
    install_requires = [
        'numpy',
        'cython',
        'psutil'
    ]
)
# print(p.args)
