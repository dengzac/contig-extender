#!/usr/bin/python
from setuptools import setup, find_packages
from Cython.Compiler import Options
import subprocess
Options.embed = 'main'
from Cython.Build import cythonize
import numpy as np
ext = cythonize('extender/extender.pyx')

flags = subprocess.run(['python3-config', '--cflags', '--ldflags'], stdout=subprocess.PIPE).stdout.decode()
flags = flags.strip().replace('\n', ' -I' + np.get_include() + ' extender/extender.c ') + ' -fPIC -o extender/extend'
# print(flags)
p  =subprocess.run(['gcc ' + flags ], shell=True)

setup(
    name = "contig-extender",
    version="0.0.1",
    ext_modules = ext, 
    include_dirs = [np.get_include()],
    install_requires = [
        'numpy',
        'cython'
    ]
)
# print(p.args)
