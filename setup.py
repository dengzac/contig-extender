#!/usr/bin/python
from setuptools import setup, find_packages, Extension
from Cython.Compiler import Options
import subprocess
Options.embed = 'main'
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np
import sysconfig, os
class NoSuffixBuilder(build_ext):
    def get_ext_filename(self, ext_name):
        filename = super().get_ext_filename(ext_name)
        suffix = sysconfig.get_config_var('EXT_SUFFIX')
        ext = os.path.splitext(filename)[1]
        return filename.replace(suffix, "") + ext
# flags = subprocess.run(['python3-config', '--cflags', '--ldflags'], stdout=subprocess.PIPE).stdout.decode()
# flags = flags.strip().replace('\n', ' -I' + np.get_include() + ' extender/extender.c ') + ' -fPIC -o extender/extend'
# print(flags)
# p  =subprocess.run(['gcc ' + flags ], shell=True)

setup(
    name = "contig-extender",
    version="0.0.1",
    ext_modules = [Extension("extender", ["extender/extender.pyx"])],
    include_dirs = [np.get_include()],
    cmdclass={'build_ext': NoSuffixBuilder},
    scripts=['extender/extender_wrapper.py'],
    install_requires = [
        'numpy',
        'cython',
        'psutil'
    ]
)
# print(p.args)
