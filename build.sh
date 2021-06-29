set -e
CFLAGS=-O0 python3 setup.py build_ext --inplace
#cp build/lib.linux-x86_64-3.6/extender.cpython-36m-x86_64-linux-gnu.so .
mv extender.so extender/extender.so

#python3 -m PyInstaller extender_wrapper.spec

echo "Build successful: output dist/extender_wrapper"
