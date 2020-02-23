set -e
CFLAGS=-O0 python setup.py build_ext --inplace
#cp build/lib.linux-x86_64-3.6/extender.cpython-36m-x86_64-linux-gnu.so .
#pyinstaller extender_wrapper.spec

mv extender.cpython-36m-x86_64-linux-gnu.so extender/extender.cpython-36m-x86_64-linux-gnu.so
echo "Build successful: output dist/extender_wrapper"
