del /f "extender\extender.c"
python setup.py build_ext --inplace
pyinstaller extender_wrapper.spec
echo "Build successful: output dist/extender_wrapper"
