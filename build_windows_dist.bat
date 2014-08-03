
setlocal EnableDelayedExpansion

call setenv.cmd /x64 /release

set DISTUTILS_USE_SDK=1

python setup.py build_ext --inplace

python setup.py bdist_wheel
