To Build:
--------
Let BUILD be the directory and SRC be the source directory.
Its best if BUILD is not same as SRC

Adjust following variables in $SRC/standalone.cmake
Default are given below. Cmake unfortunately is unable to find this or
more likely I am not skilled enough to make it find these automatically.

set(NUMPY_SRC /usr/local/lib/python2.7/site-packages/numpy/f2py/src/)
set(NUMPY_INCLUDE /usr/local/lib/python2.7/site-packages/numpy/core/include/)
set(PYTHON_INCLUDE /usr/local/include/python2.7 )

Then do following commands:

mkdir BUILD;
export CC=mpicc;export CXX=mpicxx;export FC=mpif90 (or equivalent in other shells)
cd BUILD;
cmake -DCOMPILER_FAMILY=<gnu,intel,portland> -DINSTALL_DIR=<installation_directory> -DMACHINE=<yellowstone,norgay,linux> SRC
make -j4
make install 

Here is the example I used from my csh_history
 2225  mkdir BUILD
 2226  cd BUILD/
 2227  setenv CC=mpicc;setenv CXX=mpicxx;setenv FC=mpif90;
 2228  cmake -DCOMPILER_FAMILY=intel -DINSTALL_DIR=/home/jay/code/cart -DMACHINE=yellowstone..
 2229  make -j4
 2230  make install

Note:
If COMPILER_FAMILY is not specified it is assumed as "gnu"
If INSTALL_DIR is not specified it is assumed as BUILD/cart
IF MACHINE is not specified it is assumed as yellowstone (for now - this will hopefully be automated shortly)