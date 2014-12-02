#
# build the python shared object now
#
SET(CART_SRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(CART_INTERFACEFILE ${CART_SRC}/cartInterface.f90)
SET(CART_F2PYFILE ${CART_BUILD}/interface/cart.pyf)
SET(CART_FWRAPFILE ${CART_BUILD}/interface/cart-f2pywrappers2.f90)
SET(CART_CWRAPFILE ${CART_BUILD}/interface/cartmodule.c)
SET(NUMPY_SRCFILE ${NUMPY_SRC}/fortranobject.c)
#
# All of these set in standalone.cmake automatically.  PYTHON_INCLUDE_DIRS is
# found through FIND_PACKAGE(PythonLibs).
#
include_directories(AFTER ${NUMPY_SRC} ${NUMPY_INCLUDE} ${PYTHON_INCLUDE})
#
# Include the build directory for this source directory.
#
include_directories(BEFORE ${CART_BUILD}/interface)
#
# add custom commands to generate the .pyf file from cartInterface.f90
#
ADD_CUSTOM_COMMAND(OUTPUT ${CART_F2PYFILE} 
                   COMMAND ${F2PYBIN} -m cart --overwrite-signature -h ${CART_F2PYFILE} ${CART_INTERFACEFILE}
                   MAIN_DEPENDENCY ${CART_INTERFACEFILE})
#
# add custom command to build the C and f90 source files that
# come out of the f2py wrapper
#
ADD_CUSTOM_COMMAND(OUTPUT ${CART_FWRAPFILE} ${CART_CWRAPFILE}
                   COMMAND ${F2PYBIN} ARGS ${CART_F2PYFILE} 
                   MAIN_DEPENDENCY ${CART_F2PYFILE})
#
# build the cart.so shared object
#
ADD_LIBRARY(cart SHARED ${CART_FWRAPFILE} ${CART_CWRAPFILE} 
                          ${NUMPY_SRCFILE})
SET_TARGET_PROPERTIES(cart PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${CART_BUILD}/python)
#
# Make shared object link against Python and libCART.a
#
TARGET_LINK_LIBRARIES(cart CART ${PYTHON_LIB_FLAGS})
#
# Sometimes linker has issues finding Fortran runtimes unless Fortran compiler
# invokes the linker directly, so make it call the linker instead of C.
#
SET_TARGET_PROPERTIES(cart PROPERTIES LINKER_LANGUAGE Fortran)
#
# Make the format cart.so for the output library.  On OS X we want to force
# the so extension instead of dylib.
#
SET_TARGET_PROPERTIES(cart PROPERTIES PREFIX "")
SET_TARGET_PROPERTIES(cart PROPERTIES SUFFIX ".so")
