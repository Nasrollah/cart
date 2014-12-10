# If not part of a larger project, use this file for basic compile settings.
# The most fundamental flags and definitions are included below.

# If unspecified, use the Intel compiler.
IF(NOT DEFINED COMPILER_FAMILY)
    SET(COMPILER_FAMILY "gnu")
ENDIF()

IF(COMPILER_FAMILY STREQUAL "intel")
    message("Using Intel compiler family.")
    SET(CMAKE_C_FLAGS "-fPIC -O2")
    SET(CMAKE_Fortran_FLAGS "-CB -traceback -r8 -fPIC -convert big_endian -fp-model source -vec-report0")
    #SET(CMAKE_Fortran_FLAGS "-r8 -fPIC -convert big_endian -CB -traceback")
    #-check uninit -traceback -fpe0 -vec-report0 -check pointers -check uninit")
    #SET(CMAKE_Fortran_FLAGS "-r8 -fPIC -DINTEL -O2 -fpp -convert big_endian -align sequence -heap-arrays 10 -vec-report0 -fpe0 -traceback  -check bounds -check pointers -check uninit -warn interfaces -gen-interfaces nosource")
ELSEIF(COMPILER_FAMILY STREQUAL "portland")
    message("Using Portland compiler family.")
    SET(CMAKE_C_FLAGS "-O3 -fPIC")
    SET(CMAKE_Fortran_FLAGS "-r8 -O3 -fPIC -byteswapio")
ELSEIF(COMPILER_FAMILY STREQUAL "gnu")
    message("Using GNU compiler family.")
    # Setting big_endian doesn't affect shared objects for GNU, make sure youi
    # use a CONVERT flag in any OPEN you perform in Fortran.
    SET(CMAKE_C_FLAGS "-fPIC -O3")
    SET(CMAKE_Fortran_FLAGS "-fPIC -fdefault-real-8 -O3")
ELSE()
    message( FATAL_ERROR "Unsupported compiler set, exiting.")
ENDIF()

# Find Python, and then find necessary f2py and NumPy executables, libraries,
# and include directories.
FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(PythonInterp REQUIRED)

SET(PYVER ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

#STRING(REGEX REPLACE "bin.*" "bin/f2py"
#                             F2PYBIN ${PYTHON_EXECUTABLE})
#STRING(REGEX REPLACE "bin.*" "lib/python${PYVER}/site-packages/numpy/f2py/src"
#                             NUMPY_SRC ${PYTHON_EXECUTABLE})
set(F2PYBIN f2py)
set(NUMPY_SRC      /glade/apps/opt/numpy/1.7.1/intel-autodispatch/13.0.1/lib/python2.7/site-packages/numpy/f2py/src/)
set(NUMPY_INCLUDE  /glade/apps/opt/numpy/1.7.1/intel-autodispatch/13.0.1/lib/python2.7/site-packages/numpy/core/include/)
set(PYTHON_INCLUDE /glade/apps/opt/python/2.7.5/gnu-westmere/4.7.2/include/python2.7 )
#STRING(REGEX REPLACE "bin.*"
#                     "lib/python${PYVER}/site-packages/numpy/core/include"
#                     NUMPY_INCLUDE ${PYTHON_EXECUTABLE})

#set(PYVERSION 2.7.5)
#set(PTOOLSRTE /share/apps/hiarms/ptoolsrte-0.5.1rc4)
#set (PYTHON_INCLUDE ${PTOOLSRTE}/packages/Python-${PYVERSION}/include/python2.7)


message("F2PYBIN: ${F2PYBIN}")
message("NUMPY_SRC: ${NUMPY_SRC}")
message("NUMPY_INCLUDE: ${NUMPY_INCLUDE}")
message("PYTHON_INCLUDE: ${PYTHON_INCLUDE}")
