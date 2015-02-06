#
# Contains a CMakeLists.txt file that should work for standalone CART.
# Jay Sitaraman 11/22/2014
#
PROJECT(CART)
ENABLE_LANGUAGE(Fortran C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#
# Include a file that replaces functionality of Make.defs.
#
INCLUDE("standalone.cmake")
#
# Each directory includes a CMake include file that lists
# any source file which should be compiled in the directory.
# Dependencies between Fortran modules are computed automatically,
# so there is no need to build one directory before another.
#
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/lhs/lhs.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/rhs/rhs.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/preprocessing/preprocessing.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/postprocessing/postprocessing.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/bc/bc.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/analytics/analytics.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/partition/partition.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/interface/interface.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/driver/driver.cmake")
#
# Set the top level source directory for CART.
#
SET(CART_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}) 
#
# Set the build directory - where object files are place when compiled.
#
SET(CART_BUILD ${CMAKE_CURRENT_BINARY_DIR})
#
# create a default installation dir if not provided
#
IF (NOT DEFINED INSTALL_DIR)
 SET(INSTALL_DIR ${CART_BUILD}/cart)
ENDIF()
#
# Fortran module files will be built here, so make sure preprocessor can see
# them during compile.
#
INCLUDE_DIRECTORIES(BEFORE ${CART_BUILD})
#
# CART source files without driver but 
# includes interface
#
SET(CART_SRCS ${LHSSRCS} ${RHSSRCS} ${PRESRCS} ${POSTSRCS} ${BCSRCS} 
                ${ANALYTICS_SRCS} ${INTERFACESRCS} ${PARTITIONSRCS})
#
# Build the basic libCART.a one expects.
#
ADD_LIBRARY(CART STATIC ${CART_SRCS})
SET_TARGET_PROPERTIES(CART PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY ${CART_BUILD}/lib)
#
# Build stand alone executable linking the library
#
ADD_EXECUTABLE(cart.exe ${DRIVERSRCS})
TARGET_LINK_LIBRARIES(cart.exe CART)
#
# build the verificatino tests
#
ADD_EXECUTABLE(nscheck.exe tests/src/nscheck.f90)
TARGET_LINK_LIBRARIES(nscheck.exe CART)
ADD_EXECUTABLE(linearization_test.exe tests/src/linearization_test.f90)
TARGET_LINK_LIBRARIES(linearization_test.exe CART)
#
# include interface directory to 
# build the python shared object now
#
#ADD_SUBDIRECTORY(interface)
#
# install all the required files in 
# the installation directory
#
INSTALL(FILES ${CART_SOURCE}/README.md
  DESTINATION ${INSTALL_DIR}/README)
#INSTALL(FILES ${CART_SOURCE}/interface/python/cartInterface.py 
#  DESTINATION ${INSTALL_DIR}/bin)
#INSTALL(FILES ${CART_BUILD}/python/cart.so 
#  DESTINATION ${INSTALL_DIR}/bin)
INSTALL(DIRECTORY ${CART_SOURCE}/demo
  DESTINATION ${INSTALL_DIR})
INSTALL(DIRECTORY ${CART_BUILD}/lib
  DESTINATION ${INSTALL_DIR})
INSTALL(DIRECTORY ${CART_SOURCE}/tests/accuracy
	      DESTINATION ${INSTALL_DIR}/tests)
INSTALL(DIRECTORY ${CART_SOURCE}/tests/linearization
	      DESTINATION ${INSTALL_DIR}/tests)
INSTALL(PROGRAMS ${CART_BUILD}/cart.exe
  DESTINATION ${INSTALL_DIR}/bin)
INSTALL(PROGRAMS ${CART_BUILD}/nscheck.exe
  DESTINATION ${INSTALL_DIR}/bin)
INSTALL(PROGRAMS ${CART_BUILD}/linearization_test.exe
  DESTINATION ${INSTALL_DIR}/bin)
