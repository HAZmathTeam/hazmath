
# - Try to find LAPACK
#  
#  OUTPUT:
#  LAPACK_FOUND        - system has LAPACK
#  LAPACK_LIBRARIES    - libraries for LAPACK
#
#  Xiaozhe Hu
#  02/27/2013
#  Modified   2015-08-08   --ltz
#  Modified   2020-05-29   --ltz

# Check for LAPACK library
#message(STATUS "Checking for  'BLAS' library")
find_library(LAPACK_LIBRARY lapack
  HINTS ${LAPACK_DIR} ${LAPACK_DIR}/lib ${LAPACK_DIR}/lapack/lib $ENV{LAPACK_DIR}/lib $ENV{LAPACK_DIR}/LAPACK/lib
  DOC "The LAPACK library"
  )
mark_as_advanced(LAPACK_LIBRARY)

# Collect libraries
if (LAPACK_FOUND)
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
endif()

# Try compiling and running test program
if (LAPACK_LIBRARIES)
  # Set flags for building test program
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
endif(LAPACK_LIBRARIES)


# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK
  "LAPACK could not be found. Be sure to set LAPACK_DIR correctly."
LAPACK_LIBRARY)
