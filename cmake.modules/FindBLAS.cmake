
# - Try to find BLAS
#  
#  OUTPUT:
#  BLAS_FOUND        - system has BLAS
#  BLAS_LIBRARIES    - libraries for BLAS
#
#  Xiaozhe Hu
#  02/27/2013
#  Modified   2015-08-08   --ltz
#  Modified   2020-05-29   --ltz

#message(STATUS "Checking for  'BLAS' library")
# Check for BLAS library
find_library(BLAS_LIBRARY blas
  HINTS ${BLAS_DIR} ${BLAS_DIR}/lib ${BLAS_DIR}/blas/lib $ENV{BLAS_DIR}/lib $ENV{BLAS_DIR}/BLAS/lib
  DOC "The BLAS library"
  )
mark_as_advanced(BLAS_LIBRARY)

# Collect libraries
if (BLAS_FOUND)
  set(BLAS_LIBRARIES ${BLAS_LIBRARIES})
endif()

# Try compiling and running test program
if (BLAS_LIBRARIES)
  # Set flags for building test program
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
endif(BLAS_LIBRARIES)


# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLAS
  "BLAS could not be found. Be sure to set BLAS_DIR correctly."
BLAS_LIBRARY)
