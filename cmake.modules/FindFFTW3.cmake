# - Try to find FFTW3
#  
#  OUTPUT:
#  FFTW3_FOUND        - system has FFTW3
#  FFTW3_LIBRARIES    - libraries for FFTW3
#
#  Ludmil
#  20230201

# Check for FFTW3 library
message(STATUS "Checking for  'FFTW3' library")
  find_library(FFTW3_LIBRARY ${FFTW3_LIBRARY_NAME}
  )
mark_as_advanced(FFTW3_LIBRARY)

# Collect libraries
if (FFTW3_FOUND)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES})
endif()

# Try compiling and running test program
if (FFTW3_LIBRARIES)
  # Set flags for building test program
  set(CMAKE_REQUIRED_LIBRARIES ${FFTW3_LIBRARIES})
endif(FFTW3_LIBRARIES)


# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
  "FFTW3 could not be found. Be sure to set FFTW3_DIR correctly."
FFTW3_LIBRARY)
