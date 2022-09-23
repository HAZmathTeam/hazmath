
# - Try to find HDF5
#  
#  OUTPUT:
#  HDF5_FOUND        - system has HDF5
#  HDF5_LIBRARIES    - libraries for HDF5
#
#  Xiaozhe Hu
#  02/27/2013
#  Modified   2015-08-08   --ltz
#  Modified   2020-05-29   --ltz

# Check for HDF5 library
message(STATUS "Checking for  'HDF5' library")
  find_library(HDF5_LIBRARY ${HDF5_LIBRARY_NAME}
  )
  #HINTS ${HDF5_DIR} ${HDF5_DIR}/lib ${HDF5_DIR}/hdf5/serial ${HDF5_DIR}/hdf5/openmpi DOC "The HDF5 library"
mark_as_advanced(HDF5_LIBRARY)

# Collect libraries
if (HDF5_FOUND)
  set(HDF5_LIBRARIES ${HDF5_LIBRARIES})
endif()

# Try compiling and running test program
if (HDF5_LIBRARIES)
  # Set flags for building test program
  set(CMAKE_REQUIRED_LIBRARIES ${HDF5_LIBRARIES})
endif(HDF5_LIBRARIES)


# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5
  "HDF5 could not be found. Be sure to set HDF5_DIR correctly."
HDF5_LIBRARY)
