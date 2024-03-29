# - Try to find SUITESPARSE: ALL LIBRARIES. 
# 
#  OUTPUT:
#  SUITESPARSE_FOUND        - system has SUITESPARSE
#  SUITESPARSE_INCLUDE_DIRS - include directories for SUITESPARSE
#  SUITESPARSE_LIBRARIES    - libraries for SUITESPARSE
#
#  Xiaozhe Hu
#  02/27/2013
#  Modified   2015-08-08   ludmil
#  Modified   2021-07-20   ludmil

# Find packages that SUITESPARSE depends on
find_package(BLAS)
find_package(LAPACK) 
find_package(SUITESPARSE_CONFIG)
find_package(AMD)
find_package(UMFPACK)
find_package(COLAMD)
find_package(CAMD)
find_package(CCOLAMD)
#find_package(METIS)
if(NOT APPLE) 
 find_package(RT)
endif(NOT APPLE)
#find_package(BLAS)
#find_package(LAPACK)
#find_package(UMFPACK)

message(STATUS "Checking for packages in 'SUITESPARSE'")

# Set initial values
set(SUITESPARSE_INCLUDE_DIRS "")
set(SUITESPARSE_LIBRARIES "")
mark_as_advanced(SUITESPARSE_INCLUDE_DIRS)

# Collect paths
if (AMD_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${AMD_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARY} ${AMD_LIBRARIES})
endif()
if (SUITESPARSE_CONFIG_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${SUITESPARSE_CONFIG_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${SUITESPARSE_CONFIG_LIBRARIES})
endif()
if (UMFPACK_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${UMFPACK_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARY} ${UMFPACK_LIBRARIES})
endif()
if (COLAMD_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${COLAMD_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${COLAMD_LIBRARIES})
endif()
if (CAMD_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${CAMD_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${CAMD_LIBRARIES})
endif()
if (CCOLAMD_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${CCOLAMD_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${CCOLAMD_LIBRARIES})
endif()
if (RT_FOUND AND (NOT APPLE))
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${RT_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${RT_LIBRARIES})
endif()
if (METIS_FOUND)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${METIS_INCLUDE_DIRS})
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${METIS_LIBRARIES})
endif()

if (BLAS_FOUND)
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${BLAS_LIBRARIES})
endif()
if (LAPACK_FOUND)
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

#set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} "-lgfortran")

mark_as_advanced(SUITESPARSE_INCLUDE_DIRS)
mark_as_advanced(SUITESPARSE_LIBRARIES)

# Try to run a test program that uses SUITESPARSE
# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE
  "SUITESPARSE could not be found. Be sure to set SUITESPARSE_DIR correctly."
 SUITESPARSE_LIBRARIES SUITESPARSE_INCLUDE_DIRS)
