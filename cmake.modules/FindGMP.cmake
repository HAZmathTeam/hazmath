# - Try to find GMP
# 
#  OUTPUT:
#  GMP_FOUND        - system has GMP
#  GMP_INCLUDE_DIRS - include directories for GMP
#  GMP_LIBRARIES    - libraries for GMP
#
#  Xiaozhe Hu
#  02/27/2013
#  Modified   2015-08-08   --ltz
#
message(STATUS "Checking for package 'GMP'")

# Check for header file
#find_path(GMP_INCLUDE_DIRS gmp.h
#  HINTS /usr/include/CGAL
#  DOC "Directory where the GMP header is located"
#  )
#mark_as_advanced(GMP_INCLUDE_DIRS)

# Check for GMP library
find_library(GMP_LIBRARIES gmp
  DOC "The GMP library"
  )
mark_as_advanced(GMP_LIBRARIES)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  "GMP could not be found. Be sure to set SUITESPARSE_DIR correctly."
  GMP_LIBRARIES)
