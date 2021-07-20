# - Try to find CGAL: ALL LIBRARIES. 
# 
#  OUTPUT:
#  CGAL_FOUND        - system has CGAL
#  CGAL_INCLUDE_DIRS - include directories for CGAL
#  CGAL_LIBRARIES    - libraries for CGAL
#
#  Xiaozhe Hu
#  02/27/2013
#  Modified   2021-03-13   ludmil

# Find packages that CGAL depends on
#find_package(MPFR)

message(STATUS "Checking for packages in 'CGAL' and 'GMP'")


find_package(GMP) 
 
# Check for header file
find_path(CGAL_INCLUDE_DIRS CGAL/Triangulation_data_structure.h
 DOC "Directory where the CGAL headers are located"
 )
mark_as_advanced(CGAL_INCLUDE_DIRS)

# Try to run a test program that uses CGAL
# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CGAL
  "CGAL could not be found. Be sure to set CGAL_DIR correctly."
 CGAL_INCLUDE_DIRS)
