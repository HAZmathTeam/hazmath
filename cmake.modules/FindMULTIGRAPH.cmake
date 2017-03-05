# - Try to find MULTIGRAPH source: 
# 
#  OUTPUT:
#  MULTIGRAPH_FOUND        - system has MULTIGRAPH
#  MULTIGRAPH_SOURCE_DIR - source directory for MULTIGRAPH
#
#  Ludmil 2017-03-04
#
# Check for header file
###      message("ZZZZZZZZZZZZZZZZZZ: ${MULTIGRAPH_DIR} ${USE_MULTIGRAPH}  ZZZ")
find_path(MULTIGRAPH_SOURCE_DIR solver.f
 HINTS ${MULTIGRAPH_DIR}/src
 PATH_SUFFIXES multigraph 
 DOC "Directory where the MULTIGRAPH sources are located"
 )
mark_as_advanced(MULTIGRAPH_SOURCE_DIR)


# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MULTIGRAPH
  "MULTIGRAPH could not be found. Be sure to set MULTIGRAPH_DIR correctly."
  MULTIGRAPH_SOURCE_DIR)


