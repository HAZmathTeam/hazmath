cmake_minimum_required(VERSION 3.21)

project(haznicsx LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE NB_DIR
)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

if(TARGET hazmath::hazmath)
  # built together with libhazmath ("pip install ." from the repository root)
  message(STATUS "Using libhazmath from this build")
else()
  # libhazmath, installed with -DHAZMATH_RELOCATABLE_INSTALL=ON
  find_package(hazmath 0.1 CONFIG REQUIRED)
  message(STATUS "Found HAZmath at ${hazmath_DIR}")
endif()

# The preconditioner builders (create_precond_*) are in src/haznics/helper.c of
# the HAZmath source tree, which is not part of libhazmath
set(
  HAZMATH_SOURCE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  CACHE PATH
  "HAZmath source tree, for src/haznics/helper.c"
)
set(HAZNICS_HELPER "${HAZMATH_SOURCE_DIR}/src/haznics/helper.c")
set(HAZNICS_HEADER_AWK "${HAZMATH_SOURCE_DIR}/haz_shutils/mkheaders_simple.awk")
if(NOT EXISTS "${HAZNICS_HELPER}")
  message(FATAL_ERROR "${HAZNICS_HELPER} not found; set HAZMATH_SOURCE_DIR")
endif()

# Prototypes of helper.c, generated as for the SWIG module
# (haz_shutils/headmk_haznics.sh)
find_program(AWK awk REQUIRED)
execute_process(
  COMMAND "${AWK}" -v name=haznics_add.h -f "${HAZNICS_HEADER_AWK}"
  INPUT_FILE "${HAZNICS_HELPER}"
  OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/haznics_add.h"
  RESULT_VARIABLE HAZNICS_HEADER_RESULT
)
if(NOT HAZNICS_HEADER_RESULT EQUAL 0)
  message(FATAL_ERROR "Generating haznics_add.h failed")
endif()
set_property(
  DIRECTORY
  APPEND
  PROPERTY CMAKE_CONFIGURE_DEPENDS "${HAZNICS_HELPER}" "${HAZNICS_HEADER_AWK}"
)

add_library(haznicsx_helper OBJECT "${HAZNICS_HELPER}")
set_target_properties(haznicsx_helper PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(haznicsx_helper PUBLIC hazmath::hazmath)

nanobind_add_module(_haznicsx src/haznicsx.cpp)
target_include_directories(_haznicsx PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(_haznicsx PRIVATE haznicsx_helper hazmath::hazmath)
# find a shared libhazmath where it was installed. INSTALL_RPATH_USE_LINK_PATH
# is not enough: it skips directories inside the project, such as a prefix in
# the repository
get_target_property(HAZMATH_TYPE hazmath::hazmath TYPE)
if(HAZMATH_TYPE STREQUAL "SHARED_LIBRARY")
  set_target_properties(
    _haznicsx
    PROPERTIES INSTALL_RPATH "$<TARGET_FILE_DIR:hazmath::hazmath>"
  )
endif()

install(TARGETS _haznicsx LIBRARY DESTINATION haznicsx)
