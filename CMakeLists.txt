############################################################
# CMakeLists for HAZMATH
# No Fortran compiler check. 
#
#
# Modified   2015-08-08   --ltz
# 	     2017-03-05   --ltz
#	     2017-03-10   — XHu
#	     2020-05-29   --ltz
#	     2021-03-13   --ltz
#	     2021-05-05   --ltz
#	     2021-07-18   --ltz
#	     2022-09-23   --ltz
#
############################################################
cmake_minimum_required (VERSION 3.12)

# "Delete Yourself" Prevention
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(FATAL_ERROR 
    "------------------------------------------------
Source and build directories cannot be the same.  
Please REMOVE CMakeCache.txt and CMakeFiles.
TYPE 'make help' for more information. 
-------------------------------------------------"
    )
endif()

# include helper modules
include(CheckFunctionExists)
include(CheckIncludeFile)
###########################policies
if(${CMAKE_VERSION} VERSION_LESS "3.13.0") 
  message(WARNING "Setting cmake policies CMP0078 and CMP0086 to OLD")
  cmake_policy(SET CMP0078 OLD)
  cmake_policy(SET CMP0086 OLD)  
else()
  cmake_policy(SET CMP0078 NEW)
  cmake_policy(SET CMP0086 NEW)  
endif()
####################################
set(CMAKE_VERBOSE_MAKEFILE 1) 
set(GDB 1 CACHE BOOL "debugging or not")
set(OPENMP 0 CACHE BOOL "Openmp use")

# Search for C compilers in the specified order. That will determine the rest.
if(DEFINED ENV{CC}) 
  find_program(THE_C NAMES $ENV{CC} gcc mpicc gcc6 gcc-6 gcc5 gcc-mp-6 gcc-mp-5 gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 icc clang)
else(DEFINED ENV{CC}) 
  find_program(THE_C NAMES gcc mpicc gcc6 gcc-6 gcc5 gcc-mp-6 gcc-mp-5 gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 icc clang)
endif(DEFINED ENV{CC}) 

if( THE_C )
  set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C compiler" FORCE)
else( THE_C )
  message(WARNING " the requested CC is not found: ${THE_C}; continuing with the system's one" )
endif( THE_C )

# name the project at this moment and this will set the compiler id:

# The version is used by the installed CMake package (hazmathConfigVersion.cmake)
project(HAZMATH VERSION 0.1.0 LANGUAGES C )

#message(STATUS "  ${HAZMATH_SOURCE_DIR}")



#message(STATUS "C compiler ID:${CMAKE_C_COMPILER_ID} Version:${CMAKE_C_COMPILER_VERSION}" )
message(STATUS "The system is ${CMAKE_HOST_SYSTEM}; apple=${APPLE}; unix=${UNIX}")
set(REAL_C "${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}")
############################
if(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
  string(REPLACE "gcc" "g++" C_XX ${THE_C} )
  string(REPLACE "gcc" "gfortran" F_C ${THE_C} )
  find_program(THE_CXX NAMES ${C_XX})
  find_program(THE_FC NAMES ${F_C})
elseif(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "mpicc.*") 
  string(REPLACE "mpicc" "mpicxx" C_XX ${THE_C} )
  string(REPLACE "mpicc" "mpif77" F_C ${THE_C} )
  find_program(THE_CXX NAMES ${C_XX})
  find_program(THE_FC NAMES ${F_C})
elseif( ${REAL_C} MATCHES "Intel.*" AND ${THE_C} MATCHES "icc.*" )
  find_program(THE_CXX NAMES icpc)
  find_program(THE_FC NAMES ifort)
elseif( ${REAL_C} MATCHES "Clang.*")
  find_program(THE_C NAMES clang)
  find_program(THE_CXX NAMES clang++ clang)
else()       
  message(WARNING    "
	   	     **  ${THE_C} did not match any of the preset C compilers. 
	   	     ** Continuing with the default compiler: ${CMAKE_C_COMPILER}
		     " )
  set(THE_C "0")
  set(THE_CXX "0")
  set(THE_FC "0")	      
endif(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
#
if( THE_C AND THE_CXX AND THE_FC )
  set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C   compiler" FORCE) 
  set(CMAKE_CXX_COMPILER ${THE_CXX} CACHE INTERNAL   "the C++ compiler" FORCE)
  set(CMAKE_Fortran_COMPILER ${THE_FC} CACHE INTERNAL    "the F compiler" FORCE)
endif( THE_C AND THE_CXX AND THE_FC )
# END COMPILERS SET UP................ 
##################################
enable_language(CXX) 
enable_language(Fortran) 
#
# OpenMP : defined on command line in the top Makefile
#
if(USE_OPENMP)
  find_package(OpenMP)

  if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
      "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${OpenMP_C_FLAGS}")	 
    set (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
      "${CMAKE_CXX_FLAGS_${CMAKE_CBUILD_TYPE}} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
      "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${OpenMP_C_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else(OPENMP_FOUND)
    message(WARNING  "  OpenMP was requested but not supported!")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

########## Additional compiler flags (not defined by the build

if(ADD_CFLAGS)
  set (CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
    "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CFLAGS}")	 
endif(ADD_CFLAGS)
if(ADD_CXXFLAGS)
  set (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
    "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CXXFLAGS}")	 
endif(ADD_CXXFLAGS)
if(ADD_FFLAGS)
  set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
    "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_FFLAGS}")	 
endif(ADD_FFLAGS)
#
# By default install in the source tree (lib/), as "make install" expects;
# CMAKE_INSTALL_PREFIX is ignored, as it always was. With
# -DHAZMATH_RELOCATABLE_INSTALL=ON the install goes to CMAKE_INSTALL_PREFIX
# with headers and a CMake package config (find_package(hazmath)).
set(HAZMATH_RELOCATABLE_INSTALL OFF CACHE BOOL
  "install to CMAKE_INSTALL_PREFIX with headers and a CMake package config")
if(NOT HAZMATH_RELOCATABLE_INSTALL)
  set(CMAKE_INSTALL_PREFIX "${HAZMATH_SOURCE_DIR}" )
endif()
include(GNUInstallDirs)
set(HAZMATHLIB_BASE_PATH "${HAZMATH_SOURCE_DIR}" CACHE PATH "path to HAZMATHLIB")
set(HAZMATH_INSTALL TRUE)
#
if(SHARED)
  set(HAZMATH_LIBRARY_TYPE SHARED)
else(SHARED)
  set(HAZMATH_LIBRARY_TYPE STATIC)
endif(SHARED)
if(NOT HAZMATH_RELOCATABLE_INSTALL)
  # FOR MAC OS X to find shared libs in install location
  set(CMAKE_INSTALL_NAME_DIR
    ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "path for apple")
  # FOR LINUX to find shared libs in install location
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_NAME_DIR} CACHE PATH "path for Linux")
endif()
file(GLOB HAZMATH_C_SOURCES RELATIVE ${HAZMATH_SOURCE_DIR} CONFIGURE_DEPENDS
  ${HAZMATH_SOURCE_DIR}/src/assemble/*.c 
  ${HAZMATH_SOURCE_DIR}/src/fem/*.c 
  ${HAZMATH_SOURCE_DIR}/src/mesh/*.c 
  ${HAZMATH_SOURCE_DIR}/src/amr/*.c 
  ${HAZMATH_SOURCE_DIR}/src/approximation/*.c 
  ${HAZMATH_SOURCE_DIR}/src/solver/*.c
  ${HAZMATH_SOURCE_DIR}/src/timestepping/*.c 
  ${HAZMATH_SOURCE_DIR}/src/nonlinear/*.c 
  ${HAZMATH_SOURCE_DIR}/src/utilities/*.c 
  ${HAZMATH_SOURCE_DIR}/src/eigen/*.c 
  ${HAZMATH_SOURCE_DIR}/src/interfaces/*.c
  ${HAZMATH_SOURCE_DIR}/src/amg_classical/*.c
  ${HAZMATH_SOURCE_DIR}/src/graphs/*.c
  ${HAZMATH_SOURCE_DIR}/src/assemble/*.inl 
  ${HAZMATH_SOURCE_DIR}/src/fem/*.inl
  ${HAZMATH_SOURCE_DIR}/src/mesh/*.inl
  ${HAZMATH_SOURCE_DIR}/src/solver/*.inl 
  ${HAZMATH_SOURCE_DIR}/src/timestepping/*.inl 
  ${HAZMATH_SOURCE_DIR}/src/nonlinear/*.inl 
  ${HAZMATH_SOURCE_DIR}/src/utilities/*.inl)
file(GLOB HAZMATH_CXX_SOURCES RELATIVE ${HAZMATH_SOURCE_DIR} CONFIGURE_DEPENDS
  ${HAZMATH_SOURCE_DIR}/src/amr/*.cpp)
file(GLOB HAZMATH_HAZNICS_SOURCES RELATIVE ${HAZMATH_SOURCE_DIR} 
  ${HAZMATH_SOURCE_DIR}/src/haznics/*.c)
file(GLOB HAZMATH_SWIG_SOURCES RELATIVE ${HAZMATH_SOURCE_DIR} 
  ${HAZMATH_SOURCE_DIR}/src/haznics/*.i)

file(GLOB HAZMATH_F_SOURCES RELATIVE ${HAZMATH_SOURCE_DIR} 
  ${HAZMATH_SOURCE_DIR}/src/amr/*.f) 
# auto generating headers if needed: hazmath.h is generated in the build
# directory when a source file or the generator changes, and copied to
# include/hazmath.h only when its content differs, so that editing a
# function body does not recompile the whole library. The copy runs on every
# build, which also restores a deleted include/hazmath.h; "make clean" never
# deletes it from the source tree.
add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/hazmath.h"
  COMMAND ./headmk.sh "${HAZMATH_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/hazmath.h"
  DEPENDS ${HAZMATH_C_SOURCES} ${HAZMATH_CXX_SOURCES}
    "${HAZMATH_SOURCE_DIR}/haz_shutils/headmk.sh"
    "${HAZMATH_SOURCE_DIR}/haz_shutils/mkheaders.awk"
  WORKING_DIRECTORY "${HAZMATH_SOURCE_DIR}/haz_shutils"
  COMMENT
  "HAZMATH: Autogenerating header file with C functions..."
  VERBATIM)
# Ninja must know that the copy can change include/hazmath.h, or it skips
# recompiling the objects that include it. Ninja's clean then deletes the
# header too; the next build restores it.
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(HAZMATH_HEADER_BYPRODUCTS BYPRODUCTS "${HAZMATH_SOURCE_DIR}/include/hazmath.h")
endif()
add_custom_target(headers
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_BINARY_DIR}/hazmath.h" "${HAZMATH_SOURCE_DIR}/include/hazmath.h"
  ${HAZMATH_HEADER_BYPRODUCTS}
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/hazmath.h"
  VERBATIM)

##################################################################
# For BLAS: 
##################################################################
if (USE_BLAS)
  # set the path to find specific modules
  set(CMAKE_MODULE_PATH "${HAZMATH_SOURCE_DIR}/cmake.modules")
  find_package(BLAS)
  if (BLAS_FOUND)
    add_definitions("-DWITH_BLAS=1")
    if (SHARED)
      set(CMAKE_SHARED_LINKER_FLAGS ${BLAS_LIBRARIES})
    endif(SHARED)
##    message("   BLAS was found!")
  else(BLAS_FOUND)
    message(WARNING  "    BLAS was requested but not supported!")
  endif(BLAS_FOUND)
endif(USE_BLAS)

##################################################################
# For LAPACK: 
##################################################################
if (USE_LAPACK)
  # set the path to find specific modules
  set(CMAKE_MODULE_PATH "${HAZMATH_SOURCE_DIR}/cmake.modules")
  find_package(LAPACK)
  if (LAPACK_FOUND)
    add_definitions("-DWITH_LAPACK=1")
    if (SHARED)
      set(CMAKE_SHARED_LINKER_FLAGS ${LAPACK_LIBRARIES})
    endif(SHARED)
    #message("   LAPACK was found!" ${LAPACK_LIBRARIES})
  else(LAPACK_FOUND)
    message(WARNING  " LAPACK was requested but not supported!")
  endif(LAPACK_FOUND)
endif(USE_LAPACK)

##################################################################
# For SUITESPARSE: we here require all suite sparse packages
##################################################################
if (USE_SUITESPARSE)
  # set the path to find specific modules	  
  #   message(STATUS "  ${HAZMATH_SOURCE_DIR}")
  set(CMAKE_MODULE_PATH "${HAZMATH_SOURCE_DIR}/cmake.modules")

  find_package(SUITESPARSE)

  # metis is not part of SuiteSparse, so there may be also some other metis dir. 
  ##   set(METIS_DIR "${SUITESPARSE_DIR}")

  if (SUITESPARSE_FOUND)
    add_definitions("-DWITH_SUITESPARSE=1")
    include_directories(${SUITESPARSE_INCLUDE_DIRS})
    if (SHARED)
      #      	 set(CMAKE_SHARED_LINKER_FLAGS ${SUITESPARSE_DIR}/libumfpack.dylib)
      set(CMAKE_SHARED_LINKER_FLAGS ${UMFPACK_LIBRARY})
    endif(SHARED)
  else(SUITESPARSE_FOUND)
    message(WARNING  " SUITESPARSE was requested but not supported!")
  endif(SUITESPARSE_FOUND)
endif(USE_SUITESPARSE)

##################################################################
# For HDF5: 
##################################################################
if (USE_HDF5)
  # set the path to find specific modules
  set(CMAKE_MODULE_PATH "${HAZMATH_SOURCE_DIR}/cmake.modules")
  set(HDF5_DIR /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/hdf5)
  set (HDF5_LIBRARY_NAME "hdf5_serial")
  # could be: set (HDF5_LIBRARY_NAME "hdf5_openmpi") 
  find_package(HDF5 REQUIRED COMPONENTS C CXX Fortran HL) ## COMPONENTS C CXX Fortran HL)
  if (HDF5_FOUND)
    message(STATUS "Found HDF5 at ${HDF5_DIR}")
    add_definitions("-DWITH_HDF5=1")
    if (SHARED)
      set(CMAKE_SHARED_LINKER_FLAGS ${HDF5_LIBRARIES})
    endif(SHARED)
##    message("   HDF5 was found!")
  else(HDF5_FOUND)
    message(WARNING  " HDF5 was requested but not supported!")
  endif(HDF5_FOUND)
endif(USE_HDF5)

##################################################################
# For FFTW3: 
##################################################################
if (USE_FFTW3)
  # set the path to find specific modules
  set(CMAKE_MODULE_PATH "${HAZMATH_SOURCE_DIR}/cmake.modules")
  set(FFTW3_DIR /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE})
  set (FFTW3_LIBRARY_NAME "fftw3")
#  Could be: set (FFTW3_LIBRARY_NAME "fftw3_mpi")
#  Could be: set (FFTW3_LIBRARY_NAME "fftw3_omp")
#  Could be: set (FFTW3_LIBRARY_NAME "fftw3_OTHER_THINGS")
  find_package(FFTW3 REQUIRED COMPONENTS C CXX Fortran ) ## COMPONENTS C CXX Fortran HL)
  if (FFTW3_FOUND)
    message(STATUS "Found FFTW3 at ${FFTW3_DIR}")
    add_definitions("-DWITH_FFTW3=1")
    if (SHARED)
      set(CMAKE_SHARED_LINKER_FLAGS ${FFTW3_LIBRARIES})
    endif(SHARED)
    message("   FFTW3 was found!")
  else(FFTW3_FOUND)
    message(WARNING  " FFTW3 was requested but not supported!")
  endif(FFTW3_FOUND)
endif(USE_FFTW3)

##################################################################
# For HAZNICS: we here require Python and SWIG
##################################################################
if (USE_HAZNICS)
  set(CMAKE_MODULE_PATH "${HAZMATH_SOURCE_DIR}/cmake.modules")
  
  message(STATUS "Checking for package dependencies for 'HAZNICS'")
  find_package(Python 3.5 REQUIRED COMPONENTS Interpreter Development NumPy)
  ## if you are The Special uncomment the line below to have python 3.5 below:
  ## find_package(Python 3.5 EXACT REQUIRED COMPONENTS Interpreter Development NumPy)
  
  if(Python_Development_FOUND)
    message(STATUS "Python Libs found: ${Python_LIBRARIES} version=${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}")
    include_directories(${Python_INCLUDE_DIRS})
  endif(Python_Development_FOUND)
  
  if(Python_NumPy_FOUND)
    message(STATUS "Python NumPy found: ${Python_NumPy_VERSION}")
    include_directories(${Python_NumPy_INCLUDE_DIRS})
  endif(Python_NumPy_FOUND)
  
  message(STATUS "Checking for  'SWIG'")  
  find_package(SWIG 4.0 COMPONENTS python)
  if(SWIG_FOUND)
    message(STATUS "SWIG found: ${SWIG_EXECUTABLE}")
  endif()  
  # Set initial values
  set(HAZNICS_INCLUDE_DIRS ${HAZMATH_SOURCE_DIR}/include)
  set(HAZNICS_LIBRARIES "")
  set (HAZNICS_FOUND FALSE)
  if(SWIG_FOUND AND Python_FOUND AND Python_NumPy_FOUND)
    include (UseSWIG)
    set (HAZNICS_FOUND TRUE)
    set(HAZNICS_LIBRARIES ${HAZNICS_INCLUDE_DIRS} ${Python_INCLUDE_DIRS})
    set(HAZNICS_LIBRARIES ${HAZNICS_LIBRARIES} ${Python_LIBRARIES})
    set(HAZNICS_INCLUDE_DIRS ${HAZNICS_INCLUDE_DIRS} ${SWIG_INCLUDE_DIRS})
    set(HAZNICS_LIBRARIES ${HAZNICS_LIBRARY} ${SWIG_LIBRARIES})
    add_definitions("-DWITH_HAZNICS=1")
  else()
    message(WARNING " HAZNICS was requested but not supported!")
  endif()
endif(USE_HAZNICS)
##################################################################
# For Doxygen
##################################################################
# target generating the HAZMATH documentation with doxygen
if (USE_DOXYGEN)

  find_package(Doxygen)

  if(DOXYGEN_FOUND)
    if(EXISTS ${HAZMATH_SOURCE_DIR}/haz_docs/hazmath.Doxygen.cnf.in)
      configure_file(
	${HAZMATH_SOURCE_DIR}/haz_docs/hazmath.Doxygen.cnf.in
	${CMAKE_CURRENT_BINARY_DIR}/hazmath.Doxygen.cnf @ONLY
	)
      set(DOXY_EXEC "${DOXYGEN_EXECUTABLE}")

      #if(DOXYWIZARD)
      #  find_program(WIZARD doxywizard)
      #  if(APPLE AND (NOT WIZARD) )
      #    find_program(WIZARD 
      #      /Applications/Doxygen.app/Contents/MacOS/Doxywizard)
      #   endif()
      #  if(WIZARD)
      #      set(DOXY_EXEC "${WIZARD}")
      #   endif()
      #endif(DOXYWIZARD)

      add_custom_target(docs ${DOXY_EXEC}
 	${CMAKE_CURRENT_BINARY_DIR}/hazmath.Doxygen.cnf
	WORKING_DIRECTORY 
	"${CMAKE_CURRENT_BINARY_DIR}"
	COMMENT 
	"Generating HAZMATH documentation by Doxygen" 
	VERBATIM)
    else(EXISTS ${HAZMATH_SOURCE_DIR}/haz_docs/hazmath.Doxygen.cnf.in)
      message(WARNING " Doxygen configuration file cannot be found!")
    endif(EXISTS ${HAZMATH_SOURCE_DIR}/haz_docs/hazmath.Doxygen.cnf.in)
  endif(DOXYGEN_FOUND)

endif(USE_DOXYGEN)
###############end documentation with doxygen

#################################################################
# For MATLAB: if want to use Matlab interface
#################################################################
if (USE_MATLAB)

  find_package(MATLAB)

  if (MATLAB_FOUND)
    add_definitions("-DWITH_MATLAB=1")
    include_directories(${Matlab_INCLUDE_DIRS})
  else(MATLAB_FOUND)
    message(WARNING  "MATLAB was requested but not supported!")	
  endif(MATLAB_FOUND)

endif(USE_MATLAB)

# Add include directories.
include_directories(${HAZMATHLIB_BASE_PATH}/include)
# 
# Recursively look for CMakeLists.txt in subdirs.
#    add_subdirectory("examples")

###MISMATCH CHECK:
if(NOT (CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))
  message("
**  CC and CXX mismatch. C compiler=\"${CMAKE_C_COMPILER_ID}\" ; C++ compiler=\"${CMAKE_CXX_COMPILER_ID}\"
"
    ) 
endif(NOT (CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))

##
# Build libhazmath
# depends also on the source, header and inl files 
add_library(hazmath ${HAZMATH_LIBRARY_TYPE} ${HAZMATH_C_SOURCES} ${HAZMATH_INL_SOURCES} ${HAZMATH_F_SOURCES}  ${HAZMATH_HEADERS})
add_dependencies(hazmath headers)
# PIC so that a static libhazmath can be linked into Python extensions
set_target_properties(hazmath PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(hazmath PUBLIC
  $<BUILD_INTERFACE:${HAZMATH_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hazmath>)
# the C sources use libm, so consumers of a static libhazmath need it too
if(NOT MSVC)
  target_link_libraries(hazmath PUBLIC m)
endif()

# Add links to libraries for external modules (needed for shared libraries on MAC)
# XLIBS_TO_EXPORT holds the same libraries, in the same order, as the paths
# that were found; a static libhazmath exports them to its consumers.
set(XLIBS_TO_LINK "")
set(XLIBS_TO_EXPORT "")
# blas 
if (BLAS_FOUND)
  LIST(APPEND XLIBS_TO_LINK "blas")
  LIST(APPEND XLIBS_TO_EXPORT ${BLAS_LIBRARY})
endif(BLAS_FOUND)
# lapack
if (LAPACK_FOUND)
  LIST(APPEND XLIBS_TO_LINK "lapack")
  LIST(APPEND XLIBS_TO_EXPORT ${LAPACK_LIBRARY})
endif(LAPACK_FOUND)
# suitesparse
if (SUITESPARSE_FOUND)
  LIST(APPEND XLIBS_TO_LINK ${SUITESPARSE_LIBRARIES})
  LIST(APPEND XLIBS_TO_EXPORT ${SUITESPARSE_LIBRARIES})
endif(SUITESPARSE_FOUND)
# hdf5
if (HDF5_FOUND)
  LIST(APPEND XLIBS_TO_LINK ${HDF5_LIBRARIES})
  LIST(APPEND XLIBS_TO_EXPORT ${HDF5_LIBRARIES})
endif(HDF5_FOUND)
# fftw3
if (FFTW3_FOUND)
  LIST(APPEND XLIBS_TO_LINK ${FFTW3_LIBRARIES})
  LIST(APPEND XLIBS_TO_EXPORT ${FFTW3_LIBRARIES})
endif(FFTW3_FOUND)
# matlab (no MATLAB library is called "matlab", so it is not exported)
if (MATLAB_FOUND)
  LIST(APPEND XLIBS_TO_LINK "matlab")
endif(MATLAB_FOUND)
message(STATUS " ADDITIONAL_LIBS_TO_LINK=${XLIBS_TO_LINK}" )
# transitive within this build, as before
target_link_libraries(hazmath PUBLIC "$<BUILD_INTERFACE:${XLIBS_TO_LINK}>")
# An installed shared libhazmath records its libraries itself. A static one
# exports them with the paths found here; hazmathConfig.cmake checks that
# they still exist.
if(NOT SHARED)
  target_link_libraries(hazmath PUBLIC "$<INSTALL_INTERFACE:${XLIBS_TO_EXPORT}>")
endif()

# install libhazmath
# The in-tree library and the haznics module have the configured install
# directory compiled in (rpath), so "cmake --install --prefix" cannot move them.
if(NOT HAZMATH_RELOCATABLE_INSTALL OR USE_HAZNICS)
  install(CODE "
    if(NOT CMAKE_INSTALL_PREFIX STREQUAL \"${CMAKE_INSTALL_PREFIX}\")
      message(FATAL_ERROR \"HAZmath was configured to install into ${CMAKE_INSTALL_PREFIX}; \"
        \"reconfigure with -DHAZMATH_RELOCATABLE_INSTALL=ON \"
        \"-DCMAKE_INSTALL_PREFIX=\${CMAKE_INSTALL_PREFIX} instead of passing --prefix.\")
    endif()")
endif()
if(NOT HAZMATH_RELOCATABLE_INSTALL)
  install(TARGETS hazmath
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)
else()
  # relocatable install: library, headers and CMake package config
  include(CMakePackageConfigHelpers)
  # no absolute paths baked into the installed library
  # (external libraries outside the system directories are found through
  # INSTALL_RPATH_USE_LINK_PATH)
  set_target_properties(hazmath PROPERTIES
    INSTALL_NAME_DIR "@rpath"
    INSTALL_RPATH ""
    INSTALL_RPATH_USE_LINK_PATH TRUE)
  install(TARGETS hazmath EXPORT hazmathTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  # hazmath.h and the headers it includes (haznics_add.h and mex.h are
  # only included for the SWIG and MATLAB builds)
  set(HAZMATH_PUBLIC_HEADERS hazmath.h macro.h mesh.h amr.h sparse.h dense.h
    vec.h fem.h solver.h nonlinear.h timestep.h param.h eigen.h graphs.h)
  list(TRANSFORM HAZMATH_PUBLIC_HEADERS PREPEND "${HAZMATH_SOURCE_DIR}/include/")
  install(FILES ${HAZMATH_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hazmath)
  set(HAZMATH_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/hazmath)
  install(EXPORT hazmathTargets
    NAMESPACE hazmath::
    DESTINATION ${HAZMATH_CMAKE_DIR})
  configure_package_config_file(
    ${HAZMATH_SOURCE_DIR}/cmake/hazmathConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/hazmathConfig.cmake
    INSTALL_DESTINATION ${HAZMATH_CMAKE_DIR})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/hazmathConfigVersion.cmake
    VERSION ${HAZMATH_VERSION}
    COMPATIBILITY SameMinorVersion)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hazmathConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hazmathConfigVersion.cmake
    DESTINATION ${HAZMATH_CMAKE_DIR})

endif()

if (USE_HAZNICS)
  if(HAZNICS_FOUND)
    #    add_library(haznics ${HAZMATH_LIBRARY_TYPE} ${HAZMATH_HAZNICS_SOURCES} ${HAZMATH_HAZNICS_HEADERS}) 
    set (CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
      "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} -fPIC")	 
    set (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
      "${CMAKE_CXX_FLAGS_${CMAKE_CBUILD_TYPE}} -fPIC")
    set (UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
    set (UseSWIG_MODULE_VERSION "2")
    set(CMAKE_SWIG_OUTDIR "${HAZMATH_SOURCE_DIR}/swig_files") # (swig -outdir option).
    message(STATUS
      " INFO: SWIG Lib Dir=${CMAKE_SWIG_OUTDIR}")
    set(SWIG_OUTFILE_DIR "${CMAKE_SWIG_OUTDIR}")
    swig_add_library(haznics
                 LANGUAGE python
                 TYPE SHARED
		 OUTPUT_DIR  "${CMAKE_SWIG_OUTDIR}"
		 OUTFILE_DIR "${CMAKE_SWIG_OUTDIR}"
                 SOURCES "${HAZMATH_SWIG_SOURCES}")
    add_custom_target(haznics_headers ALL ./headmk_haznics.sh "${HAZMATH_SOURCE_DIR}"
      DEPENDS ${HAZMATH_HAZNICS_SOURCES}
      WORKING_DIRECTORY "${HAZMATH_SOURCE_DIR}/haz_shutils"
      COMMENT 
      "HAZMATH: Autogenerating header file C functions for haznics/swig..."
      VERBATIM)
    add_library(haznics_helper  "${HAZMATH_HAZNICS_SOURCES}")
    add_dependencies(headers haznics_headers)
    # both compile code that includes the generated hazmath.h; with the
    # Makefile and VS generators SWIG itself runs in haznics_swig_compilation
    add_dependencies(haznics_helper headers)
    add_dependencies(haznics headers)
    if(TARGET haznics_swig_compilation)
      add_dependencies(haznics_swig_compilation headers)
    endif()
#    add_dependencies(haznics_helper haznics_headers)
    set_property(TARGET haznics PROPERTY SWIG_INCLUDE_DIRECTORIES ${HAZMATH_SOURCE_DIR}/include)
    # UseSWIG already passes -python; -py3 was removed in SWIG 4.1
    if(SWIG_VERSION VERSION_LESS "4.1")
      set_property(TARGET haznics PROPERTY SWIG_COMPILE_OPTIONS -py3 )
    endif()
    set(LIBS_TO_LINK "haznics_helper")
    LIST(APPEND LIBS_TO_LINK "hazmath")
     if (LAPACK_FOUND)
      LIST(APPEND LIBS_TO_LINK "lapack")
    endif(LAPACK_FOUND)
    if (SUITESPARSE_FOUND)
      LIST(APPEND LIBS_TO_LINK ${SUITESPARSE_LIBRARIES})
    endif(SUITESPARSE_FOUND)
     if (HDF5_FOUND)
      LIST(APPEND LIBS_TO_LINK ${HDF5_LIBRARY_NAME})
    endif(HDF5_FOUND)
     if (FFTW3_FOUND)
      LIST(APPEND LIBS_TO_LINK ${FFTW3_LIBRARY_NAME})
    endif(FFTW3_FOUND)
    LIST(APPEND LIBS_TO_LINK ${Python_LIBRARIES})
    message(STATUS " ADDITIONAL_LIBS_TO_LINK=${LIBS_TO_LINK}" )
    target_link_libraries(haznics "${LIBS_TO_LINK}")
#MAKE_SWIG_FLAGS can be used to add special flags to all swig calls.
#SWIG_MODULE_haznics_EXTRA_DEPS property

    if(HAZMATH_RELOCATABLE_INSTALL)
      # the module stays in the source tree; point it at the installed libhazmath
      set_target_properties(haznics PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
    endif()
    install(TARGETS haznics 
      LIBRARY DESTINATION "${CMAKE_SWIG_OUTDIR}"
      RUNTIME DESTINATION "${CMAKE_SWIG_OUTDIR}"
      ARCHIVE DESTINATION "${CMAKE_SWIG_OUTDIR}")
  else(HAZNICS_FOUND)
    message(WARNING  "HAZNICS not found or not requested!")
  endif(HAZNICS_FOUND)
endif(USE_HAZNICS)

