############################################################
# CMakeLists for HAZMAT (modified fom FASP)
# No Fortran compiler check. 
#
#
# Modified   2015-08-08   --ltz
############################################################
cmake_minimum_required (VERSION 2.8)
# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_VERBOSE_MAKEFILE 1) 
set(GDB 1 CACHE BOOL "debugging or not")
set(OPENMP 0 CACHE BOOL "Openmp use")
set(USE_MUMPS 0 CACHE BOOL "MUMPS use")

# Search for C compilers in the specified order. That will determine the rest.
if(DEFINED ENV{CC}) 
	find_program(THE_C NAMES $ENV{CC} gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
else(DEFINED ENV{CC}) 
	find_program(THE_C NAMES gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
endif(DEFINED ENV{CC}) 

if( THE_C )
    set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C compiler" FORCE)
else( THE_C )
    message("Warning: the requested CC is not found: ${THE_C}; continuing with the system's one" )
endif( THE_C )

# name the project at this moment and this will set the compiler id:

project(HAZMAT C )

message(STATUS " XXXXXXXXX ${HAZMAT_SOURCE_DIR}")


message(STATUS "C compiler ID:${CMAKE_C_COMPILER_ID} Version:${CMAKE_C_COMPILER_VERSION}" )
message(STATUS "The system is ${CMAKE_HOST_SYSTEM}; apple=${APPLE}; unix=${UNIX}")
set(REAL_C "${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}")
############################
if(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
     		string(REPLACE "gcc" "g++" C_XX ${THE_C} )
		find_program(THE_CXX NAMES ${C_XX})
       elseif( ${REAL_C} MATCHES "Intel.*" AND ${THE_C} MATCHES "icc.*" )
		find_program(THE_CXX NAMES icpc)
       elseif( ${REAL_C} MATCHES "Clang.*")
	 find_program(THE_C NAMES clang)
	find_program(THE_CXX NAMES clang++ clang)
	else()       
	   message("
	   	     ** WARNING: ${THE_C} did not match any of the preset C compilers. 
	   	     ** Continuing with the default compiler: ${CMAKE_C_COMPILER}
		     " )
   	   set(THE_C "0")
   	   set(THE_CXX "0")	      
    endif(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
#
	if( THE_C AND THE_CXX )
	     	 set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C   compiler" FORCE) 
	    	 set(CMAKE_CXX_COMPILER ${THE_CXX} CACHE INTERNAL   "the C++ compiler" FORCE)
	endif( THE_C AND THE_CXX )
# END COMPILERS SET UP................ 
##################################
	enable_language(CXX) 
#####	enable_language(Fortran) 

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
          message(" WARNING: OpenMP was requested but not supported!")
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
#
    set(CMAKE_INSTALL_PREFIX "${HAZMAT_SOURCE_DIR}" )
    set(HAZMATLIB_BASE_PATH "${HAZMAT_SOURCE_DIR}" CACHE PATH "path to HAZMATLIB")
    set(HAZMAT_INSTALL TRUE)
#
	if(SHARED) 
	   set(HAZMAT_LIBRARY_TYPE SHARED)
        else(SHARED)
	   set(HAZMAT_LIBRARY_TYPE STATIC)
	endif(SHARED)	   
# FOR MAC OS X to find shared libs in install location
        set(CMAKE_INSTALL_NAME_DIR 
		${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "path for apple")
# FOR LINUX to find shared libs in install location
        set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_NAME_DIR} CACHE PATH "path for Linux")
file(GLOB HAZMAT_C_SOURCES RELATIVE ${HAZMAT_SOURCE_DIR} 
	  		    ${HAZMAT_SOURCE_DIR}/src/assemble/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/fem/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/grid/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/solver/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/utilities/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/assemble/*.inl 
	  		    ${HAZMAT_SOURCE_DIR}/src/fem/*.inl
	  		    ${HAZMAT_SOURCE_DIR}/src/grid/*.inl
	  		    ${HAZMAT_SOURCE_DIR}/src/solver/*.inl 
	  		    ${HAZMAT_SOURCE_DIR}/src/utilities/*.inl)

# auto generating headers if needed
 add_custom_target(headers ./headmk.sh "${HAZMAT_SOURCE_DIR}"
		   DEPENDS ${HAZMAT_C_SOURCES}
                   WORKING_DIRECTORY "${HAZMAT_SOURCE_DIR}/haz_shutils"
                   COMMENT 
		   "HAZMAT: Autogenerating header file with C functions..."
                   VERBATIM)
##################################################################
# For SUITESPARSE: we here require all suite sparse packages
##################################################################
if (USE_SUITESPARSE)
    # set the path to find specific modules	  
#   message(STATUS "  ${HAZMAT_SOURCE_DIR}")
   set(CMAKE_MODULE_PATH "${HAZMAT_SOURCE_DIR}/cmake.modules")

   find_package(SUITESPARSE)

# metis is not part of SuiteSparse, so there may be also some other metis dir. 
##   set(METIS_DIR "${SUITESPARSE_DIR}")


   if (SUITESPARSE_FOUND)
      add_definitions("-DWITH_SUITESPARSE=1")
      include_directories(${SUITESPARSE_INCLUDE_DIRS})
   else(SUITESPARSE_FOUND)
      message("   WARNING: SUITESPARSE was requested but not supported!")
     endif(SUITESPARSE_FOUND)
endif(USE_SUITESPARSE)

# Add include directories.
    include_directories(${HAZMATLIB_BASE_PATH}/include)
# 
# Recursively look for CMakeLists.txt in subdirs.
    add_subdirectory("examples")

###MISMATCH CHECK:
if(NOT (CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))
  message("
** WARNING: CC and CXX mismatch. C compiler=\"${CMAKE_C_COMPILER_ID}\" ; C++ compiler=\"${CMAKE_CXX_COMPILER_ID}\"
"
) 
endif(NOT (CMAKE_CXX_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))

##
# Build libhazmat
# depends also on the source, header and inl files 
add_library(hazmat ${HAZMAT_LIBRARY_TYPE} ${HAZMAT_C_SOURCES} ${HAZMAT_INL_SOURCES} ${HAZMAT_HEADERS}) 

# install libhazmat
    install(TARGETS hazmat
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)
