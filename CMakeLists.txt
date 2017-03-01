############################################################
# CMakeLists for HAZMAT (modified fom FASP)
# No Fortran compiler check. 
#
#
# Modified   2015-08-08   --ltz
#TODO:should know to do "make headers" but it does not...
############################################################
cmake_minimum_required (VERSION 2.8)

# Stupidity Prevention
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
	message(FATAL_ERROR 
"------------------------------------------------
Source and build directories cannot be the same.  
Please REMOVE CMakeCache.txt and CMakeFiles.
TYPE 'make help' for more information. 
-------------------------------------------------"
)
endif()

# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_VERBOSE_MAKEFILE 1) 
set(GDB 1 CACHE BOOL "debugging or not")
set(OPENMP 0 CACHE BOOL "Openmp use")
set(USE_MUMPS 0 CACHE BOOL "MUMPS use")

# Search for C compilers in the specified order. That will determine the rest.
if(DEFINED ENV{CC}) 
	find_program(THE_C NAMES $ENV{CC} gcc-6 gcc-mp-5 gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
else(DEFINED ENV{CC}) 
	find_program(THE_C NAMES gcc-6 gcc-mp-5 gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
endif(DEFINED ENV{CC}) 

if( THE_C )
    set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C compiler" FORCE)
else( THE_C )
    message("Warning: the requested CC is not found: ${THE_C}; continuing with the system's one" )
endif( THE_C )

# name the project at this moment and this will set the compiler id:

project(HAZMAT C )

#message(STATUS " XXXXXXXXX ${HAZMAT_SOURCE_DIR}")


#message(STATUS "C compiler ID:${CMAKE_C_COMPILER_ID} Version:${CMAKE_C_COMPILER_VERSION}" )
message(STATUS "The system is ${CMAKE_HOST_SYSTEM}; apple=${APPLE}; unix=${UNIX}")
set(REAL_C "${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}")
############################
if(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
     		string(REPLACE "gcc" "g++" C_XX ${THE_C} )
		string(REPLACE "gcc" "gfortran" F_C ${THE_C} )
		find_program(THE_CXX NAMES ${C_XX})
		find_program(THE_FC NAMES ${F_C})
       elseif( ${REAL_C} MATCHES "Intel.*" AND ${THE_C} MATCHES "icc.*" )
		find_program(THE_CXX NAMES icpc)
		find_program(THE_FC NAMES ifort)
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
	if(ADD_FFLAGS)
          set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
	  	  "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_FFLAGS}")	 
	endif(ADD_FFLAGS)
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
			    ${HAZMAT_SOURCE_DIR}/src/timestepping/*.c 
			    ${HAZMAT_SOURCE_DIR}/src/nonlinear/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/utilities/*.c 
	  		    ${HAZMAT_SOURCE_DIR}/src/assemble/*.inl 
	  		    ${HAZMAT_SOURCE_DIR}/src/fem/*.inl
	  		    ${HAZMAT_SOURCE_DIR}/src/grid/*.inl
	  		    ${HAZMAT_SOURCE_DIR}/src/solver/*.inl 
			    ${HAZMAT_SOURCE_DIR}/src/timestepping/*.inl 
			    ${HAZMAT_SOURCE_DIR}/src/nonlinear/*.inl 
	  		    ${HAZMAT_SOURCE_DIR}/src/utilities/*.inl)

file(GLOB HAZMAT_F_SOURCES RELATIVE ${HAZMAT_SOURCE_DIR}
	  		    ${HAZMAT_SOURCE_DIR}/src/solver/*.f)

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

#################################################################
# For MATLAB: if want to use Matlab interface
#################################################################
if (USE_MATLAB)

   find_package(MATLAB)

   if (MATLAB_FOUND)
      add_definitions("-DWITH_MATLAB=1")
      include_directories(${Matlab_INCLUDE_DIRS})
   else(MATLAB_FOUND)
      message("    WARNING: MATLAB was requested but not supported!")	
   endif(MATLAB_FOUND)

endif(USE_MATLAB)

# Add include directories.
    include_directories(${HAZMATLIB_BASE_PATH}/include)
# 
# Recursively look for CMakeLists.txt in subdirs.
#    add_subdirectory("examples")

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
add_library(hazmat ${HAZMAT_LIBRARY_TYPE} ${HAZMAT_C_SOURCES} ${HAZMAT_INL_SOURCES} ${HAZMAT_F_SOURCES}  ${HAZMAT_HEADERS}) 

# install libhazmat
    install(TARGETS hazmat
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)
