# -*- mode: cmake -*-
# MADS
#
# USAGE: cmake CMakeLists.txt -DCMAKE_CXX_COMPILER=/opt/local/bin/c++ -DCMAKE_C_COMPILER=/opt/local/bin/gcc -DCMAKE_LINKER=/opt/local/bin/gcc
#
cmake_minimum_required (VERSION 2.8)

# enable_testing()

project (mads)

set(MADS_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/util")
set(CMAKE_MODULE_PATH  ${MADS_MODULE_PATH})

execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE VER_N)
string(REPLACE "\n" "" VER ${VER_N})
execute_process(COMMAND scripts/check_git_status OUTPUT_VARIABLE GIT_STATUS_N)
string(REPLACE "\n" "" GIT_STATUS ${GIT_STATUS_N})
message( STATUS "Git version ${VER}${GIT_STATUS}" )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "debug"
	AND NOT cmake_build_type_tolower STREQUAL "release"
	AND NOT cmake_build_type_tolower STREQUAL "lib")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, Lib (case-insensitive).")
endif()

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(PkgConfig REQUIRED)
FIND_PACKAGE(Git REQUIRED)
find_package(GLIB2 REQUIRED)
find_package(GSL REQUIRED)
find_package(LibYAML REQUIRED)
find_package(LibXML2 REQUIRED)
find_package(LibMathEval REQUIRED)

if (GSL_FOUND)
    add_definitions(${GSL_CXX_FLAGS})
    set(MADS_LIBS ${MADS_LIBS} gsl gslcblas lapack blas stdc++ "-L/opt/local/lib")
endif (GSL_FOUND)

set(CMAKE_CXX_FLAGS "-Wall -O1 -Winit-self" )
set(CMAKE_C_FLAGS "-Wall -O1 -Winit-self" )

include(CheckLibraryExists)
check_library_exists(dl dlopen "" LIBDL_FOUND)

# include(CheckTagLibFileName)
add_subdirectory(src)
add_subdirectory(wells)

if(APPLE)
    set(MADS_OS "MacOSX")
elseif(CYGWIN)
    set(MADS_OS "Cygwin")
else(APPLE)
    set(MADS_OS "${CMAKE_SYSTEM_NAME}")
endif(APPLE)

if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "-Wall -Winit-self" CACHE STRING "Default flags" FORCE)
    SET(CMAKE_CXX_FLAGS_LIB "-O1 -fPIC" CACHE STRING "Lib flags" FORCE)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Debug flags" FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE "-O1" CACHE STRING "Release flags" FORCE)
endif(CMAKE_COMPILER_IS_GNUCXX)

# Enable ExternalProject CMake module
include(ExternalProject)
 
# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/tpls)

# --- clean
add_custom_target(clean-cmake-files
   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
)

set(cmake_generated ${CMAKE_BINARY_DIR}/CMakeCache.txt
                    ${CMAKE_BINARY_DIR}/cmake_install.cmake  
                    ${CMAKE_BINARY_DIR}/Makefile
                    ${CMAKE_BINARY_DIR}/CMakeFiles
)

foreach(file ${cmake_generated})
  if (EXISTS ${file})
     file(REMOVE_RECURSE ${file})
  endif()
endforeach(file)

add_custom_target(clean-all
   COMMAND ${CMAKE_BUILD_TOOL} clean
   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
)
