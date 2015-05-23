# -*- mode: cmake -*-
# MADS
#
# USAGE:   cmake CMakeLists.txt -DCMAKE_CXX_COMPILER=<gnu_c++_compiler> -DCMAKE_C_COMPILER=<gnu_c_compiler>   -DCMAKE_LINKER=<gnu_c_compiler>
# EXAMPLE: cmake CMakeLists.txt -DCMAKE_CXX_COMPILER=/opt/local/bin/c++ -DCMAKE_C_COMPILER=/opt/local/bin/gcc -DCMAKE_LINKER=/opt/local/bin/gcc

cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

# enable_testing()

message(STATUS "-------------------------------------------------------" )
message(STATUS "MADS: Model Analyses & Decision Support (v.1.1.14) 2015")
message(STATUS "-------------------------------------------------------" )
project (mads)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
	message(FATAL_ERROR "GNU compiler version must be at least 4.8!")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message(FATAL_ERROR "You are using Clang compiler! GNU compiler is required!")
else()
    message(WARNING "You are using an unsupported compiler! GNU compiler is required!.")
endif()

set(MADS_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${MADS_MODULE_PATH})
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${PROJECT_SOURCE_DIR}/tpls/lib)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/tpls/include)
set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${PROJECT_SOURCE_DIR}/tpls/lib)
message(STATUS "Library path: ${LD_LIBRARY_PATH}" )

include(PrintVariable)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "debug"
	AND NOT cmake_build_type_tolower STREQUAL "release"
	AND NOT cmake_build_type_tolower STREQUAL "lib")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, Lib (case-insensitive).")
endif()

set(USE_GSL TRUE)
set(USE_LAPACK TRUE)
set(USE_YAML TRUE)
set(USE_XML TRUE)
set(USE_MATHEVAL TRUE)
option(USE_HTTPS "DOWNLOAD USING HTTPS." FALSE)
option(USE_GSL "USE_GSL." TRUE)
option(USE_LAPACK "USE_LAPACK." TRUE)
option(USE_YAML "USE_YAML." TRUE)
option(USE_XML "USE_XML." TRUE)
option(USE_MATHEVAL "USE_MATHEVAL." TRUE)
include(TPLVersions)
find_package(Git 1.8)
find_package(PkgConfig 0.23)
find_package(Atlas 3.10)
find_package(LAPACK)
find_package(BLAS)
find_package(gfortran)
find_package(GLIB2 ${GLIB2_VERSION})
find_package(GSL ${GSL_VERSION})
find_package(LibYAML ${LIBYAML_VERSION})
find_package(LibXML2 ${LIBXML2_VERSION})
find_package(LibMathEval ${LIBMATHEVAL_VERSION}) 
if(NOT (GIT_FOUND))
	message(WARNING "\nGIT is missing! Please install it!\nCheckout http://git-scm.com")
endif()
if(NOT (PKGCONFIG_FOUND))
	message(FATAL_ERROR "\nPkgConfig is missing! Please install it!
Checkout the MADS 'readme' file.
MAC OS X: sudo port install gsl glib2 libyaml libxml2 libmatheval atlas pkgconfig git tcsh gmake
RHEL/CentOS/Fedora: sudo yum install gsl glib2 libyaml libxml2 libmatheval atlas pkgconfig git tcsh gmake
Ubuntu/Debian: sudo apt-get install libgsl0-dev libglib2.0 libyaml-dev libxml2 libmatheval-dev liblapack-dev gfortran guile flex pkg-config git tcsh gmake")
endif()

add_definitions(${GSL_CXX_FLAGS})
# set(MADS_LIBS ${MADS_LIBS} lapack blas atlas cblas stdc++ gfortran refblas "-L/opt/local/lib" "-L/Users/monty/lib")
set(MADS_LIBS ${MADS_LIBS} stdc++ gsl gslcblas gfortran -fopenmp "-L${LD_LIBRARY_PATH}" )
# set(MADS_LIBS ${MADS_LIBS} stdc++ gfortran "-L/opt/local/lib")

set(CMAKE_CXX_FLAGS "-Wall -O3 -Winit-self -fopenmp" )
set(CMAKE_C_FLAGS "-Wall -O3 -Winit-self -fopenmp" )

include(CheckLibraryExists)
check_library_exists(dl dlopen "" LIBDL_FOUND)

# include(CheckTagLibFileName)
include(TPLtest)
set(TPL_MISSING TRUE)
if(TPL_MISSING)
	add_subdirectory(tpls)
	set(MADS_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
	include_directories("${PROJECT_SOURCE_DIR}/tpls/include")
endif()

message(STATUS "-------------------------------------------------------" )
message(STATUS "MADS: Model Analyses & Decision Support (v.1.1.14) 2015")
message(STATUS "-------------------------------------------------------" )
if(GIT_FOUND)
	execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE VER_N)
	string(REPLACE "\n" "" VER ${VER_N})
	execute_process(COMMAND scripts/check_git_status OUTPUT_VARIABLE GIT_STATUS_N)
	string(COMPARE EQUAL "${GIT_STATUS}" "" result)
	if( NOT result )
    	string(REPLACE "\n" "" GIT_STATUS ${GIT_STATUS_N})
	endif()
	message( STATUS "MADS git repo version ${VER}${GIT_STATUS}" )
endif()

add_subdirectory(src)
add_subdirectory(wells)

if(APPLE)
    set(MADS_OS "MacOSX")
elseif(CYGWIN)
    set(MADS_OS "Cygwin")
else(APPLE)
    set(MADS_OS "${CMAKE_SYSTEM_NAME}")
endif(APPLE)
message( STATUS "Operating system: ${MADS_OS}" )
message(STATUS "-------------------------------------------------------" )

if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "-Wall -Winit-self" CACHE STRING "Default flags" FORCE)
    SET(CMAKE_CXX_FLAGS_LIB "-O3 -fPIC" CACHE STRING "Lib flags" FORCE)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Debug flags" FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "Release flags" FORCE)
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
   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake)
