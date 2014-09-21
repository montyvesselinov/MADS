# -*- mode: cmake -*-
# MADS TPLs
cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

project(TPL)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH 
  ${CMAKE_MODULE_PATH} 
  ${TPL_SOURCE_DIR}/cmake)

# set up directory structure
set(CMAKE_INSTALL_PREFIX ${TPL_INSTALL_PREFIX})
set(TPL_CMAKE_BUILD_DIR "${TPL_SOURCE_DIR}/bin")
set(TPL_CMAKE_TEMPLATE_DIR "${TPL_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH  ${TPL_CMAKE_TEMPLATE_DIR})

OPTION(USE_HTTPS "DOWNLOAD USING HTTPS." FALSE)

set(TPL_BUILD_TYPE_DFLT "Debug")
if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE ${TPL_BUILD_TYPE_DFLT})
endif()
message(STATUS "TPL build type: ${CMAKE_BUILD_TYPE}")

set(TPL_DOWNLOAD_DIR_DFLT ${TPL_BINARY_DIR}/Downloads)
if(NOT TPL_DOWNLOAD_DIR)
  set(TPL_DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR_DFLT})
endif()
message(STATUS "TPL Download directory: ${TPL_DOWNLOAD_DIR}") 

set(TPL_INSTALL_PREFIX_DFLT "${TPL_BINARY_DIR}")
if ( NOT TPL_INSTALL_PREFIX )
  set(TPL_INSTALL_PREFIX ${TPL_INSTALL_PREFIX_DFLT})
endif()
message(STATUS "TPL Install directory: ${TPL_INSTALL_PREFIX}")

include(ExternalProject)
include(DefineExternalProjectArgs)
include(TPLVersions)

if(USE_GSL AND NOT GSL_FOUND)
	include(BuildGSL)
endif()

if(USE_LAPACK AND ( NOT LAPACK_FOUND OR NOT BLAS_FOUND OR NOT ATLAS_FOUND ) )
	include(BuildLAPACK)
endif()

if(USE_LIBXML2 AND NOT LIBXML2_FOUND)
	include(BuildLibXML2)
endif()

if(USE_LIBMATHEVAL AND NOT LIBMATHEVAL_FOUND)
	include(BuildLIBMATHEVAL) 
endif()
