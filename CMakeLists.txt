# -*- mode: cmake -*-
# MADS
cmake_minimum_required (VERSION 2.8)
enable_testing()
project (mads)
set(MADS_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH  ${MADS_MODULE_PATH})

find_package(LAPACK)
find_package(BLAS)
find_package(PkgConfig)

# check to use if using GSL
if (USE_GSL)
    add_definitions(-DUSE_GSL)
    set(TPL_LIBS ${TPL_LIBS}
	gsl
	gslcblas)
endif (USE_GSL)

find_package(LAPACK)
find_package(BLAS)
find_package(PkgConfig)

set(CMAKE_CXX_FLAGS "-g -Wall")

include(CheckLibraryExists)
check_library_exists(dl dlopen "" LIBDL_FOUND)


include(CheckTagLibFileName)
add_subdirectory (lm)
add_subdirectory (lhs)
add_subdirectory (misc)
add_subdirectory (mprun)
add_subdirectory (ds)
add_subdirectory (sa)
add_subdirectory (bayes)

set( objsmads ./mads.c ./mads_io.c ./mads_io_external.c ./mads_func.c ./mads_mem.c ./mads_info.c lm/opt_lm_mon.c lm/opt_lm_gsl.c lm/lu.c lm/opt_lm_ch.c misc/test_problems.c misc/anasol_contamination.c misc/io.c lhs/lhs.c )
set( objspso pso/pso-tribes-lm.c pso/Standard_PSO_2006.c pso/mopso.c )
set( objsa sa/abagus.c sa/postpua.c sa/global.c sa/do_miser.c )
set( objds ds/infogap.c ds/glue.c )
set( objsmpun mprun/mprun.c mprun/mprun_io.c )
set( objskdtree misc/kdtree-0.5.5/kdtree.c )
set( objslevmar misc/levmar-2.5/lm_m.c misc/levmar-2.5/Axb.c misc/levmar-2.5/misc.c misc/levmar-2.5/lmlec.c misc/levmar-2.5/lmbc.c misc/levmar-2.5/lmblec.c misc/levmar-2.5/lmbleic.c )
set( objsastable misc/astable/astable.c misc/astable/interpolation.c misc/astable/pqueue.c )
set( objsbayes bayes/dream.c )
set( objwells wells/wells.c )

add_executable (mads ${objsmads} ${objspso} ${objsa} ${objds} ${objsmpun} ${objskdtree} ${objslevmar} ${objsastable} ${objsbayes} ${objwells} )

# Enable ExternalProject CMake module
include(ExternalProject)
 
# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/ThirdParty)
 
ExternalProject_Add(
    libxml2
    DEPENDS zlib
    GIT_REPOSITORY git://git.gnome.org/libxml2
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ./autogen.sh --prefix=${CMAKE_BINARY_DIR} --without-python --disable-shared --without-zlib
    BUILD_IN_SOURCE 1
)

# Add gmock
ExternalProject_Add(
    googlemock
    SVN_REPOSITORY http://googlemock.googlecode.com/svn/trunk/
    TIMEOUT 30
    PATCH_COMMAND svn patch ${CMAKE_SOURCE_DIR}/gtest.patch ${CMAKE_BINARY_DIR}/ThirdParty/src/googlemock/gtest
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
    -Dgtest_force_shared_crt=ON
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)
 
# Specify include dir for googlemock and googletest
ExternalProject_Get_Property(googlemock source_dir)
include_directories(${source_dir}/include)
include_directories(${source_dir}/gtest/include)
 
if(MSVC_VERSION EQUAL 1700)
    add_definitions(-D_VARIADIC_MAX=10)
endif()
 
# Add test executable target
add_executable(MainTest ${PROJECT_SOURCE_DIR}/src/main.cpp)
 
# Create dependency of MainTest on googlemock
add_dependencies(MainTest googlemock)
 
# Specify MainTest's link libraries
ExternalProject_Get_Property(googlemock binary_dir)
target_link_libraries(MainTest
debug ${binary_dir}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${CMAKE_FIND_LIBRARY_SUFFIXES}
optimized ${binary_dir}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${CMAKE_FIND_LIBRARY_SUFFIXES})

add_revision(liblas
  GIT_REPOSITORY "https://github.com/libLAS/libLAS.git"
  GIT_TAG 1.7.0
)

add_revision(libxml2
  GIT_REPOSITORY "https://git.gnome.org/browse/libxml2"
  GIT_TAG v2.9.1
)

find_package(LibXml2 REQUIRED)
macro_log_feature( LIBXML2_FOUND "libxml2" "LibXML2 is an XML parser required by mp3tunes." "http://www.xmlsoft.org" FALSE "" "" )
find_package(GLIB2)
macro_log_feature( GLIB2_FOUND "glib2" "Required by libgpod and mp3tunes" "http://www.gtk.org" FALSE "2.x" "")

include(CheckLibraryExists)
check_library_exists(dl dlopen "" LIBDL_FOUND)

include(CheckTagLibFileName)

add_subdirectory (lm)
add_subdirectory (lhs)
add_subdirectory (misc)
add_subdirectory (mprun)
add_subdirectory (ds)
add_subdirectory (sa)
add_subdirectory (bayes)

set( objsmads ./mads.c ./mads_io.c ./mads_io_external.c ./mads_func.c ./mads_mem.c ./mads_info.c lm/opt_lm_mon.c lm/opt_lm_gsl.c lm/lu.c lm/opt_lm_ch.c misc/test_problems.c misc/anasol_contamination.c misc/io.c lhs/lhs.c )
set( objspso pso/pso-tribes-lm.c pso/Standard_PSO_2006.c pso/mopso.c )
set( objsa sa/abagus.c sa/postpua.c sa/global.c sa/do_miser.c )
set( objds ds/infogap.c ds/glue.c )
set( objsmpun mprun/mprun.c mprun/mprun_io.c )
set( objskdtree misc/kdtree-0.5.5/kdtree.c )
set( objslevmar misc/levmar-2.5/lm_m.c misc/levmar-2.5/Axb.c misc/levmar-2.5/misc.c misc/levmar-2.5/lmlec.c misc/levmar-2.5/lmbc.c misc/levmar-2.5/lmblec.c misc/levmar-2.5/lmbleic.c )
set( objsastable misc/astable/astable.c misc/astable/interpolation.c misc/astable/pqueue.c )
set( objsbayes bayes/dream.c )
set( objwells wells/wells.c )

add_executable (mads ${objsmads} ${objspso} ${objsa} ${objds} ${objsmpun} ${objskdtree} ${objslevmar} ${objsastable} ${objsbayes} ${objwells} )

# Enable ExternalProject CMake module
include(ExternalProject)
 
# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/ThirdParty)
 
# Add libXML2
ExternalProject_Add(
    libxml2
    DEPENDS zlib
    GIT_REPOSITORY git://git.gnome.org/libxml2
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ./autogen.sh --prefix=${CMAKE_BINARY_DIR} --without-python --disable-shared --without-zlib
    BUILD_IN_SOURCE 1
)

# Add gmock
ExternalProject_Add(
    googlemock
    SVN_REPOSITORY http://googlemock.googlecode.com/svn/trunk/
    TIMEOUT 30
    PATCH_COMMAND svn patch ${CMAKE_SOURCE_DIR}/gtest.patch ${CMAKE_BINARY_DIR}/ThirdParty/src/googlemock/gtest
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
    -Dgtest_force_shared_crt=ON
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)
 
# Specify include dir for googlemock and googletest
ExternalProject_Get_Property(googlemock source_dir)
include_directories(${source_dir}/include)
include_directories(${source_dir}/gtest/include)
 
if(MSVC_VERSION EQUAL 1700)
    add_definitions(-D_VARIADIC_MAX=10)
endif()
 
# Add test executable target
add_executable(MainTest ${PROJECT_SOURCE_DIR}/src/main.cpp)
 
# Create dependency of MainTest on googlemock
add_dependencies(MainTest googlemock)
 
# Specify MainTest's link libraries
ExternalProject_Get_Property(googlemock binary_dir)
target_link_libraries(MainTest
debug ${binary_dir}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${CMAKE_FIND_LIBRARY_SUFFIXES}
optimized ${binary_dir}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${CMAKE_FIND_LIBRARY_SUFFIXES})

add_revision(libxml2
  GIT_REPOSITORY "https://git.gnome.org/browse/libxml2"
  GIT_TAG v2.9.1
)

find_package(LibXml2 REQUIRED)
macro_log_feature( LIBXML2_FOUND "libxml2" "LibXML2 is an XML parser required by mp3tunes." "http://www.xmlsoft.org" FALSE "" "" )
find_package(GLIB2)
macro_log_feature( GLIB2_FOUND "glib2" "Required by libgpod and mp3tunes" "http://www.gtk.org" FALSE "2.x" "")
