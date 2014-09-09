# GNU Scientific Library GSL
#
# http://www.gnu.org/software/gsl/  and
# http://gnuwin32.sourceforge.net/packages/gsl.htm
#
# It defines the following variables:
#  GSL_FOUND - system has GSL lib
#  GSL_INCLUDE_DIRS - where to find headers
#  GSL_LIBRARIES - full path to the libraries
#  GSL_LIBRARY_DIRS, the directory where the PLplot library is found.
 
#  CMAKE_GSL_CXX_FLAGS  = Unix compiler flags for GSL, essentially "`gsl-config --cxxflags`"
#  GSL_LINK_DIRECTORIES = link directories, useful for rpath on Unix
#  GSL_EXE_LINKER_FLAGS = rpath on Unix
 
set( GSL_FOUND OFF )
set( GSL_CBLAS_FOUND OFF )
set( GSL_FIND_QUIETLY OFF )
 
# Windows, but not for Cygwin and MSys where gsl-config is available
if( WIN32 AND NOT CYGWIN AND NOT MSYS )
  # look for headers
  find_path( GSL_INCLUDE_DIR
    NAMES gsl/gsl_cdf.h gsl/gsl_randist.h gsl/gsl
    )
  if( GSL_INCLUDE_DIR )
    # look for gsl library
    find_library( GSL_LIBRARY
      NAMES gsl
    )
    if( GSL_LIBRARY )
      set( GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )
      get_filename_component( GSL_LIBRARY_DIRS ${GSL_LIBRARY} PATH )
      set( GSL_FOUND ON )
    endif( GSL_LIBRARY )
 
    # look for gsl cblas library
    find_library( GSL_CBLAS_LIBRARY
        NAMES gslcblas
      )
    if( GSL_CBLAS_LIBRARY )
      set( GSL_CBLAS_FOUND ON )
    endif( GSL_CBLAS_LIBRARY )
 
    set( GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )
  endif( GSL_INCLUDE_DIR )
 
  mark_as_advanced(
    GSL_INCLUDE_DIR
    GSL_LIBRARY
    GSL_CBLAS_LIBRARY
  )
else( WIN32 AND NOT CYGWIN AND NOT MSYS )
  if( UNIX OR MSYS )
    if (GSL_CONFIG_EXECUTABLE)
      # the config executable is defined in CMakeList.txt
    else (GSL_CONFIG_EXECUTABLE)
      find_program( GSL_CONFIG_EXECUTABLE gsl-config
        /usr/bin
        /usr/local/bin
      )
    endif ()

    message(${GSL_CONFIG_EXECUTABLE})   

    if( GSL_CONFIG_EXECUTABLE )
      set( GSL_FOUND ON )
 
      # run the gsl-config program to get cxxflags
      execute_process(
        COMMAND sh "${GSL_CONFIG_EXECUTABLE}" --cflags
        OUTPUT_VARIABLE GSL_C_FLAGS
        RESULT_VARIABLE RET
        ERROR_QUIET
        )
      if( RET EQUAL 0 )
        string( STRIP "${GSL_C_FLAGS}" GSL_C_FLAGS )
        separate_arguments( GSL_C_FLAGS )
 
        # parse definitions from cflags; drop -D* from CFLAGS
        string( REGEX MATCHALL "-D[^;]+"
          GSL_DEFINITIONS  "${GSL_C_FLAGS}" )
        string( REGEX REPLACE "-D[^;]+;" ""
          GSL_C_FLAGS "${GSL_C_FLAGS}" )
 
        # parse include dirs from cflags; drop -I prefix
        string( REGEX MATCHALL "-I[^;]+"
          GSL_INCLUDE_DIRS "${GSL_C_FLAGS}" )
        string( REPLACE "-I" ""
          GSL_INCLUDE_DIRS "${GSL_INCLUDE_DIRS}")
        string( REGEX REPLACE "-I[^;]+;" ""
          GSL_C_FLAGS "${GSL_C_FLAGS}")
 
        message("GSL_DEFINITIONS=${GSL_DEFINITIONS}")
        message("GSL_INCLUDE_DIRS=${GSL_INCLUDE_DIRS}")
        message("GSL_C_FLAGS=${GSL_C_FLAGS}")
	set(GSL_CXX_FLAGS ${GSL_C_FLAGS})
        message("GSL_CXX_FLAGS=${GSL_CXX_FLAGS}")
      else( RET EQUAL 0 )
        set( GSL_FOUND FALSE )
      endif( RET EQUAL 0 )
 
      # run the gsl-config program to get the libs
      execute_process(
        COMMAND sh "${GSL_CONFIG_EXECUTABLE}" --libs
        OUTPUT_VARIABLE GSL_LIBRARIES
        RESULT_VARIABLE RET
        ERROR_QUIET
        )
      if( RET EQUAL 0 )
        string(STRIP "${GSL_LIBRARIES}" GSL_LIBRARIES )
        separate_arguments( GSL_LIBRARIES )
 
        # extract linkdirs (-L) for rpath (i.e., LINK_DIRECTORIES)
        string( REGEX MATCHALL "-L[^;]+"
          GSL_LIBRARY_DIRS "${GSL_LIBRARIES}" )
        string( REPLACE "-L" ""
          GSL_LIBRARY_DIRS "${GSL_LIBRARY_DIRS}" )
      else( RET EQUAL 0 )
        set( GSL_FOUND FALSE )
      endif( RET EQUAL 0 )
 
      MARK_AS_ADVANCED(
	CSL_CXXFLAGS
        GSL_C_FLAGS
      )
      message( STATUS "Using GSL from ${GSL_PREFIX}" )
    else( GSL_CONFIG_EXECUTABLE )
      message( STATUS "FindGSL: gsl-config not found.")
    endif( GSL_CONFIG_EXECUTABLE )
  endif( UNIX OR MSYS )
endif( WIN32 AND NOT CYGWIN AND NOT MSYS )
 
if( GSL_FOUND )
    message( STATUS "GSL found" )
  if( NOT GSL_FIND_QUIETLY )
    message( STATUS "FindGSL: Found both GSL headers and library" )
  endif( NOT GSL_FIND_QUIETLY )
else( GSL_FOUND )
  if( GSL_FIND_REQUIRED )
    message( FATAL_ERROR "FindGSL: Could not find GSL headers or library" )
  endif( GSL_FIND_REQUIRED )
endif( GSL_FOUND )
