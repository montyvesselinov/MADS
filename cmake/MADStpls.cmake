# -*- mode: cmake -*-
# MADS TPL (Third Party Library) Definitions

# Required Libraries

##############################################################################
# GSL: GNU Scientific Library
##############################################################################
if (USE_GSL)
   if (NOT DEFINED GSL_CONFIG_EXECUTABLE)
      if (DEFINED TPL_INSTALL_PREFIX)
      	 set(GSL_CONFIG_EXECUTABLE "${TPL_INSTALL_PREFIX}/bin/gsl-config")
      endif (DEFINED TPL_INSTALL_PREFIX)
   endif(NOT DEFINED GSL_CONFIG_EXECUTABLE)
   FIND_PACKAGE(GSL REQUIRED)
   if (NOT GSL_FOUND)
       MESSAGE("Could not find GSL! MADS will rebuild GSL ...")
   else (NOT GSL_FOUND)    
      MESSAGE("\nFound GSL!  Here are the details: ")
      MESSAGE("   CMAKE_GSL_CXX_FLAGS = ${CMAKE_GSL_CXX_FLAGS}")
      MESSAGE("   GSL_INCLUDE_DIRS = ${GSL_INCLUDE_DIRS}")
      MESSAGE("   GSL_LINK_DIRECTORIES = ${GSL_LIBRARY_DIRS}")
      MESSAGE("   GSL_LIBRARIES = ${GSL_LIBRARIES}")
      MESSAGE("   GSL_EXE_LINKER_FLAGS  = ${GSL_EXE_LINKER_FLAGS}")
      include_directories(${GSL_INCLUDE_DIRS})
      link_directories(${GSL_LIBRARY_DIRS}) 
   endif (NOT GSL_FOUND)  
endif (USE_GSL)
