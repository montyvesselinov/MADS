# -*- mode: cmake -*-
# MADS Build Options
enable_language(C)
enable_language(CXX)

# check to use if using GSL
# DEFAULT: USE_GSL = FALSE
OPTION(USE_GSL "USE_GSL." FALSE)
if (USE_GSL)
    add_definitions(-DUSE_GSL)
endif (USE_GSL)

# Build static libraries
# DEFAULT: BUILD_STATIC_EXECUTABLES = TRUE
OPTION(BUILD_STATIC_EXECUTABLES "STATIC" ON)
if(PREFER_STATIC_LIBRARIES)
    # Prefer static libraries, but don't require that everything must be static. 
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
endif(PREFER_STATIC_LIBRARIES)

if(BUILD_STATIC_EXECUTABLES)
    set(CMAKE_EXE_LINKER_FLAGS -Bstatic)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(BUILD_STATIC_EXECUTABLES)

OPTION(ENABLE_DBC "Enable DBC checking" OFF)
if (ENABLE_DBC)
    message("Enabling DBC checks")
    add_definitions(-DENABLE_DBC)
endif()
