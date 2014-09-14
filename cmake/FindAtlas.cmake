SET(ATLAS_POSSIBLE_INCLUDE_PATHS
/usr/include
/usr/include/atlas
/usr/local/include
$ENV{ATLAS_DIR}
$ENV{ATLAS_DIR}/include
)
# Ubuntu's package management does not handle blas elegantly, causing
# many possible locations here.
SET(ATLAS_POSSIBLE_LIBRARY_PATHS
/usr/lib/libatlas-corei7sse3
/usr/lib/atlas-amd64sse3
/usr/lib/atlas-base
/usr/lib/sse2
/usr/lib/sse
/usr/local/lib/sse2
/usr/local/lib/sse
/usr/lib
/usr/local/lib
$ENV{ATLAS_DIR}
$ENV{ATLAS_DIR}/lib
)
FIND_PATH(ATLAS_CBLAS_INCLUDE_DIR NAMES cblas.h PATHS ${ATLAS_POSSIBLE_INCLUDE_PATHS})
FIND_PATH(ATLAS_CLAPACK_INCLUDE_DIR NAMES clapack.h PATHS ${ATLAS_POSSIBLE_INCLUDE_PATHS})
FIND_LIBRARY(ATLAS_CBLAS_LIBRARY NAMES ptcblas_r ptcblas cblas_r cblas PATHS ${ATLAS_POSSIBLE_LIBRARY_PATHS})
FIND_LIBRARY(ATLAS_ATLAS_LIBRARY NAMES atlas_r atlas PATHS ${ATLAS_POSSIBLE_LIBRARY_PATHS})
FIND_LIBRARY(ATLAS_LAPACK_ATLAS_LIBRARY NAMES alapack_r alapack lapack_atlas PATHS ${ATLAS_POSSIBLE_LIBRARY_PATHS})
SET(ATLAS_FOUND ON)
FOREACH(INCDIR ATLAS_CBLAS_INCLUDE_DIR ATLAS_CLAPACK_INCLUDE_DIR)
IF(${INCDIR})
SET(ATLAS_INCLUDE_DIR ${ATLAS_INCLUDE_DIR} ${${INCDIR}} )
ELSE(${INCDIR})
    MESSAGE(STATUS "${INCDIR} not found turning off ATLAS_FOUND")
SET(OPENCV_FOUND OFF)
ENDIF (${INCDIR})
ENDFOREACH(INCDIR)
# could make LAPACK_ATLAS optional
FOREACH(LIBNAME ATLAS_LAPACK_ATLAS_LIBRARY ATLAS_CBLAS_LIBRARY ATLAS_ATLAS_LIBRARY)
IF(${LIBNAME})
SET(ATLAS_LIBRARIES ${ATLAS_LIBRARIES} ${${LIBNAME}} )
ELSE(${LIBNAME})
    MESSAGE(STATUS "${LIBNAME} not found turning off ATLAS_FOUND")
SET(ATLAS_FOUND OFF)
ENDIF (${LIBNAME})
ENDFOREACH(LIBNAME)
IF (ATLAS_FOUND)
IF (NOT Atlas_FIND_QUIETLY)
MESSAGE(STATUS "Found Atlas libraries: ${ATLAS_LIBRARIES}")
MESSAGE(STATUS "Found Atlas include: ${ATLAS_INCLUDE_DIR}")
ENDIF (NOT Atlas_FIND_QUIETLY)
ELSE (ATLAS_FOUND)
IF (Atlas_FIND_REQUIRED)
MESSAGE(FATAL_ERROR "Could not find Atlas")
ENDIF (Atlas_FIND_REQUIRED)
ENDIF (ATLAS_FOUND)
MARK_AS_ADVANCED(
ATLAS_INCLUDE_DIR
ATLAS_CBLAS_INCLUDE_DIR
ATLAS_CLAPACK_INCLUDE_DIR
ATLAS_LIBRARIES
ATLAS_CBLAS_LIBRARY
ATLAS_ATLAS_LIBRARY
ATLAS_LAPACK_ATLAS_LIBRARY
)
set(Atlas_FOUND ${ATLAS_FOUND})
