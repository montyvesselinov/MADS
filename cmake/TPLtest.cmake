#  -*- mode: cmake -*-

set(TPL_MISSING FALSE)

if(USE_GSL AND NOT GSL_FOUND)
set(TPL_MISSING TRUE)
endif()

if(USE_LAPACK AND NOT LAPACK_FOUND)
set(TPL_MISSING TRUE)
endif()

if(USE_BLAS AND NOT BLAS_FOUND)
set(TPL_MISSING TRUE)
endif()

if(USE_ATLAS AND NOT ATLAS_FOUND)
set(TPL_MISSING TRUE)
endif()

if(USE_LIBXML2 AND NOT LIBXML2_FOUND)
set(TPL_MISSING TRUE)
endif()

if(USE_MATHEVAL AND NOT MATHEVAL_FOUND)
set(TPL_MISSING TRUE)
endif()

if(TPL_MISSING)
    message(STATUS "TPL's are missing!" )
else()
    message(STATUS "TPL's are OK!" )
endif()