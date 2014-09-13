FIND_PATH(Matheval_INCLUDEDIR
    NAMES matheval.h HINTS
    )

SET(Matheval_INCLUDE_DIRS ${Matheval_INCLUDEDIR})

FIND_LIBRARY( Matheval_LIBRARIES
    NAMES matheval HINTS
    ${Matheval_LIBDIR}
    ${Matheval_LIBRARY_DIRS}
    )

INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Matheval DEFAULT_MSG Matheval_LIBRARIES Matheval_INCLUDE_DIRS)
MARK_AS_ADVANCED(Matheval_INCLUDE_DIRS Matheval_LIBRARIES)
