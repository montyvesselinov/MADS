find_path(LibMathEval_INCLUDEDIR
    NAMES matheval.h HINTS
    )

set(LibMathEval_INCLUDE_DIRS ${Matheval_INCLUDEDIR})

find_library( LibMathEval_LIBRARIES
    NAMES matheval HINTS
    ${Matheval_LIBDIR}
    ${Matheval_LIBRARY_DIRS}
    )

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( LibMathEval DEFAULT_MSG LibMathEval_LIBRARIES LibMathEval_INCLUDE_DIRS)
mark_as_advanced(LibMathEval_INCLUDE_DIRS LibMathEval_LIBRARIES)
