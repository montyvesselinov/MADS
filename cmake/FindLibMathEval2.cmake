find_path(LibMathEval_INCLUDE
    NAMES matheval.h HINTS
	$ENV{MATHEVALROOT}/src/
	/opt/local/include
	/usr/local/include
	/usr/include
    )

find_library( LibMathEval_LIB
    NAMES matheval HINTS
	/opt/local/lib
	/usr/local/lib
	/usr/lib
    ${LibMatheval_LIBDIR}
    ${LibMatheval_LIBRARY_DIRS}
	PATH ${LibMathEval_LIBRARIES}
    )

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( LibMathEval DEFAULT_MSG LibMathEval_LIB LibMathEval_INCLUDE)
MARK_AS_ADVANCED( LibMathEval_INCLUDE_DIR LibMathEval_LIBRARIES )
