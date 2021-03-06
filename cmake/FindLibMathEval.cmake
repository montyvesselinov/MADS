PKG_CHECK_MODULES(PC_LibMathEval libmatheval)
SET(LibMathEval_DEFINITIONS ${PC_LibMathEval_CFLAGS_OTHER})

FIND_PATH(LibMathEval_INCLUDE_DIR NAMES matheval.h
   HINTS
   ${PC_LibMathEval_INCLUDEDIR}
   ${PC_LibMathEval_INCLUDE_DIRS}
   PATH_SUFFIXES matheval
   )

FIND_LIBRARY(LibMathEval_LIBRARIES NAMES matheval
   HINTS
   ${PC_LibMathEval_LIBDIR}
   ${PC_LibMathEval_LIBRARY_DIRS}
   )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibMathEval DEFAULT_MSG LibMathEval_LIBRARIES)
MARK_AS_ADVANCED(LibMathEval_INCLUDE_DIR LibMathEval_LIBRARIES LibMathEval_XMLLINT_EXECUTABLE)
