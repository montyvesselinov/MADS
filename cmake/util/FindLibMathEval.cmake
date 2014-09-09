# use pkg-config to get the directories and then use these values
# in the FIND_PATH() and FIND_LIBRARY() calls
  find_path(Matheval_INCLUDEDIR
    NAMES matheval.h HINTS
	   )
  set (Matheval_INCLUDE_DIRS ${Matheval_INCLUDEDIR})

  find_library( Matheval_LIBRARIES
    NAMES matheval HINTS
    ${Matheval_LIBDIR}
    ${Matheval_LIBRARY_DIRS}
	      )

# handle the QUIETLY and REQUIRED arguments and set Matheval_FOUND to
# TRUE if 
# all listed variables are TRUE
  include( FindPackageHandleStandardArgs )
  find_package_handle_standard_args( Matheval DEFAULT_MSG
          Matheval_LIBRARIES Matheval_INCLUDE_DIRS)

# mark_as_advanced(Matheval_INCLUDE_DIRS Matheval_LIBRARIES)
