# - Find libyaml libraries
# This module finds libyaml if it is installed and determines where the
# include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#
#  LIBYAML_FOUND               - have the libyaml libs been found
#  LIBYAML_LIBRARIES           - path to the libyaml library
#  LIBYAML_INCLUDE_DIRS        - path to where zmq.h is found

FIND_LIBRARY(LIBYAML_LIBRARY
  NAMES libyaml yaml
  PATHS
    ${LIBYAML_LIBRARIES}
)

FIND_PATH(LIBYAML_INCLUDE_DIR
  NAMES yaml.h
  PATHS
#   ${libyaml_FRAMEWORK_INCLUDES}
    ${LIBYAML_INCLUDE_DIRS}
    ${LIBYAML_INCLUDE_DIR}
)

MARK_AS_ADVANCED(
  LIBYAML_LIBRARY
  LIBYAML_INCLUDE_DIR
)
SET(LIBYAML_INCLUDE_DIRS "${LIBYAML_INCLUDE_DIR}")
SET(LIBYAML_LIBRARIES "${LIBYAML_LIBRARY}")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibYAML DEFAULT_MSG LIBYAML_LIBRARIES LIBYAML_INCLUDE_DIRS)
