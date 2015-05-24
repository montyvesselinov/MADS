#  -*- mode: cmake -*-
#function(SetURL PACKAGE ARCHIVE_FILE URL )
#if( EXISTS ${TPL_REPO_DIR}/${ARCHIVE_FILE} )
#    message(STATUS "${PACKAGE} is available in the TPL's repository for MADS (${TPL_REPO_DIR}/${ARCHIVE_FILE})" )
#	set(${URL} ${TPL_REPO_DIR}/${ARCHIVE_FILE} PARENT_SCOPE)
#elseif( EXISTS ${TPL_DOWNLOAD_DIR}/${ARCHIVE_FILE} )
#    message(STATUS "${PACKAGE} already downloaded (${TPL_DOWNLOAD_DIR}/${ARCHIVE_FILE})" )
#	set(${URL} ${${TPL_DOWNLOAD_DIR}/${ARCHIVE_FILE}} PARENT_SCOPE)
#else()
#    message(STATUS "${PACKAGE} will be downloaded (${URL})" )
#endif()
#endfunction(SetURL)

message(STATUS "Installing LibMathEval (${LIBMATHEVAL_VERSION})")
SetURL( LibMathEval ${LIBMATHEVAL_ARCHIVE_FILE} LIBMATHEVAL_URL )
message(STATUS "Installing LibMathEval using ${LIBMATHEVAL_URL}")

#function(twice varName)
#set(${varName} ${${varName}}${${varName}} PARENT_SCOPE)
#MESSAGE(STATUS "IN ${varName}")
#endfunction()
#
#SET(arg "foo")
#twice(arg)
#MESSAGE(STATUS ${arg})
#SET(URL2 "poop")
#twice(URL2)
#MESSAGE(STATUS ${URL2})

ExternalProject_Add(
	libmatheval
	URL           ${LIBMATHEVAL_URL}
	DOWNLOAD_DIR  ${TPL_DOWNLOAD_DIR}
	PREFIX        libmatheval-${LIBMATHEVAL_VERSION}
	CONFIGURE_COMMAND autoreconf -fi && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX} WORKING_DIRECTORY ${SOURCE_DIR}
	BUILD_COMMAND     make
	BUILD_IN_SOURCE   1
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND  make install
)
