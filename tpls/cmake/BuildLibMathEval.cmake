#  -*- mode: cmake -*-
message(STATUS "Installing LibMathEvaL (${LIBMATHEVAL_VERSION})")
if( EXISTS ${TPL_DOWNLOAD_DIR}/${LIBMATHEVAL_ARCHIVE_FILE} )
	message(STATUS "LibMathEvaL already downloaded" )
else()
	message(STATUS "LibMathEvaL will be downloaded" )
	ExternalProject_Add(
		libmathevaldownload
		DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
		URL          ${LIBMATHEVAL_URL}
	)
endif()

ExternalProject_Add(
	libmatheval
	SOURCE_DIR    ${TPL_DOWNLOAD_DIR}
	DOWNLOAD_NAME ${LIBMATHEVAL_ARCHIVE_FILE}
	CONFIGURE_COMMAND autoreconf -fi && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX} WORKING_DIRECTORY ${SOURCE_DIR}
	BUILD_COMMAND     make
	BUILD_IN_SOURCE   1
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND  make install
)
