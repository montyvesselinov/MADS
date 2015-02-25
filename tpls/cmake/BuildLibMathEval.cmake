#  -*- mode: cmake -*-
message(STATUS "Installing LibMathEvaL (${LIBMATHEVAL_VERSION})")
ExternalProject_Add(
	libmatheval
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
	URL          ${LIBMATHEVAL_URL}
	CONFIGURE_COMMAND autoreconf -fi && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX} WORKING_DIRECTORY ${SOURCE_DIR}
	BUILD_COMMAND     make
	BUILD_IN_SOURCE   1
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND  make install
)
