#  -*- mode: cmake -*-
message(STATUS "Installing GSL (${GSL_VERSION})")
SetURL( GSL ${GSL_ARCHIVE_FILE} GSL_URL )
message(STATUS "Installing GSL using ${GSL_URL}")

ExternalProject_Add(
	gsl
	BUILD_IN_SOURCE   1
	URL ${GSL_URL}
	PREFIX gsl-${GSL_VERSION}
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
	CONFIGURE_COMMAND ./autogen.sh && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	BUILD_COMMAND     make
	INSTALL_COMMAND   make install
)
