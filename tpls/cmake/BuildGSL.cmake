#  -*- mode: cmake -*-
message(STATUS "Installing GSL (${GSL_VERSION})")
ExternalProject_Add(
	gsl
	BUILD_IN_SOURCE   1
	URL ${GSL_URL}
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
	CONFIGURE_COMMAND ./autogen.sh && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	BUILD_COMMAND     make
	INSTALL_COMMAND   make install
)
