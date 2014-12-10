#  -*- mode: cmake -*-
message(STATUS "Installing GLIB2 (${GLIB2_VERSION})")
ExternalProject_Add(
	glib2
	BUILD_IN_SOURCE   1
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
	URL          ${GLIB2_URL}                      # URL may be a web site OR a local file
	CONFIGURE_COMMAND ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	BUILD_COMMAND     make
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
	INSTALL_COMMAND   make install
)
