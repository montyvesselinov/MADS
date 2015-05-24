#  -*- mode: cmake -*-
message(STATUS "Installing GLib2 (${GLIB2_VERSION})")
SetURL( GLib2 ${GLIB2_ARCHIVE_FILE} GLIB2_URL )
message(STATUS "Installing GLib2 using ${GLIB2_URL}")

ExternalProject_Add(
	glib2
	BUILD_IN_SOURCE   1
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
	URL          ${GLIB2_URL}
	PREFIX       glib-${GLIB2_VERSION}
	CONFIGURE_COMMAND ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	BUILD_COMMAND     make
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND   make install
)
