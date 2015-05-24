message(STATUS "Installing LibYAML (${LIBYAML_VERSION})")
SetURL( LibYAML $:{LIBYAML_ARCHIVE_FILE} LIBYAML_URL )
message(STATUS "Installing LibYAML using ${LIBYAML_URL}")

ExternalProject_Add(
	libyaml
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
	URL          ${LIBYAML_URL}
	PREFIX       libyaml-${LIBYAML_VERSION}
	CONFIGURE_COMMAND ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX} WORKING_DIRECTORY ${SOURCE_DIR}
	BUILD_COMMAND     make
	BUILD_IN_SOURCE   1
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND  make install
)
#  -*- mode: cmake -*-
