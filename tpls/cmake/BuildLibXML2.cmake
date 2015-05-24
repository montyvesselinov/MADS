#  -*- mode: cmake -*-
message(STATUS "Installing LibXML2 (${LIBXML2_VERSION})")
SetURL( LibXML2 ${LIBXML2_ARCHIVE_FILE} LIBXML2_URL )
message(STATUS "Installing LibXML2 using ${LIBXML2_URL}")

ExternalProject_Add(
	libxml2
	URL           ${LIBXML2_URL}
    DOWNLOAD_DIR  ${TPL_DOWNLOAD_DIR}
	PREFIX        libxml-${LIBXML2_VERSION}
	CONFIGURE_COMMAND ./autogen.sh --prefix=${CMAKE_BINARY_DIR}/tpls --without-python
	BUILD_IN_SOURCE 1
)
