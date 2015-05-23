#  -*- mode: cmake -*-
message(STATUS "Installing GSL (${GSL_VERSION})")
if( EXISTS ${TPL_REPO_DIR}/${GSL_ARCHIVE_FILE} )
	message(STATUS "GSL is available in the TPL's repository for MADS (${TPL_REPO_DIR}/${GSL_ARCHIVE_FILE})" )
	set(GSL_URL ${TPL_REPO_DIR}/${GSL_ARCHIVE_FILE})
elseif( EXISTS ${TPL_DOWNLOAD_DIR}/${GSL_ARCHIVE_FILE} )
	message(STATUS "GSL already downloaded (${TPL_DOWNLOAD_DIR}/${GSL_ARCHIVE_FILE})" )
	set(GSL_URL ${TPL_DOWNLOAD_DIR}/${GSL_ARCHIVE_FILE})
else()
	message(STATUS "GSL will be downloaded" )
endif()

ExternalProject_Add(
	gsl
	BUILD_IN_SOURCE   1
	URL ${GSL_URL}
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
	CONFIGURE_COMMAND ./autogen.sh && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	BUILD_COMMAND     make
	INSTALL_COMMAND   make install
)
