#  -*- mode: cmake -*-
message(STATUS "Installing GLIB2 (${GLIB2_VERSION})")
define_external_project_args(GLIB2 TARGET glib2)
set(SOURCE_DIR ${GLIB2_source_dir})
ExternalProject_Add(${GLIB2_BUILD_TARGET}
	DEPENDS   ${GLIB2_PACKAGE_DEPENDS}             # Package dependency target
	TMP_DIR   ${GLIB2_tmp_dir}                     # Temporary files directory
	STAMP_DIR ${GLIB2_stamp_dir}                   # Timestamp and log directory
	# -- Download and URL definitions
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
	URL          ${GLIB2_URL}                      # URL may be a web site OR a local file
	URL_MD5      ${GLIB2_MD5_SUM}                  # md5sum of the archive file
	# -- Configure
	SOURCE_DIR       ${GLIB2_source_dir}           # Source directory
	CONFIGURE_COMMAND ./autogen.sh && ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	# -- Build
	BINARY_DIR        ${GLIB2_source_dir}           # Build directory 
	BUILD_COMMAND     make prefix=${TPL_INSTALL_PREFIX}
	BUILD_IN_SOURCE   ${GLIB2_BUILD_IN_SOURCE}     # Flag for in source builds
	# -- Install
	INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
	INSTALL_COMMAND   make install prefix=${TPL_INSTALL_PREFIX}
	# -- Output control
	${GLIB2_logging_args}
)
