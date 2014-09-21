#  -*- mode: cmake -*-
message(STATUS "Installing GSL (${GSL_VERSION})")
define_external_project_args(GSL TARGET gsl)
set(SOURCE_DIR ${GSL_source_dir})
ExternalProject_Add(${GSL_BUILD_TARGET}
	DEPENDS          ${GSL_PACKAGE_DEPENDS}  # Package dependency target
	TMP_DIR          ${GSL_tmp_dir}
	STAMP_DIR        ${GSL_stamp_dir}
	# -- Download and URL definitions
	DOWNLOAD_DIR      ${TPL_DOWNLOAD_DIR}
	URL               ${GSL_URL}
	URL_MD5           ${GSL_MD5_SUM}
	# -- Configure
	SOURCE_DIR        ${GSL_source_dir}
	CONFIGURE_COMMAND configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX}
	# -- Build
	BINARY_DIR        ${GSL_source_dir}	
	BUILD_COMMAND     make
	BUILD_IN_SOURCE   ${GSL_BUILD_IN_SOURCE}
	# -- Install
	INSTALL_DIR       ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND   make install
	# -- Output control
	#${GSL_logging_args}
)
