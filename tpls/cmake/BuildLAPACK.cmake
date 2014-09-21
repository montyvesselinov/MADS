#  -*- mode: cmake -*-
message(STATUS "Installing LAPACK (${LAPACK_VERSION})")
define_external_project_args(LAPACK TARGET lapack)
set(SOURCE_DIR ${LAPACK_source_dir})
ExternalProject_Add(${LAPACK_BUILD_TARGET}
	DEPENDS          ${LAPACK_PACKAGE_DEPENDS}  # Package dependency target
	TMP_DIR          ${LAPACK_tmp_dir}
	STAMP_DIR        ${LAPACK_stamp_dir}
	# -- Download and URL definitions
	DOWNLOAD_DIR      ${TPL_DOWNLOAD_DIR}
	URL               ${LAPACK_URL}
	URL_MD5           ${LAPACK_MD5_SUM}
	# -- Configure
	SOURCE_DIR        ${LAPACK_source_dir}
	CONFIGURE_COMMAND cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=${TPL_INSTALL_PREFIX}
	# -- Build
	BINARY_DIR        ${LAPACK_source_dir}
	BUILD_COMMAND     make prefix=${TPL_INSTALL_PREFIX}
	BUILD_IN_SOURCE   ${LAPACK_BUILD_IN_SOURCE}
	# -- Install
	INSTALL_DIR       ${TPL_INSTALL_PREFIX}
	INSTALL_COMMAND   make install prefix=${TPL_INSTALL_PREFIX}
	# -- Output control
	#${LAPACK_logging_args}
)
