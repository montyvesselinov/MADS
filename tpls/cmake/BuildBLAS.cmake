#  -*- mode: cmake -*-
message(STATUS "Installing BLAS (${BLAS_VERSION})")
define_external_project_args(BLAS TARGET blas)
set(SOURCE_DIR ${BLAS_source_dir})
ExternalProject_Add(${BLAS_BUILD_TARGET}
	DEPENDS          ${BLAS_PACKAGE_DEPENDS}  # Package dependency target
	TMP_DIR          ${BLAS_tmp_dir}
	STAMP_DIR        ${BLAS_stamp_dir}
	# -- Download and URL definitions
	SVN_REPOSITORY    ${BLAS_SVN}
	URL_MD5           ${BLAS_MD5_SUM}
	# -- Configure
	SOURCE_DIR        ${BLAS_source_dir}
    CONFIGURE_COMMAND cd ${BLAS_source_dir} && cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=${TPL_INSTALL_PREFIX}
	# -- Build
	BINARY_DIR        ${BLAS_build_dir}
    BUILD_COMMAND     cd ${BLAS_source_dir} && make prefix=${TPL_INSTALL_PREFIX}
	BUILD_IN_SOURCE   ${BLAS_BUILD_IN_SOURCE}
	# -- Install
	INSTALL_DIR       ${TPL_INSTALL_PREFIX}
    INSTALL_COMMAND   cd ${BLAS_source_dir} && make install prefix=${TPL_INSTALL_PREFIX}
	# -- Output control
	#${BLAS_logging_args}
)
