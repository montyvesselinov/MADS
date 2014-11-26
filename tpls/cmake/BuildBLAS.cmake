#  -*- mode: cmake -*-
message(STATUS "Installing BLAS (${BLAS_VERSION})")
ExternalProject_Add(
	blas
	BUILD_IN_SOURCE   1
	SVN_REPOSITORY    ${BLAS_SVN}
    CONFIGURE_COMMAND cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=${TPL_INSTALL_PREFIX}
    BUILD_COMMAND     make prefix=${TPL_INSTALL_PREFIX}
	INSTALL_DIR       ${TPL_INSTALL_PREFIX}
    INSTALL_COMMAND   make install prefix=${TPL_INSTALL_PREFIX}
)
