#  -*- mode: cmake -*-
message(STATUS "Installing LAPACK (${LAPACK_VERSION})")
ExternalProject_Add(
	lapack
	BUILD_IN_SOURCE   1
	DOWNLOAD_DIR      ${TPL_DOWNLOAD_DIR}
	URL               ${LAPACK_URL}
    CONFIGURE_COMMAND cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=${TPL_INSTALL_PREFIX}
    BUILD_COMMAND     make prefix=${TPL_INSTALL_PREFIX}
	INSTALL_DIR       ${TPL_INSTALL_PREFIX}
    INSTALL_COMMAND   make install prefix=${TPL_INSTALL_PREFIX}
)
