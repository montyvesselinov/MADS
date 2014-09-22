#  -*- mode: cmake -*-
message(STATUS "Installing LibXL2 (${LIBXML2_VERSION})")
IF(${CMAKE_BUILD_TYPE} STREQUAL Release)
	SET(EXT_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
	SET(EXT_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
ELSE()
	SET(EXT_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
	SET(EXT_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
ENDIF()
define_external_project_args(LIBXML2 TARGET libxml2)
ExternalProject_Add(libxml2
        DEPENDS          ${LIBXML2_PACKAGE_DEPENDS}  # Package dependency target
        TMP_DIR          ${LIBXML2_tmp_dir}
        STAMP_DIR        ${LIBXML2_stamp_dir}
        # -- Download and URL definitions
	GIT_REPOSITORY git://git.gnome.org/libxml2
        DOWNLOAD_DIR      ${TPL_DOWNLOAD_DIR}
        URL               ${LIBXML2_URL}
        URL_MD5           ${LIBXML2_MD5_SUM}
        # -- Configure
        SOURCE_DIR        ${LIBXML2_source_dir}
        CONFIGURE_COMMAND ./autogen.sh CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX} --without-python
        # -- Build
        BINARY_DIR        ${LIBXML2_source_dir}
        BUILD_COMMAND     make
        BUILD_IN_SOURCE   ${LIBXML2_BUILD_IN_SOURCE}
        # -- Install
        INSTALL_DIR       ${TPL_INSTALL_PREFIX}
        INSTALL_COMMAND   make install
        # -- Output control
        #${LIBXML2_logging_args}
)
