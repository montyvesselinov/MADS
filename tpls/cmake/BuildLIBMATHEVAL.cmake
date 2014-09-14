#  -*- mode: cmake -*-
message(STATUS "Installing LibMathEvaL (${LIBMATHEVAL_VERSION})")
define_external_project_args(LIBMATHEVAL TARGET libmatheval)
set(SOURCE_DIR ${LIBMATHEVAL_source_dir})
ExternalProject_Add(${LIBMATHEVAL_BUILD_TARGET}
                    DEPENDS   ${LIBMATHEVAL_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${LIBMATHEVAL_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${LIBMATHEVAL_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${LIBMATHEVAL_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${LIBMATHEVAL_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${LIBMATHEVAL_source_dir}           # Source directory
		    CONFIGURE_COMMAND ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=${TPL_INSTALL_PREFIX} WORKING_DIRECTORY ${SOURCE_DIR}
                    # -- Build
                    BINARY_DIR        ${LIBMATHEVAL_build_dir}           # Build directory 
		    BUILD_COMMAND     make WORKING_DIRECTORY ${SOURCE_DIR}
                    BUILD_IN_SOURCE   ${LIBMATHEVAL_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
		    INSTALL_COMMAND  make install WORKING_DIRECTORY ${SOURCE_DIR}
                    # -- Output control
                    ${LIBMATHEVAL_logging_args}
		   )
