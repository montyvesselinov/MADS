#  -*- mode: cmake -*-
message(STATUS "Installing LibMathEvaL (${LIBMATHEVAL_VERSION})")
define_external_project_args(LIBMATHEVAL TARGET libmatheval)
set(SOURCE_DIR ${LIBMATHEVAL_source_dir})
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/configure-step.cmake.in ${LIBMATHEVAL_prefix_dir}/libmatheval-configure-step.cmake @ONLY)
set(LIBMATHEVAL_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${LIBMATHEVAL_prefix_dir}/libmatheval-configure-step.cmake)
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/build-step.cmake.in ${LIBMATHEVAL_prefix_dir}/libmatheval-build-step.cmake @ONLY)
set(LIBMATHEVAL_BUILD_COMMAND ${CMAKE_COMMAND} -P ${LIBMATHEVAL_prefix_dir}/libmatheval-build-step.cmake)     
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/install-step.cmake.in ${LIBMATHEVAL_prefix_dir}/libmatheval-install-step.cmake @ONLY)
set(LIBMATHEVAL_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${LIBMATHEVAL_prefix_dir}/libmatheval-build-step.cmake)     
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
                    CONFIGURE_COMMAND ${LIBMATHEVAL_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${LIBMATHEVAL_build_dir}           # Build directory 
                    BUILD_COMMAND     ${LIBMATHEVAL_BUILD_COMMAND}       # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${LIBMATHEVAL_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    INSTALL_COMMAND  ${LIBMATHEVAL_INSTALL_COMMAND}
                    # -- Output control
                    ${LIBMATHEVAL_logging_args}
		   )
