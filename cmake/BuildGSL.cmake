#  -*- mode: cmake -*-
message(STATUS "Installing GSL (${GSL_VERSION})")
define_external_project_args(GSL TARGET gsl)
set(SOURCE_DIR ${GSL_source_dir})
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/configure-step.cmake.in ${GSL_prefix_dir}/gsl-configure-step.cmake @ONLY)
set(GSL_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${GSL_prefix_dir}/gsl-configure-step.cmake)
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/build-step.cmake.in ${GSL_prefix_dir}/gsl-build-step.cmake @ONLY)
set(GSL_BUILD_COMMAND ${CMAKE_COMMAND} -P ${GSL_prefix_dir}/gsl-build-step.cmake)     
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/install-step.cmake.in ${GSL_prefix_dir}/gsl-install-step.cmake @ONLY)
set(GSL_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${GSL_prefix_dir}/gsl-install-step.cmake)     
ExternalProject_Add(${GSL_BUILD_TARGET}
                    DEPENDS   ${GSL_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${GSL_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${GSL_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${GSL_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${GSL_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${GSL_source_dir}           # Source directory
                    CONFIGURE_COMMAND ${GSL_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${GSL_build_dir}           # Build directory 
                    BUILD_COMMAND     ${GSL_BUILD_COMMAND}       # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${GSL_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    INSTALL_COMMAND  ${GSL_INSTALL_COMMAND}
                    # -- Output control
                    ${GSL_logging_args}
		   )
