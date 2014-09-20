#  -*- mode: cmake -*-
message(STATUS "Installing GLIB2 (${GLIB2_VERSION})")
define_external_project_args(GLIB2 TARGET glib2)
set(SOURCE_DIR ${GLIB2_source_dir})
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/configure-step.cmake.in ${GLIB2_prefix_dir}/glib2-configure-step.cmake @ONLY)
set(GLIB2_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${GLIB2_prefix_dir}/glib2-configure-step.cmake)
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/build-step.cmake.in ${GLIB2_prefix_dir}/glib2-build-step.cmake @ONLY)
set(GLIB2_BUILD_COMMAND ${CMAKE_COMMAND} -P ${GLIB2_prefix_dir}/glib2-build-step.cmake)     
configure_file(${TPL_CMAKE_TEMPLATE_DIR}/install-step.cmake.in ${GLIB2_prefix_dir}/glib2-install-step.cmake @ONLY)
set(GLIB2_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${GLIB2_prefix_dir}/glib2-install-step.cmake)     
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
                    CONFIGURE_COMMAND ${GLIB2_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${GLIB2_build_dir}           # Build directory 
                    BUILD_COMMAND     ${GLIB2_BUILD_COMMAND}       # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${GLIB2_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    INSTALL_COMMAND  ${GLIB2_INSTALL_COMMAND}
                    # -- Output control
                    ${GLIB2_logging_args}
		   )
