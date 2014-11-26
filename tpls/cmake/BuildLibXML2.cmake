#  -*- mode: cmake -*-
message(STATUS "Installing LibXL2 (${LIBXML2_VERSION})")
ExternalProject_Add(
	libxml2
	GIT_REPOSITORY git://git.gnome.org/libxml2
	CONFIGURE_COMMAND ./autogen.sh --prefix=${CMAKE_BINARY_DIR}/tpls --without-python
	BUILD_IN_SOURCE 1
)
