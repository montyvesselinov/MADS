macro(build_gsl install_prefix staging_prefix)
    IF(${CMAKE_BUILD_TYPE} STREQUAL Release)
	SET(EXT_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
	SET(EXT_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    ELSE()
	SET(EXT_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
	SET(EXT_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
    ENDIF()

    ExternalProject_Add(libxml2
	DEPENDS zlib
	GIT_REPOSITORY git://git.gnome.org/libxml2
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ./autogen.sh --prefix=${CMAKE_BINARY_DIR} --without-python --disable-shared --without-zlib
	BUILD_IN_SOURCE 1
    )
endmacro(build_gsl)