macro(build_netcdf install_prefix staging_prefix)
    ExternalProject_Add(
	netcdf
	URL "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.2.tar.gz"
	BUILD_IN_SOURCE 1
	CONFIGURE_COMMAND ./configure --prefix=${install_prefix} --with-pic --disable-doxygen --disable-hdf4 --disable-netcdf-4 --disable-shared --disable-dap --libdir=${install_prefix}/lib${LIB_SUFFIX} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} "CXXFLAGS=${EXT_CXX_FLAGS}" "CFLAGS=${EXT_C_FLAGS}"
	BUILD_COMMAND make
	INSTALL_COMMAND make DESTDIR=${staging_prefix} install
    )

    SET(NETCDF_LIBRARY ${staging_prefix}/${install_prefix}/lib${LIB_SUFFIX}/libnetcdf.a )
    SET(NETCDF_INCLUDE_DIR ${staging_prefix}/${install_prefix}/include )
    SET(NETCDF_FOUND ON)
endmacro(build_netcdf)
