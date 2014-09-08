# -*- mode: cmake -*-
# WELLS
cmake_minimum_required (VERSION 2.8)

enable_language(C)

set( src_wells wells.c )

add_executable ( wells ${src_wells} )
target_link_libraries( wells m )
