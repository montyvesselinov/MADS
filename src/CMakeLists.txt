# -*- mode: cmake -*-
# MADS SRC
cmake_minimum_required (VERSION 2.8)

include_directories(${CMAKEDEMO_SOURCE_DIR}/src)
link_diretories(${CMAKEDEMO_BINARY_DIR}/src)

enable_language(C)
enable_language(CXX)

add_subdirectory (lm)
add_subdirectory (lhs)
add_subdirectory (misc)
add_subdirectory (mprun)
add_subdirectory (ds)
add_subdirectory (sa)
add_subdirectory (bayes)

set( src_mads ./mads.c ./mads_io.c ./mads_io_external.c ./mads_func.c ./mads_mem.c ./mads_info.c lm/opt_lm_mon.c lm/opt_lm_gsl.c lm/lu.c lm/opt_lm_ch.c misc/test_problems.c misc/anasol_con
set( src_pso pso/pso-tribes-lm.c pso/Standard_PSO_2006.c pso/mopso.c )
set( src_sa sa/abagus.c sa/postpua.c sa/global.c sa/do_miser.c )
set( src_ds ds/infogap.c ds/glue.c )
set( src_mpun mprun/mprun.c mprun/mprun_io.c )
set( src_kdtree misc/kdtree-0.5.5/kdtree.c )
set( src_levmar misc/levmar-2.5/lm_m.c misc/levmar-2.5/Axb.c misc/levmar-2.5/misc.c misc/levmar-2.5/lmlec.c misc/levmar-2.5/lmbc.c misc/levmar-2.5/lmblec.c misc/levmar-2.5/lmbleic.c )
set( src_astable misc/astable/astable.c misc/astable/interpolation.c misc/astable/pqueue.c )
set( src_bayes bayes/dream.c )
set( sorces ${src_mads} ${src_pso} ${src_a} ${objds} ${src_mpun} ${src_kdtree} ${src_levmar} ${src_astable} ${src_bayes} )

add_executable ( mads ${sources} )
target_link_libraries(mads m) #link the math library

foreach (file ${SRC_FILES})
    get_filename_component(EXE_NAME ${file} NAME_WE)
        get_filename_component(SRC_NAME ${file} NAME)
	  message(${SRC_NAME})
	    add_executable(${EXE_NAME} ${SRC_NAME})
	      target_link_libraries(${EXE_NAME} gslib)
	        install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE_NAME} DESTINATION bin)
	    endforeach()
