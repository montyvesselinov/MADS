# -*- mode: cmake -*-
# MADS SRC
cmake_minimum_required (VERSION 3.0 FATAL_ERROR)

add_custom_target( gitversion
    DEPENDS mads_gitversion.c )

add_custom_command(TARGET gitversion
    COMMAND echo 'const char *gitversion = \"${VER}${GIT_STATUS}\"\;' > mads_gitversion.c 
	DEPENDS ../.git/HEAD ../.git/index )

set( src_mads ./mads.c ./mads_io.c ./mads_io_external.c ./mads_func.c
	./mads_mem.c ./mads_info.c ./mads_gitversion.c lm/opt_lm_mon.c lm/opt_lm_gsl.c lm/lu.c lm/opt_lm_ch.c misc/test_problems.c misc/anasol_contamination.c misc/io.c lhs/lhs.c )
set( src_pso pso/pso-tribes-lm.c pso/Standard_PSO_2006.c pso/mopso.c )
set( src_sa sa/abagus.c sa/postpua.c sa/global.c sa/do_miser.c )
set( src_ds ds/infogap.c ds/glue.c )
set( src_mpun mprun/mprun.c mprun/mprun_io.c )
set( src_kdtree misc/kdtree-0.5.5/kdtree.c )
set( src_levmar misc/levmar-2.5/lm_m.c misc/levmar-2.5/Axb.c misc/levmar-2.5/misc.c misc/levmar-2.5/lmlec.c misc/levmar-2.5/lmbc.c misc/levmar-2.5/lmblec.c misc/levmar-2.5/lmbleic.c )
set( src_astable misc/astable/astable.c misc/astable/interpolation.c misc/astable/pqueue.c )
set( src_bayes bayes/dream.cpp )
set( sources ${src_mads} ${src_pso} ${src_sa} ${src_ds} ${src_mpun} ${src_kdtree} ${src_levmar} ${src_astable} ${src_bayes} )

add_executable ( mads ${sources} )
add_dependencies( mads gitversion lapack)
target_link_libraries( mads ${MADS_LIBS} )
