PROG = mads
CMP = ./compare-results
# CMP = cp -f
# CC = g++ 
ifeq ($(OSTYPE),linux)
        DIRS = -I/home/monty/local/include -L/home/monty/local/lib
	LG = -lgfortran -lmatheval -Wl,--rpath -Wl,/home/monty/local/lib
	CC = g++
endif
ifeq ($(OSTYPE),cygwin)
        CC = gcc
endif
ifeq ($(OSTYPE),darwin)
        DIRS = -I/Users/monty/include -I/opt/local/include -L/Users/monty/lib 
	LG = -lgfortran -lmatheval
	CC = gcc
endif
# CFLAGS = -Wall -g $(DIRS) # debug
CFLAGS = -Wall $(DIRS) # release
LDLIBS = -lgsl -lm -lgslcblas -llapack -lblas $(LG) $(DIRS)
OBJSMADS = mads.o mads_io.o mads_io_external.o mads_func.o mads_mem.o mads_info.o lm/opt_lm_mon.o lm/opt_lm_gsl.o lm/lu.o lm/opt_lm_ch.o misc/test_problems.o misc/anasol_contamination.o misc/io.o lhs/lhs.o 
OBJSPSO = pso/pso-tribes-lm.o pso/Standard_PSO_2006.o pso/mopso.o abagus/abagus.o
OBJSMPUN = mprun/mprun.o mprun/mprun_io.o
OBJSKDTREE = abagus/kdtree-0.5.5/kdtree.o
OBJSLEVMAR = levmar-2.5/lm_m.o levmar-2.5/Axb.o levmar-2.5/misc.o levmar-2.5/lmlec.o levmar-2.5/lmbc.o levmar-2.5/lmblec.o levmar-2.5/lmbleic.o 
OBJSLEVMARSTYLE = levmar-2.5/lm_m.o levmar-2.5/lm_core_m.o levmar-2.5/Axb.o levmar-2.5/misc.o levmar-2.5/lmlec.o levmar-2.5/lmbc.o levmar-2.5/lmblec.o levmar-2.5/lmbleic.o 
SOURCE = $(OBJSMADS:%.o=%.c) $(OBJSPSO:%.o=%.c) $(OBJSMPUN:%.o=%.c) $(OBJSLEVMAR:%.o=%.c) $(OBJSKDTREE:%.o=%.c)
SOURCESTYLE = $(OBJSMADS:%.o=%.c) $(OBJSPSO:%.o=%.c) $(OBJSMPUN:%.o=%.c) $(OBJSLEVMARSTYLE:%.o=%.c) $(OBJSKDTREE:%.o=%.c)

all: $(PROG)

$(PROG): $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSLEVMAR) $(OBJSKDTREE)

mads.o: mads.c mads.h levmar-2.5/levmar.h
mads_io.o: mads_io.c mads.h
mads_io_external.o: mads_io_external.c mads.h
mads_func.o: mads_func.c mads.h
mads_mem.o: mads_mem.c
mads_info.o: mads_info.c
mprun/mprun.o: mprun/mprun.c mads.h
mprun/mprun_io.o: mprun/mprun_io.c mads.h
lm/opt_lm_mon.o: lm/opt_lm_mon.c mads.h
lm/opt_lm_gsl.o: lm/opt_lm_gsl.c mads.h
lm/opt_lm_ch.o: lm/opt_lm_gsl.c mads.h
lhs/lhs.o: lhs/lhs.c
misc/test_problems.o: misc/test_problems.c mads.h
misc/anasol_contamination.o: misc/anasol_contamination.c mads.h
pso/pso-tribes-lm.o: pso/pso-tribes-lm.c pso/pso.h mads.h
pso/Standard_PSO_2006.o: pso/Standard_PSO_2006.c mads.h
pso/mopso.o: pso/mopso.c pso/mopso.h
abagus/abagus.o: abagus/abagus.c mads.h abagus/kdtree-0.5.5/kdtree.o
abagus/kdtree-0.5.5/kdtree.o: abagus/kdtree-0.5.5/kdtree.c abagus/kdtree-0.5.5/kdtree.h
levmar-2.5/lm_m.o: levmar-2.5/lm_m.c levmar-2.5/lm_core_m.c levmar-2.5/levmar.h levmar-2.5/misc.h levmar-2.5/compiler.h mads.h
levmar-2.5/Axb.o: levmar-2.5/Axb.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/misc.o: levmar-2.5/misc.c levmar-2.5/misc_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/lmlec.o: levmar-2.5/lmlec.c levmar-2.5/lmlec_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/lmbc.o: levmar-2.5/lmbc.c levmar-2.5/lmbc_core.c levmar-2.5/levmar.h levmar-2.5/misc.h levmar-2.5/compiler.h
levmar-2.5/lmblec.o: levmar-2.5/lmblec.c levmar-2.5/lmblec_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/lmbleic.o: levmar-2.5/lmbleic.c levmar-2.5/lmbleic_core.c levmar-2.5/levmar.h levmar-2.5/misc.h

clean:
	rm -f $(PROG) $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSLEVMAR) $(OBJSKDTREE)


examples:
	@echo "**************************************************************************************"
	@echo "Example 1: Internal Rosenbrock Problem "
	@echo "**************************************************************************************"
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1
	@echo "**************************************************************************************"
	@echo "Example 1: DONE"
	@echo "**************************************************************************************"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo "Example 2: Internal contamination Problem "
	@echo "**************************************************************************************"
	mads example/contamination/s01 ldebug
	@echo "Example problem example/wells/w01 ..."
	cd example/wells; ../../mads w01 lmeigen
	@echo "**************************************************************************************"
	@echo "Example 2: DONE"
	@echo "**************************************************************************************"

verify: verify-internal verify-multistart1 verify-contaminant verify-multistart2 verify-external verify-parallel
	@echo VERIFICATION DONE

verify-internal:
	@echo "**************************************************************************************"
	@echo " Internal test problems "
	@echo "**************************************************************************************"
	@echo "TEST 1: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 1.1: Levenberg-Marquardt ... "
	rm -f example/rosenbrock/a01.mads_output example/rosenbrock/a01.results
	cd example/rosenbrock; ../../mads a01 test=3 opt=lm lmeigen igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.mads_output example/rosenbrock/a01.mads_output-lm-correct
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-lm-correct
	@echo ""
	@echo "TEST 1.2: Particle-Swarm ..."
	rm -f example/rosenbrock/a01.results
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-pso-correct
	@echo ""
	@echo "TEST 1.3: TRIBES ..."
	rm -f example/rosenbrock/a01.results
	cd example/rosenbrock; ../../mads a01 test=3 opt=tribes igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-tribes-correct
	@echo ""
	@echo "TEST 1.4: SQUADS ..."
	rm -f example/rosenbrock/a01.results
	cd example/rosenbrock; ../../mads a01 test=3 opt=squads igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-squads-correct
	@echo "**************************************************************************************"
	@echo "TEST 1: DONE"
	@echo ""
	@echo ""

verify-multistart1:
	@echo "**************************************************************************************"
	@echo " Multistart (paranoid) problems "
	@echo "**************************************************************************************"
	@echo "TEST 2: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 2.1: Levenberg-Marquardt ... "
	rm -f example/rosenbrock/p01lm.mads_output example/rosenbrock/p01lm.results
	cd example/rosenbrock; ../../mads p01lm test=3 dim=3 opt=lm igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet > /dev/null
	@$(CMP) example/rosenbrock/p01lm.mads_output example/rosenbrock/p01lm.mads_output-multistart-lm-correct
	@$(CMP) example/rosenbrock/p01lm.results example/rosenbrock/p01lm.results-multistart-lm-correct
	@echo ""
	@echo "TEST 2.2: SQUADS ..."
	rm -f example/rosenbrock/p01squads.mads_output example/rosenbrock/p01squads.results
	cd example/rosenbrock; ../../mads p01squads test=3 dim=3 opt=squads igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet > /dev/null
	@$(CMP) example/rosenbrock/p01squads.mads_output example/rosenbrock/p01squads.mads_output-multistart-squads-correct
	@$(CMP) example/rosenbrock/p01squads.results example/rosenbrock/p01squads.results-multistart-squads-correct
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo ""
	@echo ""

verify-contaminant:
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "TEST 3: Internal contaminant transport problems using different dispersivities ..."
	@echo "TEST 3.1: Problem example/contamination/s01 with independent dispersivities ..."
	rm -f example/contamination/s01.mads_output example/contamination/s01.results
	mads example/contamination/s01 lmeigen > /dev/null
	@$(CMP) example/contamination/s01.mads_output example/contamination/s01.mads_output-correct
	@$(CMP) example/contamination/s01.results example/contamination/s01.results-correct
	@echo ""
	@echo "TEST 3.2: Problem example/contamination/s01 with tied dispersivities ..."
	rm -f example/contamination/s01-tied_dispersivities.results
	mads example/contamination/s01-tied_dispersivities > /dev/null
	@$(CMP) example/contamination/s01-tied_dispersivities.results example/contamination/s01-tied_dispersivities.results-correct
	@echo ""
	@echo "TEST 3.3: Problem example/contamination/s01 with scaled dispersivities ..."
	rm -f example/contamination/s01-scaled_dispersivities.results
	mads example/contamination/s01-scaled_dispersivities > /dev/null
	@$(CMP) example/contamination/s01-scaled_dispersivities.results example/contamination/s01-scaled_dispersivities.results-correct
	@echo ""
	@echo "TEST 3.4: Problem example/contamination/s01 with scaled and tied dispersivities ..."
	rm -f example/contamination/s01-scaled+tied_dispersivities.results
	mads example/contamination/s01-scaled+tied_dispersivities > /dev/null
	@$(CMP) example/contamination/s01-scaled+tied_dispersivities.results example/contamination/s01-scaled+tied_dispersivities.results-correct
	@echo ""
	@echo "TEST 3.5: Problem example/contamination/s02 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f example/contamination/s02tied.results
	mads example/contamination/s02tied > /dev/null
	@$(CMP) example/contamination/s02tied.results example/contamination/s02tied.results-correct
	@echo ""
	@echo "TEST 3.6: Problem example/contamination/s02 with regularization terms for optimized model parameters ..."
	rm -f example/contamination/s02regul.results
	mads example/contamination/s02regul > /dev/null
	@$(CMP) example/contamination/s02regul.results example/contamination/s02regul.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3: DONE"
	@echo ""
	@echo ""

verify-multistart2:
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "TEST 4: Internal contaminant transport problems using different optimization techniques ..."
	@echo "TEST 4.1: Problem example/contamination/s01 IGRND ..."
	rm -f example/contamination/s01-igrnd.results example/contamination/s01-igrnd.igrnd.results
	mads example/contamination/s01-igrnd seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-igrnd.results example/contamination/s01-igrnd.results-correct
	@$(CMP) example/contamination/s01-igrnd.igrnd.results example/contamination/s01-igrnd.igrnd.results-correct
	@echo ""
	@echo "TEST 4.2: Problem example/contamination/s01 PPSD ..."
	rm -f example/contamination/s01-ppsd.mads_output example/contamination/s01-ppsd.ppsd.results
	mads example/contamination/s01-ppsd seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-ppsd.mads_output example/contamination/s01-ppsd.mads_output-correct
	@$(CMP) example/contamination/s01-ppsd.ppsd.results example/contamination/s01-ppsd.ppsd.results-correct
	@echo ""
	@echo "TEST 4.3: Problem example/contamination/s01 IGPD ..."
	rm -f example/contamination/s01-igpd.results example/contamination/s01-igpd.igpd.results
	mads example/contamination/s01-igpd seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-igpd.results example/contamination/s01-igpd.results-correct
	@$(CMP) example/contamination/s01-igpd.igpd.results example/contamination/s01-igpd.igpd.results-correct
	@echo ""
	@echo "TEST 4.4: Problem example/contamination/s01 Multi-Start LM  ..."
	rm -f example/contamination/s01-mslm.results
	mads example/contamination/s01-mslm seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-mslm.results example/contamination/s01-mslm.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 4: DONE"
	@echo ""
	@echo ""

verify-external:
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "TEST 5: Problem example/wells/w01 ..."
	rm -f example/wells/w01.mads_output example/wells/w01.results
	cd example/wells; ../../mads w01 lmeigen > /dev/null
	@$(CMP) example/wells/w01.mads_output example/wells/w01.mads_output-correct
	@$(CMP) example/wells/w01.results example/wells/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 5: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "TEST 6: Problem example/wells-short/w01 using different instruction formats ..."
	@echo "TEST 6.1: Instruction file example/wells-short/w01-v1.inst ..."
	rm -f example/wells-short/w01.results
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo ""
	@echo "TEST 6.2: Instruction file example/wells-short/w01-v2.inst ..."
	rm -f example/wells-short/w01.results
	cd example/wells-short; ln -sf w01-v2.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo ""
	@echo "TEST 6.3: Instruction file example/wells-short/w01-v3.inst ..."
	rm -f example/wells-short/w01.results
	cd example/wells-short; ln -sf w01-v3.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo ""
	@echo "TEST 6.4: Instruction file example/wells-short/w01-v4.inst ..."
	rm -f example/wells-short/w01.results 
	cd example/wells-short; ln -sf w01-v4.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 6: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "TEST 7: Problem example/wells-short/w01 using coupled parameters and regularization terms ..."
	@echo "TEST 7.1: Problem example/wells-short/w01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f example/wells-short/w01tied.results
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01tied > /dev/null
	@$(CMP) example/wells-short/w01tied.results example/wells-short/w01tied.results-correct
	@echo ""
	@echo "TEST 7.2: Problem example/wells-short/w01 with regularization terms for optimized model parameters ..."
	rm -f example/wells-short/w01regul.results
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01regul > /dev/null
	@$(CMP) example/wells-short/w01regul.results example/wells-short/w01regul.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 7: DONE"
	@echo ""
	@echo ""

verify-parallel:
	@echo "**************************************************************************************"
	@echo " Parallel execution of external problems "
	@echo "**************************************************************************************"
	@echo "TEST 8: Parallel execution of example/wells-short/w01parallel ..."
	@echo "TEST 8.1: Initial parallel execution of example/wells-short/w01parallel ..."
	rm -f example/wells-short/w01parallel.results example/wells-short/w01parallel.restart_info example/wells-short/w01parallel.restart_*.zip
	cd example/wells-short; ../../mads w01parallel np=2 eval=10 restart=0 > /dev/null
	@$(CMP) example/wells-short/w01parallel.results example/wells-short/w01parallel.results-correct
	@echo ""
	@echo "TEST 8.2: Rerun using saved results from prior parallel execution of example/wells-short/w01parallel ..."
	rm -f example/wells/w01parallel.results
	cd example/wells-short; ../../mads w01parallel np=2 eval=10 > /dev/null
	@$(CMP) example/wells-short/w01parallel.results example/wells-short/w01parallel.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 8: DONE"

clean-example:
	rm -f example/*/*.mads_output_* example/*/*.ppsd_*.results example/*/*.igpd_*.results example/*/*.igrnd_*.results example/*/*.restart_*.zip example/*/*.restart_info example/*/*.running
	rm -fR example/wells-short_w01parallel*

astyle:
	astyle $(SOURCESTYLE)
	rm -f $(SOURCESTYLE:%c=%c.orig)

tar:
	rm -f mads.tgz
	tar -cvzf mads.tgz `hg st -c -m | awk '{print $$2}'` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'` .hg

tarf:
	tar -cvzf mads.tgz `hg st -c -m | awk '{print $$2}'` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'`
