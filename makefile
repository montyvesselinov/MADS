PROG = mads
CC = gcc 
# CC = g++ 
## rainier -L/usr/lib/gcc/x86_64-redhat-linux/4.1.1
ifeq ($(OSTYPE),linux)
        DIRS = -I/home/monty/local/include -L/home/monty/local/lib
	LG = -lgfortran
else
        # DIRS = -I/opt/local/include/ -L/opt/local/lib
        DIRS = -I/Users/monty/include -L/Users/monty/lib -I/opt/local/include/ -L/opt/local/lib
endif
CFLAGS = -g -Wall $(DIRS)
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
levmar-2.5/Axb.o: levmar-2.5/Axb.c levmar-2.5/Axb_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/misc.o: levmar-2.5/misc.c levmar-2.5/misc_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/lmlec.o: levmar-2.5/lmlec.c levmar-2.5/lmlec_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/lmbc.o: levmar-2.5/lmbc.c levmar-2.5/lmbc_core.c levmar-2.5/levmar.h levmar-2.5/misc.h levmar-2.5/compiler.h
levmar-2.5/lmblec.o: levmar-2.5/lmblec.c levmar-2.5/lmblec_core.c levmar-2.5/levmar.h levmar-2.5/misc.h
levmar-2.5/lmbleic.o: levmar-2.5/lmbleic.c levmar-2.5/lmbleic_core.c levmar-2.5/levmar.h levmar-2.5/misc.h

clean:
	rm -f $(PROG) $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSLEVMAR) $(OBJSKDTREE)


examples:
	@echo "**************************************************************************************"
	@echo "**************************************************************************************"
	@echo "Example 1"
	@echo "Example problem rosenbrock ... "
	@echo "**************************************************************************************"
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1
	@echo "**************************************************************************************"
	@echo "Example 1: DONE"
	@echo "**************************************************************************************"
	@echo ""
	@echo "**************************************************************************************"
	@echo "**************************************************************************************"
	@echo "Example 2"
	mads example/contamination/s01 ldebug
	@echo "Example problem example/wells/w01 ..."
	cd example/wells; ../../mads w01 ldebug
	@echo "**************************************************************************************"
	@echo "Example 2: DONE"
	@echo "**************************************************************************************"

verify:
	#seed=2096575428
	#seed=1977879092
	@echo "**************************************************************************************"
	@echo "TEST 1: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 1.1: Levenberg-Marquardt ... "
	cd example/rosenbrock; ../../mads a01 test=3 opt=lm igrnd real=1 seed=2096575428 > /dev/null
	@./compare-results example/rosenbrock/a01.results example/rosenbrock/a01.results-lm-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.2: Particle-Swarm ..."
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1 seed=2096575428 > /dev/null
	@./compare-results example/rosenbrock/a01.results example/rosenbrock/a01.results-pso-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.3: TRIBES ..."
	cd example/rosenbrock; ../../mads a01 test=3 opt=tribes igrnd real=1 seed=2096575428 > /dev/null
	@./compare-results example/rosenbrock/a01.results example/rosenbrock/a01.results-tribes-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.4: SQUADS ..."
	cd example/rosenbrock; ../../mads a01 test=3 opt=squads igrnd real=1 seed=2096575428 > /dev/null
	@./compare-results example/rosenbrock/a01.results example/rosenbrock/a01.results-squads-correct
	@echo "**************************************************************************************"
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 1.0b: Rosenbrock problem using different multistart (paranoid) optimization techniques ..."
	@echo "TEST 1.1b: Levenberg-Marquardt ... "
	cd example/rosenbrock; ../../mads p01lm test=3 dim=3 opt=lm igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 > p01lm.outp-multistart-lm
	@./compare-results example/rosenbrock/p01lm.outp-multistart-lm example/rosenbrock/p01lm.outp-lm-correct
	@./compare-results example/rosenbrock/p01lm.results example/rosenbrock/p01lm.results-multistart-lm-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.2b: SQUADS ..."
	cd example/rosenbrock; ../../mads p01squads test=3 dim=3 opt=squads igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 > p01squads.outp-multistart-squads
	@./compare-results example/rosenbrock/p01squads.outp-multistart-squads example/rosenbrock/p01squads.outp-squads-correct
	@./compare-results example/rosenbrock/p01squads.results example/rosenbrock/p01squads.results-multistart-squads-correct
	@echo "**************************************************************************************"
	@echo "TEST 1: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 2: Problem example/contamination/s01 ..."
	mads example/contamination/s01 sindx=0.01 > /dev/null
	@./compare-results example/contamination/s01.results example/contamination/s01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 2: Problem example/wells/w01 ..."
	cd example/wells; ../../mads w01 > /dev/null
	@./compare-results example/wells/w01.results example/wells/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 3: Problem example/wells-short/w01 using different instruction formats ..."
	@echo "TEST 3.1: Instruction file example/wells-short/w01-v1.inst ..."
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01 > /dev/null
	@./compare-results example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3.2: Instruction file example/wells-short/w01-v2.inst ..."
	cd example/wells-short; ln -sf w01-v2.inst w01.inst; ../../mads w01 > /dev/null
	@./compare-results example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3.3: Instruction file example/wells-short/w01-v3.inst ..."
	cd example/wells-short; ln -sf w01-v3.inst w01.inst; ../../mads w01 > /dev/null
	@./compare-results example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3.4: Instruction file example/wells-short/w01-v4.inst ..."
	cd example/wells-short; ln -sf w01-v4.inst w01.inst; ../../mads w01 > /dev/null
	@./compare-results example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3: DONE"

create-verify:
	#seed=2096575428
	#seed=1977879092
	@echo "**************************************************************************************"
	@echo "TEST 1.0a: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 1.1a: Levenberg-Marquardt ... "
	cd example/rosenbrock; ../../mads a01 test=3 opt=lm igrnd real=1 seed=2096575428 > /dev/null
	@cp example/rosenbrock/a01.results example/rosenbrock/a01.results-lm-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.2a: Particle-Swarm ..."
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1 seed=2096575428 > /dev/null
	@cp example/rosenbrock/a01.results example/rosenbrock/a01.results-pso-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.3a: TRIBES ..."
	cd example/rosenbrock; ../../mads a01 test=3 opt=tribes igrnd real=1 seed=2096575428 > /dev/null
	@cp example/rosenbrock/a01.results example/rosenbrock/a01.results-tribes-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.4a: SQUADS ..."
	cd example/rosenbrock; ../../mads a01 test=3 opt=squads igrnd real=1 seed=2096575428 > /dev/null
	@cp example/rosenbrock/a01.results example/rosenbrock/a01.results-squads-correct
	@echo "**************************************************************************************"
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 1.0b: Rosenbrock problem using different multistart (paranoid) optimization techniques ..."
	@echo "TEST 1.1b: Levenberg-Marquardt ... "
	cd example/rosenbrock; ../../mads p01lm test=3 dim=3 opt=lm igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 > p01lm.outp-multistart-lm
	@cp example/rosenbrock/p01lm.outp-multistart-lm example/rosenbrock/p01lm.outp-lm-correct
	@cp example/rosenbrock/p01lm.results example/rosenbrock/p01lm.results-multistart-lm-correct
	@echo "**************************************************************************************"
	@echo "TEST 1.2b: SQUADS ..."
	cd example/rosenbrock; ../../mads p01squads test=3 dim=3 opt=squads igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 > p01squads.outp-multistart-squads
	@cp example/rosenbrock/p01squads.outp-multistart-squads example/rosenbrock/p01squads.outp-squads-correct
	@cp example/rosenbrock/p01squads.results example/rosenbrock/p01squads.results-multistart-squads-correct
	@echo "**************************************************************************************"
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 1: DONE"
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 2: Problem example/contamination/s01 ..."
	mads example/contamination/s01 sindx=0.01 > /dev/null
	@cp example/contamination/s01.results example/contamination/s01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 2: Problem example/wells/w01 ..."
	cd example/wells; ../../mads w01 > /dev/null
	@cp example/wells/w01.results example/wells/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 3: Problem example/wells-short/w01 using different instruction formats ..."
	@echo "TEST 3.1: Instruction file example/wells-short/w01-v1.inst ..."
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01 > /dev/null
	@cp example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3.2: Instruction file example/wells-short/w01-v2.inst ..."
	cd example/wells-short; ln -sf w01-v2.inst w01.inst; ../../mads w01 > /dev/null
	@./compare-results example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3.3: Instruction file example/wells-short/w01-v3.inst ..."
	cd example/wells-short; ln -sf w01-v3.inst w01.inst; ../../mads w01 > /dev/null
	@cp example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3.4: Instruction file example/wells-short/w01-v4.inst ..."
	cd example/wells-short; ln -sf w01-v4.inst w01.inst; ../../mads w01 > /dev/null
	@cp example/wells-short/w01.results example/wells-short/w01.results-correct
	@echo "**************************************************************************************"
	@echo "TEST 3: DONE"

astyle:
	astyle $(SOURCESTYLE)
	rm -f $(SOURCESTYLE:%c=%c.orig)

tar:
	tar -cvzf mads.tgz `hg st -c -m | awk '{print $$2}'` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'` .hg

tarf:
	tar -cvzf mads.tgz `hg st -c -m | awk '{print $$2}'` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'`
