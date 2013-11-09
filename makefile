PROG = mads
CMP = ./compare-results # MADS testing
# CMP = cp -f # Save current results for future testing DANGEROUS!

# MathEval required to evaluate expression for tied parameters and regularization terms
ifndef MATHEVAL
MATHEVAL = true
endif
# Support for YAML input files is optional 
ifndef YAML
YAML = false
endif
OS = $(shell uname -s)
ND = $(shell uname -n)
# Compilation setup
$(info OS type -- $(OS))
$(info Machine -- $(ND))
CC = gcc
CFLAGS = -Wall -O1 -Winit-self
LDLIBS = -lgsl -llapack -lstdc++
ifeq ($(OS),Linux)
# Linux
LDLIBS += -lgfortran 
$(info LINUX)
ifeq ($(ND),aquifer.lanl.gov)
$(info Machine -- AQUIFER)
YAML = true
CFLAGS += -I/home/monty/local/include
LDLIBS += -L/home/monty/local/lib -lgslcblas -lgfortran -Wl,--rpath -Wl,/home/monty/local/lib 
endif
ifeq ($(ND),madsmax)
$(info Machine -- MadsMax)
LDLIBS += -lgslcblas -lm -llapack -lblas
YAML = true
endif
ifeq ($(ND),well.lanl.gov)
$(info Machine -- WELL)
CFLAGS += -I/home/monty/local/include
LDLIBS += -L/home/monty/local/lib -L/usr/local/lib -lgslcblas -lm -llapack -lblas
YAML = true
endif
else
ifeq ($(OS),Cygwin)
# Cygwin
$(info CYGWIN)
else
ifeq ($(OS),Darwin)
# Mac
$(info MAC OS X)
CFLAGS += -I/opt/local/include
LDLIBS += -lgfortran -latlas -L/opt/local/lib
YAML = true
ifeq ($(ND),pn1246281)
LDLIBS += -lblas
endif
ifeq ($(ND),dazed.local)
$(info Machine -- Dazed)
CFLAGS += -I/Users/monty/include
LDLIBS += -lrefblas -lcblas -L/Users/monty/lib
endif
else
$(error UNKNOWN OS type -- $(OS)!)
endif
endif
endif
# MADS files
OBJSMADS = mads.o mads_io.o mads_io_external.o mads_func.o mads_mem.o mads_info.o lm/opt_lm_mon.o lm/opt_lm_gsl.o lm/lu.o lm/opt_lm_ch.o misc/test_problems.o misc/anasol_contamination.o misc/io.o lhs/lhs.o bayes/dream.o
OBJSPSO = pso/pso-tribes-lm.o pso/Standard_PSO_2006.o pso/mopso.o
OBJSA = sa/abagus.o sa/postpua.o sa/global.o sa/do_miser.o
OBJDS = ds/infogap.o ds/glue.o
OBJSMPUN = mprun/mprun.o mprun/mprun_io.o
OBJSKDTREE = misc/kdtree-0.5.5/kdtree.o
OBJSLEVMAR = misc/levmar-2.5/lm_m.o misc/levmar-2.5/Axb.o misc/levmar-2.5/misc.o misc/levmar-2.5/lmlec.o misc/levmar-2.5/lmbc.o misc/levmar-2.5/lmblec.o misc/levmar-2.5/lmbleic.o 
OBJSlEVMARSTYLE = misc/levmar-2.5/lm_m.o misc/levmar-2.5/lm_core_m.o misc/levmar-2.5/Axb.o misc/levmar-2.5/misc.o misc/levmar-2.5/lmlec.o misc/levmar-2.5/lmbc.o misc/levmar-2.5/lmblec.o misc/levmar-2.5/lmbleic.o 
OBJSASTABLE = misc/astable/astable.o misc/astable/interpolation.o misc/astable/pqueue.o
OBJSBAYES = bayes/dream.o

ifeq ($(YAML),true)
$(info YAML Support included)
OBJSMADS += mads_io_yaml.o
CFLAGS += -DYAML `pkg-config --cflags glib-2.0`
LDLIBS += -lyaml `pkg-config --libs glib-2.0`
endif

ifeq ($(MATHEVAL),true)
$(info MathEval Support included)
CFLAGS += -DMATHEVAL
LDLIBS += -lmatheval
endif

SOURCE = $(OBJSMADS:%.o=%.c) $(OBJSPSO:%.o=%.c) $(OBJSMPUN:%.o=%.c) $(OBJSA:%.o=%.c) $(OBJDS:%.o=%.c) $(OBJSLEVMAR:%.o=%.c) $(OBJSKDTREE:%.o=%.c) $(OBJSASTABLE:%.o=%.c) $(OBJSBAYES:%.o=%.cpp)
SOURCESTYLE = $(OBJSMADS:%.o=%.c) $(OBJSPSO:%.o=%.c) $(OBJSMPUN:%.o=%.c) $(OBJSA:%.o=%.c) $(OBJDS:%.o=%.c) $(OBJSLEVMARSTYLE:%.o=%.c) $(OBJSKDTREE:%.o=%.c) $(OBJSASTABLE:%.o=%.c) $(OBJSBAYES:%.o=%.cpp)

all: $(PROG)

release: $(PROG)

debug: CFLAGS += -g
debug: $(PROG)

$(PROG): $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSA) $(OBJDS) $(OBJSLEVMAR) $(OBJSKDTREE) $(OBJSASTABLE) $(OBJSBAYES)

clean:
	rm -f $(PROG) $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSA) $(OBJDS) $(OBJSLEVMAR) $(OBJSKDTREE) $(OBJSASTABLE) $(OBJSBAYES)


mads.o: mads.c mads.h misc/levmar-2.5/levmar.h
mads_io.o: mads_io.c mads.h
mads_io_yaml.o: mads_io_yaml.c mads.h
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
pso/pso-tribes-lm.o: pso/pso-tribes-lm.c pso/pso.h mads.h
pso/Standard_PSO_2006.o: pso/Standard_PSO_2006.c mads.h
pso/mopso.o: pso/mopso.c pso/mopso.h
sa/abagus.o: sa/abagus.c mads.h misc/kdtree-0.5.5/kdtree.o
sa/postpua.o: sa/postpua.c mads.h
sa/global.o: sa/global.c mads.h sa/do_miser.o
sa/do_miser.o: sa/do_miser.c sa/do_miser.h
ds/infogap.o: ds/infogap.c mads.h
ds/glue.o: ds/glue.c mads.h
misc/anasol_contamination.o: misc/anasol_contamination.c mads.h
misc/test_problems.o: misc/test_problems.c mads.h
misc/kdtree-0.5.5/kdtree.o: misc/kdtree-0.5.5/kdtree.c misc/kdtree-0.5.5/kdtree.h
misc/levmar-2.5/lm_m.o: misc/levmar-2.5/lm_m.c misc/levmar-2.5/lm_core_m.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h misc/levmar-2.5/compiler.h mads.h
misc/levmar-2.5/Axb.o: misc/levmar-2.5/Axb.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h
misc/levmar-2.5/misc.o: misc/levmar-2.5/misc.c misc/levmar-2.5/misc_core.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h
misc/levmar-2.5/lmlec.o: misc/levmar-2.5/lmlec.c misc/levmar-2.5/lmlec_core.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h
misc/levmar-2.5/lmbc.o: misc/levmar-2.5/lmbc.c misc/levmar-2.5/lmbc_core.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h misc/levmar-2.5/compiler.h
misc/levmar-2.5/lmblec.o: misc/levmar-2.5/lmblec.c misc/levmar-2.5/lmblec_core.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h
misc/levmar-2.5/lmbleic.o: misc/levmar-2.5/lmbleic.c misc/levmar-2.5/lmbleic_core.c misc/levmar-2.5/levmar.h misc/levmar-2.5/misc.h
misc/astable/astable.o: misc/astable/astable.c misc/astable/astable.h
misc/astable/interpolation.o: misc/astable/astable.c misc/astable/astable.h misc/astable/interpolation.c misc/astable/pqueue.c misc/astable/pqueue.h
misc/astable/pqueue.o: misc/astable/pqueue.c misc/astable/pqueue.h
bayes/dream.o: bayes/dream.cpp bayes/dream.h mads.h misc/astable/interpolation.o
	g++ $(CFLAGS) -c -o bayes/dream.o bayes/dream.cpp -lmisc/astable/interpolation.o

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

verify: verify-internal verify-multistart1 verify-contaminant verify-multistart2 verify-external verify-parallel verify-forward verify-sa
	@echo VERIFICATION DONE


verify-internal:
	@echo "**************************************************************************************"
	@echo " Internal test problems "
	@echo "**************************************************************************************"
	@echo "TEST 1: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 1.1: Levenberg-Marquardt ... "
	rm -f example/rosenbrock/a01.mads_output example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; ../../mads a01 test=3 opt=lm lmeigen igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.mads_output example/rosenbrock/a01.mads_output-lm-$(OS)-correct
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-lm-$(OS)-correct
	@echo ""
	@echo "TEST 1.2: Particle-Swarm ..."
	rm -f example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-pso-$(OS)-correct
	@echo ""
	@echo "TEST 1.3: TRIBES ..."
	rm -f example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; ../../mads a01 test=3 opt=tribes igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-tribes-$(OS)-correct
	@echo ""
	@echo "TEST 1.4: SQUADS ..."
	rm -f example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; ../../mads a01 test=3 opt=squads igrnd real=1 seed=2096575428 > /dev/null
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-squads-$(OS)-correct
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
	rm -f example/rosenbrock/p01lm.mads_output example/rosenbrock/p01lm.results example/rosenbrock/p01lm.running
	cd example/rosenbrock; ../../mads p01lm test=3 dim=3 opt=lm igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet > /dev/null
	@$(CMP) example/rosenbrock/p01lm.mads_output example/rosenbrock/p01lm.mads_output-multistart-lm-$(OS)-correct
	@$(CMP) example/rosenbrock/p01lm.results example/rosenbrock/p01lm.results-multistart-lm-$(OS)-correct
	@echo ""
	@echo "TEST 2.2: SQUADS ..."
	rm -f example/rosenbrock/p01squads.mads_output example/rosenbrock/p01squads.results example/rosenbrock/p01squads.running
	cd example/rosenbrock; ../../mads p01squads test=3 dim=3 opt=squads igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet > /dev/null
	@$(CMP) example/rosenbrock/p01squads.mads_output example/rosenbrock/p01squads.mads_output-multistart-squads-$(OS)-correct
	@$(CMP) example/rosenbrock/p01squads.results example/rosenbrock/p01squads.results-multistart-squads-$(OS)-correct
	@echo ""
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo ""
	@echo ""

verify-contaminant:
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "TEST 3: Internal contaminant transport problems using different dispersivities ..."
	@echo "TEST 3.1.a: Problem example/contamination/s01 with independent dispersivities (MADS text input format) ..."
	rm -f example/contamination/s01.mads_output example/contamination/s01.results example/contamination/s01.running
	./mads example/contamination/s01 obs_int=2 lmeigen > /dev/null
	@$(CMP) example/contamination/s01.mads_output example/contamination/s01.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01.results example/contamination/s01.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.1.b: Problem example/contamination/s01 with independent dispersivities (YAML input format) ..."
	rm -f example/contamination/s01_yaml.mads_output example/contamination/s01_yaml.results example/contamination/s01_yaml.running
	./mads example/contamination/s01_yaml obs_int=2 lmeigen > /dev/null
	@$(CMP) example/contamination/s01_yaml.mads_output example/contamination/s01_yaml.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01_yaml.results example/contamination/s01_yaml.results-$(OS)-correct
	@$(CMP) example/contamination/s01_yaml.results example/contamination/s01.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.2: Problem example/contamination/s01 with tied dispersivities ..."
	rm -f example/contamination/s01-tied_dispersivities.results example/contamination/s01-tied_dispersivities.running
	./mads example/contamination/s01-tied_dispersivities obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-tied_dispersivities.results example/contamination/s01-tied_dispersivities.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.3: Problem example/contamination/s01 with scaled dispersivities ..."
	rm -f example/contamination/s01-scaled_dispersivities.results example/contamination/s01-scaled_dispersivities.running
	./mads example/contamination/s01-scaled_dispersivities obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-scaled_dispersivities.results example/contamination/s01-scaled_dispersivities.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.4: Problem example/contamination/s01 with scaled and tied dispersivities ..."
	rm -f example/contamination/s01-scaled+tied_dispersivities.results example/contamination/s01-scaled+tied_dispersivities.running
	./mads example/contamination/s01-scaled+tied_dispersivities obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-scaled+tied_dispersivities.results example/contamination/s01-scaled+tied_dispersivities.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.5: Problem example/contamination/s01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f example/contamination/s01-tied.results example/contamination/s01-tied.running
	./mads example/contamination/s01-tied obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-tied.results example/contamination/s01-tied.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.6: Problem example/contamination/s01_yaml with coupled (tied) parameters based on mathematical expressions (YAML input format) ..."
	rm -f example/contamination/s01-tied_yaml.results example/contamination/s01-tied_yaml.running
	./mads example/contamination/s01-tied_yaml obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-tied_yaml.mads_output example/contamination/s01-tied_yaml.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-tied_yaml.results example/contamination/s01-tied.results-$(OS)-correct
	@$(CMP) example/contamination/s01-tied_yaml.results example/contamination/s01-tied_yaml.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.7: Problem example/contamination/s01 with regularization terms for optimized model parameters ..."
	rm -f example/contamination/s01-regul.results example/contamination/s01-regul.running
	./mads example/contamination/s01-regul obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-regul.mads_output example/contamination/s01-regul.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-regul.results example/contamination/s01-regul.results-$(OS)-correct
	@echo ""
	@echo "TEST 3.8: Problem example/contamination/s01_yaml with regularization terms for optimized model parameters (YAML input format) ..."
	rm -f example/contamination/s01-regul_yaml.results example/contamination/s01-regul_yaml.running
	./mads example/contamination/s01-regul_yaml obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-regul_yaml.mads_output example/contamination/s01-regul_yaml.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-regul_yaml.results example/contamination/s01-regul.results-$(OS)-correct
	@$(CMP) example/contamination/s01-regul_yaml.results example/contamination/s01-regul_yaml.results-$(OS)-correct
	@echo "**************************************************************************************"
	@echo ""
	@echo "TEST 3.9: Problem example/contamination/s01-multi_source with multiple source and tied model parameters (YAML input format) ..."
	rm -f example/contamination/s01-multi_source.results example/contamination/s01-multi_source.running
	./mads example/contamination/s01-multi_source obs_int=2 > /dev/null
	@$(CMP) example/contamination/s01-multi_source.mads_output example/contamination/s01-multi_source.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-multi_source.results example/contamination/s01-multi_source.results-$(OS)-correct
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
	rm -f example/contamination/s01-igrnd.results example/contamination/s01-igrnd.igrnd.results example/contamination/s01-igrnd.running
	./mads example/contamination/s01-igrnd seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-igrnd.results example/contamination/s01-igrnd.results-$(OS)-correct
	@$(CMP) example/contamination/s01-igrnd.igrnd.results example/contamination/s01-igrnd.igrnd.results-$(OS)-correct
	@echo ""
	@echo "TEST 4.2: Problem example/contamination/s01 PPSD ..."
	rm -f example/contamination/s01-ppsd.mads_output example/contamination/s01-ppsd.ppsd.results example/contamination/s01-ppsd.running
	./mads example/contamination/s01-ppsd seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-ppsd.mads_output example/contamination/s01-ppsd.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-ppsd.ppsd.results example/contamination/s01-ppsd.ppsd.results-$(OS)-correct
	@echo ""
	@echo "TEST 4.3: Problem example/contamination/s01 IGPD ..."
	rm -f example/contamination/s01-igpd.results example/contamination/s01-igpd.igpd.results example/contamination/s01-igpd.running
	./mads example/contamination/s01-igpd seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-igpd.results example/contamination/s01-igpd.results-$(OS)-correct
	@$(CMP) example/contamination/s01-igpd.igpd.results example/contamination/s01-igpd.igpd.results-$(OS)-correct
	@echo ""
	@echo "TEST 4.4: Problem example/contamination/s01 Multi-Start LM  ..."
	rm -f example/contamination/s01-mslm.results example/contamination/s01-mslm.running
	./mads example/contamination/s01-mslm seed=2096575428 > /dev/null
	@$(CMP) example/contamination/s01-mslm.results example/contamination/s01-mslm.results-$(OS)-correct
	@echo "**************************************************************************************"
	@echo "TEST 4: DONE"
	@echo ""
	@echo ""

verify-external:
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "TEST 5: Problem example/wells/w01 ..."
	rm -f example/wells/w01.mads_output example/wells/w01.results example/wells/w01.running
	cd example/wells; ../../mads w01 lmeigen > /dev/null
	@$(CMP) example/wells/w01.mads_output example/wells/w01.mads_output-$(OS)-correct
	@$(CMP) example/wells/w01.results example/wells/w01.results-$(OS)-correct
	@echo "**************************************************************************************"
	@echo "TEST 5: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "TEST 6: Problem example/wells-short/w01 using different instruction formats ..."
	@echo "TEST 6.1: Instruction file example/wells-short/w01-v1.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo ""
	@echo "TEST 6.2: Instruction file example/wells-short/w01-v2.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v2.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo ""
	@echo "TEST 6.3: Instruction file example/wells-short/w01-v3.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v3.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo ""
	@echo "TEST 6.4: Instruction file example/wells-short/w01-v4.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v4.inst w01.inst; ../../mads w01 > /dev/null
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo "**************************************************************************************"
	@echo "TEST 6: DONE"
	@echo ""
	@echo ""
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "TEST 7: Problem example/wells-short/w01 using coupled parameters and regularization terms ..."
	@echo "TEST 7.1: Problem example/wells-short/w01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f example/wells-short/w01tied.results example/wells-short/w01tied.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01tied > /dev/null
	@$(CMP) example/wells-short/w01tied.results example/wells-short/w01tied.results-$(OS)-correct
	@echo ""
	@echo "TEST 7.2: Problem example/wells-short/w01 with regularization terms for optimized model parameters (MADS text input format) ..."
	rm -f example/wells-short/w01regul.results example/wells-short/w01regul.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01regul > /dev/null
	@$(CMP) example/wells-short/w01regul.results example/wells-short/w01regul.results-$(OS)-correct
	@echo ""
	@echo "TEST 7.3: Problem example/wells-short/w01 with regularization terms for optimized model parameters (YAML input format) ..."
	rm -f example/wells-short/w01regul.results example/wells-short/w01regul.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; ../../mads w01regul_yaml > /dev/null
	@$(CMP) example/wells-short/w01regul_yaml.results example/wells-short/w01regul_yaml.results-$(OS)-correct
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
	rm -f example/wells-short/w01parallel.results example/wells-short/w01parallel.restart_info example/wells-short/w01parallel.restart_*.zip example/wells-short/w01parallel.running
	cd example/wells-short; ../../mads w01parallel np=2 eval=10 restart=0 > /dev/null
	@$(CMP) example/wells-short/w01parallel.results example/wells-short/w01parallel.results-$(OS)-correct
	@echo ""
	@echo "TEST 8.2: Rerun using saved results from prior parallel execution of example/wells-short/w01parallel ..."
	rm -f example/wells/w01parallel.results example/wells/w01parallel.running
	cd example/wells-short; ../../mads w01parallel np=2 eval=10 > /dev/null
	@$(CMP) example/wells-short/w01parallel.results example/wells-short/w01parallel.results-$(OS)-correct
	@echo "**************************************************************************************"
	@echo "TEST 8: DONE"
	@echo ""
	@echo ""

verify-forward:
	@echo "**************************************************************************************"
	@echo " Forward analytical contaminant modeling "
	@echo "**************************************************************************************"
	@echo "TEST 9: Analytical contaminant concentrations with various sources ..."
	@echo "TEST 9.1: Point source ... "
	rm -f example/forward/a01.mads_output example/forward/a01.results example/forward/a01.running
	cd example/forward; ../../mads a01 > /dev/null
	@$(CMP) example/forward/a01.mads_output example/forward/a01.mads_output-$(OS)-correct
	@$(CMP) example/forward/a01.results example/forward/a01.results-$(OS)-correct
	@echo ""
	@echo "TEST 9.2: Rectangular source ... "
	rm -f example/forward/a02.mads_output example/forward/a02.results example/forward/a02.running
	cd example/forward; ../../mads a02 > /dev/null
	@$(CMP) example/forward/a02.mads_output example/forward/a02.mads_output-$(OS)-correct
	@$(CMP) example/forward/a02.results example/forward/a02.results-$(OS)-correct
	@echo ""
	@echo "TEST 9.3: Planar gaussian source ... "
	rm -f example/forward/a03.mads_output example/forward/a03.results example/forward/a03.running
	cd example/forward; ../../mads a03 > /dev/null
	@$(CMP) example/forward/a03.mads_output example/forward/a03.mads_output-$(OS)-correct
	@$(CMP) example/forward/a03.results example/forward/a03.results-$(OS)-correct
	@echo ""
	@echo "TEST 9.4: Gaussian source ... "
	rm -f example/forward/a04.mads_output example/forward/a04.results example/forward/a04.running
	cd example/forward; ../../mads a04 > /dev/null
	@$(CMP) example/forward/a04.mads_output example/forward/a04.mads_output-$(OS)-correct
	@$(CMP) example/forward/a04.results example/forward/a04.results-$(OS)-correct
	@echo "" 
	@echo "TEST 9.5: Box source ... "
	rm -f example/forward/a05.mads_output example/forward/a05.results example/forward/a05.running
	cd example/forward; ../../mads a05 > /dev/null
	@$(CMP) example/forward/a05.mads_output example/forward/a05.mads_output-$(OS)-correct
	@$(CMP) example/forward/a05.results example/forward/a05.results-$(OS)-correct
	@echo "" 
	@echo "TEST 9.6: Box source | Levy ... "
	rm -f example/forward/a06_yaml.mads_output example/forward/a06_yaml.results example/forward/a06_yaml.running
	cd example/forward; ../../mads a06_yaml > /dev/null
	@$(CMP) example/forward/a06_yaml.mads_output example/forward/a06_yaml.mads_output-$(OS)-correct
	@$(CMP) example/forward/a06_yaml.results example/forward/a06_yaml.results-$(OS)-correct
	@echo "" 
	@echo "TEST 9.7: Box source | Symmetric Levy ... "
	rm -f example/forward/a07.mads_output example/forward/a07.results example/forward/a07.running
	cd example/forward; ../../mads a07 > /dev/null
	@$(CMP) example/forward/a07.mads_output example/forward/a07.mads_output-$(OS)-correct
	@$(CMP) example/forward/a07.results example/forward/a07.results-$(OS)-correct
	@echo "**************************************************************************************"
	@echo "TEST 9: DONE"
	@echo ""
	@echo ""

verify-sa:
	@echo "**************************************************************************************"
	@echo " Sensitivity analyses "
	@echo "**************************************************************************************"
	@echo "TEST 10: Global and local sensitivity analyses ..."
	@echo "TEST 10.1: Sobol analysis ... "
	rm -f example/sa/a01.mads_output example/sa/a01.sobol_sens_index example/sa/a01.sobol_sens_total example/sa/a01.running
	cd example/sa; ../../mads a01 test=111 gsens dim=8 real=10000 smp=lhs pardomain=0.5 seed=1517604820 > /dev/null
	@$(CMP) example/sa/a01.sobol_sens_index example/sa/a01.sobol_sens_index-$(OS)-correct
	@$(CMP) example/sa/a01.sobol_sens_total example/sa/a01.sobol_sens_total-$(OS)-correct
	@echo ""
	@echo "TEST 10: DONE"
	@echo ""

compare-os:
	./compare-results-os Linux Darwin

clean-example:
	rm -f example/*/*.mads_output_* example/*/*.ppsd_*.results example/*/*.igpd_*.results example/*/*.igrnd_*.results example/*/*.restart_*.zip example/*/*.restart_info example/*/*.running example/*/*-rerun.mads example/*/*-error.mads
	rm -fR example/wells-short_w01parallel*
	rm -f *.mads_output* *.running *.cmdline *.cmdline_hist

astyle:
	astyle $(SOURCESTYLE)
	rm -f $(SOURCESTYLE:%c=%c.orig)

tar:
	rm -f mads.git.tgz
	tar -cvzf mads.git.tgz `git ls-files` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'` .git

tarf:
	tar -cvzf mads.tgz `git ls-files` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'`
