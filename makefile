# MADS: Model Analyses & Decision Support (v.1.1.14) 2013
#
# Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
# Dan O'Malley, omalled@lanl.gov
#
# http://mads.lanl.gov
# http://www.ees.lanl.gov/staff/monty/codes/mads
# http://gitlab.com/monty/mads
#
# Licensing: GPLv3: http:#www.gnu.org/licenses/gpl-3.0.html
#
# LA-CC-10-055; LA-CC-11-035
#
# Copyright 2011.  Los Alamos National Security, LLC.  All rights reserved.
# This material was produced under U.S. Government contract DE-AC52-06NA25396 for
# Los Alamos National Laboratory, which is operated by Los Alamos National Security, LLC for
# the U.S. Department of Energy. The Government is granted for itself and others acting on its
# behalf a paid-up, nonexclusive, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, and perform publicly and display publicly. Beginning five (5) years after
# --------------- March 11, 2011, -------------------------------------------------------------------
# subject to additional five-year worldwide renewals, the Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
# material to reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY, LLC,
# NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR
# PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

MADS = mads
WELLS = wells/wells
CMP = ./compare-results # MADS testing
OUTPUT = > /dev/null
# OUTPUT = mads-debug-output
# OUTPUT =
DBG = valgrind -v --read-var-info=yes --tool=memcheck --leak-check=yes --leak-check=full --show-reachable=yes --track-origins=yes
DBG = gdb --args
DBG =
# CMP = cp -f # Save current results for future testing DANGEROUS!

# MathEval required to evaluate expression for tied parameters and regularization terms
ifndef MATHEVAL
MATHEVAL = true
endif
# Support for YAML input files 
ifndef YAML
YAML = true
endif
# Support for XML input files
ifdef XML
XML = true
endif
OS = $(shell uname -s)
ND = $(shell uname -n)
# Compilation setup
$(info MADS computationlal framework)
$(info -----------------------------)
$(info OS type -- $(OS))
$(info Machine -- $(ND))
CC = gcc
CFLAGS = -Wall -O1 -Winit-self
LDLIBS = -lgsl -llapack -lstdc++
ifeq ($(OS),Linux)
# Linux
LDLIBS += -lgfortran -lgslcblas -lm -lblas
$(info LINUX)
ifeq ($(ND),aquifer.lanl.gov)
$(info Machine -- AQUIFER)
CFLAGS += -I/home/monty/local/include-aquifer
LDLIBS = -lgsl -llapack -lstdc++ -L/home/monty/local/lib -lgslcblas -lgfortran -Wl,--rpath -Wl,/home/monty/local/lib 
endif
ifeq ($(ND),madsmax)
$(info Machine -- MadsMax)
CFLAGS += -Wno-unused-result
endif
ifeq ($(ND),well.lanl.gov)
$(info Machine -- WELL)
CFLAGS += -I/home/monty/local/include
LDLIBS += -L/home/monty/local/lib -L/usr/local/lib
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
LDLIBS += -lgfortran -lblas -L/opt/local/lib
ifeq ($(ND),bored.lanl.gov)
LDLIBS += -latlas
endif
ifeq ($(ND),pn1246281)
LDLIBS += -latlas
endif
ifeq ($(ND),macmonty.lanl.gov)
LDLIBS += -latlas
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
OBJSMADS = ./mads.o ./mads_io.o ./mads_io_external.o ./mads_func.o ./mads_mem.o ./mads_info.o lm/opt_lm_mon.o lm/opt_lm_gsl.o lm/lu.o lm/opt_lm_ch.o misc/test_problems.o misc/anasol_contamination.o misc/io.o lhs/lhs.o
OBJSPSO = pso/pso-tribes-lm.o pso/Standard_PSO_2006.o pso/mopso.o
OBJSA = sa/abagus.o sa/postpua.o sa/global.o sa/do_miser.o
OBJDS = ds/infogap.o ds/glue.o
OBJSMPUN = mprun/mprun.o mprun/mprun_io.o
OBJSKDTREE = misc/kdtree-0.5.5/kdtree.o
OBJSLEVMAR = misc/levmar-2.5/lm_m.o misc/levmar-2.5/Axb.o misc/levmar-2.5/misc.o misc/levmar-2.5/lmlec.o misc/levmar-2.5/lmbc.o misc/levmar-2.5/lmblec.o misc/levmar-2.5/lmbleic.o 
OBJSlEVMARSTYLE = misc/levmar-2.5/lm_m.o misc/levmar-2.5/lm_core_m.o misc/levmar-2.5/Axb.o misc/levmar-2.5/misc.o misc/levmar-2.5/lmlec.o misc/levmar-2.5/lmbc.o misc/levmar-2.5/lmblec.o misc/levmar-2.5/lmbleic.o 
OBJSASTABLE = misc/astable/astable.o misc/astable/interpolation.o misc/astable/pqueue.o
OBJSBAYES = bayes/dream.o
OBJWELLS = wells/wells.o

ifeq ($(YAML),true)
$(info YAML Support included)
OBJSMADS += ./mads_io_yaml.o
CFLAGS += -DMADS_YAML `pkg-config --cflags glib-2.0`
LDLIBS += -lyaml `pkg-config --libs glib-2.0`
endif

ifeq ($(XML),true)
    $(info XML Support included)
    OBJSMADS += ./mads_io_xml.o
    CFLAGS += -DMADS_XML
    LDLIBS += `xml2-config --libs`
endif

ifeq ($(MATHEVAL),true)
$(info MathEval Support included)
CFLAGS += -DMATHEVAL
LDLIBS += -lmatheval
endif

SOURCE = $(OBJSMADS:%.o=%.c) $(OBJSPSO:%.o=%.c) $(OBJSMPUN:%.o=%.c) $(OBJSA:%.o=%.c) $(OBJDS:%.o=%.c) $(OBJSLEVMAR:%.o=%.c) $(OBJSKDTREE:%.o=%.c) $(OBJSASTABLE:%.o=%.c) $(OBJSBAYES:%.o=%.cpp)
SOURCESTYLE = $(OBJSMADS:%.o=%.c) $(OBJSPSO:%.o=%.c) $(OBJSMPUN:%.o=%.c) $(OBJSA:%.o=%.c) $(OBJDS:%.o=%.c) $(OBJSLEVMARSTYLE:%.o=%.c) $(OBJSKDTREE:%.o=%.c) $(OBJSASTABLE:%.o=%.c) $(OBJSBAYES:%.o=%.cpp)
SOURCESTYLEDEL = $(OBJSMADS:%.o=%.c.orig) $(OBJSPSO:%.o=%.c.orig) $(OBJSMPUN:%.o=%.c.orig) $(OBJSA:%.o=%.c.orig) $(OBJDS:%.o=%.c.orig) $(OBJSLEVMARSTYLE:%.o=%.c.orig) $(OBJSKDTREE:%.o=%.c.orig) $(OBJSASTABLE:%.o=%.c.orig) $(OBJSBAYES:%.o=%.cpp.orig)

all: $(MADS) $(WELLS)

release: $(MADS)

debug: CFLAGS += -g
debug: $(MADS)

$(MADS): $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSA) $(OBJDS) $(OBJSLEVMAR) $(OBJSKDTREE) $(OBJSASTABLE) $(OBJSBAYES)

$(WELLS): $(OBJWELLS)

clean:
	rm -f $(MADS) $(WELLS) $(OBJWELLS) $(OBJSMADS) $(OBJSPSO) $(OBJSMPUN) $(OBJSA) $(OBJDS) $(OBJSLEVMAR) $(OBJSKDTREE) $(OBJSASTABLE) $(OBJSBAYES)


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
wells/wells.o: wells/wells.c wells/design.h wells/wells.h
bayes/dream.o: bayes/dream.cpp bayes/dream.h mads.h misc/astable/interpolation.o
	g++ $(CFLAGS) -c -o bayes/dream.o bayes/dream.cpp -lmisc/astable/interpolation.o

$(info -----------------------------)
$(info MADS testing and verification)

## Colordefinition
NO_COLOR    = \x1b[0m
OK_COLOR    = \x1b[32;01m
WARN_COLOR  = \x1b[33;01m
ERROR_COLOR = \x1b[31;01m

examples:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "Example 1: Internal Rosenbrock Problem"
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	cd example/rosenbrock; ../../mads a01 test=3 opt=pso igrnd real=1
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "Example 1: DONE"
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "Example 2: Internal contamination Problem "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	mads example/contamination/s01 ldebug
	@echo "Example problem example/wells/w01 ..."
	cd example/wells; ../../mads w01 lmeigen
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "Example 2: DONE"
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"

verify: verify-internal verify-multistart1 verify-contaminant verify-multistart2 verify-external verify-external-short verify-parallel verify-forward verify-sa
	@echo "$(OK_COLOR)"
	@echo VERIFICATION DONE
	@echo "$(NO_COLOR)"


verify-internal:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Internal test problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 1: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 1.1: Levenberg-Marquardt ... "
	rm -f example/rosenbrock/a01.mads_output example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; $(DBG) ../../mads a01 test=3 opt=lm lmeigen igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/rosenbrock/a01.mads_output example/rosenbrock/a01.mads_output-lm-$(OS)-correct
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-lm-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 1.2: Particle-Swarm ..."
	rm -f example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; $(DBG) ../../mads a01 test=3 opt=pso igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-pso-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 1.3: TRIBES ..."
	rm -f example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; $(DBG) ../../mads a01 test=3 opt=tribes igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-tribes-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 1.4: SQUADS ..."
	rm -f example/rosenbrock/a01.results example/rosenbrock/a01.running
	cd example/rosenbrock; $(DBG) ../../mads a01 test=3 opt=squads igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/rosenbrock/a01.results example/rosenbrock/a01.results-squads-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 1: DONE"

verify-multistart1:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Multistart (paranoid) problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 2: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 2.1: Levenberg-Marquardt ... "
	rm -f example/rosenbrock/p01lm.mads_output example/rosenbrock/p01lm.results example/rosenbrock/p01lm.running
	cd example/rosenbrock; $(DBG) ../../mads p01lm test=3 dim=3 opt=lm igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/rosenbrock/p01lm.mads_output example/rosenbrock/p01lm.mads_output-multistart-lm-$(OS)-correct
	@$(CMP) example/rosenbrock/p01lm.results example/rosenbrock/p01lm.results-multistart-lm-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 2.2: SQUADS ..."
	rm -f example/rosenbrock/p01squads.mads_output example/rosenbrock/p01squads.results example/rosenbrock/p01squads.running
	cd example/rosenbrock; $(DBG) ../../mads p01squads test=3 dim=3 opt=squads igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/rosenbrock/p01squads.mads_output example/rosenbrock/p01squads.mads_output-multistart-squads-$(OS)-correct
	@$(CMP) example/rosenbrock/p01squads.results example/rosenbrock/p01squads.results-multistart-squads-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo "$(NO_COLOR)"

verify-contaminant:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 3: Internal contaminant transport problems using different dispersivities ..."
	@echo "TEST 3.1.a: Problem example/contamination/s01 with independent dispersivities (MADS text input format) ..."
	rm -f example/contamination/s01.mads_output example/contamination/s01.results example/contamination/s01.running
	./mads example/contamination/s01 obs_int=2 lmeigen $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01.mads_output example/contamination/s01.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01.results example/contamination/s01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.1.b: Problem example/contamination/s01 with independent dispersivities (YAML input format) ..."
	rm -f example/contamination/s01_yaml.mads_output example/contamination/s01_yaml.results example/contamination/s01_yaml.running
	./mads example/contamination/s01_yaml obs_int=2 lmeigen $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01_yaml.mads_output example/contamination/s01_yaml.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01_yaml.results example/contamination/s01_yaml.results-$(OS)-correct
	@$(CMP) example/contamination/s01_yaml.results example/contamination/s01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.2: Problem example/contamination/s01 with tied dispersivities ..."
	rm -f example/contamination/s01-tied_dispersivities.results example/contamination/s01-tied_dispersivities.running
	./mads example/contamination/s01-tied_dispersivities obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-tied_dispersivities.results example/contamination/s01-tied_dispersivities.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.3: Problem example/contamination/s01 with scaled dispersivities ..."
	rm -f example/contamination/s01-scaled_dispersivities.results example/contamination/s01-scaled_dispersivities.running
	./mads example/contamination/s01-scaled_dispersivities obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-scaled_dispersivities.results example/contamination/s01-scaled_dispersivities.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.4: Problem example/contamination/s01 with scaled and tied dispersivities ..."
	rm -f example/contamination/s01-scaled+tied_dispersivities.results example/contamination/s01-scaled+tied_dispersivities.running
	$(DBG) ./mads example/contamination/s01-scaled+tied_dispersivities obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-scaled+tied_dispersivities.results example/contamination/s01-scaled+tied_dispersivities.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.5: Problem example/contamination/s01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f example/contamination/s01-tied.results example/contamination/s01-tied.running
	$(DBG) ./mads example/contamination/s01-tied obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-tied.results example/contamination/s01-tied.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.6: Problem example/contamination/s01_yaml with coupled (tied) parameters based on mathematical expressions (YAML input format) ..."
	rm -f example/contamination/s01-tied_yaml.results example/contamination/s01-tied_yaml.running
	$(DBG) ./mads example/contamination/s01-tied_yaml obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-tied_yaml.mads_output example/contamination/s01-tied_yaml.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-tied_yaml.results example/contamination/s01-tied.results-$(OS)-correct
	@$(CMP) example/contamination/s01-tied_yaml.results example/contamination/s01-tied_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.7: Problem example/contamination/s01 with regularization terms for optimized model parameters ..."
	rm -f example/contamination/s01-regul.results example/contamination/s01-regul.running
	$(DBG) ./mads example/contamination/s01-regul obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-regul.mads_output example/contamination/s01-regul.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-regul.results example/contamination/s01-regul.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.8: Problem example/contamination/s01_yaml with regularization terms for optimized model parameters (YAML input format) ..."
	rm -f example/contamination/s01-regul_yaml.results example/contamination/s01-regul_yaml.running
	$(DBG) ./mads example/contamination/s01-regul_yaml obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-regul_yaml.mads_output example/contamination/s01-regul_yaml.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-regul_yaml.results example/contamination/s01-regul.results-$(OS)-correct
	@$(CMP) example/contamination/s01-regul_yaml.results example/contamination/s01-regul_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	rm -f example/contamination/s01-multi_source.results example/contamination/s01-multi_source.running
	$(DBG) ./mads example/contamination/s01-multi_source obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-multi_source.mads_output example/contamination/s01-multi_source.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-multi_source.results example/contamination/s01-multi_source.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 3: DONE"
	@echo "$(NO_COLOR)"
	@echo "$(NO_COLOR)"

verify-multistart2:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 4: Internal contaminant transport problems using different optimization techniques ..."
	@echo "TEST 4.1: Problem example/contamination/s01 IGRND ..."
	rm -f example/contamination/s01-igrnd.results example/contamination/s01-igrnd.igrnd.results example/contamination/s01-igrnd.running
	$(DBG) ./mads example/contamination/s01-igrnd seed=2096575428 $(OUTPUT)
	@$(CMP) example/contamination/s01-igrnd.results example/contamination/s01-igrnd.results-$(OS)-correct
	@$(CMP) example/contamination/s01-igrnd.igrnd.results example/contamination/s01-igrnd.igrnd.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 4.2: Problem example/contamination/s01 PPSD ..."
	rm -f example/contamination/s01-ppsd.mads_output example/contamination/s01-ppsd.ppsd.results example/contamination/s01-ppsd.running
	$(DBG) ./mads example/contamination/s01-ppsd seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-ppsd.mads_output example/contamination/s01-ppsd.mads_output-$(OS)-correct
	@$(CMP) example/contamination/s01-ppsd.ppsd.results example/contamination/s01-ppsd.ppsd.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 4.3: Problem example/contamination/s01 IGPD ..."
	rm -f example/contamination/s01-igpd.results example/contamination/s01-igpd.igpd.results example/contamination/s01-igpd.running
	$(DBG) ./mads example/contamination/s01-igpd seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/contamination/s01-igpd.results example/contamination/s01-igpd.results-$(OS)-correct
	@$(CMP) example/contamination/s01-igpd.igpd.results example/contamination/s01-igpd.igpd.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 4.4: Problem example/contamination/s01 Multi-Start LM  ..."
	rm -f example/contamination/s01-mslm.results example/contamination/s01-mslm.running
	$(DBG) ./mads example/contamination/s01-mslm seed=2096575428 $(OUTPUT)
	@$(CMP) example/contamination/s01-mslm.results example/contamination/s01-mslm.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 4: DONE"
	@echo "$(NO_COLOR)"
	@echo "$(NO_COLOR)"

verify-external:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 5: Problem example/wells/w01 ..."
	rm -f example/wells/w01.mads_output example/wells/w01.results example/wells/w01.running
	cd example/wells; $(DBG) ../../mads w01 lmeigen $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells/w01.mads_output example/wells/w01.mads_output-$(OS)-correct
	@$(CMP) example/wells/w01.results example/wells/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 5: DONE"
	@echo "$(NO_COLOR)"

verify-external-short:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 6: Problem example/wells-short/w01 using different instruction formats ..."
	@echo "TEST 6.1: Instruction file example/wells-short/w01-v1.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../mads w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.2: Instruction file example/wells-short/w01-v2.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v2.inst w01.inst; $(DBG) ../../mads w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.3: Instruction file example/wells-short/w01-v3.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v3.inst w01.inst; $(DBG) ../../mads w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.4: Instruction file example/wells-short/w01-v4.inst ..."
	rm -f example/wells-short/w01.results example/wells-short/w01.running
	cd example/wells-short; ln -sf w01-v4.inst w01.inst; $(DBG) ../../mads w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01.results example/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.5: Problem example/wells-short/w01 (YAML input format) ..."
	rm -f example/wells-short/w01_yaml.results example/wells-short/w01_yaml.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../mads w01_yaml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01_yaml.results example/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.6: Problem example/wells-short/w01 (XML input format) ..."
	rm -f example/wells-short/w01_xml.results example/wells-short/w01_xml.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../mads w01_xml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01_xml.results example/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 6: DONE"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 7: Problem example/wells-short/w01 using coupled parameters and regularization terms ..."
	@echo "TEST 7.1: Problem example/wells-short/w01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f example/wells-short/w01tied.results example/wells-short/w01tied.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../mads w01tied $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01tied.results example/wells-short/w01tied.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 7.2: Problem example/wells-short/w01 with regularization terms for optimized model parameters (MADS text input format) ..."
	rm -f example/wells-short/w01regul.results example/wells-short/w01regul.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../mads w01regul $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01regul.results example/wells-short/w01regul.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 7.3: Problem example/wells-short/w01 with regularization terms for optimized model parameters (YAML input format) ..."
	rm -f example/wells-short/w01regul_yaml.results example/wells-short/w01regul_yaml.running
	cd example/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../mads w01regul_yaml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01regul_yaml.results example/wells-short/w01regul_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 7: DONE"
	@echo "$(NO_COLOR)"

verify-parallel:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Parallel execution of external problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 8: Parallel execution of example/wells-short/w01parallel ..."
	@echo "TEST 8.1: Initial parallel execution of example/wells-short/w01parallel ..."
	rm -f example/wells-short/w01parallel.results example/wells-short/w01parallel.restart_info example/wells-short/w01parallel.restart_*.zip example/wells-short/w01parallel.running
	cd example/wells-short; $(DBG) ../../mads w01parallel np=2 eval=10 restart=0 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01parallel.results example/wells-short/w01parallel.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 8.2: Rerun using saved results from prior parallel execution of example/wells-short/w01parallel ..."
	rm -f example/wells/w01parallel.results example/wells/w01parallel.running
	cd example/wells-short; $(DBG) ../../mads w01parallel np=2 eval=10 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01parallel.results example/wells-short/w01parallel.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 8: DONE"
	@echo "$(NO_COLOR)"

verify-forward:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Forward analytical contaminant modeling "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 9: Analytical contaminant concentrations with various sources ..."
	@echo "TEST 9.1: Point source ... "
	rm -f example/forward/a01.mads_output example/forward/a01.results example/forward/a01.running
	cd example/forward; $(DBG) ../../mads a01 $(OUTPUT)
	@$(CMP) example/forward/a01.mads_output example/forward/a01.mads_output-$(OS)-correct
	@$(CMP) example/forward/a01.results example/forward/a01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 9.2: Rectangular source ... "
	rm -f example/forward/a02.mads_output example/forward/a02.results example/forward/a02.running
	cd example/forward; $(DBG) ../../mads a02 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a02.mads_output example/forward/a02.mads_output-$(OS)-correct
	@$(CMP) example/forward/a02.results example/forward/a02.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 9.3: Planar gaussian source ... "
	rm -f example/forward/a03.mads_output example/forward/a03.results example/forward/a03.running
	cd example/forward; $(DBG) ../../mads a03 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a03.mads_output example/forward/a03.mads_output-$(OS)-correct
	@$(CMP) example/forward/a03.results example/forward/a03.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 9.4: Gaussian source ... "
	rm -f example/forward/a04.mads_output example/forward/a04.results example/forward/a04.running
	cd example/forward; $(DBG) ../../mads a04 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a04.mads_output example/forward/a04.mads_output-$(OS)-correct
	@$(CMP) example/forward/a04.results example/forward/a04.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.5: Box source ... "
	rm -f example/forward/a05.mads_output example/forward/a05.results example/forward/a05.running
	cd example/forward; $(DBG) ../../mads a05 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a05.mads_output example/forward/a05.mads_output-$(OS)-correct
	@$(CMP) example/forward/a05.results example/forward/a05.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.6: Box source | Levy ... "
	rm -f example/forward/a06_yaml.mads_output example/forward/a06_yaml.results example/forward/a06_yaml.running
	cd example/forward; $(DBG) ../../mads a06_yaml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a06_yaml.mads_output example/forward/a06_yaml.mads_output-$(OS)-correct
	@$(CMP) example/forward/a06_yaml.results example/forward/a06_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.7: Box source | Symmetric Levy ... "
	rm -f example/forward/a07.mads_output example/forward/a07.results example/forward/a07.running
	cd example/forward; $(DBG) ../../mads a07 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a07.mads_output example/forward/a07.mads_output-$(OS)-correct
	@$(CMP) example/forward/a07.results example/forward/a07.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.8: Box source | Levy, alpha=2 ... "
	rm -f example/forward/a08.mads_output example/forward/a08.results example/forward/a08.running
	cd example/forward; $(DBG) ../../mads a08 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a08.mads_output example/forward/a08.mads_output-$(OS)-correct
	@$(CMP) example/forward/a08.results example/forward/a08.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.9: Box source | Symmetric Levy, alpha=2 ... "
	rm -f example/forward/a09.mads_output example/forward/a09.results example/forward/a09.running
	cd example/forward; $(DBG) ../../mads a09 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/forward/a09.mads_output example/forward/a09.mads_output-$(OS)-correct
	@$(CMP) example/forward/a09.results example/forward/a09.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 9: DONE"
	@echo "$(NO_COLOR)"

verify-sa:
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Sensitivity analyses "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 10: Global and local sensitivity analyses ..."
	@echo "TEST 10.1: Sobol internal analysis ... "
	rm -f example/sa/a01.mads_output example/sa/a01.sobol_sens_index example/sa/a01.sobol_sens_total example/sa/a01.running
	cd example/sa; $(DBG) ../../mads a01 test=111 gsens dim=8 real=10000 smp=lhs pardomain=0.5 seed=1517604820 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/sa/a01.sobol_sens_index example/sa/a01.sobol_sens_index-$(OS)-correct
	@$(CMP) example/sa/a01.sobol_sens_total example/sa/a01.sobol_sens_total-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 10.2: Sobol external analysis  ..."
	rm -f example/wells-short/w01.sobol_sens_index example/wells-short/w01.sobol_sens_total example/wells-short/w01.running
	cd example/wells-short; $(DBG) ../../mads w01 gsens real=100 seed=1066732675 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01.sobol_sens_index example/wells-short/w01.sobol_sens_index-$(OS)-correct
	@$(CMP) example/wells-short/w01.sobol_sens_total example/wells-short/w01.sobol_sens_total-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 10.2: Sobol external analysis  ..."
	rm -f example/wells-short/w01.sobol_sens_index example/wells-short/w01.sobol_sens_total example/wells-short/w01.running example/wells-short/w01.restart*
	cd example/wells-short; $(DBG) ../../mads w01 gsens real=100 seed=1066732675 np=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) example/wells-short/w01.sobol_sens_index example/wells-short/w01.sobol_sens_index-$(OS)-correct
	@$(CMP) example/wells-short/w01.sobol_sens_total example/wells-short/w01.sobol_sens_total-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 10: DONE"
	@echo "$(NO_COLOR)"

compare-os:
	./compare-results-os Linux Darwin

clean-example:
	rm -f example/*/*.mads_output_* example/*/*.ppsd_*.results example/*/*.igpd_*.results example/*/*.igrnd_*.results example/*/*.restart_*.zip example/*/*.restart_info example/*/*.running example/*/*-rerun.mads example/*/*-error.mads
	rm -fR example/wells-short_w01_*
	rm -fR example/wells-short_w01parallel*
	rm -f *.mads_output* *.running *.cmdline *.cmdline_hist

astyle:
	astyle $(SOURCESTYLE)
	rm -f $(SOURCESTYLEDEL)

tar:
	rm -f mads.git.tgz
	tar -cvzf mads.git.tgz `git ls-files` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'` .git

tarf:
	tar -cvzf mads.tgz `git ls-files` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'`
