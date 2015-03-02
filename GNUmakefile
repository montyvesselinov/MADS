# MADS Model Analyses & Decision Support (v.1.1.14) 2013
#
# Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
# Dan O'Malley, omalled@lanl.gov
#
# http:./$(MADS).lanl.gov
# http://www.ees.lanl.gov/staff/monty/code./$(MADS)
# http://gitlab.com/mont./$(MADS)
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
# others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in thisDS
# material to reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY, LLC,
# NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR
# PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

## Color definitions
NO_COLOR    = \033[0m
OK_COLOR    = \033[0;32;40m
WARN_COLOR  = \033[0;33;40m
ERROR_COLOR = \033[0;31;40m
CCWARN=$(echo -e "$(WARN_COLOR)")
CCERROR=$(echo -e "$(ERROR_COLOR)")
CCNOCOL=$(echo -e "$(NO_COLOR)")
PATHPATERN="(/[^/]*)+:[0-9]+"
# 2>&1 | sed -E -e "/[Ee]rror[: ]/ s%$PATHPATERN%$CCERROR&$CCNOCOL%g" -e "/[Ww]arning[: ]/ s%$PATHPATERN%$CCWARN&$CCNOCOL%g"

SRC = ./src
OBJ = ./obj
BIN = ./bin
OBJ_DIR ?= $(OBJ)/Release
ifeq ($(MAKECMDGOALS),debug)
    OBJ_DIR = $(OBJ)/Debug
else
    ifeq ($(MAKECMDGOALS),lib)
	OBJ_DIR = $(OBJ)/Lib
    endif
endif
EXAMPLES = ./examples
MADS = $(BIN)/Release/mads
MADS_DEBUG = $(BIN)/Debug/mads
MADS_LIB = $(BIN)/Lib/libmads.so.1.1.14
WELLS = $(BIN)/wells

CMP = ./scripts/compare-results # MADS testing
# CMP = cp -f # Save current results for future testing DANGEROUS!

OUTPUT = > /dev/null
# OUTPUT = mads-debug-output
# OUTPUT =

#CFLAGS += -fsanitize=address -fno-omit-frame-pointer
#LDLIBS += -fsanitize=address -fno-omit-frame-pointer
DBG = valgrind -v --read-var-info=yes --tool=memcheck --leak-check=yes --leak-check=full --show-reachable=yes --track-origins=yes
DBG = gdb --args
DBG =

VER = $(shell git rev-parse --short HEAD)
GIT_STATUS = $(shell scripts/check_git_status)

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

# OS type
OS = $(shell uname -s)
# MACHINE name
ND = $(shell uname -n)

# Compilation setup
$(info MADS computationlal framework)
$(info Version -- $(VER)$(GIT_STATUS))
$(info ----------------------------------------------------------------------)
$(info OS type -- $(OS))
$(info Machine -- $(ND))
$(info Target  -- $(MAKECMDGOALS) $(OBJ_DIR))

CC = gcc
CXX = g++
CFLAGS = -Wall -Winit-self
LDLIBS = -lgsl -llapack -lstdc++
SONAME = soname

ifeq ($(OS),Linux)
# Linux
LDLIBS += -lgslcblas -lm -lblas
$(info LINUX)
ifeq ($(ND),aquifer.lanl.gov)
$(info Machine -- AQUIFER)
CFLAGS += -I/home/monty/local/include-aquifer
LDLIBS = -lgfortran -lgsl -llapack -lstdc++ -L/home/monty/local/lib -lgslcblas -lgfortran -Wl,--rpath,/home/monty/local/lib 
endif
ifeq ($(ND),madsmax)
$(info Machine -- MadsMax)
CFLAGS += -Wno-unused-result
endif
ifeq ($(ND),well.lanl.gov)
$(info Machine -- WELL)
CFLAGS += -I/home/monty/local/include
LDLIBS += -L/home/monty/local/lib -L/usr/local/lib -Wl,--rpath,/home/monty/local/lib
endif
ifeq ($(ND),pi-fe1.lanl.gov)
$(info Machine -- turquoise)
CFLAGS += -I/users/vvv/mads/repo-scp/tpls/include
LDLIBS += -L/users/vvv/mads/repo-scp/tpls/lib -Wl,--rpath,/users/vvv/mads/repo-scp/tpls/lib
endif
else #----------------------------------------------------
ifeq ($(OS),Cygwin)
# Cygwin
$(info CYGWIN)
else #----------------------------------------------------
ifeq ($(OS),Darwin)
# Mac
$(info MAC OS X)
SONAME = install_name
CFLAGS += -I/opt/local/include
LDLIBS += -lgfortran -lblas -L/opt/local/lib
ifeq ($(ND),bored.lanl.gov)
LDLIBS += 
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
LDLIBS += -lrefblas -lcblas -L/Users/monty/lib -Wl,--rpath,/Users/monty/lib
endif
else #----------------------------------------------------
$(error UNKNOWN OS type -- $(OS)!)
endif
endif
endif

# MADS files
SRC_MADS = $(SRC)/mads.c $(SRC)/mads_io.c $(SRC)/mads_io_external.c $(SRC)/mads_func.c $(SRC)/mads_mem.c $(SRC)/mads_info.c $(SRC)/lm/opt_lm_mon.c $(SRC)/lm/opt_lm_gsl.c $(SRC)/lm/lu.c $(SRC)/lm/opt_lm_ch.c $(SRC)/misc/test_problems.c $(SRC)/misc/anasol_contamination.c $(SRC)/misc/io.c $(SRC)/lhs/lhs.c $(SRC)/mads_gitversion.c
SRC_MADSSTYLE = $(SRC)/mads.c $(SRC)/mads_io.c $(SRC)/mads_io_external.c $(SRC)/mads_func.c $(SRC)/mads_mem.c $(SRC)/mads_info.c $(SRC)/lm/opt_lm_mon.c $(SRC)/lm/opt_lm_gsl.c $(SRC)/lm/lu.c $(SRC)/lm/opt_lm_ch.c $(SRC)/misc/test_problems.c $(SRC)/misc/anasol_contamination.c $(SRC)/misc/io.c $(SRC)/lhs/lhs.c
SRC_PSO = $(SRC)/pso/pso-tribes-lm.c $(SRC)/pso/Standard_PSO_2006.c $(SRC)/pso/mopso.c
SRC_SA = $(SRC)/sa/abagus.c $(SRC)/sa/postpua.c $(SRC)/sa/global.c $(SRC)/sa/do_miser.c
SRC_DS = $(SRC)/ds/infogap.c $(SRC)/ds/glue.c
SRC_MPUN = $(SRC)/mprun/mprun.c $(SRC)/mprun/mprun_io.c
SRC_KDTREE = $(SRC)/misc/kdtree-0.5.5/kdtree.c
SRC_LEVMAR = $(SRC)/misc/levmar-2.5/lm_m.c $(SRC)/misc/levmar-2.5/Axb.c $(SRC)/misc/levmar-2.5/misc.c $(SRC)/misc/levmar-2.5/lmlec.c $(SRC)/misc/levmar-2.5/lmbc.c $(SRC)/misc/levmar-2.5/lmblec.c $(SRC)/misc/levmar-2.5/lmbleic.c 
SRC_LEVMARSTYLE = $(SRC)/misc/levmar-2.5/lm_m.c $(SRC)/misc/levmar-2.5/lm_core_m.c $(SRC)/misc/levmar-2.5/Axb.c $(SRC)/misc/levmar-2.5/misc.c $(SRC)/misc/levmar-2.5/lmlec.c $(SRC)/misc/levmar-2.5/lmbc.c $(SRC)/misc/levmar-2.5/lmblec.c $(SRC)/misc/levmar-2.5/lmbleic.c 
SRC_ASTABLE = $(SRC)/misc/astable/astable.c $(SRC)/misc/astable/interpolation.c $(SRC)/misc/astable/pqueue.c
SRC_BAYES = $(SRC)/bayes/dream.cpp
SRC_WELLS = wells/wells.c

ifeq ($(YAML),true)
    $(info YAML Support included)
    SRC_MADS += $(SRC)/mads_io_yaml.c
    CFLAGS += -DMADS_YAML `pkg-config --cflags glib-2.0`
    LDLIBS += -lyaml `pkg-config --libs glib-2.0`
endif

ifeq ($(XML),true)
    $(info XML Support included)
    SRC_MADS += $(SRC)/mads_io_xml.c
    CFLAGS += -DMADS_XML
    LDLIBS += `xml2-config --libs`
endif

ifeq ($(MATHEVAL),true)
    $(info MathEval Support included)
    CFLAGS += -DMATHEVAL
    LDLIBS += -lmatheval
endif
CXXFLAGS = $(CFLAGS)
$(info ----------------------------------------------------------------------)

SOURCES = $(SRC_MADS) $(SRC_PSO) $(SRC_MPUN) $(SRC_SA) $(SRC_DS) $(SRC_LEVMAR) $(SRC_KDTREE) $(SRC_ASTABLE) $(SRC_BAYES)
SOURCESTYLE = $(SRC_MADSSTYLE) $(SRC_PSO) $(SRC_MPUN) $(SRC_SA) $(SRC_DS) $(SRC_LEVMARSTYLE) $(SRC_KDTREE) $(SRC_ASTABLE) $(SRC_BAYES)
SOURCESTYLEDEL = $(addsuffix .orig,$(SOURCESTYLE))
OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))
OBJECTS_RELEASE = $(patsubst $(SRC)/%,$(OBJ)/Release/%,$(OBJECTS))
OBJECTS_DEBUG = $(patsubst $(SRC)/%,$(OBJ)/Debug/%,$(OBJECTS))
OBJECTS_LIB =   $(patsubst $(SRC)/%,$(OBJ)/Lib/%,$(OBJECTS))
OBJ_WELLS = $(addsuffix .o,$(basename $(SRC_WELLS)))

all: mads wells

mads: CFLAGS += -O3
mads: release-start $(MADS)
	ln -sf ${MADS} .
	@echo "$(OK_COLOR)"
	@echo "MADS Release Version built!"
	@echo "$(NO_COLOR)"

release-start:
	@echo "$(OK_COLOR)"
	@echo "MADS Release Version ..."
	@echo "$(NO_COLOR)"

wells: $(WELLS)
	@echo "$(OK_COLOR)"
	@echo "WELLS built!"
	@echo "$(NO_COLOR)"

release: release-start $(MADS)
	ln -sf ${MADS} .
	@echo "$(OK_COLOR)"
	@echo "MADS Release Version built!"
	@echo "Execute ${MADS}"
	@echo "$(NO_COLOR)"

debug: CFLAGS += -g
debug: debug-start $(MADS_DEBUG)
	@echo "$(OK_COLOR)"
	@echo "MADS Debug Version built!"
	@echo "Execute ${MADS_DEBUG}"
	@echo "$(NO_COLOR)"

debug-start:
	@echo "$(OK_COLOR)"
	@echo "MADS Debug Version ..."
	@echo "$(NO_COLOR)"

lib: CFLAGS += -fPIC -O3
lib: lib-start $(MADS_LIB) lib-install
	@echo "$(OK_COLOR)"
	@echo "MADS Shared Library built!"
	@echo "$(NO_COLOR)"

lib-start:
	@echo "$(OK_COLOR)"
	@echo "MADS Library ..."
	@echo "$(NO_COLOR)"

lib-install:
	@echo "$(OK_COLOR)"
	@echo "Install MADS Library ..."
	@echo "$(NO_COLOR)"
	cp $(MADS_LIB) /opt/local/lib/libmads.so.1.0
	ln -sf /opt/local/lib/libmads.so.1.0 /opt/local/lib/libmads.so.1
	ln -sf /opt/local/lib/libmads.so.1.0 /opt/local/lib/libmads.dylib
	ln -sf /opt/local/lib/libmads.so.1 /opt/local/lib/libmads.so

$(MADS): $(OBJECTS_RELEASE)
	@echo "$(OK_COLOR)"
	@echo "Building MADS Release Version ..."
	@echo "$(NO_COLOR)"
	@mkdir -p $(BIN)
	@mkdir -p $(OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(OBJECTS_RELEASE) $(LDLIBS) -o $@

$(MADS_DEBUG): $(OBJECTS_DEBUG)
	@echo "$(OK_COLOR)"
	@echo "Building MADS Debug Version ..."
	@echo "$(NO_COLOR)"
	@mkdir -p $(BIN)
	@mkdir -p $(OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(OBJECTS_DEBUG) $(LDLIBS) -o $@

$(MADS_LIB): $(OBJECTS_LIB)
	@echo "$(OK_COLOR)"
	@echo "Building MADS Shared Library ..."
	@echo "$(NO_COLOR)"
	@mkdir -p $(BIN)
	@mkdir -p $(OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(LDLIBS) $(OBJECTS_LIB) -shared -Wl,-$(SONAME),libmads.so.1 -o $@

$(WELLS): $(OBJ_WELLS)
	@mkdir -p $(BIN)
	$(CC) $< -o $@ $(LDFLAGS) -lm

$(OBJ_WELLS): $(SRC_WELLS)
	@echo "$(OK_COLOR)"
	@echo "Building WELLS ..."
	@echo "$(NO_COLOR)"
	$(CC) $(CFLAGS) -c $< -o $@

clean clean-release:
	rm -f $(MADS) $(OBJECTS_RELEASE)

clean-debug:
	rm -f $(MADS_DEBUG) $(OBJECTS_DEBUG)

clean-lib:
	rm -f $(MADS_LIB) $(OBJECTS_LIB)

clean-wells:
	rm -f $(WELLS) $(OBJ_WELLS)

clean-old-mads-setup:
	rm -f *.o
	rm -fR bayes ds sa misc lhs lm example pso mprun

clean-all: clean-release clean-debug clean-lib clean-wells clean-examples
	rm -fR $(OBJ)/*
	rm -fR $(BIN)/*

# $(SRC)/mads_gitversion.c: .git/HEAD .git/COMMIT_EDITMSG
$(SRC)/mads_gitversion.c: .git/HEAD
	@echo "const char *gitversion = \"$(VER)$(GIT_STATUS)\";" > $@

$(OBJ_DIR)/%.o: $(SRC)/%.c
	@mkdir -p $(OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/mads.o: $(SRC)/mads.c $(SRC)/mads.h $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/mads_gitversion.c
$(OBJ_DIR)/mads_io.o: $(SRC)/mads_io.c $(SRC)/mads.h
$(OBJ_DIR)/mads_io_yaml.o: $(SRC)/mads_io_yaml.c $(SRC)/mads.h
$(OBJ_DIR)/mads_io_external.o: $(SRC)/mads_io_external.c $(SRC)/mads.h
$(OBJ_DIR)/mads_func.o: $(SRC)/mads_func.c $(SRC)/mads.h
$(OBJ_DIR)/mads_mem.o: $(SRC)/mads_mem.c
$(OBJ_DIR)/mads_info.o: $(SRC)/mads_info.c
$(OBJ_DIR)/mads_gitversion.o: $(SRC)/mads_gitversion.c
$(OBJ_DIR)/mprun/mprun.o: $(SRC)/mprun/mprun.c $(SRC)/mads.h
$(OBJ_DIR)/mprun/mprun_io.o: $(SRC)/mprun/mprun_io.c $(SRC)/mads.h
$(OBJ_DIR)/lm/opt_lm_mon.o: $(SRC)/lm/opt_lm_mon.c $(SRC)/mads.h
$(OBJ_DIR)/lm/opt_lm_gsl.o: $(SRC)/lm/opt_lm_gsl.c $(SRC)/mads.h
$(OBJ_DIR)/lm/opt_lm_ch.o: $(SRC)/lm/opt_lm_gsl.c $(SRC)/mads.h
$(OBJ_DIR)/lm/lu.o: $(SRC)/lm/lu.c
$(OBJ_DIR)/lhs/lhs.o: $(SRC)/lhs/lhs.c
$(OBJ_DIR)/pso/pso-tribes-lm.o: $(SRC)/pso/pso-tribes-lm.c $(SRC)/pso/pso.h $(SRC)/mads.h
$(OBJ_DIR)/pso/Standard_PSO_2006.o: $(SRC)/pso/Standard_PSO_2006.c $(SRC)/mads.h
$(OBJ_DIR)/pso/mopso.o: $(SRC)/pso/mopso.c $(SRC)/pso/mopso.h
$(OBJ_DIR)/sa/abagus.o: $(SRC)/sa/abagus.c $(SRC)/mads.h $(SRC)/misc/kdtree-0.5.5/kdtree.o
$(OBJ_DIR)/sa/postpua.o: $(SRC)/sa/postpua.c $(SRC)/mads.h
$(OBJ_DIR)/sa/global.o: $(SRC)/sa/global.c $(SRC)/mads.h $(SRC)/sa/do_miser.o
$(OBJ_DIR)/sa/do_miser.o: $(SRC)/sa/do_miser.c $(SRC)/sa/do_miser.h
$(OBJ_DIR)/ds/infogap.o: $(SRC)/ds/infogap.c $(SRC)/mads.h
$(OBJ_DIR)/ds/glue.o: $(SRC)/ds/glue.c $(SRC)/mads.h
$(OBJ_DIR)/misc/io.o: $(SRC)/misc/io.c
$(OBJ_DIR)/misc/anasol_contamination.o: $(SRC)/misc/anasol_contamination.c $(SRC)/mads.h
$(OBJ_DIR)/misc/test_problems.o: $(SRC)/misc/test_problems.c $(SRC)/mads.h
$(OBJ_DIR)/misc/kdtree-0.5.5/kdtree.o: $(SRC)/misc/kdtree-0.5.5/kdtree.c $(SRC)/misc/kdtree-0.5.5/kdtree.h
$(OBJ_DIR)/misc/levmar-2.5/lm_m.o: $(SRC)/misc/levmar-2.5/lm_m.c $(SRC)/misc/levmar-2.5/lm_core_m.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h $(SRC)/misc/levmar-2.5/compiler.h $(SRC)/mads.h
$(OBJ_DIR)/misc/levmar-2.5/Axb.o: $(SRC)/misc/levmar-2.5/Axb.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h
$(OBJ_DIR)/misc/levmar-2.5/misc.o: $(SRC)/misc/levmar-2.5/misc.c $(SRC)/misc/levmar-2.5/misc_core.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h
$(OBJ_DIR)/misc/levmar-2.5/lmlec.o: $(SRC)/misc/levmar-2.5/lmlec.c $(SRC)/misc/levmar-2.5/lmlec_core.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h
$(OBJ_DIR)/misc/levmar-2.5/lmbc.o: $(SRC)/misc/levmar-2.5/lmbc.c $(SRC)/misc/levmar-2.5/lmbc_core.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h $(SRC)/misc/levmar-2.5/compiler.h
$(OBJ_DIR)/misc/levmar-2.5/lmblec.o: $(SRC)/misc/levmar-2.5/lmblec.c $(SRC)/misc/levmar-2.5/lmblec_core.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h
$(OBJ_DIR)/misc/levmar-2.5/lmbleic.o: $(SRC)/misc/levmar-2.5/lmbleic.c $(SRC)/misc/levmar-2.5/lmbleic_core.c $(SRC)/misc/levmar-2.5/levmar.h $(SRC)/misc/levmar-2.5/misc.h
$(OBJ_DIR)/misc/astable/astable.o: $(SRC)/misc/astable/astable.c $(SRC)/misc/astable/astable.h
$(OBJ_DIR)/misc/astable/interpolation.o: $(SRC)/misc/astable/astable.c $(SRC)/misc/astable/astable.h $(SRC)/misc/astable/interpolation.c $(SRC)/misc/astable/pqueue.c $(SRC)/misc/astable/pqueue.h
$(OBJ_DIR)/misc/astable/pqueue.o: $(SRC)/misc/astable/pqueue.c $(SRC)/misc/astable/pqueue.h
$(OBJ_DIR)/wells/wells.o: $(SRC)/wells/wells.c $(SRC)/wells/design.h $(SRC)/wells/wells.h
$(OBJ_DIR)/bayes/dream.o: $(SRC)/bayes/dream.cpp $(SRC)/bayes/dream.h $(SRC)/mads.h $(SRC)/misc/astable/interpolation.o
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/bayes/dream.o $(SRC)/bayes/dream.cpp -l$(OBJ_DIR)/misc/astable/interpolation.o

examples: mads wells
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "Example 1: Internal Rosenbrock Problem"
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	cd $(EXAMPLES)/rosenbrock; ../../$(MADS) a01 test=3 opt=pso igrnd real=1
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
	mads $(EXAMPLES)/contamination/s01 ldebug
	@echo "Example problem $(EXAMPLES)/wells/w01 ..."
	cd $(EXAMPLES)/wells; ../../$(MADS) w01 lmeigen
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "Example 2: DONE"
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"

verify: mads wells verify-start verify-internal verify-multistart1 verify-contaminant verify-multistart2 verify-external verify-external-short verify-external-short2 verify-parallel verify-forward verify-sa
	@echo "$(OK_COLOR)"
	@echo VERIFICATION DONE
	@echo "$(NO_COLOR)"

verify-start:
	@echo "$(OK_COLOR)"
	@echo "$(OK_COLOR)VERIFICATION STARTING ..."
	@echo "$(NO_COLOR)"

verify-internal: mads
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Internal test problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 1: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 1.1: Levenberg-Marquardt ... "
	rm -f $(EXAMPLES)/rosenbrock/a01.mads_output $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.running
	cd $(EXAMPLES)/rosenbrock; $(DBG) ../../$(MADS) a01 test=3 opt=lm lmeigen igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/rosenbrock/a01.mads_output $(EXAMPLES)/rosenbrock/a01.mads_output-lm-$(OS)-correct
	@$(CMP) $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.results-lm-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 1.2: Particle-Swarm ..."
	rm -f $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.running
	cd $(EXAMPLES)/rosenbrock; $(DBG) ../../$(MADS) a01 test=3 opt=pso igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.results-pso-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 1.3: TRIBES ..."
	rm -f $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.running
	cd $(EXAMPLES)/rosenbrock; $(DBG) ../../$(MADS) a01 test=3 opt=tribes igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.results-tribes-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 1.4: SQUADS ..."
	rm -f $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.running
	cd $(EXAMPLES)/rosenbrock; $(DBG) ../../$(MADS) a01 test=3 opt=squads igrnd real=1 seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/rosenbrock/a01.results $(EXAMPLES)/rosenbrock/a01.results-squads-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 1: DONE"
	@echo "$(NO_COLOR)"

verify-multistart1: mads
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Multistart (paranoid) problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 2: Rosenbrock problem using different optimization techniques ..."
	@echo "TEST 2.1: Levenberg-Marquardt ... "
	rm -f $(EXAMPLES)/rosenbrock/p01lm.mads_output $(EXAMPLES)/rosenbrock/p01lm.results $(EXAMPLES)/rosenbrock/p01lm.running
	cd $(EXAMPLES)/rosenbrock; $(DBG) ../../$(MADS) p01lm test=3 dim=3 opt=lm igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/rosenbrock/p01lm.mads_output $(EXAMPLES)/rosenbrock/p01lm.mads_output-multistart-lm-$(OS)-correct
	@$(CMP) $(EXAMPLES)/rosenbrock/p01lm.results $(EXAMPLES)/rosenbrock/p01lm.results-multistart-lm-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 2.2: SQUADS ..."
	rm -f $(EXAMPLES)/rosenbrock/p01squads.mads_output $(EXAMPLES)/rosenbrock/p01squads.results $(EXAMPLES)/rosenbrock/p01squads.running
	cd $(EXAMPLES)/rosenbrock; $(DBG) ../../$(MADS) p01squads test=3 dim=3 opt=squads igrnd real=10 eval=1000 cutoff=1e-3 seed=2096575428 quiet $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/rosenbrock/p01squads.mads_output $(EXAMPLES)/rosenbrock/p01squads.mads_output-multistart-squads-$(OS)-correct
	@$(CMP) $(EXAMPLES)/rosenbrock/p01squads.results $(EXAMPLES)/rosenbrock/p01squads.results-multistart-squads-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 2: DONE"
	@echo "$(NO_COLOR)"

verify-contaminant: mads
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 3: Internal contaminant transport problems using different dispersivities ..."
	@echo "TEST 3.1.a: Problem $(EXAMPLES)/contamination/s01 with independent dispersivities (MADS text input format) ..."
	rm -f $(EXAMPLES)/contamination/s01.mads_output $(EXAMPLES)/contamination/s01.results $(EXAMPLES)/contamination/s01.running
	./$(MADS) $(EXAMPLES)/contamination/s01 obs_int=2 lmeigen $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01.mads_output $(EXAMPLES)/contamination/s01.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01.results $(EXAMPLES)/contamination/s01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.1.b: Problem $(EXAMPLES)/contamination/s01 with independent dispersivities (YAML input format) ..."
	rm -f $(EXAMPLES)/contamination/s01_yaml.mads_output $(EXAMPLES)/contamination/s01_yaml.results $(EXAMPLES)/contamination/s01_yaml.running
	./$(MADS) $(EXAMPLES)/contamination/s01_yaml obs_int=2 lmeigen $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01_yaml.mads_output $(EXAMPLES)/contamination/s01_yaml.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01_yaml.results $(EXAMPLES)/contamination/s01_yaml.results-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01_yaml.results $(EXAMPLES)/contamination/s01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.2: Problem $(EXAMPLES)/contamination/s01 with tied dispersivities ..."
	rm -f $(EXAMPLES)/contamination/s01-tied_dispersivities.results $(EXAMPLES)/contamination/s01-tied_dispersivities.running
	./$(MADS) $(EXAMPLES)/contamination/s01-tied_dispersivities obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-tied_dispersivities.results $(EXAMPLES)/contamination/s01-tied_dispersivities.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.3: Problem $(EXAMPLES)/contamination/s01 with scaled dispersivities ..."
	rm -f $(EXAMPLES)/contamination/s01-scaled_dispersivities.results $(EXAMPLES)/contamination/s01-scaled_dispersivities.running
	./$(MADS) $(EXAMPLES)/contamination/s01-scaled_dispersivities obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-scaled_dispersivities.results $(EXAMPLES)/contamination/s01-scaled_dispersivities.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.4: Problem $(EXAMPLES)/contamination/s01 with scaled and tied dispersivities ..."
	rm -f $(EXAMPLES)/contamination/s01-scaled+tied_dispersivities.results $(EXAMPLES)/contamination/s01-scaled+tied_dispersivities.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-scaled+tied_dispersivities obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-scaled+tied_dispersivities.results $(EXAMPLES)/contamination/s01-scaled+tied_dispersivities.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.5: Problem $(EXAMPLES)/contamination/s01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f $(EXAMPLES)/contamination/s01-tied.results $(EXAMPLES)/contamination/s01-tied.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-tied obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-tied.results $(EXAMPLES)/contamination/s01-tied.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.6: Problem $(EXAMPLES)/contamination/s01_yaml with coupled (tied) parameters based on mathematical expressions (YAML input format) ..."
	rm -f $(EXAMPLES)/contamination/s01-tied_yaml.results $(EXAMPLES)/contamination/s01-tied_yaml.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-tied_yaml obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-tied_yaml.mads_output $(EXAMPLES)/contamination/s01-tied_yaml.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-tied_yaml.results $(EXAMPLES)/contamination/s01-tied.results-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-tied_yaml.results $(EXAMPLES)/contamination/s01-tied_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.7: Problem $(EXAMPLES)/contamination/s01 with regularization terms for optimized model parameters ..."
	rm -f $(EXAMPLES)/contamination/s01-regul.results $(EXAMPLES)/contamination/s01-regul.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-regul obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-regul.mads_output $(EXAMPLES)/contamination/s01-regul.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-regul.results $(EXAMPLES)/contamination/s01-regul.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 3.8: Problem $(EXAMPLES)/contamination/s01_yaml with regularization terms for optimized model parameters (YAML input format) ..."
	rm -f $(EXAMPLES)/contamination/s01-regul_yaml.results $(EXAMPLES)/contamination/s01-regul_yaml.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-regul_yaml obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-regul_yaml.mads_output $(EXAMPLES)/contamination/s01-regul_yaml.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-regul_yaml.results $(EXAMPLES)/contamination/s01-regul.results-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-regul_yaml.results $(EXAMPLES)/contamination/s01-regul_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	rm -f $(EXAMPLES)/contamination/s01-multi_source.results $(EXAMPLES)/contamination/s01-multi_source.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-multi_source obs_int=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-multi_source.mads_output $(EXAMPLES)/contamination/s01-multi_source.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-multi_source.results $(EXAMPLES)/contamination/s01-multi_source.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 3: DONE"
	@echo "$(NO_COLOR)"

verify-multistart2: mads
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Internal contaminant transport problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 4: Internal contaminant transport problems using different optimization techniques ..."
	@echo "TEST 4.1: Problem $(EXAMPLES)/contamination/s01 IGRND ..."
	rm -f $(EXAMPLES)/contamination/s01-igrnd.results $(EXAMPLES)/contamination/s01-igrnd.igrnd.results $(EXAMPLES)/contamination/s01-igrnd.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-igrnd seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-igrnd.results $(EXAMPLES)/contamination/s01-igrnd.results-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-igrnd.igrnd.results $(EXAMPLES)/contamination/s01-igrnd.igrnd.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 4.2: Problem $(EXAMPLES)/contamination/s01 PPSD ..."
	rm -f $(EXAMPLES)/contamination/s01-ppsd.mads_output $(EXAMPLES)/contamination/s01-ppsd.ppsd.results $(EXAMPLES)/contamination/s01-ppsd.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-ppsd seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-ppsd.mads_output $(EXAMPLES)/contamination/s01-ppsd.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-ppsd.ppsd.results $(EXAMPLES)/contamination/s01-ppsd.ppsd.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 4.3: Problem $(EXAMPLES)/contamination/s01 IGPD ..."
	rm -f $(EXAMPLES)/contamination/s01-igpd.results $(EXAMPLES)/contamination/s01-igpd.igpd.results $(EXAMPLES)/contamination/s01-igpd.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-igpd seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-igpd.results $(EXAMPLES)/contamination/s01-igpd.results-$(OS)-correct
	@$(CMP) $(EXAMPLES)/contamination/s01-igpd.igpd.results $(EXAMPLES)/contamination/s01-igpd.igpd.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 4.4: Problem $(EXAMPLES)/contamination/s01 Multi-Start LM  ..."
	rm -f $(EXAMPLES)/contamination/s01-mslm.results $(EXAMPLES)/contamination/s01-mslm.running
	$(DBG) ./$(MADS) $(EXAMPLES)/contamination/s01-mslm seed=2096575428 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/contamination/s01-mslm.results $(EXAMPLES)/contamination/s01-mslm.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 4: DONE"
	@echo "$(NO_COLOR)"

verify-external: mads wells
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 5: Problem $(EXAMPLES)/wells/w01 ..."
	rm -f $(EXAMPLES)/wells/w01.mads_output $(EXAMPLES)/wells/w01.results $(EXAMPLES)/wells/w01.running
	cd $(EXAMPLES)/wells; $(DBG) ../../$(MADS) w01 lmeigen $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells/w01.mads_output $(EXAMPLES)/wells/w01.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/wells/w01.results $(EXAMPLES)/wells/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 5: DONE"
	@echo "$(NO_COLOR)"

verify-external-short: mads wells
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 6: Problem $(EXAMPLES)/wells-short/w01 using different instruction formats ..."
	@echo "TEST 6.1: Instruction file $(EXAMPLES)/wells-short/w01-v1.inst ..."
	rm -f $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../$(MADS) w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.2: Instruction file $(EXAMPLES)/wells-short/w01-v2.inst ..."
	rm -f $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v2.inst w01.inst; $(DBG) ../../$(MADS) w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.3: Instruction file $(EXAMPLES)/wells-short/w01-v3.inst ..."
	rm -f $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v3.inst w01.inst; $(DBG) ../../$(MADS) w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.4: Instruction file $(EXAMPLES)/wells-short/w01-v4.inst ..."
	rm -f $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v4.inst w01.inst; $(DBG) ../../$(MADS) w01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01.results $(EXAMPLES)/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.5: Problem $(EXAMPLES)/wells-short/w01 (YAML input format) ..."
	rm -f $(EXAMPLES)/wells-short/w01_yaml.results $(EXAMPLES)/wells-short/w01_yaml.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../$(MADS) w01_yaml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01_yaml.results $(EXAMPLES)/wells-short/w01_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 6.6: Problem $(EXAMPLES)/wells-short/w01 (XML input format) ..."
	rm -f $(EXAMPLES)/wells-short/w01_xml.results $(EXAMPLES)/wells-short/w01_xml.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../$(MADS) w01_xml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01_xml.results $(EXAMPLES)/wells-short/w01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "TEST 6: DONE"
	@echo "**************************************************************************************"
	@echo "$(OK_COLOR)"

verify-external-short2: mads wells
	@echo "**************************************************************************************"
	@echo " External problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 7: Problem $(EXAMPLES)/wells-short/w01 using coupled parameters and regularization terms ..."
	@echo "TEST 7.1: Problem $(EXAMPLES)/wells-short/w01 with coupled (tied) parameters based on mathematical expressions  ..."
	rm -f $(EXAMPLES)/wells-short/w01tied.results $(EXAMPLES)/wells-short/w01tied.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../$(MADS) w01tied $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01tied.results $(EXAMPLES)/wells-short/w01tied.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 7.2: Problem $(EXAMPLES)/wells-short/w01 with regularization terms for optimized model parameters (MADS text input format) ..."
	rm -f $(EXAMPLES)/wells-short/w01regul.results $(EXAMPLES)/wells-short/w01regul.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../$(MADS) w01regul $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01regul.results $(EXAMPLES)/wells-short/w01regul.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 7.3: Problem $(EXAMPLES)/wells-short/w01 with regularization terms for optimized model parameters (YAML input format) ..."
	rm -f $(EXAMPLES)/wells-short/w01regul_yaml.results $(EXAMPLES)/wells-short/w01regul_yaml.running
	cd $(EXAMPLES)/wells-short; ln -sf w01-v1.inst w01.inst; $(DBG) ../../$(MADS) w01regul_yaml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01regul_yaml.results $(EXAMPLES)/wells-short/w01regul_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 7: DONE"
	@echo "$(NO_COLOR)"

verify-parallel: mads wells
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Parallel execution of external problems "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 8: Parallel execution of $(EXAMPLES)/wells-short/w01parallel ..."
	@echo "TEST 8.1: Initial parallel execution of $(EXAMPLES)/wells-short/w01parallel ..."
	rm -f $(EXAMPLES)/wells-short/w01parallel.results $(EXAMPLES)/wells-short/w01parallel.restart_info $(EXAMPLES)/wells-short/w01parallel.restart_*.zip $(EXAMPLES)/wells-short/w01parallel.running
	cd $(EXAMPLES)/wells-short; $(DBG) ../../$(MADS) w01parallel np=2 eval=10 restart=0 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01parallel.results $(EXAMPLES)/wells-short/w01parallel.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 8.2: Rerun using saved results from prior parallel execution of $(EXAMPLES)/wells-short/w01parallel ..."
	rm -f $(EXAMPLES)/wells/w01parallel.results $(EXAMPLES)/wells/w01parallel.running
	cd $(EXAMPLES)/wells-short; $(DBG) ../../$(MADS) w01parallel np=2 eval=10 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01parallel.results $(EXAMPLES)/wells-short/w01parallel.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 8: DONE"
	@echo "$(NO_COLOR)"

verify-forward: mads
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Forward analytical contaminant modeling "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 9: Analytical contaminant concentrations with various sources ..."
	@echo "TEST 9.1: Point source ... "
	rm -f $(EXAMPLES)/forward/a01.mads_output $(EXAMPLES)/forward/a01.results $(EXAMPLES)/forward/a01.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a01 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a01.mads_output $(EXAMPLES)/forward/a01.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a01.results $(EXAMPLES)/forward/a01.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 9.2: Rectangular source ... "
	rm -f $(EXAMPLES)/forward/a02.mads_output $(EXAMPLES)/forward/a02.results $(EXAMPLES)/forward/a02.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a02 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a02.mads_output $(EXAMPLES)/forward/a02.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a02.results $(EXAMPLES)/forward/a02.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 9.3: Planar gaussian source ... "
	rm -f $(EXAMPLES)/forward/a03.mads_output $(EXAMPLES)/forward/a03.results $(EXAMPLES)/forward/a03.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a03 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a03.mads_output $(EXAMPLES)/forward/a03.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a03.results $(EXAMPLES)/forward/a03.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 9.4: Gaussian source ... "
	rm -f $(EXAMPLES)/forward/a04.mads_output $(EXAMPLES)/forward/a04.results $(EXAMPLES)/forward/a04.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a04 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a04.mads_output $(EXAMPLES)/forward/a04.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a04.results $(EXAMPLES)/forward/a04.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.5: Box source ... "
	rm -f $(EXAMPLES)/forward/a05.mads_output $(EXAMPLES)/forward/a05.results $(EXAMPLES)/forward/a05.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a05 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a05.mads_output $(EXAMPLES)/forward/a05.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a05.results $(EXAMPLES)/forward/a05.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.6: Box source | Levy ... "
	rm -f $(EXAMPLES)/forward/a06_yaml.mads_output $(EXAMPLES)/forward/a06_yaml.results $(EXAMPLES)/forward/a06_yaml.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a06_yaml $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a06_yaml.mads_output $(EXAMPLES)/forward/a06_yaml.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a06_yaml.results $(EXAMPLES)/forward/a06_yaml.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.7: Box source | Symmetric Levy ... "
	rm -f $(EXAMPLES)/forward/a07.mads_output $(EXAMPLES)/forward/a07.results $(EXAMPLES)/forward/a07.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a07 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a07.mads_output $(EXAMPLES)/forward/a07.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a07.results $(EXAMPLES)/forward/a07.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.8: Box source | Levy, alpha=2 ... "
	rm -f $(EXAMPLES)/forward/a08.mads_output $(EXAMPLES)/forward/a08.results $(EXAMPLES)/forward/a08.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a08 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a08.mads_output $(EXAMPLES)/forward/a08.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a08.results $(EXAMPLES)/forward/a08.results-$(OS)-correct
	@echo "$(NO_COLOR)" 
	@echo "TEST 9.9: Box source | Symmetric Levy, alpha=2 ... "
	rm -f $(EXAMPLES)/forward/a09.mads_output $(EXAMPLES)/forward/a09.results $(EXAMPLES)/forward/a09.running
	cd $(EXAMPLES)/forward; $(DBG) ../../$(MADS) a09 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/forward/a09.mads_output $(EXAMPLES)/forward/a09.mads_output-$(OS)-correct
	@$(CMP) $(EXAMPLES)/forward/a09.results $(EXAMPLES)/forward/a09.results-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "TEST 9: DONE"
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"

verify-sa: mads wells
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo " Sensitivity analyses "
	@echo "**************************************************************************************"
	@echo "$(NO_COLOR)"
	@echo "TEST 10: Global and local sensitivity analyses ..."
	@echo "TEST 10.1: Sobol internal analysis ... "
	rm -f $(EXAMPLES)/sa/a01.mads_output $(EXAMPLES)/sa/a01.sobol_sens_index $(EXAMPLES)/sa/a01.sobol_sens_total $(EXAMPLES)/sa/a01.running
	cd $(EXAMPLES)/sa; $(DBG) ../../$(MADS) a01 test=111 gsens dim=8 real=10000 smp=lhs pardomain=0.5 seed=1517604820 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/sa/a01.sobol_sens_index $(EXAMPLES)/sa/a01.sobol_sens_index-$(OS)-correct
	@$(CMP) $(EXAMPLES)/sa/a01.sobol_sens_total $(EXAMPLES)/sa/a01.sobol_sens_total-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 10.2: Sobol external analysis  ..."
	rm -f $(EXAMPLES)/wells-short/w01.sobol_sens_index $(EXAMPLES)/wells-short/w01.sobol_sens_total $(EXAMPLES)/wells-short/w01.running
	cd $(EXAMPLES)/wells-short; $(DBG) ../../$(MADS) w01 gsens real=100 seed=1066732675 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01.sobol_sens_index $(EXAMPLES)/wells-short/w01.sobol_sens_index-$(OS)-correct
	@$(CMP) $(EXAMPLES)/wells-short/w01.sobol_sens_total $(EXAMPLES)/wells-short/w01.sobol_sens_total-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "TEST 10.2: Sobol external analysis  ..."
	rm -f $(EXAMPLES)/wells-short/w01.sobol_sens_index $(EXAMPLES)/wells-short/w01.sobol_sens_total $(EXAMPLES)/wells-short/w01.running $(EXAMPLES)/wells-short/w01.restart*
	cd $(EXAMPLES)/wells-short; $(DBG) ../../$(MADS) w01 gsens real=100 seed=1066732675 np=2 $(OUTPUT)
	@echo "$(ERROR_COLOR)"
	@$(CMP) $(EXAMPLES)/wells-short/w01.sobol_sens_index $(EXAMPLES)/wells-short/w01.sobol_sens_index-$(OS)-correct
	@$(CMP) $(EXAMPLES)/wells-short/w01.sobol_sens_total $(EXAMPLES)/wells-short/w01.sobol_sens_total-$(OS)-correct
	@echo "$(NO_COLOR)"
	@echo "$(OK_COLOR)"
	@echo "**************************************************************************************"
	@echo "TEST 10: DONE"
	@echo "$(NO_COLOR)"

compare-os:
	./compare-results-os Linux Darwin

clean-examples:
	find . -name "*.mads_output_*" -print0 | xargs -0 rm
	rm -f $(EXAMPLES)/*/*.ppsd_*.results $(EXAMPLES)/*/*.igpd_*.results $(EXAMPLES)/*/*igrnd-0000* $(EXAMPLES)/*/*.igrnd_*.results $(EXAMPLES)/*/*.restart_*.zip $(EXAMPLES)/*/*.restart_info $(EXAMPLES)/*/*.running $(EXAMPLES)/*/*-rerun.mads $(EXAMPLES)/*/*-error.mads
	rm -fR $(EXAMPLES)/wells-short_w01_*
	rm -fR $(EXAMPLES)/wells-short_w01parallel*
	rm -f *.mads_output* *.running *.cmdline *.cmdline_hist

astyle:
	astyle $(SOURCESTYLE)
	rm -f $(SOURCESTYLEDEL)

tar: astyle
	rm -f mads.git.tgz
	tar -cvzf mads.git.tgz `git ls-files` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'` .git

tarf: astyle
	tar -cvzf mads.tgz `git ls-files` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'`
