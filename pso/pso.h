// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
//
// http://www.ees.lanl.gov/staff/monty/codes/mads
//
// LA-CC-10-055; LA-CC-11-035
//
// Copyright 2011.  Los Alamos National Security, LLC.  All rights reserved.
// This material was produced under U.S. Government contract DE-AC52-06NA25396 for
// Los Alamos National Laboratory, which is operated by Los Alamos National Security, LLC for
// the U.S. Department of Energy. The Government is granted for itself and others acting on its
// behalf a paid-up, nonexclusive, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, and perform publicly and display publicly. Beginning five (5) years after
// --------------- March 11, 2011, -------------------------------------------------------------------
// subject to additional five-year worldwide renewals, the Government is granted for itself and
// others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
// material to reproduce, prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//
// NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY, LLC,
// NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
// RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR
// PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

// Structures
struct objfunc
{
	int size;
	double *f;
};
struct position
{
	int size;
	double *x;
	struct objfunc f;
};
struct particle
{
	int id;
	int strategy;
	int status;
	struct position x;
	struct position xBest;
	struct position xLast;
	struct objfunc fBestPrev;
};
struct tribe
{
	int size; // tribe size = number of particles within the tribe
	int best; // best particle
	int status;
	struct particle *part;
	struct objfunc fBestPrev;
};
struct swarm
{
	int size; // swarm size = number of tribes
	int status;
	int tr_best;
	struct tribe *trib;
	struct position best;
	struct objfunc fBestPrev;
	struct objfunc fBestStag;
};
struct problem
{
	int init;
	int maxEval;
	int repeat;
	int D;
	float lmfactor;
	float *min;
	float *max;
	float *dx;
	int nPhi;
	int *code;
	float *objective;
	int *valSize;
	float **val;
	float *ival;
	struct objfunc maxError;
	struct position pos_success;
	int success;
};
