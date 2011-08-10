// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
//
// http://www.ees.lanl.gov/staff/monty/codes/mads
//
// LA-CC-10-055; LA-CC-11-035
//
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
