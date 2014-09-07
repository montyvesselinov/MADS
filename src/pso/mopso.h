// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
// http://www.ees.lanl.gov/staff/monty/codes/mads
//
// LA-CC-10-055; LA-CC-11-035
//
// Based on TRIBES-D developed by Maurice Clerc (see below)
//
/*
     TRIBES-D, a fully adaptive parameter-free particle swarm optimiser
		 for real heterogeneous problems
                             -------------------
    last update           : 2008-02-01
    email                 : Maurice.Clerc@WriteMe.com

 ---- Why?
 In the complete TRIBES, there is an assumption
 saying that it is possible to define an Euclidean distance between
 two positions in the search space. In particular some points are
 chosen inside hyperspheres.

 However, for a lot of real problems this is meaningless.
 For example, if one variable x1 is a cost, and another one x2 a weight,
 defining a distance by sqrt(x1*x1 + x2*x2) is completely arbitrary.
 Well, defining _any_ distance may be arbitrary.

 TRIBES-D is "TRIBES without distance".
 Or, more precisely, without multi-dimensional distance.
 However this is not an algorithm for combinatorial problems.
 We still need to be able to compute the distance between
 two points along one dimension, typically |x-y|.

 Of course, a drawback is that it is not very good on some artificial
 problems, particularly when the search space is a hypercube, and when
 the maximum number of fitness evaluations is small.
 I have added a few such problems from the CEC 2005 benchmark
 so that you can nevertheless try TRIBES-D on them.

 Conversely, when the search space is not a hypercube, but a
 hyperparallelepid, and moreover when some dimensions are discrete,
 and the other ones continuous (heterogeneous problems)
 TRIBES-D may be pretty good.

 Note that in this version, the acceptable values for a variable
 can be not only an interval, but any given list of values.

 Also, compared to previous TRIBES versions,
 it is far better for multiobjective problems.

 ---- Principles
- generate a swarm (typically one tribe one particle)
- at each time step
	- the informer of the "current" particle is the best particle
		of its tribe (the shaman)
	- if the particle is the shaman, select an informer at random amongst
			the other shamans (mono-objective case)
			or in the archive (multiobjective case).
			Note that the probability distribution to do that is
			not necessarily uniform
	- apply a strategy of move depending on the recent past

	From time to time (the delay may be different for each tribe):
		- check whether the tribe is bad or not
		- if bad, increase its size by generating a particle
		- if good and enough particles, remove the worst particle

	From time to time:
		- check whether the swarm is good or not
		- if bad, add a new tribe
		- if good, and if enough tribes, remove the worst tribe

---- About the strategies
	Depending on its recent past a particle may be good, neutral, or bad.
	So there are three strategies, one for each case.
	However, in order to add a kind of "natural selection" amongst strategies,
	a good particle keeps the same strategy (with a probability 0.5),
	no matter which one it is.

---- Tricks

- for a bad tribe there are two ways to increase diversity
	(see swarmAdapt()):
	- generate a completely new particle
	- modify the best position (the memory) of a existent one.
		However this is done for just one dimension
		(the one with the smallest discrepancy over the whole tribe)

---- About multiobjective problems
In TRIBES-D, you can ask to solve simultaneously several problems.
First approach (when the codes functions are positive):
		the algorithm just tries to find the best compromise
		I am not sure it is very useful, but who knows?
		And anyway it costs nothing.

The second approach (negative code functions) is multiobjective.
	The algorithm is looking for a Pareto front.
	It keeps up to date an archive of non dominated positions
	(a classical method).
	What is less clasical is that it keeps it over several runs
	(which can therefore be seen as manual"restarts")
	So, if the maximum number of fitness evaluations is FEmax,
	you may try different strategies. For example just one run
	with FEmax, or 10 runs with FEmax/10. Actually the second method
	is usually better, although the algorithm is already pretty good
	when launching just one run.
	See fArchive.txt for the final result.

	For multibojective problems, the algorithm adaptively uses
	two kinds of comparisons between positions:
	either classical dominance or an extension of the DWA
	(Dynamic Weighted Aggregation) method.

	It also makes use of the Crowding Distance method, but not for
	all particles, only for the best one of each tribe.
	Unfortunately these method needs a parameter, in order to select
	the "guide". I tried here to replace this user defined parameter by
	an adaptive one, depending on the number of tribes
	(see archiveCrowDistSelect()).
	This is not completely satisfying, though.

---- About randomness
The standard rand() function in C is quite bad. That is why I use
for years a better one (KISS).
Results are significantly different. Not necessarily better, for
the bad randomness of C implies a "granularity", and thanks to this
artifact, the solution may be found faster.

---- About a few useless little things
You may have noted that a few variables in some structures
are not used (for example fPrevBest).

Similarly some features are noted "EXPERIMENT", or simply "commented"
There are here just for tests.

---- How to use TRIBES-D?
Data are read from the file problem.txt

In order to define your own problem, you have to modify this file and
positionEval().

Have fun, and keep me posted!

*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define RAND_MAX_KISS ((unsigned long) 4294967295)

#define DMax 12 // Max. number of dimensions 
#define fMax 4+1 // Maximum number of objective functions Note: number of "real" objective functions +1
#define rMax 500 // Maximum number of runs
#define initOptionsNb 5 // Number of options for initialisation of a particle
#define partMax	20 //Maximum number of particles in a tribe
#define tribMax	40 //Maximum number of tribes
#define valMax 11 // Maximum number of acceptable values on a dimension
#define almostZero 0.000000000000001
#define almostInfinite 100000000
#define archiveMax 100 // Maximum number of archived positions in the case of multiobjective optimization

// Structures
struct fitness {
	int size;
	double f[fMax]; };
struct problem {
	int D; 
	float min[DMax]; 
	float max[DMax]; 
	float dx[DMax];
	int fNb; 
	int code[fMax];
	float objective[fMax];
	struct fitness errorMax;
	float evalMax;
	int repeat;
	int valSize[DMax];
	float val[DMax][valMax]; };
struct position	{
	int size;
	double x[DMax];
	struct fitness f; };
struct particle	{
	int label;
	struct position x;
	struct position xBest;
	struct position xPrev;
	int strategy;
	struct fitness fBestPrev;
	int status; };
struct tribe {
	int size; // tribe size = number of particles within the tribe
	struct particle part[partMax];
	int best;
	struct fitness fBestPrev;
	int status; }; 
struct swarm {
	int size; // NOTE: swarm size = number of tribes
	struct tribe trib[tribMax];
	struct position best;
	struct fitness fBestPrev;
	int status; 
	struct fitness fBestStag; };	
// Specific to multiobjective
struct archived { struct position x;double crowD; };	
struct distRank { double dist; int rank; };

//----------------------------------------- Subroutines
unsigned long rand_kiss();
void seed_rand_kiss(unsigned long seed);
double alea(double a, double b);
double aleaGauss(double mean, double std_dev);
int aleaInteger(int a, int b);
void archiveCrowDist();
struct position archiveCrowDistSelect(int size);
void archiveDisplay();
struct fitness archiveFitnessVar();
void archiveLocalSearch(struct problem pb);
void archiveSave(struct archived archiv[],int archiveNb,FILE *fArchive);
double archiveSpread();
int compareAdaptF(struct fitness f1[],int run,int fCompare);
static int compareCrowD (void const *a, void const *b); // For qsort
static int compareDR (void const *a,void const *b); // For qsort
static int compareFn (void const *a, void const *b); // For qsort
static int compareXY (void const *a, void const *b); // For qsort
struct position constraint(struct problem pb,struct position pos);
int fitnessCompare(struct fitness f1,struct fitness f2,int compareType);
void fitnessDisplay(struct fitness f);
double fitnessDist(struct fitness f1, struct fitness f2);
double fitnessTot(struct fitness f,int type);
double granul(double value,double granul);
double maxXY(double x, double y);
double minXY(double x,double y);
struct particle particleInit(struct problem pb, int option, struct position guide1,struct position guide2, struct swarm S);
void positionArchive(struct position pos);
void positionDisplay(struct position pos);
struct fitness positionEval(struct problem pb, struct position x);
void positionSave(FILE*fRun, struct position pos,int run);
struct position positionUpdate(struct problem pb,struct particle par,struct particle informer);
void problemDisplay(struct problem pb);																	
struct problem problemRead(FILE *fProblem);	
struct swarm problemSolve(struct problem pb,int compareType,int run);
int sign(double x);
struct swarm swarmAdapt(struct problem pb,struct swarm S,int compareType);
void swarmDisplay(struct swarm S);
struct swarm swarmInit(struct problem pb,int compareType);
struct swarm swarmLocalSearch(struct problem pb,struct swarm S);
struct swarm swarmMove(struct problem pb,struct swarm S,int compareType,int run);
int swarmTotPart(struct swarm S);
int	tribeBest(struct tribe T,int compareType);
void tribeDisplay(struct tribe T);
struct tribe tribeInit(struct problem pb, int partNb,int compareType, struct swarm S);
int tribeVarianceMin(struct tribe T);
void wFAdapt(int fNb,int run);

//----------------------------------------- Global variables
struct archived archiv[archiveMax+1];
int arch;
int archiveNb;
struct fitness archiveVar;
double epsilon_vector[fMax]; // For epsilon-dominance
float eval; // Number of fitness evaluations
struct fitness f1[rMax];
int fCompare;
int fn;
double fitMax[fMax]; // Maximum fitness value found during the process
double fitMin[fMax]; // Minimum fitness value found during the process
int iter;
int iterLocalSearchNb;
int iterSwarmAdapt; // Number of iterations between two swarm adaptations
int iterSwarmStag;
int iterTribeAdapt[tribMax]; // The same, but for each tribe
int label; // label (integer number) of the last generated particle
int multiObj; // Flag for multiobjective problem
float o[DMax]; // Offset, in particular for CEC 2005 benchmark
int overSizeSwarm;// Nb of times the swarm tends to generate too many tribes
int overSizeTribe;// Nb of times a tribe tends to generate too many particles
int restart;
int restartNb;
int verbose; // Read on the problem file. The higher it is, the more verbose is the program
double wF[fMax+1]; // Dynamic penalties
