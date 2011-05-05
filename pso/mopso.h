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
