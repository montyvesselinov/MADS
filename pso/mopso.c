// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
//
// http://www.ees.lanl.gov/staff/monty/codes/mads
//
// LA-CC-10-055; LA-CC-11-035

#include "../mads.h"
// #include "mopso.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define RAND_MAX_KISS ((unsigned long) 4294967295)

#define TRUE 1
#define FALSE 0
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
struct fitness
{
	int size;
	double f[fMax];
};
struct problem
{
	int D;
	int init;
	float ival[DMax];
	float min[DMax];
	float max[DMax];
	float dx[DMax];
	int fNb;
	int code[fMax];
	float objective[fMax];
	struct fitness errorMax;
	int evalMax;
	int repeat;
	int valSize[DMax];
	float val[DMax][valMax];
};
struct position
{
	int size;
	double x[DMax];
	struct fitness f;
};
struct particle
{
	int label;
	struct position x;
	struct position xBest;
	struct position xPrev;
	int strategy;
	struct fitness fBestPrev;
	int status;
};
struct tribe
{
	int size; // tribe size = number of particles within the tribe
	struct particle part[partMax];
	int best;
	struct fitness fBestPrev;
	int status;
};
struct swarm
{
	int size; // NOTE: swarm size = number of tribes
	struct tribe trib[tribMax];
	struct position best;
	struct fitness fBestPrev;
	int status;
	struct fitness fBestStag;
};
// Specific to multiobjective
void problem_init( struct opt_data *op, struct problem *pr );
struct archived { struct position x; double crowD; };
struct distRank { double dist; int rank; };

//----------------------------------------- Subroutines
int mopso( struct opt_data *op );
static unsigned long rand_kiss();
static void seed_rand_kiss( unsigned long seed );
static double random_double( double a, double b );
static double random_gaussian( double mean, double std_dev );
static int random_int( int a, int b );
static void archive_crowding_dist();
static struct position archiveCrowDistSelect( int size );
static void archive_print();
static struct fitness archiveFitnessVar();
static void archive_local_search( struct problem pb );
static void archive_save( struct archived archiv[], int archiveNb, FILE *fArchive );
static double archive_spread();
static int compare_adapt_f( struct fitness f1[], int run, int fCompare );
static int compare_crowding_dist( void const *a, void const *b ); // For qsort
static int compare_dist_rank( void const *a, void const *b ); // For qsort
static int compare_fit( void const *a, void const *b ); // For qsort
static int compare_double( void const *a, void const *b ); // For qsort
static struct position particle_check( struct problem pb, struct position pos );
static int compare_particle_fit( struct fitness f1, struct fitness f2, int compare_type );
static void fitness_print( struct fitness f );
static double fitness_dist( struct fitness f1, struct fitness f2 );
static double fitness_total( struct fitness f, int type );
static double granularity( double value, double granul );
static double maxXY( double x, double y );
static double minXY( double x, double y );
static struct particle particle_init( struct problem pb, int option, struct position guide1, struct position guide2, struct swarm S );
static void position_archive( struct position pos );
static void position_print( struct position pos );
static struct fitness position_eval( struct problem pb, struct position x );
static void position_save( FILE *fRun, struct position pos, int run );
static struct position position_update( struct problem pb, struct particle par, struct particle informer );
static void problem_print( struct problem pb );
// static struct problem problemRead( FILE *fProblem );
static struct swarm pso_solver( struct problem pb, int compare_type, int run );
static int sign( double x );
static struct swarm swarm_adapt( struct problem pb, struct swarm S, int compare_type );
static void swarm_print( struct swarm S );
static struct swarm swarm_init( struct problem pb, int compare_type );
static struct swarm swarm_local_search( struct problem pb, struct swarm S );
static struct swarm swarm_move( struct problem pb, struct swarm S, int compare_type, int run );
static int swarm_particle_count( struct swarm S );
static int	tribe_best_particle( struct tribe T, int compare_type );
static void tribe_print( struct tribe T );
static struct tribe tribe_init( struct problem pb, int partNb, int compare_type, struct swarm S );
static int tribe_varmin_dimension( struct tribe T );
static void modify_weights( int fNb, int run );
double irand();
//----------------------------------------- Global variables
int lmo_count;
int lmo_flag;
struct archived archiv[archiveMax + 1];
int arch;
int archiveNb;
struct fitness archiveVar;
double epsilon_vector[fMax]; // For epsilon-dominance
int eval; // Number of functional evaluations
struct fitness f1[rMax];
int fCompare;
int fn;
double phi_max[fMax]; // Maximum fitness value found during the process
double phi_min[fMax]; // Minimum fitness value found during the process
int iter;
int iterLocalSearchNb;
int iterSwarmAdapt; // Number of iterations between two swarm adaptations
int iterSwarmStag;
int iterTribeAdapt[tribMax]; // The same, but for each tribe
int label = 0; // label (integer number) of the last generated particle
int multiObj; // Flag for multiobjective problem
float o[DMax]; // Offset, in particular for CEC 2005 benchmark
int overSizeSwarm;// Nb of times the swarm tends to generate too many tribes
int overSizeTribe;// Nb of times a tribe tends to generate too many particles
int restart;
int restartNb;
int debug_level; // Read on the problem file. The higher it is, the more debug_level is the program
double wF[fMax + 1]; // Dynamic penalties

static struct position position_lm( struct opt_data *op, struct problem pb, struct position pos );
void Transform( double *x, void *data, double *f );
void DeTransform( double *x, void *data, double *f );
int get_seed( );
int optimize_lm( struct opt_data *op );
static struct problem problemset( struct opt_data *op );
struct opt_data *gop;
double *res;

//--------Kiss variables
static unsigned int kiss_x = 1;
static unsigned int kiss_y = 2;
static unsigned int kiss_z = 4;
static unsigned int kiss_w = 8;
static unsigned int kiss_carry = 0;
static unsigned int kiss_k;
static unsigned int kiss_m;
//--------Internal random generator
int *irand_seed;

int mopso( struct opt_data *op )
{
	FILE *fRun; // To save the run
	FILE *fArchive; // To save the Pareto front
	struct position bestBest;
	struct problem pb;
	struct swarm S;
	double successRate[fMax], errorMean[fMax], errorTot;
	char filename[80];
	int debug, i, n, r, eval_total;
	gop = op;
	irand_seed = &op->cd->seed;
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); }
	else { seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); if( op->cd->pdebug ) tprintf( "Current seed: %d\n", op->cd->seed ); }
	overSizeSwarm = 0;
	overSizeTribe = 0;
	pb.fNb = 1;
	if( ( res = ( double * ) malloc( op->od->nObs * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); exit( 1 ); }
	if( op->cd->pdebug )
	{
		sprintf( filename, "%s.runs", op->root );
		fRun = fopen( filename, "w" );
		sprintf( filename, "%s.archive", op->root );
		fArchive = fopen( filename, "w" );
	}
	pb = problemset( op );
	problem_print( pb );
	errorTot = 0;
	eval_total = 0;
	archiveNb = 0;
	fCompare = 0; // Kind of comparison, to begin (see fitnessCompare() )
	bestBest.size = pb.fNb; //Prepare final result
	for( n = 0; n < pb.fNb; n++ )
	{
		bestBest.f.f[n] = almostInfinite;
		successRate[n] = 0;
		epsilon_vector[n] = 0.0;
	}
	if( multiObj )
	{
		// Prepare epsilon-dominance
		// as this setting is outside of the loop on runs
		// each run will take advantage of the previous ones
		for( n = 0; n < pb.fNb; n++ ) epsilon_vector[n] = 1e30; // to set epsilon-dominance, set epsilon[n] to huge number
		// for(n=pb.fNb;n<pb.fNb;n++) epsilon_vector[n]=almostZero; // to eliminate epsilon-dominance, set epsilon[n] to zero
	}
	else
		for( n = 0; n < pb.fNb; n++ ) epsilon_vector[n] = 0.0;
	for( n = 0; n < pb.fNb; n++ )
	{
		phi_min[n] = 1e30; // Very big value to begin
		phi_max[n] = 0; // Very low value to begin
	}
	lmo_count = restart = restartNb = 0;
	tprintf( "Particle-Swarm Optimization started ...\n" );
	lmo_flag = 0;
	if( strstr( op->cd->opt_method, "lm" ) != NULL ) lmo_flag = 1;
	for( r = 0; r < pb.repeat; r++ )
	{
		if( restart == 0 )
		{
			eval = 0; // Init nb of fitness evaluations
			if( pb.repeat > 1 ) tprintf( "\nRUN %i fCompare %i\n", r + 1, fCompare );
		}
		else restart = 0;
		S = pso_solver( pb, fCompare, r );
		if( lmo_flag )
		{
			fitness_print( S.best.f );
			tprintf( "COUPLED CALIBRATION: Particle-Swarm Optimization followed by Levenberg-Marquardt Optimization\n" );
			S.best = position_lm( gop, pb, S.best );
			lmo_count++;
			tprintf( "phi %g\n", gop->phi );
			fitness_print( S.best.f );
		}
		f1[r] = S.best.f;
		if( r < pb.repeat - 1 )
			fCompare = compare_adapt_f( f1, r, fCompare ); // Define the next comparison method
		if( restart == 1 )
		{
			restartNb = restartNb + 1;
			tprintf( "Restart at eval %d\n", eval );
			r--;
		}
		else
		{
			if( restartNb > 0 ) tprintf( "Restarts = %i\n", restartNb );
			tprintf( "Final " ); swarm_print( S );
			tprintf( "Final result after %i time steps and %d evaluations:\n", iter, eval );
			if( debug_level && pb.repeat > 1 )
			{
				debug = op->cd->debug; op->cd->debug = 3;
				func_global( S.best.x, gop, res );
				op->cd->debug = debug;
			}
			fitness_print( S.best.f );
			if( debug_level ) position_print( S.best );
			if( op->cd->pdebug ) position_save( fRun, S.best, r ); // Save the result
			for( n = 0; n < pb.fNb; n++ )
			{
				errorMean[n] += S.best.f.f[n];
				if( S.best.f.f[n] < pb.errorMax.f[n] ) successRate[n]++;
			}
			if( compare_particle_fit( S.best.f, bestBest.f, 2 ) == 1 ) bestBest = S.best; // Save the best solution
			eval_total += eval;
		}
	}
	tprintf( "Particle-Swarm Optimization completed successfully!\n" );
	if( pb.repeat > 1 )
		for( n = 0; n < pb.fNb; n++ )
		{
			errorMean[n] /= pb.repeat;
			successRate[n] /= pb.repeat;
		}
	if( debug_level && pb.repeat > 1 )
	{
		tprintf( "\nMean errors:" );
		for( n = 0; n < pb.fNb; n++ ) tprintf( " %g", errorMean[n] );
		tprintf( "\nSuccess rate: " );
		for( n = 0; n < pb.fNb; n++ ) tprintf( " %f", successRate[n] );
		tprintf( "\nBest result:\n" ); // position_print(bestBest);
	}
	debug = op->cd->debug; op->cd->debug = 3;
	func_global( bestBest.x, gop, res );
	op->cd->debug = debug;
	fitness_print( bestBest.f );
	tprintf( "Total number of functional evaluations %d\n", eval_total );
	tprintf( "Total number of Leverberg-Marquardt optimizations %d\n", lmo_count );
	if( pb.repeat > 1 ) tprintf( "Number of Particle-Swarm Optimizations %d\n", pb.repeat );
	if( pb.repeat > 1 ) tprintf( "Average number of evaluations for each Particle-Swarm Optimization %g\n", ( ( double ) eval_total / pb.repeat ) );
	if( overSizeTribe > 0 ) tprintf( "WARNING: Tribe size has been constrained %i times\n", overSizeTribe );
	if( overSizeSwarm > 0 ) tprintf( "WARNING: Swarm size has been constrained %i times\n", overSizeSwarm );
	if( multiObj ) archive_save( archiv, archiveNb, fArchive );
	DeTransform( bestBest.x, op, bestBest.x );
	for( i = 0; i < op->pd->nOptParam; i++ )
		op->pd->var[op->pd->var_index[i]] = bestBest.x[i];
	op->cd->standalone = 1;
	if( op->cd->pdebug ) fclose( fRun );
	if( op->cd->pdebug ) fclose( fArchive );
	free( res );
	return EXIT_SUCCESS;
}

static struct fitness position_eval( struct problem pb, struct position x )
{
	struct fitness fit;
	double f;
	int i;

	eval++;
	fit.size = pb.fNb;
	func_global( x.x, gop, res ); // evaluation ... either internal of external
	f = 0;
	for( i = 0; i < gop->od->nObs; i++ ) f += res[i] *res[i];
	fit.f[0] = fabs( f - pb.objective[0] );
	for( i = 0; i < pb.fNb; i++ ) // Save the min and the max fitness ever found
	{
		phi_min[i] = minXY( phi_min[i], fit.f[i] );
		phi_max[i] = maxXY( phi_max[i], fit.f[i] );
	}
	if( multiObj && archiveNb > 0 ) // Additional "fitness" for multiobjective
	{
		archiv[archiveNb].x.f = fit; // Temporary put in the archive
		archiveNb++;
		fit.f[pb.fNb] = archive_spread();
		archiveNb--; // Virtually removed from the archive
	}
	if( debug_level > 2 )
	{
		tprintf( "position evaluation:" );
		fitness_print( fit );
		if( debug_level > 3 ) position_print( x );
	}
	return fit;
}

static struct problem problemset( struct opt_data *op )
{
	int d;
	struct problem pb;
	double *opt_var, *tr_var;

	pb.D = op->pd->nOptParam;
	if( pb.D > DMax )
	{
		fprintf( stderr, "\nToo high dimension %i", pb.D );
		fprintf( stderr, "\nYou may increase DMax (%i)\n", DMax );
		exit( 1 );
	}
	if( ( opt_var = ( double * ) malloc( pb.D *sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( tr_var = ( double * ) malloc( pb.D *sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	for( d = 0; d < pb.D; d++ )
		opt_var[d] = op->pd->var[ op->pd->var_index[d] ];
	Transform( opt_var, op, tr_var );
	pb.init = 1;
	if( op->cd->sintrans == 0 )
		for( d = 0; d < pb.D; d++ ) // typically used for test problems
		{
			pb.ival[d] = ( float ) tr_var[d];
			pb.min[d] = op->pd->var_min[op->pd->var_index[d]];
			pb.max[d] = op->pd->var_max[op->pd->var_index[d]];
			pb.dx[d] = op->pd->var_dx[op->pd->var_index[d]]; // defines discretization
			pb.valSize[d] = 0;
		}
	else
		for( d = 0; d < pb.D; d++ ) // typically used for actual problems
		{
			pb.ival[d] = ( float ) tr_var[d];
			pb.min[d] = -M_PI / 2;
			pb.max[d] = M_PI / 2;
			pb.dx[d] = 0; // defines discretization
			pb.valSize[d] = 0;
		}
	pb.fNb = 1;
	if( pb.fNb > fMax - 1 )
	{
		fprintf( stderr, "Too many functions: %i (max %i) \n", pb.fNb, fMax );
		exit( 1 );
	}
	pb.code[0] = -1;
	pb.objective[0] = 0;
	pb.errorMax.f[0] = op->cd->phi_cutoff;
	pb.errorMax.size = pb.fNb;
	pb.evalMax = op->cd->maxeval;
	if( op->cd->nretries == 0 ) pb.repeat = 1;
	else pb.repeat = op->cd->nretries;
	op->cd->standalone = 0;

	if( pb.repeat > rMax )
	{
		fprintf( stderr, "Too many runs (max %i) \n", rMax );
		exit( 1 );
	}
	debug_level = op->cd->pdebug;
	multiObj = FALSE;
	return pb;
}

static double random_double( double a, double b )
{
	double 	r;
	// Normally, RAND_MAX = 32767 = 2^15-1
	// Not very good. Replaced it by KISS
	// r = ( double )rand_kiss() / RAND_MAX_KISS;
	// r = ( double ) rand() / RAND_MAX;
	r = irand();
	r = a + r * ( b - a );
	return r; /* Random real  number between a and b */
}

static double random_gaussian( double mean, double stddev )
{
	double v1, v2, w;
	do
	{
		v1 = 2.0 * random_double( 0, 1 ) - 1; // polar form of the Box-Muller transformation to obtain
		v2 = 2.0 * random_double( 0, 1 ) - 1; // a pseudo random number from a Gaussian distribution
		w = v1 * v1 + v2 * v2;
	}
	while( w >= 1.0 );
	return( v1 * sqrt( -2.0 * log( w ) / w ) * stddev + mean );
}

static int random_int( int min, int max )
{
	double 	r;
	if( min >= max ) return min;
	r = random_double( min, max );
	return ( int )floor( r + 0.5 ); /* Random integer number between min and max */
}

static struct position archiveCrowDistSelect( int size )
{
	// Choose a random position, according to a non uniform distribution
	// of the crowding distance
	// the bigger the CD, the higher the probability
	// Also the bigger the number of tribes (i.e. "size")
	// the "sharper" the distribution

	int n;
	double pr;
	archive_crowding_dist(); // Compute the Crowding Distances
	qsort( archiv, archiveNb, sizeof( archiv[0] ), compare_crowding_dist ); // Sort the archive by increasing Crowding Distance
	pr = 2 * log( 1 + size ); // Choose at random according to a non uniform distribution:
	n = pow( random_double( 0, pow( archiveNb - 1, pr ) ), 1. / pr );
	return archiv[n].x;
}

static void archive_crowding_dist() // Compute the crowding distances in archive[n] (global variable)
{
	double dist, max, min;
	int fNb, n;
	fNb = archiv[0].x.f.size;
	for( n = 0; n < archiveNb; n++ ) archiv[n].crowD = 1; // Initialise
	for( fn = 0; fn < fNb; fn++ )
	{
		qsort( archiv, archiveNb, sizeof( archiv[0] ), compare_fit ); // Sort archive according to f[fn]  NOTE: fn is a global variable
		for( n = 0; n < archiveNb - 1; n++ ) // For each position find the distance to the nearest one
		{
			if( n == 0 ) min = phi_min[fn]; else min = archiv[n - 1].x.f.f[fn];
			if( n == archiveNb - 2 ) max = phi_max[fn]; else max = archiv[n + 1].x.f.f[fn];
			dist = max - min;
			if( dist < epsilon_vector[fn] ) epsilon_vector[fn] = dist; // Prepare epsilon-dominance
			archiv[n].crowD = archiv[n].crowD * dist; // Contribution to the hypervolume
		}
	}
}

static void archive_print()
{
	int n;
	tprintf( "Archive %i positions\n", archiveNb );
	for( n = 0; n < archiveNb; n++ )
	{
		position_print( archiv[n].x );
		tprintf( " crowd %g", archiv[n].crowD );
	}
	tprintf( "\n" );
}

static struct fitness archiveFitnessVar()
{
	// Variance of the archive for each fitness
	// More precisely return for each fitness g(var)
	// where g is a decreasing funtion of var:
	// The smaller the better
	double mean, var, z;
	int m, n;
	struct fitness varF;

	varF.size = archiv[0].x.f.size;

	for( n = 0; n < varF.size; n++ )
	{
		mean = 0;
		for( m = 0; m < archiveNb; m++ )
			mean += archiv[m].x.f.f[n];
		mean /= archiveNb;
		var = 0;
		for( m = 0; m < archiveNb; m++ )
		{
			z = archiv[m].x.f.f[n] - mean;
			var += z * z;
		}
		varF.f[n] = var / archiveNb;
	}
	for( n = 0; n < varF.size; n++ )
		varF.f[n] = 1. / ( 1 + varF.f[n] ); // Transform into something to minimise
	return varF;
}

static void archive_local_search( struct problem pb )
{
	/* Work with the global variables
	 archiv[]
	 archiveNb
	 arch

	 Let fNb be the dimension of the fitness space.
	 Loop:
	 Define a simplex of dimension fNb
	 Find a point "inside". Check if it is a good one
	 or
	 Find a point "outside". Check if it is a good one
	 */
	struct distRank dR[archiveMax - 1];
	struct position simplex[fMax + 1];
	struct position xIn, xOut;
	int fNb, m, n, r, d, out;
	fNb = archiv[0].x.f.size; // Dimension of the fitness space
	if( archiveNb < fNb + 1 ) return;
	if( archiveNb < archiveMax ) return;
	if( iterLocalSearchNb < archiveNb ) return;
	tprintf( "Iter %i Eval %d: Local search from archive\n", iter, eval );
	iterLocalSearchNb = 0;
	xIn.size = archiv[0].x.size;
	xOut.size = xIn.size;
	out = random_int( 0, 1 );
	for( n = 0; n < archiveNb - fNb; n++ )
	{
		// Define a simplex
		simplex[0] = archiv[n].x; // First element
		m = 0;
		for( r = 0; r < archiveNb; r++ )
		{
			if( r == n ) continue;
			dR[m].dist = fitness_dist( archiv[n].x.f, archiv[m].x.f ); // Compute the distances to the others
			dR[m].rank = r;
			m = m + 1;
		}
		// Find the fNb nearest ones in the archive in order to complete the simplex
		fNb = archiv[0].x.f.size;
		qsort( dR, archiveNb - 1, sizeof( dR[0] ), compare_dist_rank );
		for( m = 0; m < fNb; m++ )
			simplex[m + 1] = archiv[dR[m].rank].x;
		// Define a new point
		//out=aleaInteger(0,1); // TO TRY
		//out=1;
		//		tprintf("\n out %i",out);
		switch( out )
		{
			case 0:	// Inside the simplex
				for( d = 0; d < xIn.size; d++ )
				{
					xIn.x[d] = 0;
					for( m = 0; m < fNb + 1; m++ )
						xIn.x[d] += simplex[m].x[d];
				}
				for( d = 0; d < xIn.size; d++ )
					xIn.x[d] /= ( fNb + 1 );
				xIn = particle_check( pb, xIn );
				xIn.f = position_eval( pb, xIn );
				position_archive( xIn );
				break;
			case 1: // Outside the simplex
				for( d = 0; d < xOut.size; d++ )
				{
					xOut.x[d] = 0;
					for( m = 1; m < fNb + 1; m++ ) // Partial simplex, dimension D-1
						xOut.x[d] += simplex[m].x[d];
				}
				for( d = 0; d < xOut.size; d++ )
					xOut.x[d] /= fNb; // Gravity center
				for( d = 0; d < xOut.size; d++ )
					xOut.x[d] -= ( simplex[0].x[d] - xOut.x[d] ); // Reflection of the first vertex
				xOut = particle_check( pb, xOut );
				xOut.f = position_eval( pb, xOut );
				position_archive( xOut );
				break;
		}
		if( arch == 1 ) // If "good" position ..
		{
			// ... there is nevertheless a small probability
			// that the "curvature" is not the same for the next point.
			// and that it is interesting to try the other case
			if( random_double( 0, 1 ) < 1 / ( fNb + 1 ) ) out = 1 - out;
		}
		if( arch == 0 ) // If "bad" position ...
		{
			// ... there is a high probability that
			// the "curvature" is the same for the next point
			// and that is interesting to try the other case
			if( random_double( 0, 1 ) < fNb / ( fNb + 1 ) ) out = 1 - out;
		}
	}
}

static void archive_save( struct archived archiv[], int archiveNb, FILE *fArchive )
{
	int D, fNb, m, n;
	tprintf( "Save archive (%i positions)\n", archiveNb );
	D = archiv[0].x.size;
	fNb = archiv[0].x.f.size;
	tprintf( "%i fitnesses, and %i coordinates\n", fNb, D );
	for( n = 0; n < fNb; n++ ) fprintf( fArchive, "f%i ", n + 1 );
	for( n = 0; n < D; n++ ) fprintf( fArchive, "x%i ", n + 1 );
	fprintf( fArchive, "\n" );
	for( n = 0; n < archiveNb; n++ )
	{
		for( m = 0; m < fNb; m++ )
			fprintf( fArchive, "%g ", archiv[n].x.f.f[m] ); // Fitnesses
		for( m = 0; m < D; m++ )
			fprintf( fArchive, "%g ", archiv[n].x.x[m] ); // Position
		fprintf( fArchive, "\n" );
	}
}

static double archive_spread()
{
	double dMin[archiveMax], diversity1, diversity2, dMean1, dMean2, noSpread, noSpread2, z, d;
	int m, n;
	if( archiveNb < 2 ) return almostInfinite;
	// For each position, find the  nearest one
	// and the farest one, in the fitness space
	// Without the new position
	diversity1 = 0;
	for( n = 0; n < archiveNb - 1; n++ )
	{
		dMin[n] = almostInfinite;
		for( m = 0; m < archiveNb - 1; m++ )
		{
			if( m == n ) continue;
			d = fitness_dist( archiv[n].x.f, archiv[m].x.f );
			if( d < dMin[n] ) dMin[n] = d;
		}
		diversity1 += dMin[n];
	}
	dMean1 = diversity1 / ( archiveNb - 1 );
	// With the new position
	diversity2 = 0;
	for( n = 0; n < archiveNb; n++ )
	{
		dMin[n] = almostInfinite;
		for( m = 0; m < archiveNb; m++ )
		{
			if( m == n ) continue;
			d = fitness_dist( archiv[n].x.f, archiv[m].x.f );
			if( d < dMin[n] ) dMin[n] = d;
		}
		diversity2 += dMin[n];
	}
	dMean2 = diversity2 / archiveNb;
	noSpread2 = 0; // "noSpread" estimation
	for( n = 0; n < archiveNb; n++ )
	{
		z = dMean2 - dMin[n];
		noSpread2 += z * z;
	}
	noSpread = sqrt( noSpread2 ); // Initial value
	// Take distances (in the fitness space) into account
	// Distance between the new position and the first archived
	d = fitness_dist( archiv[0].x.f, archiv[archiveNb - 1].x.f );
	for( n = 1; n < archiveNb - 1; n++ ) // Distance min to the others
	{
		z = fitness_dist( archiv[n].x.f, archiv[archiveNb - 1].x.f );
		if( z < d ) d = z;
	}
	if( d < dMean1 ) noSpread += dMean1; // Penalty
	return noSpread;
}

static int compare_adapt_f( struct fitness f1[], int runs, int fCompare )
{
	/*
	 Modify the type of comparison, if multiobjective
	 Possible methods:
	 - adaptive weights
	 - epsilon-dominance
	 */
	double mean[fMax], var[fMax];
	double dfMax, dfMin;
	double z;
	int fSame[fMax];
	int fNb;
	int compType;
	int n;
	int same;
	int r;
	if( runs == 0 ) return fCompare; // Not possible for the first run
	fNb = f1[0].size;
	if( fNb == 1 ) return fCompare; // Nothing to do for mono-objective
	// Is the diversity big enough?
	for( n = 0; n < fNb; n++ )
	{
		mean[n] = 0;
		for( r = 0; r <= runs; r++ )
			mean[n] += f1[r].f[n];
		mean[n] = mean[n] / ( runs + 1 );
	}
	for( n = 0; n < fNb; n++ )
	{
		var[n] = 0;
		for( r = 0; r <= runs; r++ )
			var[n] = pow( f1[r].f[n] - mean[n], 2 );
		var[n] = sqrt( var[n] / ( runs + 1 ) ); // Standard deviation
	}
	for( n = 0; n < fNb; n++ )
	{
		dfMin = mean[n] - 3 * var[n];
		dfMax = mean[n] + 3 * var[n];
		if( dfMin < phi_min[n] ) {fSame[n] = 0; continue;}
		if( dfMax > phi_max[n] ) {fSame[n] = 0; continue;}
		// Fuzzy rule
		fSame[n] = 0;
		z = random_double( phi_min[n], mean[n] );
		if( z < dfMin ) fSame[n] = 1;
		z = random_double( mean[n], phi_max[n] );
		if( z > dfMax ) fSame[n] = 1;
	}
	same = 0;
	for( n = 0; n < fNb; n++ )
		if( fSame[n] == 1 ) same++;
	// Either 0 or 1
	if( same == fNb ) compType = 1 - fCompare; else compType = fCompare;
	return compType;
}

static int compare_crowding_dist( void const *a, void const *b ) // QSORT: compare crowding distances in the archive
{
	struct archived const *pa = ( struct archived const * ) a;
	struct archived const *pb = ( struct archived const * ) b;
	if( pa->crowD > pb->crowD ) return 1;
	if( pa->crowD < pb->crowD ) return -1;
	return 0;
}

static int compare_dist_rank( void const *a, void const *b ) // QSORT: compare distancences
{
	struct distRank const *pa = ( struct distRank const * ) a;
	struct distRank const *pb = ( struct distRank const * ) b;
	if( pa->dist > pb->dist ) return 1;
	if( pa->dist < pb->dist ) return -1;
	return 0;
}

static int compare_fit( void const *a, void const *b ) // QSORT: compare f[fn] in archive
{
	struct archived const *pa = ( struct archived const * ) a;
	struct archived const *pb = ( struct archived const * ) b;
	if( pa->x.f.f[fn] > pb->x.f.f[fn] ) return 1; // NOTE: needs the global variable fn
	if( pa->x.f.f[fn] < pb->x.f.f[fn] ) return -1;
	return 0;
}

static int compare_double( void const *a, void const *b ) // QSORT: compare double numbers
{
	double const *pa = ( double const * ) a;
	double const *pb = ( double const * ) b;
	if( *pa > *pb ) return 1;
	if( *pa < *pb ) return -1;
	return 0;
}

static struct position particle_check( struct problem pb, struct position pos )
{
	struct position posG;
	int rank, d, n;

	posG = pos;
	// Granularity
	for( d = 0; d < pb.D; d++ )
		if( pb.valSize[d] == 0 )
			if( pb.dx[d] > 0 )
				posG.x[d] = granularity( pos.x[d], pb.dx[d] );

	// Clamping
	for( d = 0; d < pb.D; d++ )
	{
		if( pb.valSize[d] == 0 ) // Variable Interval
		{
			if( posG.x[d] > pb.max[d] ) posG.x[d] = pb.max[d];
			else if( posG.x[d] < pb.min[d] ) posG.x[d] = pb.min[d];
		}
		else // Variable List: Find the nearest one
		{
			rank = 0;
			for( n = 1; n < pb.valSize[d]; n++ )
				if( fabs( pb.val[d][n] - posG.x[d] ) < fabs( pb.val[d][rank] - posG.x[d] ) ) rank = n;
			posG.x[d] = pb.val[d][rank];
		}
	}
	return posG;
}

static int compare_particle_fit( struct fitness f1, struct fitness f2, int compare_type )
{
	// 1 => f1 better than f2
	// 0 equivalent
	// -1 worse
	int better = 0;
	int n;
	double wF1 = 0;
	double wF2 = 0;
	int worse = 0;
	switch( compare_type )
	{
		case 1: // Using penalties
			for( n = 0; n < f1.size; n++ )
			{
				wF1 += wF[n] * f1.f[n];
				wF2 += wF[n] * f2.f[n];
			}
			if( wF1 < wF2 ) return 1;
			if( wF2 < wF1 ) return -1;
			return 0;
		default: //0  Epsilon-dominance
			for( n = 0; n < f1.size; n++ )
			{
				if( f1.f[n] < f2.f[n] + epsilon_vector[n] ) {better++; continue;}
				if( f1.f[n] > f2.f[n] - epsilon_vector[n] ) {worse++;	continue;}
			}
			if( better > 0 && worse == 0 ) return 1;
			if( better == 0 && worse > 0 ) return -1;
			if( !multiObj ) return 0;
			//Multiobjective difficult to decide. Use the "spread" criterion
			n = f1.size; // Contain the "noSpread" value that should be as small as possible
			if( f1.f[n] < f2.f[n] ) return 1;
			if( f1.f[n] > f2.f[n] ) return -1;
			return 0;
		case 2: // Pure dominance, for multiobjective
			// Only to decide if a position has to be archived or not
			for( n = 0; n < f1.size; n++ )
			{
				if( f1.f[n] < f2.f[n] ) {better++; continue;}
				if( f1.f[n] > f2.f[n] ) {worse++; continue;}
			}
			if( better > 0 && worse == 0 ) return 1;
			if( better == 0 && worse > 0 ) return -1;
			return 0;
	}
}

static void fitness_print( struct fitness f )
{
	int i;
	if( f.size == 1 )
	{
		tprintf( "Objective function: %g\n", f.f[0] );
		return;
	}
	tprintf( "Objective functions (%d): ", f.size );
	for( i = 0; i < f.size; i++ )
		tprintf( " %g", f.f[i] );
	if( multiObj ) tprintf( " noSpread %g\n", f.f[f.size] );
	else tprintf( "\n" );
}

static double fitness_dist( struct fitness f1, struct fitness f2 )
{
	double t, dist = 0;
	int i;
	for( i = 0; i < f1.size; i++ )
	{
		t = f1.f[i] - f2.f[i];
		dist += t * t;
	}
	return sqrt( dist );
}

static double fitness_total( struct fitness fit, int type ) // Total fitness (weighted if multiobjective)
{
	double error = 0;
	int i;
	if( fit.size <= 1 ) return fit.f[0];
	switch( type )
	{
		default:
		case 0:
			for( i = 0; i < fit.size; i++ )
				error += fit.f[i] * fit.f[i];
			error = sqrt( error );
			break;
		case 1: // Weighted (wF is a global variable)
			for( i = 0; i < fit.size; i++ )
				error += wF[i] * fabs( fit.f[i] );
			break;
	}
	return error;
}

static double granularity( double value, double granul ) // Modify a value according to a granularity
{
	double x;
	if( granul < almostZero ) return value;// Pseudo-continuous (depending on the machine)
	if( value >= 0 ) x = granul * floor( value / granul + 0.5 );
	else x = ( -granul * floor( -value / granul + 0.5 ) );
	return x;
}

static double maxXY( double x, double y )
{
	if( x > y ) return x; else return y;
}

static double minXY( double x, double y )
{
	if( x < y ) return x; else return y;
}

static struct particle particle_init( struct problem pb, int initOption, struct position guide1, struct position guide2, struct swarm S )
{
	struct particle part = {0};
	double mean, range, sort_vec[2 * tribMax *partMax];
	int i, ip, it, k, count, rank;

	part.x.size = pb.D; // Initial number of particles equal number of dimensions (optimized variables)
	switch( initOption )
	{
		default:
		case 0: // Random position
			for( k = 0; k < pb.D; k++ )
				part.x.x[k] = random_double( pb.min[k], pb.max[k] );
			break;
		case 1: // Guided
			for( k = 0; k < pb.D; k++ )
			{
				if( pb.valSize[k] == 0 ) // variable defined over interval
				{
					range = fabs( guide1.x[k] - guide2.x[k] );
					part.x.x[k] = random_double( guide1.x[k] - range, guide1.x[k] + range ); // TO TRY: aleaGauss instead)
				}
				else // List
				{
					i = random_int( 0, pb.valSize[k] - 1 ); // WARNING: still purely random, not very satisfying
					part.x.x[k] = pb.val[k][i];
				}
			}
			break;
		case 2: // On the corners
			for( k = 0; k < pb.D; k++ )
			{
				if( random_double( 0, 1 ) < 0.5 ) part.x.x[k] = pb.min[k];
				else part.x.x[k] = pb.max[k];
			}
			break;
		case 3: // Biggest empty hyperparallelepid (Binary Search)
			// TO TRY: for multiobjective, one may use the archive instead of S as the list of known positions
			for( k = 0; k < pb.D; k++ )
			{
				sort_vec[0] = pb.min[k];
				sort_vec[1] = pb.max[k];
				count = 2;
				for( it = 0; it < S.size; it++ )
					if( S.trib[it].size > 0 )
						for( ip = 0; ip < S.trib[it].size; ip++ )
						{
							sort_vec[count] = S.trib[it].part[ip].x.x[k]; // List of known coordinates
							sort_vec[count + 1] = S.trib[it].part[ip].xBest.x[k];
							count += 2;
						}
				if( count > 1 ) qsort( sort_vec, count, sizeof( double ), compare_double ); // Sort the list of known coordinates
				rank = 1;
				for( i = 2; i < count; i++ )
					if( sort_vec[i] - sort_vec[i - 1] > sort_vec[rank] - sort_vec[rank - 1] ) rank = i; // Find the biggest empty interval
				part.x.x[k] = random_double( sort_vec[rank - 1], sort_vec[rank] ); // Select a random position "centered" on the interval
			}
			break;
		case 4: // Centered EXPERIMENT
			for( k = 0; k < pb.D; k++ )
			{
				mean = random_double( pb.min[k], pb.max[k] );
				range = ( 1. / 3 ) * maxXY( pb.max[k] - mean, mean - pb.min[k] );
				part.x.x[k] = random_gaussian( mean, range );
			}
			break;
		case 5: // User supplied initial values
			for( k = 0; k < pb.D; k ++ )
				part.x.x[k] = pb.ival[k];
			break;
	}
	part.x = particle_check( pb, part.x ); // Take particle_checks into account
	part.x.f = position_eval( pb, part.x ); // Evaluate the position
	if( multiObj ) position_archive( part.x ); // Archive (multiobjective)
	part.xPrev = part.xBest = part.x;
	part.fBestPrev = part.x.f;
	part.label = ++label;
	part.strategy = 0; //TO TRY: aleaInteger(0,2);
	return part;
}

static void position_archive( struct position pos )
{
	// Note: global variables
	// archiv[], archiveNb
	int archiveStore;
	int i, j, cmp;
	int dominanceType = 2; // 2 => pure-dominance
	arch = 1; // Global variable
	if( archiveNb > 0 ) // It is not the first one. Check if dominated
		for( i = 0; i < archiveNb; i++ )
		{
			cmp = compare_particle_fit( archiv[i].x.f, pos.f, dominanceType );
			if( cmp == 1 )
			{
				arch = 0; // Dominated, don't keep it
				break;
			}
		}
	if( arch == 1 )
	{
		iterLocalSearchNb++; // Useful to decide
		// when to perform local search
		// Remove the dominated positions
		if( archiveNb > 1 )
			for( i = 0; i < archiveNb; i++ )
			{
				cmp = compare_particle_fit( pos.f, archiv[i].x.f, dominanceType );
				if( cmp == 1 )
				{
					if( i < archiveNb - 1 )
						for( j = i; j < archiveNb - 1; j++ )
							archiv[j] = archiv[j + 1];
					archiveNb = archiveNb - 1;
				}
			}
		// TO TRY: one may also remove one of two "too similar" positions
		if( archiveNb < archiveMax ) // Store the position
			archiveStore = archiveNb;
		else
		{
			archiveNb = archiveMax;
			archiveStore = 0;
			for( i = 1; i < archiveNb; i++ ) // Find the most "crowded" archived position
				if( archiv[i].crowD < archiv[archiveStore].crowD ) archiveStore = i;
		}
		archiv[archiveStore].x = pos;
		if( archiveNb < archiveMax ) archiveNb++;
	}
}

static void position_print( struct position pos )
{
	int d;
	tprintf( "position:" );
	for( d = 0; d < pos.size; d++ )
		tprintf( " %g", pos.x[d] );
	tprintf( "\nfitness:" );
	for( d = 0; d < pos.f.size; d++ )
		tprintf( " %g", pos.f.f[d] );
	if( multiObj ) tprintf( " noSpread %g\n", pos.f.f[pos.f.size] );
	else tprintf( "\n" );
}

static void position_save( FILE *fRun, struct position pos, int run )
{
	int d;
	// Run
	fprintf( fRun, "%i %i %d ", run, iter, eval );
	// Fitness
	for( d = 0; d < pos.f.size; d++ )
		fprintf( fRun, "%g ", pos.f.f[d] );
	// Position
	for( d = 0; d < pos.size; d++ )
		fprintf( fRun, "%g ", pos.x[d] );
	fprintf( fRun, "\n" );
}

static struct position position_update( struct problem pb, struct particle par, struct particle informer )
{
	double c1, c2;
	double c3 = 0;
	int d;
	int Dim;
	double error1 = 0;
	double error2 = 0;
	double noise;
	struct position pos, pos1, pos2;
	double dx;
	double radius;
	int strat[3];
	int type;

	double w1 = 0.74; // => gaussian fifty-fifty
	double w2 = 0.72; // < 1/(2*ln(2)) Just for tests
	double c = 1.19; // < 0.5 + ln(2) Just for tests

	// Define the three strategies
	strat[0] = 0; // For good particle
	strat[1] = 1; // For neutral particle
	strat[2] = 2; // For bad particle

	pos.size = pb.D;
	Dim = informer.xBest.size; // May be zero if there is no informer
	switch( par.strategy )
	{
		case -1: // Random
			for( d = 0; d < pb.D; d++ )
				pos.x[d] = random_double( pb.min[d], pb.max[d] ); // Try something within the range
			par.strategy = 0;
			break;
		default: // Improved Bare bones
		case 0:
			for( d = 0; d < Dim; d++ )
			{
				dx = w1 * fabs( informer.xBest.x[d] - par.xBest.x[d] );
				pos.x[d] = random_gaussian( informer.xBest.x[d], dx );
			}
			break;
		case 1: // Pivot by dimension
		case 2: // Pivot by dimension + noise
			if( multiObj ) type = 1; else type = 0;
			if( Dim > 0 ) // If there is an informer
			{
				error1 = fitness_total( par.xBest.f, type );
				error2 = fitness_total( informer.xBest.f, type );
				c3 = error1 + error2;
				if( c3 > almostZero ) c1 = error2 / c3;
				else  c1 = random_double( 0, 1 );
				c2 = 1 - c1;
				for( d = 0; d < Dim; d++ )
				{
					radius = fabs( par.xBest.x[d] - informer.xBest.x[d] );
					pos1.x[d] = random_double( par.xBest.x[d] - radius, par.xBest.x[d] + radius );
					pos2.x[d] = random_double( informer.xBest.x[d] - radius, informer.xBest.x[d] + radius );
					pos.x[d] = c1 * pos1.x[d] + c2 * pos2.x[d];
				}
			}
			else // No other informer than itself
				for( d = 0; d < Dim; d++ )
				{
					radius = maxXY( fabs( par.xBest.x[d] - pb.min[d] ), fabs( par.xBest.x[d] - pb.max[d] ) );
					dx = ( double ) 3.0 * radius;
					pos.x[d] = random_gaussian( par.xBest.x[d], dx );
				}
			if( par.strategy == 1 ) break;
			if( c3 > almostZero ) noise = random_gaussian( 0, ( error1 - error2 ) / c3 );
			else              noise = random_gaussian( 0, 1 );
			for( d = 0; d < pb.D; d++ )
				pos.x[d] = pos.x[d] + noise * ( par.xBest.x[d] - pos.x[d] ); // Add noise
			break;
		case 99: // Standard PSO for testing
			for( d = 0; d < Dim; d++ )
			{
				dx = w2 * ( par.x.x[d] - par.xPrev.x[d] );
				dx += random_double( 0, c ) * ( par.xBest.x[d] - par.x.x[d] ) + random_double( 0, c ) * ( informer.xBest.x[d] - par.x.x[d] );
				pos.x[d] = par.xPrev.x[d] + dx;
			}
			break;
	}
	pos = particle_check( pb, pos );
	pos.f = position_eval( pb, pos );
//	position_print(pos);
//	pos=position_lm(gop,pb,pos);
//	position_print(pos);
	if( multiObj ) position_archive( pos ); // Archive
	return pos;
}

static struct position position_lm( struct opt_data *op, struct problem pb, struct position pos )
{
	struct position p;
	int d;

	p = pos;
	for( d = 0; d < pb.D; d++ )
		op->pd->var[op->pd->var_index[d]] = p.x[d];
	d = gop->cd->neval;
	optimize_lm( op );
	eval += gop->cd->neval - d;
	for( d = 0; d < pb.D; d++ )
		p.x[d] = asin( sin( op->pd->var[op->pd->var_index[d]] ) );
	p.f.f[0] = op->phi;
//	position_print(p);
//	p.f=position_eval(pb,p);
	return p;
}

static void problem_print( struct problem pb )
{
	int d, n;
	tprintf( "\nNumber of parameters (parameter space dimensions) %i\n", pb.D );
	if( debug_level > 1 )
	{
		tprintf( "Search space:\n" );
		for( d = 0; d < pb.D; d++ )
		{
			if( pb.valSize[d] == 0 )
				tprintf( "[%g  %g] dx=%g init %g\n", pb.min[d], pb.max[d], pb.dx[d], pb.ival[d] );
			else
			{
				tprintf( "List: " );
				for( n = 0; n < pb.valSize[d]; n++ )
					tprintf( "%g ", pb.val[d][n] );
				tprintf( "\n" );
			}
		}
	}
	if( pb.fNb > 1 )
		for( n = 0; n < pb.fNb; n++ )
			tprintf( "Objective function #%d: Target %g Acceptable error %g\n", n, pb.objective[n], pb.errorMax.f[n] );
	else if( pb.objective[0] > 0 || pb.errorMax.f[0] > 0 )
		tprintf( "Objective function: Target %g Acceptable error %g\n", pb.objective[0], pb.errorMax.f[0] );
	tprintf( "Maximum number of functional evaluations: %d\n", pb.evalMax );
	if( pb.repeat != 1 ) tprintf( "Number of runs: %i\n", pb.repeat );
}

static struct swarm pso_solver( struct problem pb, int compare_type, int run )
{
	int n, tr, stop = 0;
	struct swarm S = {0};

	for( n = 0; n < pb.fNb; n++ ) wF[n] = 1; // Initials penalties (for multiobjective)
	iterLocalSearchNb = 0; // Prepare local search (for multiobjective)
	S = swarm_init( pb, compare_type ); //Initialisation of the swarm
	if( debug_level == 0 ) { tprintf( "Initial " ); swarm_print( S ); }
	iter = 0;
	iterSwarmAdapt = 0; // Last iteration at which the swarm has been adapted
	for( tr = 0; tr < S.size; tr++ ) iterTribeAdapt[tr] = 0;
	while( stop == 0 )
	{
		for( tr = 0; tr < S.size; tr++ )
			iterTribeAdapt[tr] = iterTribeAdapt[tr] + 1;
		iterSwarmAdapt = iterSwarmAdapt + 1;
		iterSwarmStag = iterSwarmStag + 1;
		archiveVar = archiveFitnessVar();
		S = swarm_move( pb, S, compare_type, run );
		S = swarm_adapt( pb, S, compare_type );
//		S = swarm_lm( pb, S );
		if( !multiObj && compare_particle_fit( S.best.f, pb.errorMax, 1 ) == 1 ) {stop = 1; continue;}
		if( multiObj ) archive_local_search( pb );
		else S = swarm_local_search( pb, S );
		iter++;
		//if(restart==1) {stop=1;continue;} // For future automatic restart
		if( eval >= pb.evalMax )  {stop = 1; continue;}
		if( !multiObj && compare_particle_fit( S.best.f, pb.errorMax, 1 ) == 1 ) {stop = 1; continue;}
	}
	return S;
}

static int sign( double x )
{
	if( x > 0 ) return 1;
	if( x < 0 ) return -1;
	return 0;
}

static struct swarm swarm_adapt( struct problem pb, struct swarm S0, int compare_type )
{
	int adaptSwarm;
	int add_part_count = 0;
	int add_tribe_count = 0;
	int d;
	int disturbPart;
	struct fitness f1, f2;
	int improve_count;
	int initOption;
	int linkNb;
	int pa;
	int partNb;
	int partWorst;
	double pr;
	int s;
	struct swarm St;
	int status;
	int tr;
	int tribeNb;
	int trWorst;

	St = S0;
	tribeNb = St.size;
	for( tr = 0; tr < St.size; tr++ )
	{
		iterTribeAdapt[tr] = 0;
		improve_count = 0;
		for( pa = 0; pa < St.trib[tr].size; pa++ ) // Status of the the tribe
		{
			f1 = St.trib[tr].part[pa].xPrev.f;
			f2 = St.trib[tr].part[pa].x.f;
			if( compare_particle_fit( f2, f1, compare_type ) == 1 ) // If f2 is better than f1
				improve_count++;
		}
		pr = improve_count / ( float )St.trib[tr].size; // Fuzzy rule
		if( random_double( 0, 1 ) > pr ) St.trib[tr].status = 0; // Bad tribe
		else             St.trib[tr].status = 1; // Good tribe
		switch( St.trib[tr].status )
		{
			case 0: // Bad tribe.
				/*
				 The idea is to increase diversity.
				 This is done by adding sometimes a completely new particle
				 and by disturbing another one a bit, typically along
				 just one dimension.
				 */
				disturbPart = random_double( 0, 1 ) < 1 - 1. / St.trib[tr].size;
				if( disturbPart ) // Disturb a particle (not the best one)
				{
					if( St.trib[tr].size > 1 )
					{
						do { pa = random_int( 0, St.trib[tr].size - 1 ); }
						while( pa == St.trib[tr].best );
						d = random_int( 0, pb.D - 1 );
						//d=tribeVarianceMin(St.trib[tr]); // EXPERIMENT
						St.trib[tr].part[pa].xBest.x[d] = random_double( pb.min[d], pb.max[d] );
						// Should be MODIFIED. Too costly for just one modification
						St.trib[tr].part[pa].xBest = particle_check( pb, St.trib[tr].part[pa].xBest );
						St.trib[tr].part[pa].xBest.f = position_eval( pb, St.trib[tr].part[pa].xBest );
					}
				}
				partNb = St.trib[tr].size;
				if( partNb < partMax )
				{
					if( multiObj ) initOption = 3;
					else initOption = random_int( 0, 3 );
					St.trib[tr].part[partNb] = particle_init( pb, initOption, St.best, St.trib[tr].part[St.trib[tr].best].xBest, St );
					St.trib[tr].size = partNb + 1;
					add_part_count++;  // Add a new particle
					f1 = St.trib[tr].part[partNb].xBest.f; // By extraordinary chance this particle might be the best of the tribe
					f2 = St.trib[tr].part[St.trib[tr].best].xBest.f;
					if( compare_particle_fit( f1, f2, compare_type ) == 1 )
						St.trib[tr].best = partNb;
				}
				else
				{
					overSizeTribe++;
					if( debug_level > 0 ) tprintf( "WARNING: Cannot add a particle (increase partMax = %i)\n", partMax );
				}
				break;
			case 1: // Good tribe
				if( St.trib[tr].size < 2 ) break;
				partWorst = 0;
				for( pa = 1; pa < St.trib[tr].size; pa++ ) // Find the worst particle
				{
					f1 = St.trib[tr].part[partWorst].xBest.f;
					f2 = St.trib[tr].part[pa].xBest.f;
					if( compare_particle_fit( f2, f1, compare_type ) == -1 )
						partWorst = pa;
				}
				if( partWorst == St.trib[tr].best ) break; // It might be also the best. In that case, don't remove it
				if( partWorst < St.trib[tr].size - 1 ) // Remove it from the tribe
					for( pa = partWorst; pa < St.trib[tr].size - 1; pa++ )
						St.trib[tr].part[pa] = St.trib[tr].part[pa + 1];
				St.trib[tr].size--;
				break;
		}
	}
	if( St.size > 1 )
	{
		linkNb = St.size * ( St.size - 1 ) / 2;
		adaptSwarm = iterSwarmAdapt >= linkNb / 2;
	}
	else adaptSwarm = TRUE; // Always true if there is just one tribe
	if( adaptSwarm ) // Reinitialise the previous best, and the iteration counter
	{
		St.fBestPrev = St.best.f;
		iterSwarmAdapt = 0;
		// Status of the swarm
		if( compare_particle_fit( St.best.f, St.fBestPrev, compare_type ) <= 0 ) St.status = 0; // Bad swarm
		else St.status = 1; // Good swarm
		switch( St.status )
		{
			case 0: // Bad swarm
				if( St.size < tribMax )
				{
					St.trib[St.size] = tribe_init( pb, 1, compare_type, S0 );
					f1 = St.trib[St.size].part[St.trib[St.size].best].xBest.f;
					f2 = St.best.f;
					if( compare_particle_fit( f1, f2, compare_type ) == 1 )
					{
						St.fBestPrev = St.best.f;
						St.best = St.trib[St.size].part[St.trib[St.size].best].xBest;
					}
					St.size++;
					add_tribe_count++; // Add a new tribe
				}
				else
				{
					overSizeSwarm++;
					if( debug_level > 0 ) tprintf( "WARNING: Cannot add a tribe (increase tribMax = %i)\n", tribMax );
				}
				break;
			case 1: // Good swarm
				if( St.size > 1 )
				{
					trWorst = 0;
					for( tr = 1; tr < St.size; tr++ ) // Find the worst tribe
					{
						f1 = St.trib[trWorst].part[St.trib[trWorst].best].xBest.f;
						f2 = St.trib[tr].part[St.trib[tr].best].xBest.f;
						if( compare_particle_fit( f2, f1, compare_type ) == -1 ) trWorst = tr;
					}
					if( trWorst < St.size - 1 ) // Remove the worst tribe
					{
						for( tr = trWorst; tr < St.size - 1; tr++ )
							St.trib[tr] = St.trib[tr + 1];
					}
					St.size--;
					if( debug_level > 3 ) tprintf( "remove tribe %i\n", trWorst );
				}
				break;
		}
	}
	for( tr = 0; tr < tribeNb; tr++ ) // Potentially modify particle strategy (strategies are not immediately modified for new particles
	{
		for( pa = 0; pa < St.trib[tr].size; pa++ )
		{
			if( St.trib[tr].part[pa].strategy == -1 )  continue;
			f1 = St.trib[tr].part[pa].xPrev.f;
			f2 = St.trib[tr].part[pa].x.f;
			status = 3 * compare_particle_fit( f2, f1, compare_type );
			f1 = St.trib[tr].part[pa].fBestPrev;
			f2 = St.trib[tr].part[pa].xBest.f;
			status += compare_particle_fit( f2, f1, compare_type );
			if( status <= -2 ) St.trib[tr].part[pa].status = -1;
			else
			{
				if( status >= 3 ) St.trib[tr].part[pa].status = 1;
				else St.trib[tr].part[pa].status = 0;
			}
			switch( St.trib[tr].part[pa].status )
			{
				case -1: // Bad particle
					do { s = random_int( 0, 2 ); }
					while( s == St.trib[tr].part[pa].strategy );
					St.trib[tr].part[pa].strategy = s; // Try another strategy
					break;
				case 0: // Normal particle
					St.trib[tr].part[pa].strategy = random_int( 1, 2 );
					break;
				case 1: // Good particle
					if( random_double( 0, 1 ) < 0.5 ) St.trib[tr].part[pa].strategy = 0; // Modify the strategy with a probability 0.5
					break;
			}
		}
	}
	if( debug_level )
	{
		if( add_part_count > 0 || add_tribe_count > 0 ) swarm_print( St );
		else if( debug_level > 1 ) swarm_print( St );
	}
	return St;
}

static void swarm_print( struct swarm S )
{
	int nTotPart = 0, it;
	for( it = 0; it < S.size; it++ )
		nTotPart += S.trib[it].size; // Total number of particles
	if( nTotPart > 1 ) tprintf( "Swarm: %i particles", nTotPart );
	else tprintf( "Swarm: %i particle", nTotPart );
	if( S.size > 1 ) tprintf( " %i tribes", S.size );
	else tprintf( " %i tribe", S.size );
	if( debug_level == 1 )
		for( it = 0; it < S.size; it++ )
			tprintf( " %i", S.trib[it].size );
	tprintf( "\n" );
	if( debug_level > 1 )
		for( it = 0; it < S.size; it++ )
		{
			tprintf( "tribe %i: ", it + 1 );
			tribe_print( S.trib[it] );
		}
}

static struct swarm swarm_init( struct problem pb, int compare_type )
{
	struct fitness f1, f2;
	struct swarm S = {0};
	int partNb;
	int tr;
	int tribNb;

	tribNb = 1; // Initial number of tribes
	partNb = 1; // Initial number of particles in each tribe
	for( tr = 0; tr < tribNb; tr++ )
		S.trib[tr] = tribe_init( pb, partNb, compare_type, S );
	S.size = tribNb; // Number of tribes
	// Warning: must be done AFTER tribeInit()
	// Find the best position of the swarm
	S.best = S.trib[0].part[S.trib[0].best].xBest;
	if( S.size > 1 )
		for( tr = 1; tr < S.size; tr++ )
		{
			f1 = S.trib[tr].part[S.trib[tr].best].xBest.f;
			f2 = S.best.f;
			if( compare_particle_fit( f1, f2, compare_type ) == 1 ) S.best = S.trib[tr].part[S.trib[tr].best].xBest;
		}
	S.fBestPrev = S.best.f;
	S.fBestStag = S.best.f;
	S.status = 0;
	if( debug_level > 0 ) { tprintf( "Initial " ); swarm_print( S ); }
	return S;
}

static struct swarm swarm_local_search( struct problem pb, struct swarm S )
{
	// EXPERIMENT
	struct swarm St;
	int d;
	double dist;
	struct fitness f2;
	int m, mm;
	int out;
	int pa, nTotPart = 0;
	int shaman;
	struct position simplex[fMax + 1];
	int tr;
	struct position xNew;
	//	double w2=0.74;
	double z;
	int option = 1;

	St = S;
	xNew.size = pb.D;
	switch( option )
	{
		case 0:
			for( tr = 0; tr < St.size; tr++ ) // For each tribe
			{
				if( St.trib[tr].size < pb.D + 1 ) continue; // If there is not enough particles do nothing
				if( debug_level > 0 ) tprintf( "Iter %i Eval %d:  Local search ...", iter, eval );
				mm = St.trib[tr].best; // Define a simplex
				for( m = 0; m < pb.D + 1; m++ )
				{
					mm += m; if( mm > St.trib[tr].size - 1 ) mm = 0;
					simplex[m] = St.trib[tr].part[mm].x;
				}
				// Define a new point
				//out=aleaInteger(0,1); // TO TRY
				// out=1;
				out = 0;
				switch( out )
				{
					case 0:	// Inside the simplex
						for( d = 0; d < xNew.size; d++ )
						{
							xNew.x[d] = 0;
							for( m = 0; m < pb.D + 1; m++ )
								xNew.x[d] += simplex[m].x[d] / simplex[m].f.f[0];
						}
						for( d = 0; d < xNew.size; d++ )
							xNew.x[d] /= ( pb.D + 1 );
						break;
					case 1: // Outside the simplex
						for( d = 0; d < xNew.size; d++ )
						{
							xNew.x[d] = 0;
							for( m = 1; m < pb.D + 1; m++ ) // Simplex dimension D
								xNew.x[d] += simplex[m].x[d];
						}
						for( d = 0; d < xNew.size; d++ )
							xNew.x[d] /= pb.D; // Gravity center
						for( d = 0; d < xNew.size; d++ )
							xNew.x[d] -= ( simplex[0].x[d] - xNew.x[d] ); // Reflection of the first vertex
						break;
				}
				xNew = particle_check( pb, xNew );
				xNew.f = position_eval( pb, xNew );
				f2 = St.trib[tr].part[St.trib[tr].best].xBest.f;
				if( compare_particle_fit( xNew.f, f2, 0 ) == 1 ) // Possibly update tribe best
				{
					St.trib[tr].part[St.trib[tr].best].xBest = xNew;
					St.trib[tr].part[St.trib[tr].best].x = xNew;
				}
				f2 = St.best.f;
				if( compare_particle_fit( xNew.f, f2, 0 ) == 1 ) St.best = xNew; // Possibly update swarm best
			} // end "for(tr
			break;
		case 1: // "Around" the shaman
			for( tr = 0; tr < St.size; tr++ ) // For each tribe
			{
				if( St.trib[tr].size < 2 ) continue; // Small size if(St.trib[tr].size<pb.D+1) continue;
				shaman = St.trib[tr].best;
				for( d = 0; d < pb.D; d++ )
				{
					dist = almostInfinite;
					for( pa = 0; pa < St.trib[tr].size; pa++ ) // find the nearest position
					{
						if( pa == shaman ) continue;
						z = fabs( St.trib[tr].part[pa].xBest.x[d] - St.trib[tr].part[shaman].xBest.x[d] );
						if( z < dist ) dist = z;
					}
					z = St.trib[tr].part[shaman].xBest.x[d]; // define a random intermediate coordinate
					//xNew.x[d]=z+RandomGauss(z,w2*dist);
					xNew.x[d] = z + ( 1 - 2 * random_double( 0, 1 ) ) * random_double( 0, dist );
				}
				xNew = particle_check( pb, xNew );
				xNew.f = position_eval( pb, xNew );
				if( compare_particle_fit( xNew.f, St.trib[tr].part[shaman].xBest.f, 0 ) == 1 ) // Update tribe best
				{
					St.trib[tr].part[shaman].xBest = xNew;
					St.trib[tr].part[shaman].x = xNew;
				}
				if( compare_particle_fit( xNew.f, St.best.f, 0 ) == 1 )
					St.best = xNew; // Swarm best
				nTotPart += S.trib[tr].size; // Total number of particles
			}
			break;
	}
	return St;
}

static struct swarm swarm_move( struct problem pb, struct swarm S, int compare_type, int run )
{
	struct fitness f1, f2;
	struct particle informer = {0};
	int pa;
	int sh = 0;
	int shamanNb;
	int shBest;
	int shList[tribMax];
	struct swarm St;
	int tr;

	St = S;
	modify_weights( pb.fNb, run ); // Penalties (global variable wF)
	St.fBestPrev = S.best.f; // Save the currenr best result of the whole swarm
	for( tr = 0; tr < St.size; tr++ ) // Cicle all tribes
	{
		if( debug_level > 2 )
		{
			tprintf( "Best result witin tribe #%i", tr );
			fitness_print( St.trib[tr].part[St.trib[tr].best].xBest.f );
		}
		St.trib[tr].fBestPrev = St.trib[tr].part[St.trib[tr].best].xBest.f; // Save the previous best result
		for( pa = 0; pa < St.trib[tr].size; pa++ )
		{
			St.trib[tr].part[pa].xPrev = St.trib[tr].part[pa].x; // Save previous information
			St.trib[tr].part[pa].fBestPrev = St.trib[tr].part[pa].xBest.f;
			if( pa != St.trib[tr].best ) // ... if it is not the shaman
				informer = St.trib[tr].part[St.trib[tr].best];
			else // For the shaman, it is a bit more complicated
			{
				if( multiObj ) // Select an informer in the archive
					informer.xBest = archiveCrowDistSelect( St.size );
				else // Select an informer (attractor) in the swarm
				{
					if( St.trib[tr].size == 1 ) // Just one particle => the informer is the shaman itself
						informer = St.trib[tr].part[pa];
					else // There are several tribes. Look for an external informer
					{
						shamanNb = random_int( 1, St.size ); // Number of informer shamans
						for( sh = 0; sh < shamanNb; sh++ ) // Build a random list of shamans
							shList[sh] = random_int( 0, St.size - 1 );
						shBest = shList[0];
						if( shamanNb > 1 )
							for( sh = 1; sh < shamanNb; sh++ ) // Find the best
							{
								f1 = St.trib[shBest].part[St.trib[shBest].best].xBest.f;
								f2 = St.trib[shList[sh]].part[St.trib[shList[sh]].best].xBest.f;
								if( compare_particle_fit( f2, f1, compare_type ) == 1 ) shBest = shList[sh];
							}
						f1 = St.trib[shBest].part[St.trib[shBest].best].xBest.f;
						f2 = St.trib[tr].part[pa].xBest.f;
						if( compare_particle_fit( f1, f2, compare_type ) )
							informer = St.trib[shBest].part[St.trib[shBest].best]; // The informer is another shaman
						else
							informer = St.trib[tr].part[pa]; // Not better than itself => no external informer
					}
				}
				if( debug_level > 2 )
				{
					tprintf( "particle %i is informed by %i", St.trib[tr].part[pa].label, informer.label );
					fitness_print( informer.xBest.f );
					tprintf( "\n" );
				}
			}
			St.trib[tr].part[pa].x = position_update( pb, St.trib[tr].part[pa], informer ); // Move the particle
			f1 = St.trib[tr].part[pa].x.f;
			f2 = St.trib[tr].part[pa].xBest.f;
			if( compare_particle_fit( f1, f2, compare_type ) == 1 ) // Possibly update xBest of the particle
			{
				St.trib[tr].part[pa].xBest = St.trib[tr].part[pa].x;
				if( debug_level > 2 )
				{
					tprintf( "improvement => " );
					fitness_print( f1 );
					tprintf( "\n" );
				}
			}
			f2 = St.trib[tr].part[St.trib[tr].best].xBest.f;
			if( compare_particle_fit( f1, f2, compare_type ) == 1 ) // Possibly update tribe best
			{
				St.trib[tr].best = pa;
				if( debug_level > 2 )
				{
					tprintf( "best particle is now %i", St.trib[tr].part[pa].label );
					fitness_print( f1 );
					tprintf( "\n" );
				}
			}
			f2 = St.best.f;
			if( compare_particle_fit( f1, f2, compare_type ) == 1 ) // Possibly update swarm best
				St.best = St.trib[tr].part[pa].x;
		}
	}
	return St;
}

static int swarm_particle_count( struct swarm S ) // Compute the total number of particles
{
	int partNb = 0;
	int tr;
	for( tr = 0; tr < S.size; tr++ )
		partNb += S.trib[tr].size;
	return partNb;
}

static int tribe_best_particle( struct tribe T, int compare_type )
{
	// Find the best particle (shaman)
	int ip;
	int ibest = 0;
	if( T.size > 1 )
		for( ip = 1; ip < T.size; ip++ )
			if( compare_particle_fit( T.part[ip].xBest.f, T.part[ibest].xBest.f, compare_type ) == 1 )
				ibest = ip;
	return ibest;
}

static void tribe_print( struct tribe T )
{
	int pa;
	tprintf( " %i particles\n", T.size );
	tprintf( "Labels    :" );
	for( pa = 0; pa < T.size; pa++ )
		tprintf( " %3i", T.part[pa].label );
	tprintf( "\nStrategies:" );
	for( pa = 0; pa < T.size; pa++ )
		tprintf( " %3i", T.part[pa].strategy );
	tprintf( "\n" );
}

static struct tribe tribe_init( struct problem pb, int partNb, int compare_type, struct swarm S )
{
	struct tribe newTribe;
	struct position dumm;
	int ip, init_option;

	newTribe.best = 0;
	newTribe.status = 0;
	newTribe.size = minXY( partNb, partMax );
	for( ip = 0; ip < newTribe.size; ip++ )
	{
		if( pb.init == 1 ) { pb.init = 0; init_option = 5; } // User provided initial values
		else init_option = 3;
		newTribe.part[ip] = particle_init( pb, init_option, dumm, dumm, S );
	}
	newTribe.best = tribe_best_particle( newTribe, compare_type );
	newTribe.fBestPrev = newTribe.part[newTribe.best].xBest.f;
	return newTribe;
}

static int tribe_varmin_dimension( struct tribe T )
{
	int i, ip, dim_varmin = -1;
	double mean_dim, var_dim, min_var, z;
	min_var = almostInfinite; // Arbitrary very big value
	for( i = 0; i < T.part[0].x.size; i++ )
	{
		mean_dim = 0;
		for( ip = 0; ip < T.size; ip++ )
			mean_dim += T.part[ip].xBest.x[i]; // computed on the xBest positions
		mean_dim /= T.size;
		var_dim = 0;
		for( ip = 0; ip < T.size; ip++ )
		{
			z = T.part[ip].xBest.x[i] - mean_dim;
			var_dim += z * z;
		}
		var_dim /= T.size;
		if( var_dim < min_var ) { dim_varmin = i; min_var = var_dim; }
	}
	return dim_varmin; // Return the dimension on which the variance of the tribe is minimum
}

static void modify_weights( int fNb, int run )
{
	// Modify the penalties
	// Global variables: wF[]
	// NOTE: this is the most empirical part of the algorithm
	// It may certainly be largely improved
	int deb;
	int n;
	int m;
	m = random_int( 0, 1 );
	switch( m )
	{
		case 0:
			for( n = 0; n < fNb; n++ ) wF[n] = 0;
			m = random_int( 0, fNb - 1 ); wF[m] = 1;
			break;
		case 1: // (Extended) Dynamic Weighted Aggregation
			deb = random_int( 0, fNb - 1 );
			for( n = 0; n < fNb; n++ )
			{
				m = n + deb; if( m >= fNb ) m = fNb - m;
				wF[m] = fabs( sin( 2 * M_PI * ( n + 1 ) / fNb ) );
			}
			break;
	}
}
//================================================== KISS
/*
 A good pseudo-random numbers generator

 The idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
 x(n)=a*x(n-1)+1 mod 2^32
 y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
 z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
 2^32*(2^32-1)*(2^63+2^32-1) > 2^127
 */

static void seed_rand_kiss( unsigned long seed )
{
	kiss_x = seed | 1;
	kiss_y = seed | 2;
	kiss_z = seed | 4;
	kiss_w = seed | 8;
	kiss_carry = 0;
}

static unsigned long rand_kiss()
{
	kiss_x = kiss_x * 69069 + 1;
	kiss_y ^= kiss_y << 13;
	kiss_y ^= kiss_y >> 17;
	kiss_y ^= kiss_y << 5;
	kiss_k = ( kiss_z >> 2 ) + ( kiss_w >> 3 ) + ( kiss_carry >> 2 );
	kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
	kiss_z = kiss_w;
	kiss_w = kiss_m;
	kiss_carry = kiss_k >> 30;
	return kiss_x + kiss_y + kiss_w;
}
