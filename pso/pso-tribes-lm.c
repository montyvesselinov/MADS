// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
//
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
*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "../mads.h"
#include "pso.h"

#define TRUE 1
#define FALSE 0
#define MAXPHI 4+1 // Maximum number of objective functions Note: number of "real" objective functions +1
#define MAXRUNS 500 // Maximum number of runs
#define MAXPART 100
#define MAXTRIBE 60 //Maximum number of tribes
#define MAXVALUES 11 // Maximum number of acceptable values on a dimension
#define MAXARCHIVE 100 // Maximum number of archived positions in the case of multiobjective optimization

int pso_tribes( struct opt_data *op );
// Specific to multiobjective
struct archived { struct position x; double crow_dist; };
struct distRank { double dist; int rank; };
//----------------------------------------- Subroutines
void archive_crowding_dist();
void archive_print();
void archive_local_search( struct problem *pb );
void archive_save( struct archived archiv[], int nArchive, FILE *fArchive );
double archive_spread();
struct objfunc archive_phi_var( struct problem *pb );
struct position archiveCrowDistSelect( int size );
int adapt_compare( struct objfunc phi[], int run, int compare_type );
int compare_particles( struct objfunc *phi1,  struct objfunc *phi2, int compare_type );
void   objfunc_print( struct objfunc *phi );
double objfunc_dist( struct objfunc *phi1,  struct objfunc *phi2 );
double objfunc_total( struct objfunc *phi, int type );
double granularity( double value, double granul );
void particle_init( struct problem *pb, int option, struct position *guide1, struct position *guide2, struct swarm *S, struct particle *P );
void position_archive( struct position *pos );
void position_print( struct position *pos );
void position_save( FILE *fRun, struct position *pos, int run );
void position_eval( struct problem *pb, struct position *P );
void position_check( struct problem *pb, struct position *P );
void position_lm( struct opt_data *op, struct problem *pb, struct position *P );
void position_update( struct problem *pb, struct particle *P, struct particle informer );
void problem_init( struct opt_data *op, struct problem *pr );
void problem_print( struct problem *pb );
void pso_solver( struct problem *pb, int compare_type, int run, struct swarm *S );
void swarm_adapt( struct problem *pb, struct swarm *S, int compare_type );
void swarm_print( struct swarm *S );
void swarm_init( struct problem *pb, int compare_type, struct swarm *S );
void swarm_local_search( struct problem *pb, struct swarm *S );
void swarm_move( struct problem *pb, struct swarm *S, int compare_type, int run );
int  swarm_particle_count( struct swarm *S );
int  tribe_shaman( struct tribe *T, int compare_type );
void tribe_print( struct tribe *T );
void tribe_init( struct problem *pb, int nPart, int compare_type, struct swarm *S, struct tribe *T );
int  tribe_varmin_dimension( struct tribe *T );
void modify_weights( int nPhi, int run );
void swarm_lm( struct problem *pb, struct swarm *S );
void seed_rand_kiss( unsigned int seed );
unsigned int rand_kiss();
static int compare_crowding_dist( void const *a, void const *b ); // For qsort
static int compare_dist_rank( void const *a, void const *b );
static int compare_fit( void const *a, void const *b );
static int compare_double( void const *a, void const *b ); // For qsort
double random_double( double a, double b );
double random_gaussian( double mean, double std_dev );
double maxXY( double x, double y );
double minXY( double x, double y );
int random_int( int a, int b );
int sign( double x );
void set_tribe( struct problem *pb, struct tribe *T );
void set_particle( struct problem *pb, struct particle *P );
void set_swarm( struct problem *pb, struct swarm *S );
void set_position( struct problem *pb, struct position *pos );
void set_objfunc( struct problem *pb, struct objfunc *phi );
void copy_position( struct position *ps, struct position *pd );
void copy_particle( struct particle *ps, struct particle *pd );
void free_tribe( struct tribe *T );
void free_swarm( struct swarm *S );
void free_particle( struct particle *P );
void free_position( struct position *pos );
void free_objfunc( struct objfunc *phi );
double irand();
// Elsewhere ...
void Transform( double *x, void *data, double *f );
void DeTransform( double *x, void *data, double *f );
int get_seed( );
int optimize_lm( struct opt_data *op );
float **float_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
//----------------------------------------- Global variables
int lmo_count;
int lmo_flag;
struct archived multiobj_archive[MAXARCHIVE + 1];
int arch = 0;
int nArchive = 0;
struct objfunc archiveVar;
double epsilon_vector[MAXPHI]; // For epsilon-dominance
int eval; // Number of functional evaluations
struct objfunc phi_list[MAXRUNS];
int compare_type;
int fn;
double phi_max[MAXPHI]; // Maximum objfunc value found during the process
double phi_min[MAXPHI]; // Minimum objfunc value found during the process
int iter = 0;
int nLocalSearchIter = 0;
int nSwarmAdaptIter = 0; // Number of iterations between two swarm adaptations
int count_not_good_tribes = 0;
int nTribeAdaptIter[MAXTRIBE]; // The same, but for each tribe
int id_global = 0; // id (integer number) of the last generated particle
int multiObj = FALSE; // Flag for multi-objective problem
int nExceedSizeSwarm; // Number of times the swarm tends to generate too many tribes
int nExceedSizeTribe; // Number of times a tribe tends to generate too many particles
int restart;
int nRestarts;
int debug_level;
int lm_wait = 0;
double phi_lm_wait = HUGE_VAL;
double phi_lm_init;
double phi_weights[MAXPHI + 1]; // Dynamic penalties

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

int pso_tribes( struct opt_data *op )
{
	FILE *fRun; // To save the run
	FILE *fArchive; // To save the Pareto front
	struct problem pb;
	struct swarm S = {0};
	struct position bestBest;
	double successRate[MAXPHI], errorMean[MAXPHI], errorTot;
	int debug, i, n, r, eval_total;
	char filename[255];
	pb.nPhi = 1;
	gop = op;
	irand_seed = &op->cd->seed;
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); }
	else { seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); if( op->cd->pdebug ) tprintf( "Current seed: %d\n", op->cd->seed ); }
	nExceedSizeSwarm = nExceedSizeTribe = 0;
	pb.lm_factor = op->cd->lm_factor;
	if( ( res = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); exit( 1 ); }
	if( op->cd->pdebug )
	{
		sprintf( filename, "%s.runs", op->root );
		fRun = fopen( filename, "w" );
		sprintf( filename, "%s.archive", op->root );
		fArchive = fopen( filename, "w" );
	}
	problem_init( op, &pb );
	if( op->cd->pdebug ) problem_print( &pb );
	set_swarm( &pb, &S );
	set_position( &pb, &bestBest );
	set_position( &pb, &( pb ).pos_success );
	errorTot = 0;
	eval_total = 0;
	nArchive = 0;
	compare_type = 0; // Kind of comparison, to begin (see objfuncCompare() )
	bestBest.size = multiobj_archive[0].x.size = pb.nPhi; //Prepare final result
	for( n = 0; n < pb.nPhi; n++ )
	{
		bestBest.f.f[n] = HUGE_VAL;
		successRate[n] = 0;
	}
	if( multiObj )
	{
		// Prepare epsilon-dominance
		// as this setting is outside of the loop on runs
		// each run will take advantage of the previous ones
		// to set epsilon-dominance, set epsilon[n] to huge number
		for( n = pb.nPhi; n < pb.nPhi; n++ ) epsilon_vector[n] = DBL_EPSILON;
	}
	else
		for( n = 0; n < pb.nPhi; n++ ) epsilon_vector[n] = DBL_EPSILON;
	for( n = 0; n < pb.nPhi; n++ )
	{
		// Very big value to begin
		phi_min[n] = 1e30;
		// Very low value to begin
		phi_max[n] = 0;
	}
	lmo_count = restart = nRestarts = 0;
	lmo_flag = 0;
	if( strstr( op->cd->opt_method, "lm" ) != NULL || strncmp( op->cd->opt_method, "squad", 5 ) == 0 ) { op->cd->squads = 1; lmo_flag = 1; }
	if( lmo_flag )
		tprintf( "SQUADS: Coupled Particle-Swarm and Levenberg-Marquardt Optimization ... " );
	else
		tprintf( "Particle-Swarm Optimization TRIBES ... " );
	if( op->cd->pdebug )  tprintf( "\n" );
	else fflush( stdout );
	for( i = 0; i < MAXRUNS; i++ )
		set_objfunc( &pb, &phi_list[i] );
	for( r = 0; r < pb.repeat; r++ )
	{
		if( restart == 0 )
		{
			eval = op->cd->neval;
			if( pb.repeat > 1 ) tprintf( "\nRun #%i compare_type %i\n", r + 1, compare_type );
		}
		else restart = 0; // TODO this statement does not make sense; needs work!
		pso_solver( &pb, compare_type, r, &S );
		if( op->cd->pdebug ) { tprintf( "\nFinal results: " ); swarm_print( &S ); }
		for( i = 0; i < pb.nPhi; i++ )
			phi_list[r].f[i] = S.best.f.f[i];
		if( r < pb.repeat - 1 )
			compare_type = adapt_compare( phi_list, r, compare_type ); // Define the next comparison method
		if( restart == 1 )
		{
			nRestarts++;
			tprintf( "Restart at evaluation #%d\n", eval );
			r--;
		}
		else
		{
			if( nRestarts > 0 ) tprintf( "Restarts = %i\n", nRestarts );
			op->phi = S.best.f.f[0];
			if( op->cd->pdebug > 1 )
			{
				debug = op->cd->fdebug; op->cd->fdebug = 3;
				func_global( S.best.x, op, res ); // evaluate the best result
				op->cd->fdebug = debug;
				objfunc_print( &S.best.f );
				tprintf( " current position %d: ", r ); position_print( &S.best );
			}
			if( op->cd->pdebug ) position_save( fRun, &S.best, r );
			for( n = 0; n < pb.nPhi; n++ )
			{
				errorMean[n] += S.best.f.f[n];
				if( S.best.f.f[n] < pb.maxError.f[n] ) successRate[n]++;
			}
			if( compare_particles( &S.best.f, &bestBest.f, 2 ) == 1 ) copy_position( &S.best, &bestBest );
			eval_total += eval;
		}
	}
	if( op->cd->pdebug )
	{
		if( lmo_flag ) tprintf( "\nOptimization completed after %i swarm steps, %d LM steps, %d evaluations:\n", iter, lmo_count, eval );
		else tprintf( "\nOptimization completed after %i swarm steps and %d evaluations:\n", iter, eval );
	}
	if( op->cd->pdebug > 1 )
	{
		tprintf( "\nOptimization results:\n\n" );
		if( op->cd->fdebug < 3 ) { debug = op->cd->fdebug; op->cd->fdebug = 3; }
		else debug = 0;
	}
	func_global( bestBest.x, op, res ); // evaluate the best BEST result
	if( op->cd->pdebug > 1 && debug ) op->cd->fdebug = debug;
	for( n = 0; n < pb.nPhi; n++ )
	{
		errorMean[n] /= pb.repeat;
		successRate[n] /= pb.repeat;
	}
	if( op->cd->pdebug && pb.repeat > 1 )
	{
		tprintf( "\nMean errors:" );
		for( n = 0; n < pb.nPhi; n++ ) tprintf( " %g", errorMean[n] );
		tprintf( "\nSuccess rate: " );
		for( n = 0; n < pb.nPhi; n++ ) tprintf( " %g", successRate[n] );
		tprintf( "\nBest result: " );
		objfunc_print( &bestBest.f );
	}
	if( pb.repeat > 1 )
	{
		tprintf( "Number of optimization retries %d\n", pb.repeat );
		tprintf( "Average number of evaluations per optimization %g\n", ( ( double ) eval_total / pb.repeat ) );
	}
	if( nExceedSizeTribe > 0 ) tprintf( "WARNING: Tribe size has been constrained %i times\n", nExceedSizeTribe );
	if( nExceedSizeSwarm > 0 ) tprintf( "WARNING: Swarm size has been constrained %i times\n", nExceedSizeSwarm );
	if( multiObj && op->cd->pdebug ) archive_save( multiobj_archive, nArchive, fArchive );
	DeTransform( bestBest.x, op, bestBest.x );
	for( i = 0; i < op->pd->nOptParam; i++ )
		op->pd->var[op->pd->var_index[i]] = bestBest.x[i];
	if( op->cd->pdebug ) { fclose( fRun ); fclose( fArchive ); }
	free( res );
	free_position( &bestBest );
	free_position( &( pb ).pos_success );
	free_objfunc( &pb.maxError );
	free_swarm( &S );
	free( pb.ival ); free( pb.min ); free( pb.max ); free( pb.dx ); free( pb.valSize ); free( pb.objective ); free( pb.code );
	free_matrix( ( void ** ) pb.val, pb.D );
	for( i = 0; i < MAXRUNS; i++ )
		free_objfunc( &phi_list[i] );
	op->cd->lmstandalone = 1;
	return EXIT_SUCCESS;
}

void set_swarm( struct problem *pb, struct swarm *S )
{
	int i;
	S->size = 0;
	S->status = 0;
	if( ( S->trib = ( struct tribe * ) malloc( MAXTRIBE * sizeof( struct tribe ) ) ) == NULL ) { tprintf( "No memory 1!\n" ); exit( 1 ); }
	for( i = 0; i < MAXTRIBE; i++ )
		set_tribe( pb, &S->trib[i] );
	set_position( pb, &S->best );
	set_objfunc( pb, &S->fBestPrev );
	set_objfunc( pb, &S->fBestStag );
}

void free_swarm( struct swarm *S )
{
	int i;
	for( i = 0; i < MAXTRIBE; i++ )
		free_tribe( &S->trib[i] );
	free_position( &S->best );
	free_objfunc( &S->fBestPrev );
	free_objfunc( &S->fBestStag );
}

void copy_tribe( struct tribe *ps, struct tribe *pd )
{
	int i;
	pd->size = ps->size;
	pd->best = ps->best;
	pd->status = ps->status;
	for( i = 0; i < ps->size; i++ ) copy_particle( &ps->part[i], &pd->part[i] );
	for( i = 0; i < ps->fBestPrev.size; i++ ) pd->fBestPrev.f[i] = ps->fBestPrev.f[i];
}

void set_tribe( struct problem *pb, struct tribe *T )
{
	int i;
	T->size = 0;
	T->status = 0;
	if( ( T->part = ( struct particle * ) malloc( MAXPART * sizeof( struct particle ) ) ) == NULL ) { tprintf( "No memory 2!\n" ); exit( 1 ); }
	for( i = 0; i < MAXPART; i++ )
		set_particle( pb, &T->part[i] );
	set_objfunc( pb, &T->fBestPrev );
}

void free_tribe( struct tribe *T )
{
	int i;
	for( i = 0; i < MAXPART; i++ )
		free_particle( &T->part[i] );
	free_objfunc( &T->fBestPrev );
}

void copy_particle( struct particle *ps, struct particle *pd )
{
	int i;
	pd->id = ps->id;
	pd->strategy = ps->strategy;
	pd->status = ps->status;
	copy_position( &ps->x, &pd->x );
	copy_position( &ps->xBest, &pd->xBest );
	copy_position( &ps->xLast, &pd->xLast );
	for( i = 0; i < ps->fBestPrev.size; i++ ) pd->fBestPrev.f[i] = ps->fBestPrev.f[i];
}

void free_particle( struct particle *P )
{
	free_position( &P->x );
	free_position( &P->xBest );
	free_position( &P->xLast );
	free_objfunc( &P->fBestPrev );
}

void set_particle( struct problem *pb, struct particle *P )
{
	P->strategy = -1;
	set_position( pb, &P->x );
	set_position( pb, &P->xBest );
	set_position( pb, &P->xLast );
	set_objfunc( pb, &P->fBestPrev );
}

void copy_position( struct position *ps, struct position *pd )
{
	int i;
	if( ps->size == 0 ) return;
	pd->size = ps->size;
	pd->f.size = ps->f.size;
	for( i = 0; i < ps->size; i++ ) pd->x[i] = ps->x[i];
	for( i = 0; i < ps->f.size; i++ ) pd->f.f[i] = ps->f.f[i];
}

void free_position( struct position *pos )
{
	free( pos->x );
	free_objfunc( &pos->f );
}

void set_position( struct problem *pb, struct position *pos )
{
	pos->size = ( *pb ).D;
	if( ( pos->x = ( double * ) malloc( ( *pos ).size * sizeof( double ) ) ) == NULL ) { tprintf( "No memory 3!\n" ); exit( 1 ); }
	set_objfunc( pb, &pos->f );
}

void free_objfunc( struct objfunc *phi )
{
	free( phi-> f );
}

void set_objfunc( struct problem *pb, struct objfunc *phi )
{
	phi->size = ( *pb ).nPhi;
	if( ( phi->f = ( double * ) malloc( ( *phi ).size * sizeof( double ) ) ) == NULL ) { tprintf( "No memory 4!\n" ); exit( 1 ); }
}

void position_eval( struct problem *pb, struct position *pos )
{
	double f;
	int i;
	( *pos ).f.size = ( *pb ).nPhi;
	if( gop->global_success )
		f = HUGE_VAL;
	else
	{
		func_global( ( *pos ).x, gop, res ); // evaluation ... either internal of external
		f = gop->phi;
		eval++;
	}
	( *pos ).f.f[0] = fabs( f - ( *pb ).objective[0] );
	for( i = 0; i < ( *pb ).nPhi; i++ ) // Save the min and the max objfunc ever found
	{
		phi_min[i] = minXY( phi_min[i], ( *pos ).f.f[i] );
		phi_max[i] = maxXY( phi_max[i], ( *pos ).f.f[i] );
	}
	if( multiObj && nArchive > 0 ) // Additional "objfunc" for multi-objective
	{
		for( i = 0; i < ( *pb ).nPhi; i++ )
			multiobj_archive[nArchive].x.f.f[i] = ( *pos ).f.f[i]; // Temporary put in the archive
		nArchive++;
		( *pos ).f.f[( *pb ).nPhi] = archive_spread();
		nArchive--; // Virtually removed from the archive
	}
	if( gop->global_success == 0 && ( gop->cd->check_success && gop->success ) )
	{
		if( debug_level && gop->cd->fdebug == 0 ) tprintf( "PSO Success: Predictions are within the predefined calibration bounds (within position_eval)!\n" );
		copy_position( pos, &( *pb ).pos_success );
		gop->global_success = 1;
	}
	if( debug_level > 2 )
	{
		if( debug_level == 3 ) { tprintf( "eval: " ); objfunc_print( &( *pos ).f ); }
		else if( debug_level > 3 ) { tprintf( "eval: " ); position_print( pos ); }
	}
}

void problem_init( struct opt_data *op, struct problem *pb )
{
	int d;
	double *opt_var, *tr_var;
	( *pb ).D = op->pd->nOptParam;
	if( ( opt_var = ( double * ) malloc( ( *pb ).D * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( tr_var = ( double * ) malloc( ( *pb ).D * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->ival = ( double * ) malloc( ( *pb ).D * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->min = ( double * ) malloc( ( *pb ).D * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->max = ( double * ) malloc( ( *pb ).D * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->dx = ( double * ) malloc( ( *pb ).D * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->valSize = ( int * ) malloc( ( *pb ).D * sizeof( int ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->objective = ( double * ) malloc( MAXPHI * sizeof( double ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	if( ( pb->code = ( int * ) malloc( MAXPHI * sizeof( int ) ) ) == NULL ) { tprintf( "No memory!\n" ); exit( 1 ); }
	pb->val = double_matrix( ( *pb ).D, MAXVALUES );
	set_objfunc( pb, &pb->maxError );
	for( d = 0; d < ( *pb ).D; d++ )
		opt_var[d] = op->pd->var[ op->pd->var_index[d] ];
	Transform( opt_var, op, tr_var );
	op->cd->fdebug = 10; func_global( tr_var, op, op->od->res ); op->cd->fdebug = 0;
	if( op->cd->sintrans == 0 )
		for( d = 0; d < ( *pb ).D; d++ ) // typically used for test problems
		{
			( *pb ).ival[d] = tr_var[d];
			( *pb ).min[d] = op->pd->var_min[op->pd->var_index[d]];
			( *pb ).max[d] = op->pd->var_max[op->pd->var_index[d]];
			( *pb ).dx[d] = op->pd->var_dx[op->pd->var_index[d]]; // defines discretization
			( *pb ).valSize[d] = 0;
		}
	else
		for( d = 0; d < ( *pb ).D; d++ ) // typically used for actual problems
		{
			( *pb ).ival[d] = tr_var[d];
			( *pb ).min[d] = -M_PI / 2;
			( *pb ).max[d] = M_PI / 2;
			( *pb ).dx[d] = 0; // defines discretization
			( *pb ).valSize[d] = 0;
		}
	if( ( *pb ).nPhi > MAXPHI - 1 )
	{
		tprintf( "ERROR: Too many multi-objective functions: %i (max %i) \n", ( *pb ).nPhi, MAXPHI );
		exit( 1 );
	}
	( *pb ).code[0] = -1;
	( *pb ).objective[0] = 0;
	( *pb ).maxError.f[0] = op->cd->phi_cutoff;
	( *pb ).maxError.size = ( *pb ).nPhi;
	( *pb ).maxEval = op->cd->maxeval;
	if( op->cd->nretries == 0 )( *pb ).repeat = 1;
	else( *pb ).repeat = op->cd->nretries;
	op->cd->lmstandalone = 0;
	if( ( *pb ).repeat > MAXRUNS )
	{
		tprintf( "MADS Quits: Too many retries (max %i < %i) \n", MAXRUNS, ( *pb ).repeat );
		exit( 1 );
	}
	debug_level = op->cd->pdebug;
	multiObj = FALSE;
	free( opt_var ); free( tr_var );
}

double random_double( double a, double b )
{
	double r;
//	r = ( double ) rand_kiss() / UINT32_MAX;
//	r = ( double ) rand() / RAND_MAX;
	r = irand();
	r = a + r * ( b - a );
	return r;
}

double random_gaussian( double mean, double stddev )
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

int random_int( int min, int max )
{
	double  r;
	if( min >= max ) return min;
	r = random_double( min, max );
	return( ( int ) floor( r + 0.5 ) );
}

struct position archiveCrowDistSelect( int size )
{
	// Choose a random position, according to a non uniform distribution
	// of the crowding distance
	// the bigger the CD, the higher the probability
	// Also the bigger the number of tribes (i.e. "size")
	// the "sharper" the distribution

	int n;
	double pr;
	archive_crowding_dist(); // Compute the Crowding Distances
	qsort( multiobj_archive, nArchive, sizeof( multiobj_archive[0] ), compare_crowding_dist ); // Sort the archive by increasing Crowding Distance
	pr = 2 * log( ( double ) 1 + size ); // Choose at random according to a non uniform distribution:
	n = pow( random_double( 0, pow( ( double ) nArchive - 1, pr ) ), ( double ) 1. / pr );
	return multiobj_archive[n].x;
}

void archive_crowding_dist() // Compute the crowding distances in archive[n] (global variable)
{
	double dist, max, min;
	int nPhi, n;
	nPhi = multiobj_archive[0].x.f.size;
	for( n = 0; n < nArchive; n++ ) multiobj_archive[n].crow_dist = 1;
	for( fn = 0; fn < nPhi; fn++ )
	{
		qsort( multiobj_archive, nArchive, sizeof( multiobj_archive[0] ), compare_fit ); // Sort archive according to f[fn]  NOTE: fn is a global variable
		for( n = 0; n < nArchive - 1; n++ ) // For each position find the distance to the nearest one
		{
			if( n == 0 ) min = phi_min[fn]; else min = multiobj_archive[n - 1].x.f.f[fn];
			if( n == nArchive - 2 ) max = phi_max[fn]; else max = multiobj_archive[n + 1].x.f.f[fn];
			dist = max - min;
			if( dist < epsilon_vector[fn] ) epsilon_vector[fn] = dist; // Prepare epsilon-dominance
			multiobj_archive[n].crow_dist = multiobj_archive[n].crow_dist * dist; // Contribution to the hypervolume
		}
	}
}

void archive_print()
{
	int n;
	tprintf( "Archive %i positions\n", nArchive );
	for( n = 0; n < nArchive; n++ )
	{
		position_print( &multiobj_archive[n].x );
		tprintf( " crowd %g", multiobj_archive[n].crow_dist );
	}
	tprintf( "\n" );
}

struct objfunc archive_phi_var( struct problem *pb ) // Variance of the archive for each objfunc
{
	// More precisely return for each objfunc g(var) where g is a decreasing funtion of var:
	// The smaller the better
	double mean, var, z;
	struct objfunc var_phi;
	int m, n;

	set_objfunc( pb, &var_phi );
	var_phi.size = multiobj_archive[0].x.f.size;

	for( n = 0; n < var_phi.size; n++ )
	{
		mean = 0;
		for( m = 0; m < nArchive; m++ )
			mean += multiobj_archive[m].x.f.f[n];
		mean /= nArchive;
		var = 0;
		for( m = 0; m < nArchive; m++ )
		{
			z = multiobj_archive[m].x.f.f[n] - mean;
			var += z * z;
		}
		var_phi.f[n] = var / nArchive;
	}
	for( n = 0; n < var_phi.size; n++ )
		var_phi.f[n] = 1. / ( 1 + var_phi.f[n] ); // Transform into something to minimise
	free_objfunc( &var_phi );
	return var_phi;
}

void archive_local_search( struct problem *pb )
{
	/*
	 Let nPhi be the dimension of the objfunc space.
	 Loop:
	 Define a simplex of dimension nPhi
	 Find a point "inside". Check if it is a good one
	 or
	 Find a point "outside". Check if it is a good one
	*/
	struct distRank dR[MAXARCHIVE - 1];
	struct position simplex[MAXPHI + 1];
	struct position xIn, xOut;
	int nPhi, m, n, r, d, out;
	set_position( pb, &xIn );
	set_position( pb, &xOut );
	for( n = 0; n < MAXPHI + 1; n++ )
		set_position( pb, &simplex[n] );
	nPhi = multiobj_archive[0].x.f.size; // Dimension of the objfunc space
	if( nArchive < nPhi + 1 ) return;
	if( nArchive < MAXARCHIVE ) return;
	if( nLocalSearchIter < nArchive ) return;
	tprintf( "Iter %i Eval %d: Local search from archive\n", iter, eval );
	nLocalSearchIter = 0;
	xIn.size = multiobj_archive[0].x.size;
	xOut.size = xIn.size;
	out = random_int( 0, 1 );
	for( n = 0; n < nArchive - nPhi; n++ ) // Define a simplex
	{
		simplex[0] = multiobj_archive[n].x; // First element
		m = 0;
		for( r = 0; r < nArchive; r++ )
		{
			if( r == n ) continue;
			dR[m].dist = objfunc_dist( &multiobj_archive[n].x.f, &multiobj_archive[m].x.f ); // Compute the distances to the others
			dR[m].rank = r;
			m++;;
		}
		nPhi = multiobj_archive[0].x.f.size; // Find the nPhi nearest ones in the archive in order to complete the simplex
		qsort( dR, nArchive - 1, sizeof( dR[0] ), compare_dist_rank );
		for( m = 0; m < nPhi; m++ )
			simplex[m + 1] = multiobj_archive[dR[m].rank].x;
		// Define a new point
		//out=aleaInteger(0,1); // TO TRY
		//out=1;
		switch( out )
		{
			case 0: // Inside the simplex
				for( d = 0; d < xIn.size; d++ )
				{
					xIn.x[d] = 0;
					for( m = 0; m < nPhi + 1; m++ )
						xIn.x[d] += simplex[m].x[d];
				}
				for( d = 0; d < xIn.size; d++ )
					xIn.x[d] /= ( nPhi + 1 );
				position_check( pb, &xIn );
				position_eval( pb, &xIn );
				position_archive( &xIn );
				break;
			case 1: // Outside the simplex
				for( d = 0; d < xOut.size; d++ )
				{
					xOut.x[d] = 0;
					for( m = 1; m < nPhi + 1; m++ )
						xOut.x[d] += simplex[m].x[d]; // Partial simplex, dimension D-1
				}
				for( d = 0; d < xOut.size; d++ ) xOut.x[d] /= nPhi;
				for( d = 0; d < xOut.size; d++ ) xOut.x[d] -= ( simplex[0].x[d] - xOut.x[d] );
				position_check( pb, &xOut );
				position_eval( pb, &xOut );
				position_archive( &xOut );
				break;
		}
		if( arch == 1 ) // If "good" position ..
		{
			// ... there is nevertheless a small probability
			// that the "curvature" is not the same for the next point.
			// and that it is interesting to try the other case
			if( random_double( 0, 1 ) < 1 / ( nPhi + 1 ) ) out = 1 - out;
		}
		if( arch == 0 ) // If "bad" position ...
		{
			// ... there is a high probability that
			// the "curvature" is the same for the next point
			// and that is interesting to try the other case
			if( random_double( 0, 1 ) < nPhi / ( nPhi + 1 ) ) out = 1 - out;
		}
	}
	free_position( &xIn );
	free_position( &xOut );
	for( n = 0; n < MAXPHI + 1; n++ )
		free_position( &simplex[n] );
}

void archive_save( struct archived archiv[], int nArchive, FILE *fArchive )
{
	int D, nPhi, m, n;
	tprintf( "Save archive (%i positions)\n", nArchive );
	D = archiv[0].x.size;
	nPhi = archiv[0].x.f.size;
	tprintf( "%i objfunces, and %i coordinates\n", nPhi, D );
	for( n = 0; n < nPhi; n++ ) fprintf( fArchive, "f%i ", n + 1 );
	for( n = 0; n < D; n++ ) fprintf( fArchive, "x%i ", n + 1 );
	fprintf( fArchive, "\n" );
	for( n = 0; n < nArchive; n++ )
	{
		for( m = 0; m < nPhi; m++ ) fprintf( fArchive, "%g ", archiv[n].x.f.f[m] );
		for( m = 0; m < D; m++ ) fprintf( fArchive, "%g ", archiv[n].x.x[m] );
		fprintf( fArchive, "\n" );
	}
}

double archive_spread()
{
	double dMin[MAXARCHIVE], diversity1, diversity2, dMean1, dMean2, phi_spread, z, d;
	int m, n;
	if( nArchive < 2 ) return HUGE_VAL;
	// For each position, find the  nearest one
	// and the farest one, in the objfunc space
	// Without the new position
	diversity1 = 0;
	for( n = 0; n < nArchive - 1; n++ )
	{
		dMin[n] = HUGE_VAL;
		for( m = 0; m < nArchive - 1; m++ )
		{
			if( m == n ) continue;
			d = objfunc_dist( &multiobj_archive[n].x.f, &multiobj_archive[m].x.f );
			if( d < dMin[n] ) dMin[n] = d;
		}
		diversity1 += dMin[n];
	}
	dMean1 = diversity1 / ( nArchive - 1 );
	diversity2 = 0;
	for( n = 0; n < nArchive; n++ )
	{
		dMin[n] = HUGE_VAL;
		for( m = 0; m < nArchive; m++ )
		{
			if( m == n ) continue;
			d = objfunc_dist( &multiobj_archive[n].x.f, &multiobj_archive[m].x.f );
			if( d < dMin[n] ) dMin[n] = d;
		}
		diversity2 += dMin[n];
	}
	dMean2 = diversity2 / nArchive;
	phi_spread = 0; // "phi_spread" estimation
	for( n = 0; n < nArchive; n++ )
	{
		z = dMean2 - dMin[n];
		phi_spread += z * z;
	}
	phi_spread = sqrt( phi_spread ); // Initial value
	// Take distances (in the objfunc space) into account
	d = objfunc_dist( &multiobj_archive[0].x.f, &multiobj_archive[nArchive - 1].x.f ); // Distance between the new position and the first archived
	for( n = 1; n < nArchive - 1; n++ )
	{
		z = objfunc_dist( &multiobj_archive[n].x.f, &multiobj_archive[nArchive - 1].x.f ); // Distance min to the others
		if( z < d ) d = z;
	}
	if( d < dMean1 ) phi_spread += dMean1; // Penalty
	return phi_spread;
}

int adapt_compare( struct objfunc phi1[], int runs, int compare_type )
{
	/*
	 Modify the type of comparison, if multiobjective
	 Possible methods:
	 - adaptive weights
	 - epsilon-dominance
	 */
	double mean[MAXPHI], var[MAXPHI];
	double maxdf, mindf;
	double z;
	int fSame[MAXPHI];
	int nPhi;
	int new_compare_type;
	int n;
	int same;
	int r;
	if( runs == 0 ) return compare_type; // Not possible for the first run
	nPhi = phi_list[0].size;
	if( nPhi == 1 ) return compare_type; // Nothing to do for mono-objective
	for( n = 0; n < nPhi; n++ ) // Is the diversity big enough?
	{
		mean[n] = 0;
		for( r = 0; r <= runs; r++ )
			mean[n] += phi_list[r].f[n];
		mean[n] /= ( runs + 1 );
	}
	for( n = 0; n < nPhi; n++ )
	{
		var[n] = 0;
		for( r = 0; r <= runs; r++ )
			var[n] = pow( phi_list[r].f[n] - mean[n], 2 );
		var[n] = sqrt( var[n] / ( runs + 1 ) ); // Standard deviation
	}
	for( n = 0; n < nPhi; n++ )
	{
		mindf = mean[n] - 3 * var[n];
		maxdf = mean[n] + 3 * var[n];
		if( mindf < phi_min[n] ) {fSame[n] = 0; continue;}
		if( maxdf > phi_max[n] ) {fSame[n] = 0; continue;}
		fSame[n] = 0;
		z = random_double( phi_min[n], mean[n] ); // Fuzzy rule
		if( z < mindf ) fSame[n] = 1;
		z = random_double( mean[n], phi_max[n] ); // Fuzzy rule
		if( z > maxdf ) fSame[n] = 1;
	}
	same = 0;
	for( n = 0; n < nPhi; n++ )
		if( fSame[n] == 1 ) same++;
	if( same == nPhi ) new_compare_type = 1 - compare_type; else new_compare_type = compare_type; // Either 0 or 1
	return new_compare_type;
}

static int compare_crowding_dist( void const *a, void const *b ) // QSORT: compare crowding distances in the archive
{
	struct archived const *pa = ( struct archived const * ) a;
	struct archived const *pb = ( struct archived const * ) b;
	if( pa->crow_dist > pb->crow_dist ) return 1;
	if( pa->crow_dist < pb->crow_dist ) return -1;
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

void position_check( struct problem *pb, struct position *P )
{
	int rank, d, n;
	for( d = 0; d < ( *pb ).D; d++ ) // Granularity
		if( ( *pb ).valSize[d] == 0 )
			if( ( *pb ).dx[d] > 0 )
				( *P ).x[d] = granularity( ( *P ).x[d], ( *pb ).dx[d] );
	for( d = 0; d < ( *pb ).D; d++ ) // Clamping
	{
		if( ( *pb ).valSize[d] == 0 ) // Variable Interval
		{
			if( ( *P ).x[d] > ( *pb ).max[d] )( *P ).x[d] = ( *pb ).max[d];
			else if( ( *P ).x[d] < ( *pb ).min[d] )( *P ).x[d] = ( *pb ).min[d];
		}
		else // Variable List: Find the nearest one
		{
			rank = 0;
			for( n = 1; n < ( *pb ).valSize[d]; n++ )
				if( fabs( ( *pb ).val[d][n] - ( *P ).x[d] ) < fabs( ( *pb ).val[d][rank] - ( *P ).x[d] ) ) rank = n;
			( *P ).x[d] = ( *pb ).val[d][rank];
		}
	}
}

int compare_particles( struct objfunc *phi1, struct objfunc *phi2, int compare_type )
{
	int n, better = 0, worse = 0;
	double phi_weights1 = 0, phi_weights2 = 0;
	switch( compare_type )
	{
		case 1: // Using penalties
			for( n = 0; n < ( *phi1 ).size; n++ )
			{
				phi_weights1 += phi_weights[n] * ( *phi1 ).f[n];
				phi_weights2 += phi_weights[n] * ( *phi2 ).f[n];
			}
			if( phi_weights1 < phi_weights2 - DBL_EPSILON ) return  1; // better
			if( phi_weights1 > phi_weights2 + DBL_EPSILON ) return -1; // worse
			return 0;                                                   // no change
		case 0:
		default: // Epsilon-dominance
			for( n = 0; n < ( *phi1 ).size; n++ )
			{
				if( ( *phi1 ).f[n] < ( *phi2 ).f[n] - epsilon_vector[n] ) { better++; continue; }
				if( ( *phi1 ).f[n] > ( *phi2 ).f[n] + epsilon_vector[n] ) { worse++; continue; }
			}
//			tprintf( "Compare %d %d: %f %f %f\n", better, worse, ( *phi1 ).f[0], ( *phi2 ).f[0], epsilon_vector[0] );
			if( better >  0 && worse == 0 ) return  1; // better
			if( better == 0 && worse >  0 ) return -1; // worse
			if( !multiObj ) return 0;                  // no change
			// Multi-objective difficult to decide. Use the "spread" criterion
			// Contain the "phi_spread" value that should be as small as possible
			n = ( *phi1 ).size;
			if( ( *phi1 ).f[n] < ( *phi2 ).f[n] - DBL_EPSILON ) return  1; // better
			if( ( *phi1 ).f[n] > ( *phi2 ).f[n] + DBL_EPSILON ) return -1; // worse
			return 0;                                       // no change
		case 2: // Pure dominance
			for( n = 0; n < ( *phi1 ).size; n++ )
			{
				if( ( *phi1 ).f[n] < ( *phi2 ).f[n] ) { better++; continue; }
				if( ( *phi1 ).f[n] > ( *phi2 ).f[n] ) { worse++; continue; }
			}
			if( better >  0 && worse == 0 ) return  1; // better
			if( better == 0 && worse >  0 ) return -1; // worse
			return 0;                                  // no change
		case 3: // for maxError check; OF below the cutoff value
			for( n = 0; n < ( *phi1 ).size; n++ )
			{
				if( ( *phi1 ).f[n] < ( *phi2 ).f[n] + epsilon_vector[n] ) { better++; continue; }
				if( ( *phi1 ).f[n] > ( *phi2 ).f[n] + epsilon_vector[n] ) { worse++; continue; }
			}
			if( better >  0 && worse == 0 ) return  1; // better
			if( better == 0 && worse >  0 ) return -1; // worse
			if( !multiObj ) return 0;                  // no change
			// Multi-objective difficult to decide. Use the "spread" criterion
			// Contain the "phi_spread" value that should be as small as possible
			n = ( *phi1 ).size;
			if( ( *phi1 ).f[n] < ( *phi2 ).f[n] + DBL_EPSILON ) return  1; // better
			if( ( *phi1 ).f[n] > ( *phi2 ).f[n] + DBL_EPSILON ) return -1; // worse
			return 0;                                       // no change
	}
}

void objfunc_print( struct objfunc *phi )
{
	int i;
	if( ( *phi ).size == 1 )
	{
		tprintf( "OF: %g", ( *phi ).f[0] );
		return;
	}
	tprintf( "OF (%d): ", ( *phi ).size );
	for( i = 0; i < ( *phi ).size; i++ )
		tprintf( " %g", ( *phi ).f[i] );
	if( multiObj ) tprintf( " phi_spread %g\n", ( *phi ).f[( *phi ).size] );
	else tprintf( "\n" );
}

double objfunc_dist( struct objfunc *phi1,  struct objfunc *phi2 )
{
	double t, dist = 0;
	int i;
	for( i = 0; i < ( *phi1 ).size; i++ )
	{
		t = ( *phi1 ).f[i] - ( *phi2 ).f[i];
		dist += t * t;
	}
	return sqrt( dist );
}

double objfunc_total( struct objfunc *phi, int type ) // Total objfunc (weighted if multi-objective)
{
	double error = 0;
	int i;
	if( ( *phi ).size <= 1 ) return ( *phi ).f[0];
	switch( type )
	{
		default:
		case 0:
			for( i = 0; i < ( *phi ).size; i++ )
				error += ( *phi ).f[i] * ( *phi ).f[i];
			error = sqrt( error );
			break;
		case 1:
			for( i = 0; i < ( *phi ).size; i++ )
				error += phi_weights[i] * fabs( ( *phi ).f[i] ); // Weighted (phi_weights is a global variable)
			break;
	}
	return error;
}

double granularity( double value, double granul ) // Modify a value according to a granularity
{
	double x;
	if( granul < DBL_EPSILON ) return value; // Pseudo-continuous (depending on the machine)
	if( value >= 0 ) x = granul * floor( value / granul + 0.5 );
	else x = ( -granul * floor( -value / granul + 0.5 ) );
	return x;
}

double maxXY( double x, double y )
{
	if( x > y ) return x; else return y;
}

double minXY( double x, double y )
{
	if( x < y ) return x; else return y;
}

void particle_init( struct problem *pb, int initOption, struct position *guide1, struct position *guide2, struct swarm *S, struct particle *P )
{
	double mean, range, sort_vec[2 * MAXTRIBE * MAXPART], min, max, x;
	int i, ip, it, k, count, rank, tr_best;
	( *P ).x.size = ( *pb ).D; // Initial number of particles equal number of dimensions (optimized variables)
	switch( initOption )
	{
		default:
		case 0: // Random position
			for( k = 0; k < ( *pb ).D; k++ )
				( *P ).x.x[k] = random_double( ( *pb ).min[k], ( *pb ).max[k] );
			break;
		case 1: // Guided
			for( k = 0; k < ( *pb ).D; k++ )
			{
				if( ( *pb ).valSize[k] == 0 ) // variable defined over interval
				{
					range = fabs( ( *guide1 ).x[k] - ( *guide2 ).x[k] );
					( *P ).x.x[k] = random_double( ( *guide1 ).x[k] - range, ( *guide1 ).x[k] + range ); // TO TRY: random_gauss instead)
				}
				else // List
				{
					i = random_int( 0, ( *pb ).valSize[k] - 1 ); // WARNING: still purely random, not very satisfying
					( *P ).x.x[k] = ( *pb ).val[k][i];
				}
			}
			break;
		case 2: // On the corners
			for( k = 0; k < ( *pb ).D; k++ )
			{
				if( random_double( 0, 1 ) < 0.5 )( *P ).x.x[k] = ( *pb ).min[k];
				else( *P ).x.x[k] = ( *pb ).max[k];
			}
			break;
		case 3: // Biggest empty hyper-parallelepiped (Binary Search)
			// TODO: TO TRY for multi-objective, one may use the archive instead of (*S) as the list of known positions
			for( k = 0; k < ( *pb ).D; k++ )
			{
				sort_vec[0] = ( *pb ).min[k]; // List of known coordinates
				sort_vec[1] = ( *pb ).max[k]; // List of known coordinates
				count = 2;
				for( it = 0; it < ( *S ).size; it++ )
					for( ip = 0; ip < ( *S ).trib[it].size; ip++ )
					{
						sort_vec[count++] = ( *S ).trib[it].part[ip].x.x[k]; // List of known coordinates
						sort_vec[count++] = ( *S ).trib[it].part[ip].xBest.x[k]; // List of known coordinates
					}
				if( count > 1 ) qsort( sort_vec, count, sizeof( double ), compare_double ); // Sort the list of known coordinates
				rank = 1;
				for( i = 2; i < count; i++ ) // Find the biggest empty interval
					if( sort_vec[i] - sort_vec[i - 1] > sort_vec[rank] - sort_vec[rank - 1] ) rank = i;
				( *P ).x.x[k] = random_double( sort_vec[rank - 1], sort_vec[rank] ); // Select a random position "centered" on the interval
			}
			break;
		case 4: // Centered in the model domain EXPERIMENT NOT USED AT THE MOMENT
			for( k = 0; k < ( *pb ).D; k++ )
			{
				mean = random_double( ( *pb ).min[k], ( *pb ).max[k] );
				range = ( 1. / 3 ) * maxXY( ( *pb ).max[k] - mean, mean - ( *pb ).min[k] );
				( *P ).x.x[k] = random_gaussian( mean, range );
			}
			break;
		case 5: // User supplied initial values
			for( k = 0; k < ( *pb ).D; k ++ )
				( *P ).x.x[k] = ( *pb ).ival[k];
			break;
		case 6: // Centered to all the shamans in the swarm (applied by LM speedup to reset bad shamans)
			for( k = 0; k < ( *pb ).D; k++ )
			{
				mean = 0;
				min = HUGE_VAL;
				max = 0;
				count = 0;
				for( it = 0; it < ( *S ).size; it++ )
					for( ip = 0; ip < ( *S ).trib[it].size; ip++ )
					{
						tr_best = ( *S ).trib[it].best;
						x = ( *S ).trib[it].part[tr_best].x.x[k]; // List of known coordinates
						mean += x;
						count++;
					}
				if( count > 1 ) mean /= count;
				( *P ).x.x[k] = mean;
			}
			break;
	}
	position_check( pb, &( *P ).x );
	position_eval( pb, &( *P ).x );
	if( multiObj ) position_archive( &( *P ).x ); // Archive (multiobjective)
	for( k = 0; k < ( *pb ).D; k++ )
		( *P ).xLast.x[k] = ( *P ).xBest.x[k] = ( *P ).x.x[k];
	for( k = 0; k < ( *pb ).nPhi; k++ )
		( *P ).fBestPrev.f[k] = ( *P ).xLast.f.f[k] = ( *P ).xBest.f.f[k] = ( *P ).x.f.f[k];
	( *P ).id = ++id_global;
	( *P ).strategy = 0; //TO TRY: aleaInteger(0,2);
}

void position_archive( struct position *pos )
{
	int archiveStore;
	int i, j, cmp;
	int dominanceType = 2; // 2 => pure-dominance
	arch = 1; // Global variable
	if( nArchive > 0 ) // It is not the first one. Check if dominated
		for( i = 0; i < nArchive; i++ )
		{
			cmp = compare_particles( &multiobj_archive[i].x.f, &( *pos ).f, dominanceType );
			if( cmp == 1 )
			{
				arch = 0; // Dominated, don't keep it
				break;
			}
		}
	if( arch == 1 )
	{
		nLocalSearchIter++; // Useful to decide when to perform local search
		if( nArchive > 1 ) // Remove the dominated positions
			for( i = 0; i < nArchive; i++ )
			{
				cmp = compare_particles( &( *pos ).f, &multiobj_archive[i].x.f, dominanceType );
				if( cmp == 1 )
				{
					if( i < nArchive - 1 )
						for( j = i; j < nArchive - 1; j++ )
							multiobj_archive[j] = multiobj_archive[j + 1];
					nArchive--;
				}
			}
		// TO TRY: one may also remove one of two "too similar" positions
		if( nArchive < MAXARCHIVE ) archiveStore = nArchive; // (*S)ore the position
		else
		{
			nArchive = MAXARCHIVE;
			archiveStore = 0;
			for( i = 1; i < nArchive; i++ ) // Find the most "crowded" archived position
				if( multiobj_archive[i].crow_dist < multiobj_archive[archiveStore].crow_dist ) archiveStore = i;
		}
		multiobj_archive[archiveStore].x = *pos;
		if( nArchive < MAXARCHIVE ) nArchive++;
	}
}

void position_print( struct position *pos )
{
	int d;
	tprintf( "position:" );
	for( d = 0; d < ( *pos ).size; d++ )
		tprintf( " %g", ( *pos ).x[d] );
	tprintf( " | OF:" );
	for( d = 0; d < ( *pos ).f.size; d++ )
		tprintf( " %g", ( *pos ).f.f[d] );
	if( multiObj ) tprintf( " phi_spread %g", ( *pos ).f.f[( *pos ).f.size] );
}

void position_save( FILE *fRun, struct position *pos, int run )
{
	int d;
	fprintf( fRun, "%i %i %d ", run, iter, eval );
	for( d = 0; d < ( *pos ).f.size; d++ )
		fprintf( fRun, "%g ", ( *pos ).f.f[d] ); // Fitness
	for( d = 0; d < ( *pos ).size; d++ )
		fprintf( fRun, "%g ", ( *pos ).x[d] ); // Position
	fprintf( fRun, "\n" );
}

void position_update( struct problem *pb, struct particle *P, struct particle informer )
{
	double c1, c2, c3 = 0;
	double error1 = 0, error2 = 0, noise, dx, radius, x1, x2;
	int d;
	int Dim;
	int type;
	double w1 = 0.74; // => gaussian fifty-fifty
	double w2 = 0.72; // < 1/(2*ln(2)) Just for tests
	double c = 1.19; // < 0.5 + ln(2) Just for tests
	( *P ).x.size = ( *pb ).D;
	Dim = informer.xBest.size; // May be zero if there is no informer
	switch( ( *P ).strategy )
	{
		case -1: // Random
			for( d = 0; d < ( *pb ).D; d++ )
				( *P ).x.x[d] = random_double( ( *pb ).min[d], ( *pb ).max[d] ); // Try something within the range
			( *P ).strategy = 0;
			break;
		default:
		case 0: // Improved Bare bones
			for( d = 0; d < Dim; d++ )
			{
				dx = w1 * fabs( informer.xBest.x[d] - ( *P ).xBest.x[d] );
				if( dx < DBL_EPSILON ) // informer = current particle; i.e. no informer
				{
					dx = maxXY( fabs( ( *P ).xBest.x[d] - ( *pb ).min[d] ), fabs( ( *P ).xBest.x[d] - ( *pb ).max[d] ) );
					( *P ).x.x[d] = random_gaussian( informer.xBest.x[d], dx ); // informer.xBest.x[d] = ( *P ).xBest.x[d]
				}
				else // there is an informer
					( *P ).x.x[d] = random_gaussian( informer.xBest.x[d], dx );
			}
			break;
		case 1: // Pivot by dimension
		case 2: // Pivot by dimension + noise
			if( multiObj ) type = 1; else type = 0;
			if( Dim > 0 ) // If there is an informer
			{
				error1 = objfunc_total( &( *P ).xBest.f, type );
				error2 = objfunc_total( &informer.xBest.f, type );
				c3 = error1 + error2;
				if( c3 > DBL_EPSILON ) c1 = error2 / c3;
				else                   c1 = random_double( 0, 1 );
				c2 = 1.0 - c1;
				for( d = 0; d < Dim; d++ )
				{
					radius = fabs( ( *P ).xBest.x[d] - informer.xBest.x[d] );
					x1 = random_double( ( *P ).xBest.x[d] - radius, ( *P ).xBest.x[d] + radius );
					x2 = random_double( informer.xBest.x[d] - radius, informer.xBest.x[d] + radius );
					( *P ).x.x[d] = c1 * x1 + c2 * x2;
				}
			}
			else // No other informer than itself
			{
				for( d = 0; d < ( *P ).xBest.size; d++ )
				{
					radius = maxXY( fabs( ( *P ).xBest.x[d] - ( *pb ).min[d] ), fabs( ( *P ).xBest.x[d] - ( *pb ).max[d] ) );
					dx = ( double ) 3.0 * radius;
					( *P ).x.x[d] = random_gaussian( ( *P ).xBest.x[d], dx );
				}
				break;
			}
			if( ( *P ).strategy == 1 ) break; // strategy 1
			if( c3 > DBL_EPSILON ) noise = random_gaussian( 0, ( error1 - error2 ) / c3 ); // strategy 2
			else                   noise = random_gaussian( 0, 1 );
			for( d = 0; d < ( *pb ).D; d++ )
				( *P ).x.x[d] += noise * ( ( *P ).xBest.x[d] - ( *P ).x.x[d] ); // Add noise
			break;
		case 99: // Standard PSO updating for testing
			for( d = 0; d < Dim; d++ )
			{
				dx = w2 * ( ( *P ).x.x[d] - ( *P ).xLast.x[d] );
				dx += random_double( 0, c ) * ( ( *P ).xBest.x[d] - ( *P ).x.x[d] ) + random_double( 0, c ) * ( informer.xBest.x[d] - ( *P ).x.x[d] );
				( *P ).x.x[d] = ( *P ).xLast.x[d] + dx;
			}
			break;
	}
	position_check( pb, &( *P ).x );
	position_eval( pb, &( *P ).x );
	if( multiObj ) position_archive( &( *P ).x ); // Archive
}

void position_lm( struct opt_data *op, struct problem *pb, struct position *P )
{
	int d;
//	for( d = 0; d < ( *pb ).D; d++ )
//		tprintf( "lm %g\n", ( *P ).x[d] );
	for( d = 0; d < ( *pb ).D; d++ )
		op->pd->var[op->pd->var_index[d]] = ( *P ).x[d];
//	DeTransform( ( *P ).x, op, ( *P ).x ); // there is no need to detransform
//	for( d = 0; d < ( *pb ).D; d++ )
//		tprintf( "lm %g\n", ( *P ).x[d] );
	d = op->cd->neval;
	op->phi = ( *P ).f.f[0];
	optimize_lm( op );
	eval += op->cd->neval - d; // add the number of evaluations performed within LM
	if( op->cd->sintrans )
	{
		for( d = 0; d < ( *pb ).D; d++ )
			( *P ).x[d] = asin( sin( op->pd->var[op->pd->var_index[d]] ) ); // keep the estimates within the initial range ...
		( *P ).f.f[0] = op->phi;
	}
	else
	{
		for( d = 0; d < ( *pb ).D; d++ )
			( *P ).x[d] = op->pd->var[op->pd->var_index[d]]; // transfer back the values
		if( op->cd->pardx > DBL_EPSILON )
		{
			position_check( pb, P );
			position_eval( pb, P );
		}
	}
}

void problem_print( struct problem *pb )
{
	int d, n;
	tprintf( "Number of parameters (parameter space dimensions) %i\n", ( *pb ).D );
	if( debug_level > 1 )
	{
		tprintf( "Search space:\n" );
		for( d = 0; d < ( *pb ).D; d++ )
		{
			if( ( *pb ).valSize[d] == 0 )
				tprintf( "[%g  %g] dx=%g init=%g\n", ( *pb ).min[d], ( *pb ).max[d], ( *pb ).dx[d], ( *pb ).ival[d] );
			else
			{
				tprintf( "List: " );
				for( n = 0; n < ( *pb ).valSize[d]; n++ )
					tprintf( "%g ", ( *pb ).val[d][n] );
				tprintf( "\n" );
			}
		}
	}
	if( ( *pb ).nPhi > 1 )
		for( n = 0; n < ( *pb ).nPhi; n++ )
			tprintf( "Objective function #%d: Target %g Acceptable error %g\n", n, ( *pb ).objective[n], ( *pb ).maxError.f[n] );
	else if( ( *pb ).objective[0] > 0 || ( *pb ).maxError.f[0] > 0 )
		tprintf( "Objective function: Target %g Acceptable error %g\n", ( *pb ).objective[0], ( *pb ).maxError.f[0] );
	tprintf( "Maximum number of functional evaluations: %d\n", ( *pb ).maxEval );
	if( ( *pb ).repeat != 1 ) tprintf( "Number of optimization retries: %i\n", ( *pb ).repeat );
}

void pso_solver( struct problem *pb, int compare_type, int run, struct swarm *S )
{
	int n, tr, shaman;
	double phi_current_best, min, max, r1, r2;
	for( n = 0; n < ( *pb ).nPhi; n++ ) phi_weights[n] = 1; // Initials penalties (for multiobjective)
	nLocalSearchIter = 0; // Prepare local search (for multiobjective)
	if( gop->cd->check_success ) gop->global_success = gop->success = 0;
	swarm_init( pb, compare_type, S ); //Initialization of the swarm
	if( ( gop->cd->check_success && gop->success ) || gop->global_success )
	{
		if( debug_level ) tprintf( "Success: Initial model predictions are within the predefined calibration bounds!\n" );
		return;
	}
	phi_lm_init = ( *S ).best.f.f[0]; // Store the first best phi
	if( lmo_flag )
	{
		min = HUGE_VAL; max = 0;
		for( tr = 0; tr < ( *S ).size; tr++ )
		{
			shaman = ( *S ).trib[tr].best;
			phi_current_best = ( *S ).trib[tr].part[shaman].xBest.f.f[0];
			if( phi_current_best > max ) max = phi_current_best;
			if( phi_current_best < min ) min = phi_current_best;
		}
		r1 = ( max - min ) / min;
		r2 = 1 + ( *pb ).lm_factor / 10;
		if( r1 < r2 )
		{
			lm_wait = 1;
			phi_lm_wait = min / ( *pb ).lm_factor * 2;
			if( debug_level ) tprintf( "Skip LM search because OF range of tribes' shamans is small till OF %g is less than %g (min %g max %g ratio %g < %g)!\n", min, phi_lm_wait, min, max, r1, r2 );
		}
		else if( debug_level ) tprintf( "Do LM search because OF range of tribes' shamans is sufficient (min %g max %g ratio %g >= %g)!\n", min, max, r1, r2 );
	}
	if( debug_level ) swarm_print( S );
	if( ( gop->cd->check_success && gop->success ) || gop->global_success )
	{
		if( debug_level ) tprintf( "Success: Initial model predictions are within the predefined calibration bounds!\n" );
		return;
	}
	iter = 0;
	nSwarmAdaptIter = 0; // Last iteration at which the swarm has been adapted
	for( tr = 0; tr < ( *S ).size; tr++ ) nTribeAdaptIter[tr] = 0;
	while( 1 )
	{
		for( tr = 0; tr < ( *S ).size; tr++ ) nTribeAdaptIter[tr]++;
		nSwarmAdaptIter++;
		archiveVar = archive_phi_var( pb );
		if( debug_level > 2 ) tprintf( "Swarm move ...\n" );
		swarm_move( pb, S, compare_type, run ); // SWARM MOVE
		if( debug_level > 2 )
		{
			tprintf( "MO " );
			for( tr = 0; tr < ( *S ).size; tr++ )
				tprintf( "T %d OF %g ", tr + 1, ( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f.f[0] );
			tprintf( ": Tbest %d OF %g TOF %g\n", ( *S ).tr_best + 1, ( *S ).trib[( *S ).tr_best].part[( *S ).trib[( *S ).tr_best].best].xBest.f.f[0], ( *S ).best.f.f[0] );
		}
		if( ( gop->cd->check_success && gop->success ) || gop->global_success ) break; // Success: Predictions are within the predefined calibration bounds
		if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "Maximum number of evaluations is achieved!\n" ); break; }
		if( debug_level > 2 ) tprintf( "Swarm adapt ...\n" );
		swarm_adapt( pb, S, compare_type ); // SWARM ADAPT
		if( debug_level > 2 )
		{
			tprintf( "AD " );
			for( tr = 0; tr < ( *S ).size; tr++ )
				tprintf( "T %d OF %g ", tr + 1, ( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f.f[0] );
			tprintf( ": Tbest %d OF %g TOF %g\n", ( *S ).tr_best + 1, ( *S ).trib[( *S ).tr_best].part[( *S ).trib[( *S ).tr_best].best].xBest.f.f[0], ( *S ).best.f.f[0] );
		}
		if( ( gop->cd->check_success && gop->success ) || gop->global_success ) break; // Success: Predictions are within the predefined calibration bounds
		if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "Maximum number of evaluations is achieved!\n" ); break; }
		if( gop->od->nTObs > 0 && lmo_flag )
		{
			if( debug_level > 2 ) tprintf( "LM search ...\n" );
			swarm_lm( pb, S ); // LM SEARCH
			if( debug_level > 2 )
			{
				tprintf( "LM " );
				for( tr = 0; tr < ( *S ).size; tr++ )
					tprintf( "T %d OF %g ", tr + 1, ( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f.f[0] );
				tprintf( ": Tbest %d OF %g TOF %g\n", ( *S ).tr_best + 1, ( *S ).trib[( *S ).tr_best].part[( *S ).trib[( *S ).tr_best].best].xBest.f.f[0], ( *S ).best.f.f[0] );
			}
		}
		if( !multiObj && compare_particles( &( *S ).best.f, &( *pb ).maxError, 3 ) == 1 ) { if( debug_level ) tprintf( "Success: OF is minimized below the cutoff value! (%g<%g)\n", ( *S ).best.f.f[0], ( *pb ).maxError.f[0] ); break; }
		if( ( gop->cd->check_success && gop->success ) || gop->global_success ) break; // Success: Predictions are within the predefined calibration bounds
		if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "Maximum number of evaluations is achieved!\n" ); break; }
		if( debug_level > 2 ) tprintf( "Local search ...\n" );
		if( multiObj ) archive_local_search( pb ); else swarm_local_search( pb, S ); // LOCAL SEARCH
		if( debug_level > 2 )
		{
			tprintf( "LO " );
			for( tr = 0; tr < ( *S ).size; tr++ )
				tprintf( "T %d OF %g ", tr + 1, ( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f.f[0] );
			tprintf( ": Tbest %d OF %g TOF %g\n", ( *S ).tr_best + 1, ( *S ).trib[( *S ).tr_best].part[( *S ).trib[( *S ).tr_best].best].xBest.f.f[0], ( *S ).best.f.f[0] );
		}
		//if( restart ) break; // For future automatic restart
		if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "Maximum number of evaluations is achieved!\n" ); break; }
		if( !multiObj && compare_particles( &( *S ).best.f, &( *pb ).maxError, 3 ) == 1 ) { if( debug_level ) tprintf( "Success: OF is minimized below the cutoff value! (%g<%g)\n", ( *S ).best.f.f[0], ( *pb ).maxError.f[0] ); break; }
		if( ( gop->cd->check_success && gop->success ) || gop->global_success ) break; // Success: Predictions are within the predefined calibration bounds
		iter++;
		if( debug_level > 2 ) tprintf( "\n" );
	}
	if( ( gop->cd->check_success && gop->success ) || gop->global_success )
	{
		if( debug_level ) tprintf( "PSO Success: Predictions are within the predefined calibration bounds!\n" );
		if( compare_particles( &( *S ).best.f, &( *pb ).pos_success.f, 0 ) == 1 )
			if( debug_level ) tprintf( "Parameter space location producing predictions within the predefined calibration ranges does not have the lowest OF!\n" );
		/*
				tprintf( "ggggggg\n" );
				for( n = 0; n < ( *pb ).pos_success.size; n++ )
					tprintf( "%g\n", ( *pb ).pos_success.x[n] );
				func( ( *pb ).pos_success.x, gop, res ); // evaluate the best BEST result
		*/
		copy_position( &( *pb ).pos_success, &( *S ).best );
	}
}

int sign( double x )
{
	if( x > 0 ) return 1;
	if( x < 0 ) return -1;
	return 0;
}

void swarm_adapt( struct problem *pb, struct swarm( *S ), int compare_type )
{
	double pr;
	int adaptSwarm;
	int add_part_count = 0;
	int add_tribe_count = 0;
	int d;
	int disturbPart;
	int improve_count;
	int initOption;
	int pa;
	int nPart;
	int part_worst;
	int s;
	int status;
	int i, tr;
	int nTribe;
	int tribe_worst;
	nTribe = ( *S ).size;
	for( tr = 0; tr < ( *S ).size; tr++ )
	{
		nPart = ( *S ).trib[tr].size;
		nTribeAdaptIter[tr] = 0;
		if( ( *S ).trib[tr].status == -1 )( *S ).trib[tr].status = 0; // Very bad tribe
		else
		{
			improve_count = 0;
			for( pa = 0; pa < ( *S ).trib[tr].size; pa++ ) // Status of the the tribe
			{
				if( compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).trib[tr].part[pa].xLast.f, compare_type ) == 1 )
					improve_count++;
			}
			pr = ( double ) improve_count / ( *S ).trib[tr].size; // Fuzzy rule
			if( random_double( 0, 1 ) > pr )( *S ).trib[tr].status = 0; // Bad tribe
			else( *S ).trib[tr].status = 1;                             // Good tribe
		}
		switch( ( *S ).trib[tr].status )
		{
			case 0: // Bad tribe.
				/* The idea is to increase diversity. This is done by adding sometimes a completely new particle
				and by disturbing another one a bit, typically along just one dimension. */
				disturbPart = ( double ) random_double( 0, 1 ) < ( ( double ) -1 / ( *S ).trib[tr].size + 1 );
				if( disturbPart )
					if( ( *S ).trib[tr].size > 1 )
					{
						if( debug_level > 2 ) tprintf( "Bad tribe %i: one randomly selected particle is randomly disturbed ", tr + 1 );
						do { pa = random_int( 0, ( *S ).trib[tr].size - 1 ); }
						while( pa == ( *S ).trib[tr].best ); // Do not disturb the best one ...
						d = random_int( 0, ( *pb ).D - 1 );
						// d = tribeVarianceMin(St.trib[tr]); // EXPERIMENT
						( *S ).trib[tr].part[pa].x.x[d] = random_double( ( *pb ).min[d], ( *pb ).max[d] );
						position_check( pb, &( *S ).trib[tr].part[pa].x );
						position_eval( pb, &( *S ).trib[tr].part[pa].x );
						// tprintf( "mmm %i %i\n", pa, ( *S ).trib[tr].best );
						// tribe_print( &( *S ).trib[tr] );
						if( compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, compare_type ) == 1 ) // Possibly update tribe best
						{
							( *S ).trib[tr].best = pa; // By chance this particle is the best of the tribe
							copy_position( &( *S ).trib[tr].part[pa].x, &( *S ).trib[tr].part[pa].xBest ); // copy to tribes best position
							if( debug_level > 2 )
							{
								tprintf( "best particle in tribe #%d is now #%i ", tr + 1, ( *S ).trib[tr].part[pa].id );
								objfunc_print( &( *S ).trib[tr].part[pa].x.f );
								tprintf( "\n" );
							}
						}
						else if( debug_level > 2 ) tprintf( "\n" );
						if( compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).best.f, compare_type ) == 1 ) // Possibly update swarm best
						{
							copy_position( &( *S ).trib[tr].part[pa].x, &( *S ).best ); // BEST COPY
							( *S ).tr_best = tr;
						}
					}
				if( nPart < MAXPART )
				{
					if( multiObj ) initOption = 3;
					else initOption = random_int( 0, 3 );
					if( debug_level > 2 ) tprintf( "Bad tribe %i: new particle is added ", tr + 1 );
					add_part_count++;
					particle_init( pb, initOption, &( *S ).best, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest, S, &( *S ).trib[tr].part[nPart] );
					( *S ).trib[tr].size++; // Add a new particle
					// tprintf( "mmm %i %i\n ", initOption, nPart);
					if( compare_particles( &( *S ).trib[tr].part[nPart].x.f, &( *S ).trib[tr].part[nPart].xBest.f, compare_type ) == 1 ) // Possibly update particles best
						copy_position( &( *S ).trib[tr].part[nPart].x, &( *S ).trib[tr].part[nPart].xBest ); // copy to particles best position
					if( compare_particles( &( *S ).trib[tr].part[nPart].x.f, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, compare_type ) == 1 ) // Possibly update tribe best
					{
						( *S ).trib[tr].best = nPart; // By chance this particle is the best of the tribe
						if( debug_level > 2 )
						{
							tprintf( "best particle in tribe #%d is now #%i ", tr + 1, ( *S ).trib[tr].part[nPart].id );
							objfunc_print( &( *S ).trib[tr].part[nPart].x.f );
							tprintf( "\n" );
						}
					}
					else if( debug_level > 2 ) tprintf( "\n" );
					if( compare_particles( &( *S ).trib[tr].part[nPart].x.f, &( *S ).best.f, compare_type ) == 1 ) // Possibly update swarm best
					{
						copy_position( &( *S ).trib[tr].part[nPart].x, &( *S ).best ); // By chance this particle is the best of the swarm; BEST COPY
						( *S ).tr_best = tr;
					}
					nPart++; // Change the total internal count because a new particle is added
				}
				else
				{
					nExceedSizeTribe++;
					if( debug_level > 0 ) tprintf( "\nWARNING: Cannot add a particle (increase MAXPART = %i)\n", MAXPART );
				}
				break;
			case 1: // Good tribe
				if( ( *S ).trib[tr].size < 2 ) break;
				part_worst = 0;
				for( pa = 1; pa < ( *S ).trib[tr].size; pa++ ) // Find the worst particle
				{
					if( compare_particles( &( *S ).trib[tr].part[pa].xBest.f, &( *S ).trib[tr].part[part_worst].xBest.f, compare_type ) == -1 )
						part_worst = pa;
				}
				if( part_worst == ( *S ).trib[tr].best ) break; // It might be also the best. In that case, don't remove it
				if( part_worst < ( *S ).trib[tr].size - 1 ) // Remove it from the tribe
					for( pa = part_worst; pa < ( *S ).trib[tr].size - 1; pa++ )
						copy_particle( &( *S ).trib[tr].part[pa + 1], &( *S ).trib[tr].part[pa] );
				if( debug_level > 2 ) tprintf( "Good tribe %i: worst particle is removed %i\n", tr + 1, part_worst + 1 );
				( *S ).trib[tr].size--;
				break;
		}
	}
	if( ( *S ).size > 1 ) adaptSwarm = ( nSwarmAdaptIter >= ( *S ).size * ( ( *S ).size - 1 ) / 4 );
	else                  adaptSwarm = TRUE; // Always true if there is just one tribe
	if( count_not_good_tribes > ( *S ).size / 3 ) { ( *S ).status = 0; adaptSwarm = TRUE; } // Bad swarm after nStagnantLMiterations LM attempts
	if( adaptSwarm ) // Reinitialize the previous best, and the iteration counter
	{
//		Original TRIBES algorithm is updating the previous best before the test; the swarm is always bad!
//		for( i = 0; i < ( *pb ).nPhi; i++ )
//			( *S ).fBestPrev.f[i] = ( *S ).best.f.f[i];
		nSwarmAdaptIter = 0;
		if( ( *S ).status == -1 )( *S ).status = 0; // Very bad swarm according to LM
		else
		{
			if( compare_particles( &( *S ).best.f, &( *S ).fBestPrev, compare_type ) <= 0 )( *S ).status = 0; // Bad swarm
			else( *S ).status = 1; // Good swarm
		}
		// update after the swarm is defined good or bad
		for( i = 0; i < ( *pb ).nPhi; i++ )
			( *S ).fBestPrev.f[i] = ( *S ).best.f.f[i];
		switch( ( *S ).status )
		{
			case 0: // Bad swarm
				if( ( *S ).size < MAXTRIBE )
				{
					if( debug_level > 2 ) tprintf( "Bad swarm: tribe is added\n" );
					tribe_init( pb, 1, compare_type, S, &( *S ).trib[tr] );
					tr = ( *S ).size++; // Add a new tribe
					add_tribe_count++; // Add a new tribe
					if( debug_level > 2 ) tprintf( "\n" );
					if( compare_particles( &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, &( *S ).best.f, compare_type ) == 1 )
					{
						for( i = 0; i < ( *pb ).nPhi; i++ )
							( *S ).fBestPrev.f[i] = ( *S ).best.f.f[i];
						copy_position( &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest, &( *S ).best ); // BEST COPY
						( *S ).tr_best = tr;
					}
				}
				else
				{
					nExceedSizeSwarm++;
					if( debug_level > 0 ) tprintf( "WARNING: Cannot add a tribe (increase MAXTRIBE = %i)\n", MAXTRIBE );
				}
				break;
			case 1: // Good swarm
				if( ( *S ).size > 1 )
				{
					tribe_worst = 0;
					if( debug_level > 2 ) tprintf( "Good swarm: worst tribe %i is removed\n", tribe_worst + 1 );
					for( tr = 1; tr < ( *S ).size; tr++ ) // Find the worst tribe
					{
						if( compare_particles( &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, &( *S ).trib[tribe_worst].part[( *S ).trib[tribe_worst].best].xBest.f, compare_type ) == -1 )
							tribe_worst = tr;
					}
					if( tribe_worst < ( *S ).size - 1 ) // Remove the worst tribe
						for( tr = tribe_worst; tr < ( *S ).size - 1; tr++ )
							copy_tribe( &( *S ).trib[tr + 1], &( *S ).trib[tr] );
					( *S ).size--;
				}
				break;
		}
	}
	for( tr = 0; tr < nTribe; tr++ ) // Potentially modify particle strategy (strategies are not immediately modified for new particles
	{
		for( pa = 0; pa < ( *S ).trib[tr].size; pa++ )
		{
			if( ( *S ).trib[tr].part[pa].strategy == -1 ) continue;
			status = 3 * compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).trib[tr].part[pa].xLast.f, compare_type );
			status += compare_particles( &( *S ).trib[tr].part[pa].xBest.f, &( *S ).trib[tr].part[pa].fBestPrev, compare_type );
			if( status <= -2 )( *S ).trib[tr].part[pa].status = -1;
			else
			{
				if( status >= 3 )( *S ).trib[tr].part[pa].status = 1;
				else( *S ).trib[tr].part[pa].status = 0;
			}
			switch( ( *S ).trib[tr].part[pa].status )
			{
				case -1: // Bad particle
					do { s = random_int( 0, 2 ); }
					while( s == ( *S ).trib[tr].part[pa].strategy );
					( *S ).trib[tr].part[pa].strategy = s; // Try another strategy
					break;
				case 0: // Normal particle
					( *S ).trib[tr].part[pa].strategy = random_int( 1, 2 );
					break;
				case 1: // Good particle
					if( random_double( 0, 1 ) < 0.5 )( *S ).trib[tr].part[pa].strategy = -1; // Modify the strategy with a probability 0.5
					break;
			}
		}
	}
	if( debug_level )
	{
		if( add_part_count > 0 || add_tribe_count > 0 ) swarm_print( S );
		else if( debug_level > 1 ) swarm_print( S );
	}
}

void swarm_print( struct swarm *S )
{
	int nTotPart = 0, it;
	for( it = 0; it < ( *S ).size; it++ )
		nTotPart += ( *S ).trib[it].size; // Total number of particles
	if( debug_level > 2 ) tprintf( "\n" );
	if( nTotPart > 1 ) tprintf( "Swarm: %i tribes", ( *S ).size );
	else tprintf( "Swarm: %i tribe", ( *S ).size );
	if( debug_level >= 1 )
	{
		tprintf( " [" );
		for( it = 0; it < ( *S ).size; it++ )
			tprintf( " %i", ( *S ).trib[it].size );
		tprintf( " ]" );
	}
	if( ( *S ).size > 1 ) tprintf( " | %i particles | ", nTotPart );
	else tprintf( " | %i particle | ", ( *S ).size );
	if( gop->cd->pdebug ) tprintf( "OF %g E %d (%d) S %d\n", ( *S ).best.f.f[0], eval, gop->cd->neval, gop->global_success );
	if( debug_level > 3 )
		for( it = 0; it < ( *S ).size; it++ )
		{
			tprintf( "Tribe %i: ", it + 1 );
			tribe_print( &( *S ).trib[it] );
		}
}

void swarm_init( struct problem *pb, int compare_type, struct swarm *S )
{
	int nPart, nTribe, i;
	if( gop->cd->init_particles > 0 ) nTribe = gop->cd->init_particles; // Imported number of tribe
	else { nTribe = 10 + ( int ) 2 * sqrt( ( double )( *pb ).D ); if( nTribe < ( *pb ).D ) nTribe = ( *pb ).D; } // Initial number of tribes
	if( nTribe > MAXTRIBE ) nTribe = MAXTRIBE;
	( *S ).size = nTribe; // Number of tribes
	nPart = 1; // Initial number of particles in each tribe
	( *pb ).init = 1;
	for( i = 0; i < nTribe; i++ )
		tribe_init( pb, nPart, compare_type, S, &( *S ).trib[i] );
	copy_position( &( *S ).trib[0].part[( *S ).trib[0].best].xBest, &( *S ).best ); // The best position of the swarm
	( *S ).tr_best = 0;
	if( ( *S ).size > 1 )
		for( i = 1; i < ( *S ).size; i++ )
		{
			if( compare_particles( &( *S ).trib[i].part[( *S ).trib[i].best].xBest.f, &( *S ).best.f, compare_type ) == 1 )
			{
				copy_position( &( *S ).trib[i].part[( *S ).trib[i].best].xBest, &( *S ).best ); // BEST COPY
				( *S ).tr_best = i;
			}
		}
	for( i = 0; i < ( *pb ).nPhi; i++ )
		( *S ).fBestPrev.f[i] = ( *S ).fBestStag.f[i] = ( *S ).best.f.f[i];
	( *S ).status = 0;
}

void swarm_lm( struct problem *pb, struct swarm( *S ) )
{
	struct particle best = {0};
	double phi_current_best, r1;
	int pa, nTotPart = 0;
	int shaman, count_bad_tribes = 0;
	int tr;
	if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "LM optimization cannot be performed; the maximum number of evaluations is achieved!\n" ); return; }
	if( lm_wait == 1 )
	{
		phi_current_best = ( *S ).best.f.f[0];
		if( phi_lm_wait > phi_current_best )
		{
			if( debug_level ) tprintf( "Restart LM search because OF %g is less than %g!\n", phi_current_best, phi_lm_wait );
			lm_wait = 0;
		}
		else
		{
			if( debug_level ) tprintf( "Skip LM search until OF %g is less than %g!\n", phi_current_best, phi_lm_wait );
			return;
		}
	}
	for( tr = 0; tr < ( *S ).size; tr++ )
		nTotPart += ( *S ).trib[tr].size; // Total number of particles
	set_particle( pb, &best );
	if( nTotPart > 2 * ( *pb ).D )
	{
		count_not_good_tribes = 0; // Counter for not very successful LM iterations
		for( tr = 0; tr < ( *S ).size; tr++ )
		{
			shaman = ( *S ).trib[tr].best;
			if( debug_level )
			{
				tprintf( "Shaman move using LM (particle #%d in tribe #%d consisting of %d particles) ... ", shaman + 1, tr + 1, ( *S ).trib[tr].size );
				fflush( stdout );
				if( debug_level > 1 ) { tprintf( "\nold " ); position_print( &( *S ).trib[tr].part[shaman].xBest ); }
			}
			lmo_count++;
			phi_current_best = ( *S ).trib[tr].part[shaman].xBest.f.f[0];
			position_lm( gop, pb, &( *S ).trib[tr].part[shaman].xBest );
			if( debug_level ) tprintf( " OF %g -> %g E %d", phi_current_best, ( *S ).trib[tr].part[shaman].xBest.f.f[0], eval );
			if( debug_level > 1 ) { tprintf( "new " ); position_print( &( *S ).trib[tr].part[shaman].xBest ); }
			copy_position( &( *S ).trib[tr].part[shaman].xBest, &( *S ).trib[tr].part[shaman].x );
			if( gop->cd->check_success && gop->success )
			{
				if( debug_level ) tprintf( "LM Success: Predictions are within the predefined calibration bounds!\n" );
				copy_position( &( *S ).trib[tr].part[shaman].xBest, &( *pb ).pos_success );
				gop->phi = ( *pb ).pos_success.f.f[0];
				break;
			}
//			if( !multiObj && compare_particles( &( *S ).best.f, &( *pb ).maxError, 3 ) == 1 )
			if( ( *pb ).maxError.f[0] > DBL_EPSILON && ( *S ).best.f.f[0] < ( *pb ).maxError.f[0] )
			{
				if( debug_level ) tprintf( "LM Success: OF is minimized below the cutoff value! (%.15g<%.15g)\n", ( *S ).best.f.f[0], ( *pb ).maxError.f[0] );
				return;
			}
			if( compare_particles( &( *S ).trib[tr].part[shaman].xBest.f, &( *S ).best.f, 0 ) == 1 )
			{
				count_not_good_tribes = -( *S ).size;
				copy_position( &( *S ).trib[tr].part[shaman].xBest, &( *S ).best ); // BEST COPY
				gop->phi = ( *S ).best.f.f[0];
				( *S ).tr_best = tr;
			}
			if( !multiObj && compare_particles( &( *S ).best.f, &( *pb ).maxError, 3 ) == 1 ) { if( debug_level ) tprintf( "LM Success: OF is minimized below the cutoff value! (%g<%g)\n", ( *S ).best.f.f[0], ( *pb ).maxError.f[0] ); break; }
			r1 = ( double ) 1 + ( *pb ).lm_factor / 20;
			if( ( *S ).trib[tr].part[shaman].xBest.f.f[0] > phi_current_best / r1 )
			{
				count_bad_tribes++;
				if( debug_level ) tprintf( " Bad shaman (%g > %g / %g)!\n", ( *S ).trib[tr].part[shaman].xBest.f.f[0], phi_current_best, r1 );
				// ( *S ).trib[tr].status = -1; // Bad tribe IF LABELED BAD TRIBE PERFORMANCE DIMINISHES
				particle_init( pb, 6, &( *S ).best, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest, S, &( *S ).trib[tr].part[shaman] ); // reset shaman based on current shamans
				if( debug_level > 1 ) position_print( &( *S ).trib[tr].part[shaman].x );
				for( pa = 0; pa < ( *S ).trib[tr].size; pa++ ) // Possibly update tribe best
					if( compare_particles( &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, &( *S ).trib[tr].part[pa].xBest.f, 0 ) == 1 )
					{
						copy_position( &( *S ).trib[tr].part[pa].x, &( *S ).trib[tr].part[pa].xBest );
						( *S ).trib[tr].best = pa;
						count_bad_tribes--;
					}
				if( compare_particles( &( *S ).trib[tr].part[( *S ).trib[tr].best].x.f, &( *S ).best.f, 0 ) == 1 )
				{
					copy_position( &( *S ).trib[tr].part[( *S ).trib[tr].best].x, &( *S ).best ); // BEST COPY; Possibly update swarm best
					( *S ).tr_best = tr;
				}
			}
			else
			{
				count_not_good_tribes++;
				if( debug_level ) tprintf( " Good shaman (%g <= %g / %g)!\n", ( *S ).trib[tr].part[shaman].xBest.f.f[0], phi_current_best, r1 );
			}
			if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "LM optimization cannot be performed; the maximum number of evaluations is achieved!\n" ); return; }
		}
		if( count_bad_tribes >= ( *S ).size )
		{
			lm_wait = 1;
			phi_lm_wait = ( *S ).best.f.f[0] / ( *pb ).lm_factor;
			if( debug_level ) tprintf( "Skip LM search till OF %g is less than %g!\n", ( *S ).best.f.f[0], phi_lm_wait );
			( *S ).status = -1; // Bad swarm --- LM not happy with the swarm; try to add a tribe
			nSwarmAdaptIter = ( *S ).size * ( *S ).size; // Force swarm adaption
		}
		gop->phi = ( *S ).best.f.f[0];
	}
	// else if(( *S ).size > ( *pb ).D )   // TODO EXPLORE DIFFERENT
	// else if( ( nTotPart > ( double ) 0.9 * ( *pb ).lm_factor * ( *pb ).D || ( double ) ( *pb ).lm_factor * ( *pb ).maxEval < ( double ) 2.0 * eval ) )  // EXPLORE DIFFERENT
	/*
		else if( nTotPart > ( double ) 100 * ( *pb ).D )
		{
			if( debug_level )
			{
				tprintf( "Best particle move using LM ... (particles %d > %f * dimension %d) ", nTotPart, ( *pb ).lm_factor / 2, ( *pb ).D );
				fflush( stdout );
				if( debug_level > 1 ) { tprintf( "\nold " ); position_print( &( *S ).best ); }
			}
			if( eval >= ( *pb ).maxEval ) { if( debug_level ) tprintf( "LM optimizaiton cannot be performed; the maximum number of evaluations is achieved!\n" ); return; }
			lmo_count++;
			phi_current_best = ( *S ).best.f.f[0];
			position_lm( gop, pb, &( *S ).best );
			if( debug_level ) tprintf( "OF %g -> %g\n", phi_current_best, ( *S ).best.f.f[0] );
			if( debug_level > 1 ) { tprintf( "new " ); position_print( &( *S ).best ); }
			tr = ( *S ).tr_best;
			shaman = ( *S ).trib[tr].best;
			if( debug_level > 1 ) tprintf( "BEST Shaman moved (particle #%d in tribe #%d consisting of %d particles OF %g -> %g)!\n", shaman + 1, tr + 1, ( *S ).trib[tr].size, ( *S ).trib[tr].part[shaman].xBest.f.f[0], ( *S ).best.f.f[0] );
			gop->phi = ( *S ).trib[tr].part[shaman].xBest.f.f[0] = ( *S ).best.f.f[0];
			copy_position( &( *S ).trib[tr].part[shaman].xBest, &( *S ).trib[tr].part[shaman].xLast );
			copy_position( &( *S ).best, &( *S ).trib[tr].part[shaman].xBest );
			if( gop->cd->check_success && gop->success ) { copy_position( &( *S ).best, &( *pb ).pos_success ); if( debug_level ) tprintf( "LM Success: Predictions are within the predefined calibration bounds!\n" ); return; }
			if( !multiObj && compare_particles( &( *S ).best.f, &( *pb ).maxError, 3 ) == 1 ) { if( debug_level ) tprintf( "LM Success: OF is minimized below the cutoff value! (%g<%g)\n", ( *S ).best.f.f[0], ( *pb ).maxError.f[0] ); return; }
			else
			{
				nSwarmAdaptIter = ( *S ).size * ( *S ).size;
				( *S ).status = -1; // LM not happy with the swarm; try to add a tribe
			}
		}
	*/
}

void swarm_local_search( struct problem *pb, struct swarm( *S ) ) // Does not add particles; adjusts the best ones only
{
	struct position simplex[MAXPHI + 1];
	struct position new_position;
	double dist, z;
	int m, mm, d;
	int out;
	int pa;
	int shaman;
	int tr;
	int option = 1;
//	option = aleaInteger(0,1); // TODO TO TRY
	set_position( pb, &new_position );
	new_position.size = ( *pb ).D; // new_position is not a new particle
	switch( option )
	{
		case 0: // currently case=1 is executed only; TODO try case=0
			for( m = 0; m < MAXPHI + 1; m++ ) set_position( pb, &simplex[m] );
			for( tr = 0; tr < ( *S ).size; tr++ ) // For each tribe
			{
				if( ( *S ).trib[tr].size < ( *pb ).D + 1 ) continue; // If there is not enough particles do anything
				if( debug_level > 2 ) tprintf( "Local search tribe %i ", tr + 1 );
				mm = ( *S ).trib[tr].best; // Define a simplex
				for( m = 0; m < ( *pb ).D + 1; m++ )
				{
					mm += m;
					if( mm > ( *S ).trib[tr].size - 1 ) mm = 0;
					simplex[m] = ( *S ).trib[tr].part[mm].x;
				}
				//out=aleaInteger(0,1); // TODO TO TRY
				// out=1;
				out = 0; // inside simplex
				switch( out ) // Define a new point
				{
					case 0: // Inside the simplex
						if( debug_level > 2 ) tprintf( "inside simplex " );
						for( d = 0; d < new_position.size; d++ )
						{
							new_position.x[d] = 0;
							for( m = 0; m < ( *pb ).D + 1; m++ ) new_position.x[d] += simplex[m].x[d] / simplex[m].f.f[0];
						}
						for( d = 0; d < new_position.size; d++ ) new_position.x[d] /= ( ( *pb ).D + 1 );
						break;
					case 1: // Outside the simplex
						if( debug_level > 2 ) tprintf( "outside simplex " );
						for( d = 0; d < new_position.size; d++ )
						{
							new_position.x[d] = 0;
							for( m = 1; m < ( *pb ).D + 1; m++ ) new_position.x[d] += simplex[m].x[d];
						}
						for( d = 0; d < new_position.size; d++ ) new_position.x[d] /= ( *pb ).D; // Gravity center
						for( d = 0; d < new_position.size; d++ ) new_position.x[d] -= ( simplex[0].x[d] - new_position.x[d] ); // Reflection of the first vertex
						break;
				}
				position_check( pb, &new_position );
				position_eval( pb, &new_position ); // if there is no improvement the new position will be forgotten
				if( compare_particles( &new_position.f, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, 0 ) == 1 ) // Possibly update tribe best
				{
					if( debug_level > 2 ) tprintf( " SUCCESS" );
					copy_position( &new_position, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest );
					copy_position( &new_position, &( *S ).trib[tr].part[( *S ).trib[tr].best].x );
					if( compare_particles( &new_position.f, &( *S ).best.f, 0 ) == 1 ) // Possibly update swarm best
					{
						copy_position( &new_position, &( *S ).best ); // BEST COPY
						( *S ).tr_best = tr;
					}
				}
				if( debug_level > 2 ) tprintf( "\n" );
			}
			for( m = 0; m < MAXPHI + 1; m++ ) free_position( &simplex[m] );
			break;
		case 1: // Random position between the shaman and the nearest point in the tribe
			for( tr = 0; tr < ( *S ).size; tr++ )
			{
				if( ( *S ).trib[tr].size < ( *pb ).D ) continue; // Small size
				if( debug_level > 2 ) tprintf( "Local search tribe %i near the shaman ", tr + 1 );
				shaman = ( *S ).trib[tr].best;
				for( d = 0; d < ( *pb ).D; d++ ) // find the nearest position to the shaman
				{
					dist = HUGE_VAL;
					for( pa = 0; pa < ( *S ).trib[tr].size; pa++ )
					{
						if( pa == shaman ) continue;
						z = fabs( ( *S ).trib[tr].part[pa].xBest.x[d] - ( *S ).trib[tr].part[shaman].xBest.x[d] );
						if( z < dist ) dist = z;
					}
					new_position.x[d] = ( *S ).trib[tr].part[shaman].xBest.x[d] + ( ( double ) 1.0 - 2.0 * random_double( 0, 1 ) ) * random_double( 0, dist ); // define a random intermediate coordinate
					//new_position.x[d]=z+RandomGauss(z,w2*dist);
				}
				position_check( pb, &new_position );
				position_eval( pb, &new_position );
				if( compare_particles( &new_position.f, &( *S ).trib[tr].part[shaman].xBest.f, 0 ) == 1 ) // Update bests
				{
					if( debug_level > 2 ) tprintf( " SUCCESS" );
					copy_position( &new_position, &( *S ).trib[tr].part[shaman].xBest );
					copy_position( &new_position, &( *S ).trib[tr].part[shaman].x );
					if( compare_particles( &new_position.f, &( *S ).best.f, 0 ) == 1 )
					{
						copy_position( &new_position, &( *S ).best ); // BEST COPY
						( *S ).tr_best = tr;
					}
				}
				if( debug_level > 2 ) tprintf( "\n" );
			}
			break;
	}
	free_position( &new_position );
}

void swarm_move( struct problem *pb, struct swarm( *S ), int compare_type, int run )
{
	struct particle informer = {0};
	int pa;
	int sh = 0;
	int nShaman;
	int shBest;
	int shList[MAXTRIBE];
	int tr, i;
	set_particle( pb, &informer );
	informer.xBest.size = 0;
	modify_weights( ( *pb ).nPhi, run ); // Penalties (global variable phi_weights)
	for( i = 0; i < ( *pb ).nPhi; i++ )
		( *S ).fBestPrev.f[i] = ( *S ).best.f.f[i]; // Save the current best result of the whole swarm
	for( tr = 0; tr < ( *S ).size; tr++ )
	{
		if( debug_level > 2 )
		{
			tprintf( "Tribe %i move: (before) best particle is #%i ", tr + 1, ( *S ).trib[tr].best + 1 );
			if( debug_level == 3 ) objfunc_print( &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f );
			else position_print( &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest );
			tprintf( "\n" );
		}
		for( i = 0; i < ( *pb ).nPhi; i++ )
			( *S ).trib[tr].fBestPrev.f[i] = ( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f.f[i]; // Save the previous best result for the tribe
		for( pa = 0; pa < ( *S ).trib[tr].size; pa++ )
		{
			for( i = 0; i < ( *pb ).nPhi; i++ )
				( *S ).trib[tr].part[pa].fBestPrev.f[i] = ( *S ).trib[tr].part[pa].xBest.f.f[i]; // Save the previous best result for the partiles
			copy_position( &( *S ).trib[tr].part[pa].x, &( *S ).trib[tr].part[pa].xLast ); // Save the last particle position
			if( pa != ( *S ).trib[tr].best ) // not shaman
				copy_particle( &( *S ).trib[tr].part[( *S ).trib[tr].best], &informer ); // INFORMER is the shaman
			else // shaman
			{
				if( multiObj ) // Select an informer in the archive
					informer.xBest = archiveCrowDistSelect( ( *S ).size ); // TODO NEEDS FIX!!!!!!!!!!!!!!!!
				else // Select an informer (attractor) in the swarm
				{
					if( ( *S ).trib[tr].size == 1 ) // Just one particle => the informer is the shaman itself
					{
						copy_particle( &( *S ).trib[tr].part[pa], &informer );
						( *S ).trib[tr].part[pa].strategy = -1; // the particle should inform itself ... random position
					}
					else // There are several tribes. Look for an external informer
					{
						nShaman = random_int( 1, ( *S ).size ); // Random count of informers (shamans) list
						for( sh = 0; sh < nShaman; sh++ )
							shList[sh] = random_int( 0, ( *S ).size - 1 ); // Build a random list of shamans
						shBest = shList[0];
						if( nShaman > 1 )
							for( sh = 1; sh < nShaman; sh++ ) // Find the best
								if( compare_particles( &( *S ).trib[shList[sh]].part[( *S ).trib[shList[sh]].best].xBest.f, &( *S ).trib[shBest].part[( *S ).trib[shBest].best].xBest.f, compare_type ) == 1 )
									shBest = shList[sh];
						if( compare_particles( &( *S ).trib[shBest].part[( *S ).trib[shBest].best].xBest.f, &( *S ).trib[tr].part[pa].xBest.f, compare_type ) )
							copy_particle( &( *S ).trib[shBest].part[( *S ).trib[shBest].best], &informer );
						else
							copy_particle( &( *S ).trib[tr].part[pa], &informer );
					}
				}
			}
			if( debug_level > 2 )
				tprintf( "Tribe %i move: particle #%i is informed by #%i ", tr + 1, ( *S ).trib[tr].part[pa].id, informer.id );
			if( ( *S ).trib[tr].part[pa].id == informer.id )( *S ).trib[tr].part[pa].strategy = -1; // the particle should inform itself ... random position
			position_update( pb, &( *S ).trib[tr].part[pa], informer );
			if( compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).trib[tr].part[pa].xBest.f, compare_type ) == 1 )
			{
				if( debug_level > 2 ) { tprintf( " improvement " ); objfunc_print( &( *S ).trib[tr].part[pa].xBest.f ); tprintf( " => " ); objfunc_print( &( *S ).trib[tr].part[pa].x.f ); tprintf( "\n" ); }
				copy_position( &( *S ).trib[tr].part[pa].x, &( *S ).trib[tr].part[pa].xBest );
			}
			else if( debug_level > 2 ) tprintf( "\n" );
			if( compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).trib[tr].part[( *S ).trib[tr].best].xBest.f, compare_type ) == 1 ) // Possibly update tribe best
			{
				( *S ).trib[tr].best = pa;
				if( debug_level > 2 )
				{
					tprintf( "Tribe %i move: (after) the best particle is #%i ", tr + 1, ( *S ).trib[tr].part[pa].id );
					objfunc_print( &( *S ).trib[tr].part[pa].x.f ); tprintf( "\n" );
				}
			}
			if( compare_particles( &( *S ).trib[tr].part[pa].x.f, &( *S ).best.f, compare_type ) == 1 ) // Possibly update swarm best
			{
				copy_position( &( *S ).trib[tr].part[pa].x, &( *S ).best ); // BEST COPY
				( *S ).tr_best = tr;
			}
		}
	}
	free_particle( &informer );
}

int swarm_particle_count( struct swarm( *S ) )
{
	int nPart = 0;
	int tr;
	for( tr = 0; tr < ( *S ).size; tr++ )
		nPart += ( *S ).trib[tr].size;
	return nPart;
}

int tribe_shaman( struct tribe *T, int compare_type )
{
	int ip, ibest = 0;
	if( ( *T ).size > 1 )
		for( ip = 1; ip < ( *T ).size; ip++ )
			if( compare_particles( &( *T ).part[ip].xBest.f, &( *T ).part[ibest].xBest.f, compare_type ) == 1 )
				ibest = ip;
	return ibest;
}

void tribe_print( struct tribe *T )
{
	int pa;
	tprintf( " %i particles\n", ( *T ).size );
	tprintf( "Labels    :" );
	for( pa = 0; pa < ( *T ).size; pa++ )
		tprintf( " %3i", ( *T ).part[pa].id );
	tprintf( "Labels    :" );
	for( pa = 0; pa < ( *T ).size; pa++ )
		tprintf( " %g", ( *T ).part[pa].xBest.f.f[0] );
	tprintf( "\nStrategies:" );
	for( pa = 0; pa < ( *T ).size; pa++ )
		tprintf( " %3i", ( *T ).part[pa].strategy );
	tprintf( "\n" );
}

void tribe_init( struct problem *pb, int nPart, int compare_type, struct swarm( *S ), struct tribe *T ) // S is not needed here; it is needed by particle_init
{
	int i, init_option;
	( *T ).status = 0;
	( *T ).size = minXY( nPart, MAXPART );
	for( i = 0; i < ( *T ).size; i++ )
	{
		if( ( *pb ).init == 1 ) { ( *pb ).init = 0; init_option = 5; } // User provided initial values
		else init_option = 3; // Biggest empty hyper-parallelepiped
		particle_init( pb, init_option, &( *T ).part[i].x, &( *T ).part[i].x, S, &( *T ).part[i] ); // Arguments 3 & 4 are dummy
	}
	( *T ).best = tribe_shaman( T, compare_type );
	for( i = 0; i < ( *pb ).nPhi; i++ )
		( *T ).fBestPrev.f[i] = ( *T ).part[( *T ).best].xBest.f.f[i] = ( *T ).part[( *T ).best].x.f.f[i];
}

int tribe_varmin_dimension( struct tribe *T )
{
	int i, ip, dim_varmin = -1;
	double mean_dim, var_dim, min_var, z;
	min_var = HUGE_VAL;
	for( i = 0; i < ( *T ).part[0].x.size; i++ )
	{
		mean_dim = 0;
		for( ip = 0; ip < ( *T ).size; ip++ )
			mean_dim += ( *T ).part[ip].xBest.x[i]; // computed on the xBest positions
		mean_dim /= ( *T ).size;
		var_dim = 0;
		for( ip = 0; ip < ( *T ).size; ip++ )
		{
			z = ( *T ).part[ip].xBest.x[i] - mean_dim;
			var_dim += z * z;
		}
		var_dim /= ( *T ).size;
		if( var_dim < min_var ) { dim_varmin = i; min_var = var_dim; }
	}
	return dim_varmin; // Return the dimension on which the variance of the tribe is minimum
}

void modify_weights( int nPhi, int run ) // Modify the penalties defined in the global variable phi_weights[]
{
	// NOTE: this is the most empirical part of the algorithm; It may certainly be largely improved
	int deb, m, n;
	m = random_int( 0, 1 );
	switch( m )
	{
		case 0:
			for( n = 0; n < nPhi; n++ ) phi_weights[n] = 0;
			m = random_int( 0, nPhi - 1 ); phi_weights[m] = 1;
			break;
		case 1: // (Extended) Dynamic Weighted Aggregation
			deb = random_int( 0, nPhi - 1 );
			for( n = 0; n < nPhi; n++ )
			{
				m = n + deb; if( m >= nPhi ) m = nPhi - m;
				phi_weights[m] = fabs( sin( 2 * M_PI * ( n + 1 ) / nPhi ) );
			}
			break;
	}
}

//================================================== KISS
/* A good pseudo-random numbers generator
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
2^32*(2^32-1)*(2^63+2^32-1) > 2^127*/

void seed_rand_kiss( unsigned int seed )
{
	kiss_x = seed | 1;
	kiss_y = seed | 2;
	kiss_z = seed | 4;
	kiss_w = seed | 8;
	kiss_carry = 0;
}

unsigned int rand_kiss()
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

double irand()
{
	int k;
	if( *irand_seed <= 0 ) { tprintf( "ERROR: the seed for random generator is improperly set!\n" ); exit( 1 ); }
	k = *irand_seed / 127773;
	*irand_seed = 16807 * ( *irand_seed - k * 127773 ) - k * 2836;
	if( *irand_seed < 0 ) *irand_seed += 2147483647;
	return( ( double )( *irand_seed ) * 4.656612875E-10 );
}
