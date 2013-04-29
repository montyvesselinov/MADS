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


/* Standard PSO version 2006

Motivation
Quite often some authors say they compare their PSO versions
to the "standard one" ... which is never the same!
So the idea is to define a real standard at least for one year, validated by some
researchers of the field, in particular James Kennedy and Maurice Clerc.
This PSO version does not intend to be the best one on the market (in particular there is
no adaptation of the swarm size nor of the coefficients) but simply very near of the
original version (1995) with just a few improvements based on some recent works.

So referring to "standard PSO 2006" would mean exactly this version with the default values
detailed below as,for example, referring to "standard PSO 2006 (w=0.8)" would mean almost
this version but with a non standard first cognitive/confidence coefficient.

Parameters
S := swarm size
K := maximum number of particles _informed_ by a given one
T := topology of the information links
w := first cognitive/confidence coefficient
c := second cognitive/confidence coefficient
R := random distribution of c
B := rule "to keep the particle in the box"

Equations
For each particle and each dimension
v(t+1) = w*v(t) + R(c)*(p(t)-x(t)) + R(c)*(g(t)-x(t))
x(t+1) = x(t) + v(t+1)
where
v(t) := velocity at time t
x(t) := position at time t
p(t) := best previous position of the particle
g(t) := best previous position of the informants of the particle

Default values
S = 10+2*sqrt(D) where D is the dimension of the search space
K = 3
T := randomly modified after each step if there has been no improvement
w = 1/(2*ln(2))
c = 0.5 + ln(2)
R = U(0..c), i.e. uniform distribution on [0, c]
B := set the position to the min. (max.) value and the velocity to zero

About information links topology
A lot of works have been done about this topic. The main result is there is no
"best" topology. Hence the random approach used here. Note that is does not mean that
each particle is informed by K ones: the number of particles that informs a given one
may be any value between 1 (for each particle informs itself) and S.

About initialisation
Initial positions are chosen at random inside the search space (which is
supposed to be a hyperparallelepid, and even often a hypercube), according to an uniform
distribution. This is not the best way, but the one of the original PSO.

Each initial velocity is simply defined as the difference of two random
positions. It is simple, and needs no additional parameter.
However, again, it is not the best approach. The resulting distribution is not even
uniform, as for any method that uses an uniform distribution independently for each
component. The mathematically correct approach needs to use an uniform
distribution inside a hypersphere. It is not that difficult, and indeed used in some PSO
versions, but quite different from the original one.

Some results with the standard values
You may want to recode this program in another language. Also you may want to modify it
for your own purposes. Here are some results on classical test functions to help you to
check your code.
Dimension D=30
Acceptable error eps=0
Objective value f_min=0
Number of runs n_exec_max=50
Number of evaluations for each run eval_max=30000

Problem                            Mean best value    Standard deviation
Parabola/Sphere on [-100, 100]^D        0                  0
Griewank on [-600, 600]^D               0.018              0.024
Rosenbrock/Banana on [-30, 30]^D       50.16              36.9
Rastrigin on [-5.12, 5.12]^D           48.35              14.43
Ackley on [-32, 32]^D                   1.12               0.85

Last updates
2006-02-27 Fixed a bug about minimal best value over several runs
2006-02-16 Fixed a bug (S_max for P, V, X, instead of D_max), thanks to Manfred Stickel
2006-02-16 replaced k by i x by xs (in perf()), because of possible confusion with K and X
2006-02-13  added De Jong's f4
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
//#include "alea.c"
#include "kdtree-0.5.5/kdtree.h"
#include "../mads.h"

#define	D_max 100  // Max number of dimensions of the search space
#define	S_max 100 // Max swarm size
#define R_max 200 // Max number of runs
//#define E_max 100000// Max number of evaluations


// Structures
struct velocity
{
	int size;
	double v[D_max];
};


struct position
{
	int size;
	double x[D_max];
	double f;
};

// Sub-programs
double perf_pssa( int s, int function ); // Fitness evaluation
double round_to_res( double x, double dx );
int compare_pos( double x[], double y[], int size );
void read_in( struct kdtree *kdtree, int size, double eps, int check_success, double *finv, int *ind, char *file );
void write_loc( double of, double *x, int x_size, int *ind );
void print_collected( struct kdtree *kd, double *eps, double *dmax, double *dxmin, FILE *fid ); 	// Save result
void print_particles( struct position *X, FILE *fid );
// in Standard_PSO_2006.c
double alea( double a, double b );
int alea_integer( int a, int b );
// in ../mads.c
int get_seed( );
void seed_rand_kiss( unsigned long seed );

// Global variables
int best; // Best of the best position (rank in the swarm)
int D; // Search space dimension
double E; // exp(1). Useful for some test functions
double f_min; // Objective(s) to reach
int LINKS[S_max] [S_max]; // Information links
int nb_eval; // Total number of evaluations
double pi; // Useful for some test functions
struct position P[S_max]; // Best positions found by each particle
int S; // Swarm size
struct velocity V[S_max]; // Velocities
struct position X[S_max]; // Positions
//struct position G[E_max]; // Eaten positions
double xmin[D_max], xmax[D_max]; // Intervals defining the search space
double dx[D_max]; // Resolution of search
struct opt_data *gop;
double *res;


// File(s)
FILE *f_run;
//--------Internal random generator
int *irand_seed;

// =================================================
int abagus( struct opt_data *op )
{
	double c; // Second onfidence coefficient
	int d; // Current dimension
	double eps, err; // Admissible error
	double eps_mean; // Average error
	double error; // Error for a given position
	double error_prev; // Error after previous iteration
	int eval_max; // Max number of evaluations
	double eval_mean; // Mean number of evaluations
	int function; // Code of the objective function
	int g; // Rank of the best informant
	int init_links; // Flag to (re)init or not the information links
	int i, j, k;
	int K; // Max number of particles informed by a given one
	int m;
	//double mean_best[R_max];
	double min; // Best result through several runs
	int n_exec, n_exec_max; // Nbs of executions
	int n_failure; // Number of failures
	int s; // Rank of the current particle
	//double t1, t2;
	//double variance;
	double w; // First confidence coefficient
	int f_ind = 0, f_ind_bad = 0, f_ind_old = 0; // finv index
	struct kdtree *kd; // pointer to kdtree
	struct kdtree *kdbad; // pointer to kdtree of points not to keep
	struct kdres *kdset; // nearest neighbor search results
	struct kdres *kdsetbad; // nearest neighbor search results
	double *pch; // retrieved phi from kd tree search
	double f, *finv, *finvbad; // phi and adjusted phi value
	struct position G; // Eaten positions
	double dmax = 0; // maximum dimension length
	double dxmin = 100000; // minimum parameter dx
	char filename[80]; // output filename
	int enrgy_add; // energy to add for good position
	int energy; // energy of swarm
	int energy_prev; // energy from previous iteration
	int kdsize = 0; // size of kdset
	int kdsizebad = 0; // size of bad kdset
	int old_pos = 0, old_bad_pos = 0, new_pos = 0; // number of particles visiting old and new positions each iteration
	int t_old_pos = 0, t_old_bad_pos = 0, t_new_pos = 0; // total number of particles visiting old and new positions
	//double expl_rate, old_expl_rate = 1; // rate of new to old positions
	int mflag; // indicates if particles are actually moving
	double domain_size = 1;
	double cell_size = 1;
	int nb_eval_left;
	int reinvert_flag = 0;
	FILE *f_out;
	int iter = 0;
	op->cd->compute_phi = 1;
	if( ( res = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); exit( 1 ); }
	eval_max = op->cd->maxeval; // Max number of evaluations for each run
	if( ( finv = ( double * ) malloc( eval_max * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); exit( 1 ); }
	if( ( finvbad = ( double * ) malloc( eval_max * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); exit( 1 ); }
	irand_seed = &op->cd->seed;
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); }
	else { seed_rand_kiss( op->cd->seed ); srand( op->cd->seed ); if( op->cd->pdebug ) tprintf( "Current seed: %d\n", op->cd->seed ); }
	E = exp( 1 );
	pi = acos( -1 );
	D = op->pd->nOptParam; // Search space dimension
	gop = op;
	kd = (struct kdtree *) kd_create( D ); // initialize kdtree
	kdbad = (struct kdtree *) kd_create( D ); // initialize kdtree
	energy = op->cd->energy;
	enrgy_add = energy / 10;
	// D-cube data
	for( d = 0; d < D; d++ )
	{
		xmin[d] = op->pd->var_min[op->pd->var_index[d]];
		xmax[d] = op->pd->var_max[op->pd->var_index[d]];
		dx[d] = op->pd->var_dx[op->pd->var_index[d]];
		// Find largest dimension in parameter space
		if( dmax < ( xmax[d] - xmin[d] ) ) dmax = xmax[d] - xmin[d];
		// Find smallest dx among parameters
		if( dxmin > dx[d] ) dxmin = dx[d];
		// Calculate multi-dimensional size of domain
		domain_size *= fabs( xmax[d] - xmin[d] );
		// Calculate multi-dimensional size of cell of discretized domain
		cell_size *= dx[d];
	}
	// Determine eps by max possible phi within success
	if( op->cd->check_success )
	{
		eps = 0;
		for( k = 0; k < op->od->nTObs; k++ )
		{
			i = op->od->obs_well_index[k];
			j = op->od->obs_time_index[k];
			if( op->wd->obs_log[i][j] == 0 )
			{
				if( ( op->wd->obs_target[i][j] - op->wd->obs_min[i][j] ) > ( op->wd->obs_max[i][j] - op->wd->obs_target[i][j] ) )
					err =  op->wd->obs_target[i][j] - op->wd->obs_min[i][j];
				else err = op->wd->obs_max[i][j] - op->wd->obs_target[i][j];
				if( op->cd->objfunc_type != SSR )
				{
					if( op->cd->objfunc_type == SSD0 ) err = 0;
					if( op->cd->objfunc_type == SSDA )
						err = sqrt( fabs( err ) );
				}
			}
			else
			{
				if( ( op->wd->obs_target[i][j] - op->wd->obs_min[i][j] ) > ( op->wd->obs_max[i][j] - op->wd->obs_target[i][j] ) )
					err = log10( op->wd->obs_target[i][j] ) - log10( op->wd->obs_min[i][j] );
				else err = log10( op->wd->obs_max[i][j] ) - log10( op->wd->obs_target[i][j] );
			}
			eps += pow( err * op->wd->obs_weight[i][j], 2 );
		}
		if( op->cd->pdebug ) printf( "Max OF within success: %g\n", eps );
	}
	else eps = op->cd->phi_cutoff; // If success option is not selected, use phi_cutoff
	if( eps < op->cd->phi_cutoff ) { tprintf( "phi_cutoff > Max OF within success (%lf > %lf), abagus will use phi_cutoff\n", op->cd->phi_cutoff, eps ); eps = op->cd->phi_cutoff;  } // If max eps within success is less than phi_cutoff
	f_min = 0; // Objective value
	n_exec_max = 1; // Numbers of runs
	// Read in previous results
	if( op->cd->infile[0] != 0 )
	{
		read_in( kd, D, eps, op->cd->check_success, finv, &f_ind, op->cd->infile );
		f_ind_old = f_ind;
	}
	// Write accepted locations from input file to output file
	// Open output file
	sprintf( filename, "%s.abagus", op->root ); // TODO rename pssa to abagus to avoid confusion
	if( ( f_run = fopen( filename, "w" ) ) == NULL ) { tprintf( "File %s cannot be opened to write results!\n", filename ); exit( 0 ); }
	fprintf( f_run, "Number OF parameters...\n" ); // Write header
	for( d = 0; d < D; d++ ) G.x[d] = xmin[d] + 0.5 * ( xmax[d] - xmin[d] ); // Determine center parameter space for search
	if( op->cd->infile[0] != 0 )
	{
		if( f_ind > 0 )
		{
			kdset = (struct kdres *) kd_nearest_range( kd, G.x, dmax );
			i = 0;
			while( !kd_res_end( kdset ) )
			{
				i++;
				pch = ( double * ) kd_res_item( kdset, G.x );
				//fprintf( f_run, "%d ", i + 1 );
				f = eps - ( *pch - eps );
				write_loc( f, G.x, D, &i );
				//fprintf( f_run, " %lf", f );
				//for( d = 0; d < D; d++ )
				//	fprintf( f_run, " %lf", G.x[d] );
				//fprintf( f_run, "\n" );
				kd_res_next( kdset );
			}
			tprintf( "\n%d locations accepted from %s\n\n", kd_res_size( kdset ), op->cd->infile );
			kd_res_free( kdset );
		}
		else if( op->cd->check_success ) tprintf( "\nNo solutions found in %s within target ranges\n\n", op->cd->infile );
		else tprintf( "\nNo solutions found in %s at phi cutoff = %g\n\n", op->cd->infile, eps );
	}
	if( n_exec_max > R_max ) n_exec_max = R_max;
	//-----------------------------------------------------  PARAMETERS
	S = 10 + ( int )( 2 * sqrt( D ) ); if( S > S_max ) S = S_max;
	K = 3;
	//w = 1 / ( 2 * log( 2 ) ); c = 0.5 + log( 2 );
	//w = 1.2; c = 1.7;
	w = 0.5; c = 0.7;
	if( op->cd->pdebug ) tprintf( "\n Swarm size %i", S );
	if( op->cd->pdebug ) tprintf( "\n coefficients %f %f \n", w, c );
	//----------------------------------------------------- INITIALISATION
	//t1 = clock(); // Init time
	// Initialisation of information variables
	n_exec = 0; eval_mean = 0; eps_mean = 0; n_failure = 0;
//init:
	n_exec = n_exec + 1;
	// Set first particle to IVs from mads input file
	for( d = 0; d < D; d++ )
		X[0].x[d] = op->pd->var[op->pd->var_index[d]];
	if( op->cd->pdebug )
	{
		tprintf( "\nParticle 1 position from initial values:\n" );
		for( d = 0; d < D; d++ )
			tprintf( "%g\n", X[0].x[d] );
		tprintf( "\n" );
	}
	// Randomly initialize remaining particles
	for( s = 1; s < S; s++ )   // Positions and velocities
	{
		X[s].size = D; V[s].size = D;
		for( d = 0; d < D; d++ )
		{
			X[s].x[d] = alea( xmin[d], xmax[d] );
			X[s].x[d] = round_to_res( X[s].x[d], dx[d] );
			V[s].v[d] = ( alea( xmin[d], xmax[d] ) - X[s].x[d] ) / 2; // Non uniform
			// V[s].v[d] = ( xmin[d]-xmax[d] )*(0.5-alea(0,1)); //Uniform. 2006-02-24
		}
		// Transform all parameter to have unit dx
	}
	// First evaluations
	nb_eval = 0;
	for( s = 0; s < S; s++ )
	{
		if( f_ind > 0 || f_ind_bad > 0 )
		{
			kdset = (struct kdres *) kd_nearest_range( kd, X[s].x, dxmin * 0.9 );
			kdsetbad = (struct kdres *) kd_nearest_range( kdbad, X[s].x, dxmin * 0.9 );
			if( kd_res_size( kdset ) == 0 && kd_res_size( kdsetbad ) == 0 )
				X[s].f = fabs( perf_pssa( s, function ) - f_min );
			else if( kd_res_size( kdset ) == 1 && kd_res_size( kdsetbad ) == 0 )
			{
				pch = ( double * ) kd_res_item( kdset, X[s].x );
				X[s].f = *pch;
			}
			else if( kd_res_size( kdsetbad ) == 1 && kd_res_size( kdset ) == 0 )
			{
				pch = ( double * ) kd_res_item( kdsetbad, X[s].x );
				X[s].f = *pch;
			}
			else tprintf( "Warning: Point(s) within grid resolution of proposal point!\n" );
			kd_res_free( kdset );
			kd_res_free( kdsetbad );
		}
		else X[s].f = fabs( perf_pssa( s, function ) - f_min );
		//}
		// invert OF values if below eps
		//for( s = 0; s < S; s++ )
		//{
		if( ( X[s].f < eps && ( ! op->cd->check_success ) ) || ( op->cd->check_success && op->success ) )
		{
			finv[f_ind] = eps + ( eps - X[s].f );
			kd_insert( kd, X[s].x, &finv[f_ind] );
			f_ind++;
			write_loc( X[s].f, X[s].x, D, &f_ind );
			energy += enrgy_add;
		}
		else if( kdsize == 0 && kdsizebad == 0 )
		{
			finvbad[f_ind_bad] = X[s].f;
			kd_insert( kdbad, X[s].x, &finvbad[f_ind_bad] );
			f_ind_bad++;
		}
		P[s] = X[s]; // Best position = current one
	}
	// Find the best
	best = 0;
	for( s = 1; s < S; s++ )
		if( P[s].f < P[best].f ) best = s;
	error =  P[best].f ; // Current min error
	if( n_exec == 1 ) min = error;
	error_prev = error; // Previous min error
	energy_prev = energy; // Previous energy
	init_links = 1; // So that information links will be initialized
	//---------------------------------------------- ITERATIONS
loop:
	iter++;
	new_pos = old_pos = old_bad_pos = 0;
	if( init_links == 1 )
	{
		// Who informs who, at random
		for( s = 0; s < S; s++ )
		{
			for( m = 0; m < S; m++ ) LINKS[m] [s] = 0;  // Init to "no link"
			LINKS[s] [s] = 1; // Each particle informs itself
		}
		for( m = 0; m < S; m++ )  // Other links
		{
			for( i = 0; i < K; i++ )
			{
				s = alea_integer( 0, S - 1 );
				LINKS[m] [s] = 1;
			}
		}
	}
	// The swarm MOVES
	for( s = 0; s < S; s++ )  // For each particle ...
	{
		// .. find the best informant
		g = s;
		for( m = 0; m < S; m++ )
		{
			if( LINKS[m] [s] == 1 && P[m].f < P[g].f ) g = m;
		}
		// ... compute the new velocity, and move
		mflag = 0;
		for( d = 0; d < D; d++ )
		{
			V[s].v[d] = w * V[s].v[d] + alea( 0, c ) * ( P[s].x[d] - X[s].x[d] );
			V[s].v[d] = V[s].v[d] + alea( 0, c ) * ( P[g].x[d] - X[s].x[d] );
			if( V[s].v[d] > 0.5 * dx[d] ) mflag = 1;
		}
		if( mflag == 1 ) // use velocity to move
			for( d = 0; d < D; d++ )
			{
				X[s].x[d] = X[s].x[d] + V[s].v[d];
				X[s].x[d] = round_to_res( X[s].x[d], dx[d] );
			}
		else // reinitialize particle
			for( d = 0; d < D; d++ )
			{
				X[s].x[d] = alea( xmin[d], xmax[d] );
				V[s].v[d] = ( alea( xmin[d], xmax[d] ) - X[s].x[d] ) / 2; // Non uniform
				// V[s].v[d] = ( xmin[d]-xmax[d] )*(0.5-alea(0,1)); //Uniform. 2006-02-24
			}
		// ... interval confinement (keep in the box)
		for( d = 0; d < D; d++ )
		{
			if( X[s].x[d] < xmin[d] ) { X[s].x[d] = xmin[d]; V[s].v[d] = 0; }
			if( X[s].x[d] > xmax[d] ) { X[s].x[d] = xmax[d]; V[s].v[d] = 0; }
			X[s].x[d] = round_to_res( X[s].x[d], dx[d] );
		}
		// ... evaluate the new position
		if( f_ind > 0 || f_ind_bad > 0 )
		{
			reinvert_flag = 1;
			kdset = (struct kdres *) kd_nearest_range( kd, X[s].x, dxmin * 0.9 );
			kdsetbad = (struct kdres *) kd_nearest_range( kdbad, X[s].x, dxmin * 0.9 );
			kdsize = kd_res_size( kdset );
			kdsizebad = kd_res_size( kdsetbad );
			if( kdsize == 0 && kdsizebad == 0 )
			{
				X[s].f = fabs( perf_pssa( s, function ) - f_min );
				new_pos++;
			}
			else if( kdsize == 1 && kdsizebad == 0 )
			{
				pch = ( double * ) kd_res_item( kdset, X[s].x );
				X[s].f = *pch;
				old_pos++;
				op->success = 0;
				// Check if finv was set back to f during discretization interval reduction
				if( X[s].f < eps )
				{
					reinvert_flag = 0;
					*pch = eps + ( eps - X[s].f ); // invert f again
					energy += enrgy_add;
				}
			}
			else if( kdsizebad == 1 && kdsize == 0 )
			{
				pch = ( double * ) kd_res_item( kdsetbad, X[s].x );
				X[s].f = *pch;
				old_bad_pos++;
				//op->success = 0;
			}
			else tprintf( "Warning: Multiple points within grid resolution of proposal point!\n" );
			kd_res_free( kdset );
			kd_res_free( kdsetbad );
		}
		else { X[s].f = fabs( perf_pssa( s, function ) - f_min ); new_pos++; }
		// if below eps, insert into kdtree
		//if( X[s].f < eps )
		//TODO: success checking needs to be fixed
		if( ( X[s].f < eps && ( ! op->cd->check_success ) && reinvert_flag ) || ( op->cd->check_success && op->success ) )
		{
			finv[f_ind] = eps + ( eps - X[s].f ); // invert f
			kd_insert( kd, X[s].x, &finv[f_ind] );
			f_ind++;
			write_loc( X[s].f, X[s].x, D, &f_ind );
			energy += enrgy_add;
			//X[s].f = finv[f_ind]; // set phi to inverted value
		}
		else if( kdsize == 0 && kdsizebad == 0 )
		{
			finvbad[f_ind_bad] = X[s].f;
			kd_insert( kdbad, X[s].x, &finvbad[f_ind_bad] );
			f_ind_bad++;
		}
		// ... check for the best previous position if greater than eps
		if( X[s].f < P[s].f )
		{
			P[s] = X[s];
			// ... update the best of the bests
			if( P[s].f < P[best].f ) best = s;
		}
		energy--; // Move takes energy
	}
	// If no improvement, information links will be reinitialized
	error = P[best].f;
	if( error >= error_prev ) init_links = 1;
	else init_links = 0;
	error_prev = error;
	// If particles are not exploring new locations enough, modify w and c
//	if( old_pos > 0 ) expl_rate = new_pos / ( double ) old_pos;
//	else { expl_rate = old_expl_rate = 10000; }
	/*	if( w > 0.001 && c < 3 && c > 0.001 )
		{
			if( expl_rate < 0.1 ) {w *= 1.001; c *= 1.001;}
			else if( expl_rate > 30 ) {w *= 0.999; c *= 0.999; }
			else if( expl_rate < old_expl_rate ) {w *= 1.00001; c *= 1.00001;}
			else if( expl_rate > old_expl_rate ) {w *= 0.99999; c *= 0.99999;}
		}*/
	// Reset best location if expl_rate is too low
	//if( 2 * old_pos > new_pos )
	if( energy < energy_prev )
	{
		for( s = 0; s < S; s++ ) P[s] = X[s];
		best = 0;
		for( s = 1; s < S; s++ )
			if( P[s].f < P[best].f ) best = s;
		init_links = 0;
	}
	energy_prev = energy;
	if( op->cd->pdebug > 1 ) tprintf( "evals: %d; n_found: %d; old_pos: %d; new_pos: %d; old_bad_pos: %d; energy: %d\n", nb_eval, f_ind, old_pos, new_pos, old_bad_pos, energy );
	t_new_pos += new_pos; t_old_pos += old_pos; t_old_bad_pos += old_bad_pos;
	// Check if dx can be reduced
	nb_eval_left = eval_max - nb_eval;
	if( ( nb_eval * pow( 2, D ) * 5 ) < nb_eval_left && energy <= 0 && f_ind > 0 ) // only do if an acceptable solution has been found
	{
		cell_size = 1;
		for( d = 0; d < D; d++ )
		{
			dx[d] /= 2;
			cell_size *= dx[d];
			if( dxmin > dx[d] ) dxmin = dx[d];
		}
		tprintf( "Decreasing discretization intervals by factor of 2.\n" );
		tprintf( "dx = [" );
		for( d = 0; d < D; d++ ) tprintf( " %g", dx[d] );
		tprintf( " ]\n" );
		energy = 1;
		// Reset finv to f for recollection
		for( i = 0; i < f_ind; i++ )
			finv[i] = eps + ( eps - finv[i] ); // return inverted (collected) OFs to f for recollection
		// Set first particle to IVs from mads input file, presumably obtained by optimization
		for( d = 0; d < D; d++ )
			X[0].x[d] = op->pd->var[op->pd->var_index[d]];
		energy = op->cd->energy; // Put energy back to original
	}
	// Print out collected positions
	if( op->cd->pdebug == 10 )
	{
		sprintf( filename, "%s-%08d.collected", op->root, iter );
		if( ( f_out = fopen( filename, "w" ) ) == NULL ) { tprintf( "File %s cannot be opened to write results!\n", filename ); exit( 0 ); }
		print_collected( kd, &eps, &dmax, &dxmin, f_out );
		fclose( f_out );
		sprintf( filename, "%s-%08d.particles", op->root, iter );
		if( ( f_out = fopen( filename, "w" ) ) == NULL ) { tprintf( "File %s cannot be opened to write results!\n", filename ); exit( 0 ); }
		print_particles( X, f_out );
		fclose( f_out );
	}
	// Check if finished
	if( nb_eval < eval_max && energy > 0 ) goto loop;
	//if( nb_eval < eval_max ) goto loop;
	if( error > eps ) n_failure = n_failure + 1;
	// Result display
	if( op->cd->pdebug ) tprintf( "Current seed: %d\n", op->cd->seed );
	if( op->cd->infile[0] != 0 )
	{
		tprintf( "\n%d solutions from %s\n", f_ind_old, op->cd->infile );
		tprintf( "%d new solutions\n", f_ind - f_ind_old );
	}
	kdset = kd_nearest_range( kd, X[s].x, dmax );
	kdsetbad = kd_nearest_range( kdbad, X[s].x, dmax );
	kdsize = kd_res_size( kdset );
	kdsizebad = kd_res_size( kdsetbad );
	kd_res_free( kdset );
	kd_res_free( kdsetbad );
	tprintf( "\n%d total solutions collected\n", f_ind );
	tprintf( "%d total solutions rejected\n", f_ind_bad );
	tprintf( "%d calculated solutions\n", nb_eval );
	tprintf( "%d revisits to saved solutions\n", t_old_pos );
	tprintf( "%d revisits to bad solutions\n", t_old_bad_pos );
	tprintf( "\nFraction of domain with collected solutions: %g\n", ( double ) f_ind * cell_size / domain_size );
	if( energy > 0 ) tprintf( "\nExploration may not be complete! Swarm energy still left (energy = %d)! Increase evals!\n\n", energy );
	else
	{
		tprintf( "cell_size / domain_size = %g\n", cell_size / domain_size );
		tprintf( "Swarm energy used up: energy = %d\n\n", energy );
	}
	/*	// Save result
		for( d = 0; d < D; d++ ) G.x[d] = xmin[d] + 0.5 * ( xmax[d] - xmin[d] );
		if( f_ind > 0 )
		{
			kdset = kd_nearest_range( kd, G.x, dmax );
			fprintf( f_run, "Number OF parameters...\n" );
			i = 0;
			while( !kd_res_end( kdset ) )
			{
				pch = ( double * ) kd_res_item( kdset, G.x );
				fprintf( f_run, "%d ", i + 1 );
				f = eps - ( *pch - eps );
				fprintf( f_run, " %lf", f );
				for( d = 0; d < D; d++ )
					fprintf( f_run, " %lf", G.x[d] );
				fprintf( f_run, "\n" );
				i++;
				kd_res_next( kdset );
			}
			kd_res_free( kdset );
			tprintf( "\nResults written to %s\n\n", filename );
		}
		else if( op->cd->check_success ) tprintf( "\nNo solutions found within observation ranges (success)\n\n");
		else tprintf( "\nNo solutions found at phi cutoff = %g\n\n", eps );
	*/
	kd_free( kd );
	fclose( f_run );
	/*    // Compute some statistical information
	    if ( error < min ) min = error;
	    eval_mean = eval_mean + nb_eval;
	    eps_mean = eps_mean + error;
	    mean_best[n_exec - 1] = error;

	    if ( n_exec < n_exec_max ) goto init;

	    // END. Display some statistical information
	    t2 = clock();
	    tprintf( "\n\n Total clocks %.0f", t2 - t1 );
	    eval_mean = eval_mean / ( double ) n_exec;
	    eps_mean = eps_mean / ( double ) n_exec;
	    tprintf( "\n\n Eval. (mean)= %f", eval_mean );
	    tprintf( "\n Error (mean) = %f", eps_mean );

	    // Variance
	    variance = 0;
	    for ( d = 0; d < n_exec_max; d++ ) variance = variance + ( mean_best[d] - eps_mean ) * ( mean_best[d] - eps_mean );
	    variance = sqrt( variance / n_exec_max );
	    tprintf( "\n Std. dev. %f", variance );

	    // Success rate and minimum value
	    tprintf( "\n Success rate = %.2f%%\n", 100 * ( 1 - n_failure / ( double ) n_exec ) );
	    if ( n_exec > 1 ) tprintf( "\n Best min value = %f", min );
	*/
	op->cd->compute_phi = 0;
	return 0;
}

//===========================================================
/*  double alea( double a, double b )
  { // random number (uniform distribution) in [a b]
    double r;
     r=(double) rand(); r=r/RAND_MAX;
    return a + r * ( b - a );
  }
  //===========================================================
  int alea_integer( int a, int b )
  { // Integer random number in [a b]
    int ir;
    double r;
    r = alea( 0, 1 ); ir = ( int )( a + r * ( b + 1 - a ) );
    if ( ir > b ) ir = b;
    return ir;
  }*/
//==========================================================
double round_to_res( double x, double dx )
{
	double x_rnd, resol;
	int n;
	n = ( int )( x / dx );
	resol = x - n * dx ;
	if( fabs( resol ) <= dx / 2 )
		x_rnd = x - resol;
	else if( x > 0 )
		x_rnd = x + ( dx - resol );
	else // fabs(resol) > dx && x < 0
		x_rnd = x - ( dx - fabs( resol ) );
	return x_rnd;
}
//===========================================================
int compare_pos( double x[], double y[], int size )
{
	int i;
	for( i = 0; i < size; i++ )
		if( x[i] != y[i] ) return 0;
	return 1;
}

void read_in( struct kdtree *kdtree, int size, double eps, int check_success, double *finv, int *ind, char *file )
{
	FILE *fl;
	char buf[500];
	double of, pars[size];
	int i;
	if( check_success ) tprintf( "\nAssuming samples in %s meet current criteria for success!", file );
	tprintf( "\nReading previous results from %s...", file );
	fl = fopen( file, "r" );
	if( fl == NULL ) { tprintf( "\nError opening %s\n", file ); exit( 0 ); }
	fgets( buf, sizeof buf, fl );
	while( fscanf( fl, "%*d %lf", &of ) > 0 )
	{
		if( of < eps || check_success )
		{
			finv[*ind] = eps + ( eps - of );
			for( i = 0; i < size; i++ )
				fscanf( fl, "%lf", &pars[i] );
			fscanf( fl, " \n" );
			kd_insert( kdtree, pars, &finv[*ind] );
			( *ind )++;
		}
	}
	fclose( fl );
	tprintf( "Done\n" );
}

void write_loc( double of, double *x, int x_size, int *ind )
{
	int d;
	fprintf( f_run, "%d ", ( *ind ) );
	//f = eps - ( *pch - eps );
	fprintf( f_run, " %lf", of );
	for( d = 0; d < x_size; d++ )
		fprintf( f_run, " %lf", x[d] );
	fprintf( f_run, "\n" );
	fflush( f_run );
}

void print_collected( struct kdtree *kd, double *eps, double *dmax, double *dxmin, FILE *fid ) 	// Save result
{
	int d;
	struct kdres *kdset; // nearest neighbor search results
	double *pch;
	struct position G;
	for( d = 0; d < D; d++ ) G.x[d] = xmin[d] + 0.5 * ( xmax[d] - xmin[d] );
//	if( f_ind > 0 )
//	{
	kdset = kd_nearest_range( kd, G.x, *dmax );
	fprintf( fid, "OF parameters... dxmin = %g\n", *dxmin );
	while( !kd_res_end( kdset ) )
	{
		pch = ( double * ) kd_res_item( kdset, G.x );
		//fprintf( fid, "%d ", i + 1 );
		//f = eps - ( *pch - eps );
		fprintf( fid, "%lf", *pch );
		for( d = 0; d < D; d++ )
			fprintf( fid, " %lf", G.x[d] );
		fprintf( fid, "\n" );
		//i++;
		kd_res_next( kdset );
	}
	kd_res_free( kdset );
	//printf( "\nResults written to %s\n\n", filename );
//	}
}

void print_particles( struct position *X, FILE *fid )
{
	int s, d;
	fprintf( fid, "OF parameters...\n" );
	for( s = 0; s < S; s++ )
	{
		fprintf( fid, "%lf", X[s].f );
		for( d = 0; d < D; d++ )
			fprintf( fid, " %lf", X[s].x[d] );
		fprintf( fid, "\n" );
	}
}

//===========================================================
double perf_pssa( int s, int function )
{
	// Evaluate the fitness value for the particle of rank s
	int d;
	int i, j, k;
	double f, p, xd, x1, x2;
	double sum1, sum2;
	double t0, tt, t1;
	struct position xs;
	// For Foxholes problem
	static int a[2] [25] =
	{
		{
			-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32
		},
		{
			-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32
		}
	};
	// For polynomial fitting problem
	int const M = 60;
	double py, y = -1, dx = ( double ) M;
	nb_eval = nb_eval + 1;
	xs = X[s];
	switch( function )
	{
		case 0: // Parabola (Sphere)
			f = 0;
			p = 0; // Shift
			for( d = 0; d < D; d++ )
			{
				xd = xs.x[d] - p;
				f = f + xd * xd;
			}
			break;
		case 1: // De Jong's f4
			f = 0;
			p = 0; // Shift
			for( d = 0; d < D; d++ )
			{
				xd = xs.x[d] - p;
				f = f + ( d + 1 ) * xd * xd * xd * xd;
			}
			break;
		case 2: // Griewank
			f = 0;
			p = 1;
			for( d = 0; d < D; d++ )
			{
				xd = xs.x[d];
				f = f + xd * xd;
				p = p * cos( xd / sqrt( d + 1 ) );
			}
			f = f / 4000 - p + 1;
			break;
		case 3: // Rosenbrock
			f = 0;
			t0 = xs.x[0];
			for( d = 1; d < D; d++ )
			{
				t1 = xs.x[d];
				tt = 1 - t0;
				f += tt * tt;
				tt = t1 - t0 * t0;
				f += 100 * tt * tt;
				t0 = t1;
			}
			break;
		case 4: // Step
			f = 0;
			for( d = 0; d < D; d++ ) f = f + ( int ) xs.x[d];
			break;
		case 6: //Foxholes 2D
			f = 0;
			for( j = 0; j < 25; j++ )
			{
				sum1 = 0;
				for( d = 0; d < 2; d++ )
				{
					sum1 = sum1 + pow( xs.x[d] - a[d] [j], 6 );
				}
				f = f + 1 / ( j + 1 + sum1 );
			}
			f = 1 / ( 0.002 + f );
			break;
		case 7: // Polynomial fitting problem
			// on [-100 100]^9
			f = 0;
			dx = 2 / dx;
			for( i = 0; i <= M; i++ )
			{
				py = xs.x[0];
				for( d = 1; d < D; d++ )
				{
					py = y * py + xs.x[d];
				}
				if( py < -1 || py > 1 ) f += ( 1 - py ) * ( 1 - py );
				y += dx;
			}
			py = xs.x[0];
			for( d = 1; d < D; d++ ) py = 1.2 * py + xs.x[d];
			py = py - 72.661;
			if( py < 0 ) f += py * py;
			py = xs.x[0];
			for( d = 1; d < D; d++ ) py = -1.2 * py + xs.x[d];
			py = py - 72.661;
			if( py < 0 ) f += py * py;
			break;
		case 8: // Clerc's f1, Alpine function, min 0
			f = 0;
			for( d = 0; d < D; d++ )
			{
				xd = xs.x[d];
				f += fabs( xd * sin( xd ) + 0.1 * xd );
			}
			break;
		case 9: // Rastrigin. Minimum value 0. Solution (0,0 ...0)
			k = 10;
			f = 0;
			for( d = 0; d < D; d++ )
			{
				xd = xs.x[d];
				f += xd * xd - k * cos( 2 * pi * xd );
			}
			f += D * k;
			break;
		case 10: // Ackley
			sum1 = 0;
			sum2 = 0;
			for( d = 0; d < D; d++ )
			{
				xd = xs.x[d];
				sum1 += xd * xd;
				sum2 += cos( 2 * pi * xd );
			}
			y = D;
			f = ( -20 * exp( -0.2 * sqrt( sum1 / y ) ) - exp( sum2 / y ) + 20 + E );
			break;
		case 13: // 2D Tripod function (Louis Gacogne)
			// Search [-100, 100]^2. min 0 on (0  -50)
			x1 = xs.x[0];
			x2 = xs.x[1];
			if( x2 < 0 )
			{
				f = fabs( x1 ) + fabs( x2 + 50 );
			}
			else
			{
				if( x1 < 0 )
					f = 1 + fabs( x1 + 50 ) + fabs( x2 - 50 );
				else
					f = 2 + fabs( x1 - 50 ) + fabs( x2 - 50 );
			}
			break;
		case 17: // KrishnaKumar
			f = 0;
			for( d = 0; d < D - 1; d++ )
			{
				f = f + sin( xs.x[d] + xs.x[d + 1] ) + sin( 2 * xs.x[d] * xs.x[d + 1] / 3 );
			}
			break;
		case 18: // Eason 2D (usually on [-100,100]
			// Minimum -1  on (pi,pi)
			x1 = xs.x[0]; x2 = xs.x[1];
			f = -cos( x1 ) * cos( x2 ) / exp( ( x1 - pi ) * ( x1 - pi ) + ( x2 - pi ) * ( x2 - pi ) );
			break;
		case 19: // mads contaminant transport model
			//Transform( xs.x, gop, xs.x );
			func_global( xs.x, gop, res );
			f = gop->phi;
			break;
	}
	return f;
}

