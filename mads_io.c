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

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>

#include "mads.h"

/* Functions here */
int parse_cmd( char *buf, struct calc_data *cd );
int load_problem( char *filename, int argn, char *argv[], struct opt_data *op );
int save_problem( char *filename, struct opt_data *op );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc2( char *filename, char *filename2, struct opt_data *op );
void compute_btc( char *filename, struct opt_data *op );
static char *strsave( const char *s, const char *lim );
char **shellpath( void );
void freeshellpath( char *shellpath[] );
unsigned maxpathlen( char *path[], const char *base );
void execvepath( char *path[], const char *base, char *const argv[], char *const envp[] );
int count_lines( char *filename );
int count_cols( char *filename, int row );

/* Functions elsewhere */
char **char_matrix( int maxCols, int maxRows );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
int Ftest( char *filename );
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );
int set_test_problems( struct opt_data *op );
void *malloc_check( const char *what, size_t n );

int parse_cmd( char *buf, struct calc_data *cd )
{
	int w;
	char *sep = " \t\n", *word;
	cd->opt_method = ( char * ) malloc( 50 * sizeof( char ) ); cd->opt_method[0] = 0;
	cd->smp_method = ( char * ) malloc( 50 * sizeof( char ) ); cd->smp_method[0] = 0;
	cd->paran_method = ( char * ) malloc( 50 * sizeof( char ) ); cd->paran_method[0] = 0;
	cd->infile = ( char * ) malloc( 255 * sizeof( char ) ); cd->infile[0] = 0;
	cd->restart_zip_file = ( char * ) malloc( 255 * sizeof( char ) ); cd->restart_zip_file[0] = 0;
	strcpy( cd->opt_method, "lm" );
	cd->problem_type = UNKNOWN;
	cd->calib_type = SIMPLE;
	cd->solution_type = EXTERNAL;
	cd->objfunc_type = SSR;
	cd->check_success = 0;
	cd->debug = 0;
	cd->fdebug = 0;
	cd->ldebug = 0;
	cd->leigen = 0;
	cd->pdebug = 0;
	cd->mdebug = 0;
	cd->odebug = 0;
	cd->insdebug = 0;
	cd->tpldebug = 0;
	cd->pardebug = 0;
	cd->sintrans = 1;
	cd->plogtrans = -1;
	cd->ologtrans = -1;
	cd->oweight = -1;
	cd->num_proc = -1;
	cd->restart = 1;
	cd->nreal = 0;
	cd->niter = 0;
	cd->init_particles = -1;
	cd->nretries = 0;
	cd->seed = cd->seed_init = 0;
	cd->maxeval = 5000;
	cd->paranoid = 0;
	cd->phi_cutoff = 0;
	cd->truth = 0.1;
	cd->test_func = -1;
	cd->test_func_dim = 0;
	cd->energy = 0;
	cd->ireal = 0;
	cd->sindx = 0.0000001;
	cd->lindx = 0.001;
	cd->pardx = 0;
	cd->pardomain = 100;
	cd->lmfactor = 10.0;
	for( word = strtok( buf, sep ); word; word = strtok( NULL, sep ) )
	{
		w = 0;
		if( strcasestr( word, "crea" ) ) { w = 1; cd->problem_type = CREATE; }
		if( strcasestr( word, "forw" ) ) { w = 1; cd->problem_type = FORWARD; }
		if( strcasestr( word, "cali" ) ) { w = 1; cd->problem_type = CALIBRATE; }
		if( strcasestr( word, "lsen" ) ) { w = 1; if( cd->problem_type == CALIBRATE ) cd->leigen = 1; else cd->problem_type = LOCALSENS; }
		if( strcasestr( word, "eige" ) ) { w = 1; if( cd->problem_type == CALIBRATE ) cd->leigen = 1; else cd->problem_type = EIGEN; }
		if( strcasestr( word, "mont" ) ) { w = 1; cd->problem_type = MONTECARLO; }
		if( strcasestr( word, "gsen" ) ) { w = 1; cd->problem_type = GLOBALSENS; }
		if( strcasestr( word, "abag" ) ) { w = 1; cd->problem_type = ABAGUS; }
		if( strcasestr( word, "infogap" ) ) { w = 1; cd->problem_type = INFOGAP; }
		if( strcasestr( word, "postpua" ) ) { w = 1; cd->problem_type = POSTPUA; }
		if( strcasestr( word, "sing" ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = SIMPLE; }
		if( strcasestr( word, "igpd" ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = IGPD; }
		if( strcasestr( word, "ppsd" ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = PPSD; }
		if( strcasestr( word, "igrnd" ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = IGRND; }
		if( strcasestr( word, "leig" ) ) { w = 1; cd->problem_type = CALIBRATE; cd->leigen = 1;  }
		if( strcasestr( word, "energy=" ) ) { w = 1; sscanf( word, "energy=%d", &cd->energy ); }
		if( strcasestr( word, "lmfactor=" ) ) { w = 1; sscanf( word, "lmfactor=%lf", &cd->lmfactor ); }
		if( strcasestr( word, "infile=" ) ) { w = 1; sscanf( word, "infile=%s", cd->infile ); }
		if( strcasestr( word, "real=" ) ) { w = 1; sscanf( word, "real=%d", &cd->nreal ); }
		if( strcasestr( word, "iter=" ) ) { w = 1; sscanf( word, "iter=%d", &cd->niter ); }
		if( strcasestr( word, "eval=" ) ) { w = 1; sscanf( word, "eval=%d", &cd->maxeval ); }
		if( strcasestr( word, "case=" ) ) { w = 1; sscanf( word, "case=%d", &cd->ireal ); }
		if( strcasestr( word, "retry" ) ) { w = 1; sscanf( word, "retry=%d", &cd->nretries ); }
		if( strcasestr( word, "particles" ) ) { w = 1; if( sscanf( word, "particles=%d", &cd->init_particles ) != 1 ) cd->init_particles = -1; }
		if( strcasestr( word, "opt=" ) ) { w = 1; if( cd->problem_type == UNKNOWN ) cd->problem_type = CALIBRATE; sscanf( word, "opt=%s", cd->opt_method ); }
		if( strcasestr( word, "smp=" ) ) { w = 1; sscanf( word, "smp=%s", cd->smp_method ); }
		if( strcasestr( word, "paran" ) ) { w = 1; cd->paranoid = 1; sscanf( word, "paran=%s", cd->paran_method ); }
		if( strcasestr( word, "nosin" ) ) { w = 1; cd->sintrans = 0; }
		if( strcasestr( word, "plog" ) ) { w = 1; if( sscanf( word, "plog=%d", &cd->plogtrans ) != 1 ) cd->plogtrans = 1; }
		if( strcasestr( word, "olog" ) ) { w = 1; if( sscanf( word, "olog=%d", &cd->ologtrans ) != 1 ) cd->ologtrans = 1; }
		if( strcasestr( word, "oweight" ) ) { w = 1; if( sscanf( word, "oweight=%d", &cd->oweight ) != 1 ) cd->oweight = 1; }
		if( strcasestr( word, "succ" ) ) { w = 1; cd->check_success = 1; }
		if( strcasestr( word, "cutoff=" ) ) { w = 1; sscanf( word, "cutoff=%lf", &cd->phi_cutoff ); }
		if( strcasestr( word, "truth=" ) ) { w = 1; sscanf( word, "truth=%lf", &cd->truth ); }
		if( strcasestr( word, "sindx=" ) ) { w = 1; cd->sintrans = 1; sscanf( word, "sindx=%lf", &cd->sindx ); if( cd->sindx < DBL_EPSILON ) cd->sindx = 0.0000001; }
		if( strcasestr( word, "lindx=" ) ) { w = 1; cd->sintrans = 0; sscanf( word, "lindx=%lf", &cd->lindx ); if( cd->lindx < DBL_EPSILON ) cd->lindx = 0.001; }
		if( strcasestr( word, "pardx" ) ) { w = 1; cd->sintrans = 0; sscanf( word, "pardx=%lf", &cd->pardx ); if( cd->pardx < DBL_EPSILON ) cd->pardx = 0.1; }
		if( strcasestr( word, "pardomain=" ) ) { w = 1; cd->sintrans = 0; sscanf( word, "pardomain=%lf", &cd->pardomain ); if( cd->pardomain < DBL_EPSILON ) cd->pardomain = 100; }
		if( strcasestr( word, "seed=" ) ) { w = 1; sscanf( word, "seed=%d", &cd->seed ); cd->seed_init = cd->seed; }
		if( strcasestr( word, "np" ) ) { w = 1; cd->num_proc = 0; sscanf( word, "np=%d", &cd->num_proc ); if( cd->num_proc <= 0 ) cd->num_proc = 0; }
		if( strcasestr( word, "restart" ) ) { w = 1; sscanf( word, "restart=%d", &cd->restart ); if( cd->restart < 0 || cd->restart > 1 ) cd->restart = -1; }
		if( strcasestr( word, "rstfile=" ) ) { w = 1; sscanf( word, "rstfile=%s", cd->restart_zip_file ); cd->restart = -1; }
		if( strncasecmp( word, "debug", 5 ) == 0 ) { w = 1; if( sscanf( word, "debug=%d", &cd->debug ) == 0 || cd->debug == 0 ) cd->debug = 1; } // Global debug
		if( strcasestr( word, "fdebug" ) ) { w = 1; sscanf( word, "fdebug=%d", &cd->fdebug ); if( cd->fdebug == 0 ) cd->fdebug = 1; }
		if( strcasestr( word, "ldebug" ) ) { w = 1; sscanf( word, "ldebug=%d", &cd->ldebug ); if( cd->ldebug == 0 ) cd->ldebug = 1; }
		if( strcasestr( word, "pdebug" ) ) { w = 1; sscanf( word, "pdebug=%d", &cd->pdebug ); if( cd->pdebug == 0 ) cd->pdebug = 1; }
		if( strcasestr( word, "mdebug" ) ) { w = 1; sscanf( word, "mdebug=%d", &cd->mdebug ); if( cd->mdebug == 0 ) cd->mdebug = 1; }
		if( strcasestr( word, "odebug" ) ) { w = 1; sscanf( word, "odebug=%d", &cd->odebug ); if( cd->odebug != 1 ) cd->odebug = 1; }
		if( strcasestr( word, "insdebug" ) ) { w = 1; sscanf( word, "insdebug=%d", &cd->insdebug ); if( cd->insdebug == 0 ) cd->insdebug = 1; }
		if( strcasestr( word, "tpldebug" ) ) { w = 1; sscanf( word, "tpldebug=%d", &cd->tpldebug ); if( cd->tpldebug == 0 ) cd->tpldebug = 1; }
		if( strcasestr( word, "pardebug" ) ) { w = 1; sscanf( word, "pardebug=%d", &cd->pardebug ); if( cd->pardebug == 0 ) cd->pardebug = 1; }
		if( strcasecmp( word, "ssr" ) == 0 ) { w = 1; cd->objfunc_type = SSR; }
		if( strcasecmp( word, "ssd0" ) == 0 ) { w = 1; cd->objfunc_type = SSD0; }
		if( strcasecmp( word, "ssda" ) == 0 ) { w = 1; cd->objfunc_type = SSDA; }
		if( strcasecmp( word, "ssdr" ) == 0 ) { w = 1; cd->objfunc_type = SSDR; }
		if( strcasestr( word, "test" ) ) { w = 1; cd->test_func = 1; cd->test_func_dim = 2; sscanf( word, "test=%d", &cd->test_func ); ( *cd ).solution_type = TEST; }
		if( strcasestr( word, "dim=" ) ) { w = 1; sscanf( word, "dim=%d", &cd->test_func_dim ); if( cd->test_func_dim < 2 ) cd->test_func_dim = 2; }
		if( strcasestr( word, "npar=" ) ) { w = 1; sscanf( word, "npar=%d", &cd->test_func_npar ); if( cd->test_func_npar < 2 ) cd->test_func_npar = 2; }
		if( strcasestr( word, "nobs=" ) ) { w = 1; sscanf( word, "nobs=%d", &cd->test_func_nobs ); if( cd->test_func_nobs < 2 ) cd->test_func_nobs = 2; }
		if( strcasestr( word, "poi" ) ) { w = 1; ( *cd ).solution_type = POINT; }
		if( strcasestr( word, "rec" ) ) { w = 1; if( strcasestr( word, "ver" ) )( *cd ).solution_type = PLANE3D; else( *cd ).solution_type = PLANE; }
		if( strcasestr( word, "box" ) ) { w = 1; ( *cd ).solution_type = BOX; }
		if( w == 0 ) { printf( "\nERROR: Unknown keyword \'%s\'!\nExecute 'mads' without arguments to list acceptable keywords!\n", word ); return( -1 ); }
	}
	if( cd->seed != 0 ) cd->seed *= -1; // Modify the seed to show that is imported
	if( cd->problem_type == UNKNOWN ) { cd->problem_type = CALIBRATE; cd->calib_type = SIMPLE; }
	if( ( cd->problem_type == MONTECARLO || cd->calib_type == IGRND || cd->problem_type == GLOBALSENS || cd->problem_type == ABAGUS ) && cd->nreal == 0 ) cd->nreal = 100;
	if( cd->nretries > 0 && cd->problem_type == CALIBRATE && strncasecmp( cd->opt_method, "lm", 2 ) == 0 ) { strcat( cd->opt_method, "_paran" ); cd->paranoid = 1; }
	if( cd->test_func )
	{
		if( cd->test_func < 40 ) printf( "Test Function #%d Dimensionality %d\n", cd->test_func, cd->test_func_dim );
		else
		{
			if( cd->test_func == 41 ) cd->test_func_npar = 4;
			else cd->test_func_npar = 2;
			if( cd->test_func_nobs < 2 ) cd->test_func_nobs = 100;
			printf( "Test Function #%d Parameters %d Observations %d\n", cd->test_func, cd->test_func_npar, cd->test_func_nobs );
		}
	}
	printf( "Problem type: " );
	switch( cd->problem_type )
	{
		case CREATE: printf( "create a calibration input file based on forward run (no calibration)" ); break;
		case FORWARD: printf( "forward run (no calibration)" ); break;
		case CALIBRATE: printf( "calibration" ); break;
		case LOCALSENS: printf( "sensitivity analysis" ); break;
		case EIGEN: printf( "eigen analysis" ); break;
		case MONTECARLO: printf( "monte-carlo analysis (realizations = %d)", cd->nreal ); break;
		case GLOBALSENS: printf( "global sensitivity analysis (realizations = %d)", cd->nreal ); break;
		case ABAGUS: printf( "abagus: agent-based global uncertainty and sensitivity analysis" ); break;
		case INFOGAP: printf( "Info-gap decision analysis" ); break;
		case POSTPUA: printf( "predictive uncertainty analysis of sampling results" ); break;
		default: printf( "WARNING: unknown problem type; calibration assumed" ); cd->problem_type = CALIBRATE; break;
	}
	printf( "\n" );
	if( cd->problem_type == CALIBRATE )
	{
		printf( "\nCalibration technique: " );
		switch( cd->calib_type )
		{
			case IGRND: printf( "sequential calibration using a set of random initial values (realizations = %d)", cd->nreal ); break;
			case IGPD: printf( "sequential calibration using a set discretized initial values" ); break;
			case PPSD: printf( "sequential calibration using partial parameter parameter discretization" ); break;
			case SIMPLE: printf( "single calibration using initial guesses provided in the input file" ); break;
			default: printf( "WARNING: unknown calibration type ASSUMED: single calibration using initial guesses provided in the input file" ); cd->calib_type = SIMPLE; break;
		}
		printf( "\nOptimization method: opt=%s | ", cd->opt_method );
		if( strncasecmp( cd->opt_method, "squad", 5 ) == 0 || ( strcasestr( cd->opt_method, "pso" ) && strcasestr( cd->opt_method, "lm" ) ) )
			printf( "SQUADS: Coupled Particle-Swarm and Levenberg-Marquardt optimization\n" );
		else if( strncasecmp( cd->opt_method, "lm", 2 ) == 0 ) { printf( "Levenberg-Marquardt optimization\n" ); if( cd->calib_type == SIMPLE ) cd->leigen = 1; }
		else if( strcasestr( cd->opt_method, "pso" ) || strncasecmp( cd->opt_method, "swarm", 5 ) == 0 || strncasecmp( cd->opt_method, "tribe", 5 ) == 0 )
			printf( "Particle-Swarm optimization\n" );
		else { printf( "WARNING: Unknown method! Levenberg-Marquardt optimization assumed\n" ); strcpy( cd->opt_method, "lm" ); }
		if( cd->nretries > 0 ) printf( "Number of calibration retries = %d\n", cd->nretries );
		if( cd->niter < 0 ) cd->niter = 0;
		if( cd->niter > 0 ) printf( "Number of Levenberg-Marquardt iterations = %d\n", cd->niter );
		else printf( "Number of Levenberg-Marquardt iterations = will be computed internally\n" );
		if( strcasestr( cd->opt_method, "apso" ) || strcasestr( cd->opt_method, "tribe" ) || strcasestr( cd->opt_method, "squad" ) )
		{
			if( cd->init_particles > 1 ) printf( "Number of particles = %d\n", cd->init_particles );
			if( cd->init_particles == -1 ) printf( "Number of particles = will be computed internally\n" );
		}
		if( cd->leigen == 1 ) printf( "Eigen analysis will be performed of the final optimization results\n" );
		printf( "\nGlobal termination criteria:\n" );
		printf( "1: Maximum number of evaluations = %d\n", cd->maxeval );
		printf( "2: Objective function cutoff value: " );
		if( cd->phi_cutoff <= DBL_EPSILON ) printf( "NOT implemented (ADD keyword cutoff=[value] to implement)\n" );
		else printf( "%g\n", cd->phi_cutoff );
		printf( "3: Observations within defined calibration ranges: " );
		if( cd->check_success == 1 ) printf( "implemented (keyword success)\n" );
		else printf( "NOT implemented (ADD keyword success to implement)\n" );
		printf( "Objectve function: " );
		switch( cd->objfunc_type )
		{
			case SSR: printf( "sum of squared residuals" ); break;
			case SSDR: printf( "sum of squared discrepancies and squared residuals" ); break;
			case SSDA: printf( "sum of squared discrepancies and residuals" ); break;
			case SSD0: printf( "sum of squared discrepancies" ); break;
			default: printf( "unknown value; sum of squared residuals assumed" ); cd->objfunc_type = SSR; break;
		}
		printf( "\n" );
	}
	if( cd->sintrans == 0 ) printf( "\nSin transformation of the model parameters is not applied!\n" );
	else printf( "\nSin transformation of the model parameters is applied!\n" );
	if( cd->plogtrans == 1 ) printf( "\nLog transformation enforced on all parameters!\n" );
	else if( cd->plogtrans == 0 ) printf( "\nLog transformation is not applied to any parameters!\n" );
	if( cd->ologtrans == 1 ) printf( "\nLog transformation enforced on all calibration targets!\n" );
	else if( cd->ologtrans == 0 ) printf( "\nLog transformation is not applied to any calibration targets!\n" );
	if( cd->oweight == 1 ) printf( "\nUnit residual weights for all calibration targets!\n" );
	else if( cd->oweight == 0 ) printf( "\nWARNING: Zero residual weights for all calibration targets (for testing and manupulation purposes)!\n" );
	else if( cd->oweight == 2 ) printf( "\nResidual weights reversely propotional to observations for all calibration targets!\n" );
	if( cd->problem_type == ABAGUS )
	{
		if( cd->energy == 0 ) { cd->energy = 10000; printf( "\nInitial particle energy set to default value: %d\n", cd->energy );}
		else printf( "\nInitial particle energy set to: %d\n", cd->energy );
		cd->sintrans = 0; printf( "\nsine tranformation disabled for ABAGUS runs" );
	}
	if( cd->problem_type == ABAGUS && cd->infile[0] != 0 ) { printf( "\nResults in %s to be read into kdtree\n", cd->infile );}
	if( cd->problem_type == INFOGAP && cd->infile[0] == 0 ) { printf( "\nInfile must be specified for infogap run\n" ); exit( 0 ); }
	if( cd->problem_type == POSTPUA && cd->infile[0] == 0 ) { printf( "\nInfile must be specified for postpua run\n" ); exit( 0 ); }
	if( cd->smp_method[0] != 0 )
	{
		printf( "\nSampling method: " );
		if( strncasecmp( cd->smp_method, "olhs", 4 ) == 0 ) printf( "Optimal Latin Hyper Cube (LHS) (if real <= 500 IDLHS; if real > 500 LHS)\n" );
		else if( strncasecmp( cd->smp_method, "lhs", 3 ) == 0 ) printf( "Standard Latin Hyper Cube (LHS)\n" );
		else if( strncasecmp( cd->smp_method, "idlhs", 5 ) == 0 ) printf( "Improved Distance Latin Hyper Cube (IDLHS)\n" );
		else if( strncasecmp( cd->smp_method, "random", 5 ) == 0 ) printf( "Pure random\n" );
		else { printf( "WARNING: Unknown (rnd=%s); Optimal Latin Hyper Cube selected (if real <= 500 IDLHS; if real > 500 LHS)\n", cd->smp_method ); strcpy( cd->smp_method, "olhs" ); }
	}
	if( cd->paran_method[0] != 0 )
	{
		printf( "\nParanoid Sampling method: " );
		if( strncasecmp( cd->paran_method, "olhs", 4 ) == 0 ) printf( "Optimal Latin Hyper Cube (LHS) (if real <= 500 IDLHS; if real > 500 LHS)\n" );
		else if( strncasecmp( cd->paran_method, "lhs", 3 ) == 0 ) printf( "Standard Latin Hyper Cube (LHS)\n" );
		else if( strncasecmp( cd->paran_method, "idlhs", 5 ) == 0 ) printf( "Improved Distance Latin Hyper Cube (IDLHS)\n" );
		else if( strncasecmp( cd->paran_method, "random", 5 ) == 0 ) printf( "Pure random\n" );
		else { printf( "WARNING: Unknown (rnd=%s); Optimal Latin Hyper Cube selected (if real <= 500 IDLHS; if real > 500 LHS)\n", cd->paran_method ); strcpy( cd->paran_method, "olhs" ); }
	}
	printf( "\nDebug level (general): debug=%d\n", cd->debug );
	if( cd->debug || cd->fdebug ) printf( "Debug level for functional evaluations: fdebug=%d\n", cd->fdebug );
	if( cd->debug && cd->problem_type == CALIBRATE )
	{
		printf( "Debug level for Levenberg-Marquardt optimization progress: ldebug= %d\n", cd->ldebug );
		printf( "Debug level for Particle-Swarm optimization progress: pdebug= %d\n", cd->pdebug );
		printf( "Debug level for objective function progress: odebug=%d\n", cd->odebug );
	}
	if( ( cd->debug || cd->mdebug ) && cd->problem_type != CREATE && cd->problem_type != EIGEN )
		printf( "Debug level for random sets: mdebug=%d\n", cd->mdebug );
	if( ( cd->debug || cd->pardebug ) && cd->num_proc > 1 )
		printf( "Debug level for parallel execution: pardebug=%d\n", cd->pardebug );
	if( ( cd->debug || cd->tpldebug || cd->insdebug ) && cd->solution_type == EXTERNAL )
	{
		printf( "Debug level for template file: tpldebug=%d\n", cd->tpldebug );
		printf( "Debug level for instruction file: insdebug=%d\n", cd->insdebug );
	}
	return( 1 );
}

int load_problem( char *filename, int argn, char *argv[], struct opt_data *op )
{
	FILE *infile;
//	FILE *infileb;
	double x0, y0, x, y, d, alpha, beta;
	char buf[1000], *file, **path, exec[1000];
	int  i, j, k, bad_data, nofile = 0, skip = 0;
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od;
	struct obs_data *preds;
	struct well_data *wd;
	struct grid_data *gd;
	struct extrn_data *ed;
	cd = op->cd;
	pd = op->pd;
	od = op->od;
	preds = op->preds;
	wd = op->wd;
	gd = op->gd;
	ed = op->ed;
	bad_data = 0;
	if( ( infile = fopen( filename, "r" ) ) == NULL )
	{
		printf( "File \'%s\' cannot be opened to read problem information!\n", filename );
		sprintf( filename, "%s.in", op->root );
		if( ( infile = fopen( filename, "r" ) ) == NULL )
			nofile = 1;
		else
			printf( "WARNING: File \'%s\' is opened to read problem information!\n", filename );
	}
	if( nofile == 0 )
	{
		// Read commands in the file
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ":" ); fscanf( infile, "%[^\n]s\n", buf ); fscanf( infile, "\n" );
		if( sscanf( buf, "%i", &( *pd ).nParam ) == 1 ) { buf[0] = 0; skip = 1; }
	}
	else buf[0] = 0;
	// Add commands provided as arguments
	for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
	if( parse_cmd( buf, cd ) == -1 ) { sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	// Read Solution Type
	cd->solution_id = ( char * ) malloc( 50 * sizeof( char ) ); ( *cd ).solution_id[0] = 0;
	if( nofile == 0 && skip == 0 ) { fscanf( infile, "%[^:]s", buf ); fscanf( infile, ":" ); fscanf( infile, "%s\n", ( *cd ).solution_id ); sscanf( ( *cd ).solution_id, "%d",  &( *cd ).solution_type ); }
	if( strcasestr( ( *cd ).solution_id, "ext" ) )( *cd ).solution_type = EXTERNAL;
	if( strcasestr( ( *cd ).solution_id, "poi" ) )( *cd ).solution_type = POINT;
	if( strcasestr( ( *cd ).solution_id, "rec" ) ) { if( strcasestr( ( *cd ).solution_id, "ver" ) )( *cd ).solution_type = PLANE3D; else( *cd ).solution_type = PLANE; }
	if( strcasestr( ( *cd ).solution_id, "box" ) )( *cd ).solution_type = BOX;
	if( strcasestr( ( *cd ).solution_id, "test" ) || ( *cd ).test_func >= 0 ) { ( *cd ).solution_type = TEST; od->nObs = 0; }
	if( nofile )
	{
		if( ( *cd ).solution_type != TEST || ( *cd ).problem_type == INFOGAP )
		{
			printf( "\nERROR: Input file is needed!\n\n" );
			return( -1 );
		}
		else
		{
			printf( "\nWARNING: Input file is not needed because this is a test problem!\n" );
			//	return( 1 );
		}
	}
	printf( "\nModel: " );
	switch( ( *cd ).solution_type )
	{
		case EXTERNAL: { printf( "external" ); strcpy( ( *cd ).solution_id, "external" ); break; }
		case POINT: { printf( "internal | point contaminant source" ); strcpy( ( *cd ).solution_id, "point" ); break; }
		case PLANE: { printf( "internal | rectangular contaminant source" ); strcpy( ( *cd ).solution_id, "rect" ); break; }
		case PLANE3D: { printf( "internal | rectangular contaminant source with vertical flow component" ); strcpy( ( *cd ).solution_id, "rect_vert" ); break; }
		case BOX: { printf( "internal | box contaminant source" ); strcpy( ( *cd ).solution_id, "box" ); break; }
		case TEST: { printf( "internal | test optimization problem #%d: ", ( *cd ).test_func ); set_test_problems( op ); printf( " dimensionality %d", ( *pd ).nOptParam ); strcpy( ( *cd ).solution_id, "test" ); break; }
		default: printf( "WARNING! UNDEFINED model type!" );
	}
	printf( "\n" );
	if( ( *cd ).solution_type == TEST ) return( 1 );
	// Read parameters
	if( skip == 0 ) fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &( *pd ).nParam );
	printf( "\nNumber of model parameters: %d\n", ( *pd ).nParam );
	pd->var_id = char_matrix( ( *pd ).nParam, 50 );
	pd->var = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	cd->var = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc( ( *pd ).nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc( ( *pd ).nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	( *pd ).nOptParam = ( *pd ).nFlgParam = 0;
	for( i = 0; i < ( *pd ).nParam; i++ )
	{
		fscanf( infile, "%[^:]s", pd->var_id[i] );
		if( cd->debug ) printf( "%-26s: ", pd->var_id[i] );
		fscanf( infile, ": %lf %d %d %lf %lf %lf\n", &( *pd ).var[i], &( *pd ).var_opt[i], &( *pd ).var_log[i], &( *pd ).var_dx[i], &( *pd ).var_min[i], &( *pd ).var_max[i] );
		( *cd ).var[i] = ( *pd ).var[i];
		if( cd->debug ) printf( "init %17.12g opt %1d log %1d step %7g min %7g max %7g\n", ( *pd ).var[i], ( *pd ).var_opt[i], ( *pd ).var_log[i], ( *pd ).var_dx[i], ( *pd ).var_min[i], ( *pd ).var_max[i] );
		if( ( *pd ).var_opt[i] == 1 )( *pd ).nOptParam++;
		else if( ( *pd ).var_opt[i] == 2 )
		{
			( *pd ).nFlgParam++;
			if( ( *cd ).problem_type == CALIBRATE && ( *cd ).calib_type != PPSD )
				( *pd ).nOptParam++;
		}
		if( ( *pd ).var_opt[i] >= 1 )
		{
			if( pd->var_max[i] < pd->var[i] || pd->var_min[i] > pd->var[i] )
			{
				printf( "ERROR: Parameter initial value is outside the specified min/max range! " );
				printf( "Parameter %s: %g min %g max %g\n", pd->var_id[i], pd->var[i], pd->var_min[i], pd->var_max[i] );
				bad_data = 1;
			}
			if( pd->var_max[i] <= pd->var_min[i] )
			{
				printf( "ERROR: Parameter min/max range is not correctly specified! " );
				printf( "Parameter %s: min %g max %g\n", pd->var_id[i], pd->var_min[i], pd->var_max[i] );
				bad_data = 1;
			}
			if( cd->plogtrans == 1 )( *pd ).var_log[i] = 1;
			else if( cd->plogtrans == 0 )( *pd ).var_log[i] = 0;
			if( ( *pd ).var_log[i] == 1 )
			{
				if( pd->var_min[i] < 0 || pd->var[i] < 0 )
				{
					printf( "ERROR: Parameter cannot be log transformed\n" );
					printf( "Parameter %s: min %g max %g\n", pd->var_id[i], pd->var_min[i], pd->var_max[i] );
					if( cd->plogtrans ) { pd->var_log[i] = 0; pd->var_range[i] = pd->var_max[i] - pd->var_min[i]; continue; }
					else bad_data = 1;
				}
				if( pd->var_dx[i] < 2 ) d = ( pd->var_max[i] - pd->var_min[i] ) / pd->var_dx[i];
				if( pd->var[i] < DBL_EPSILON ) pd->var[i] = DBL_EPSILON;
				if( pd->var_min[i] < DBL_EPSILON ) pd->var_min[i] = DBL_EPSILON;
				pd->var[i] = log10( pd->var[i] );
				pd->var_min[i] = log10( pd->var_min[i] );
				pd->var_max[i] = log10( pd->var_max[i] );
				if( pd->var_dx[i] < 2 ) pd->var_dx[i] = ( pd->var_max[i] - pd->var_min[i] ) / d;
				else pd->var_dx[i] = log10( pd->var_dx[i] );
			}
			pd->var_range[i] = pd->var_max[i] - pd->var_min[i];
			if( pd->var_dx[i] > DBL_EPSILON ) cd->pardx = 1; // discretizations is ON
		}
	}
	if( bad_data ) { sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if( ( *cd ).problem_type == CALIBRATE && ( *cd ).calib_type == PPSD && ( *pd ).nFlgParam == 0 )
	{
		printf( "WARNING: Partial parameter-space discretization (PPSD) is selected.\nHowever no parameters are flagged!\nSingle calibration will be performed using the initial guesses provided in the input file!\n" );
		( *cd ).calib_type = SIMPLE;
	}
	pd->var_index = ( int * ) malloc( ( *pd ).nOptParam * sizeof( int ) );
	if( cd->debug ) printf( "\n" );
	printf( "Number of optimized parameters = %d\n", ( *pd ).nOptParam );
	for( k = i = 0; i < ( *pd ).nParam; i++ )
		if( ( *pd ).var_opt[i] == 1 || ( ( *pd ).var_opt[i] > 1 && ( *cd ).calib_type != PPSD ) )
		{
			if( cd->debug ) printf( "%-26s: init %17.12g step %8.3g min %7g max %7g\n", pd->var_id[i], pd->var[i], ( *pd ).var_dx[i], ( *pd ).var_min[i], ( *pd ).var_max[i] );
			( *pd ).var_index[k++] = i;
		}
	if( cd->debug ) printf( "\n" );
	printf( "Number of flagged parameters = %d\n", ( *pd ).nFlgParam );
	for( i = 0; i < ( *pd ).nParam; i++ )
		if( ( *pd ).var_opt[i] == 2 )
			if( cd->debug ) printf( "%-26s: init %17.12g step %6g min %6g max %6g\n", pd->var_id[i], pd->var[i], ( *pd ).var_dx[i], ( *pd ).var_min[i], ( *pd ).var_max[i] );
	if( ( *pd ).nParam == ( *pd ).nOptParam && cd->debug ) printf( "\nNO fixed parameters\n" );
	else
	{
		if( cd->debug ) printf( "\n" );
		printf( "Number of fixed parameters = %d\n", ( *pd ).nParam - ( *pd ).nOptParam );
		for( i = 0; i < ( *pd ).nParam; i++ )
			if( ( *pd ).var_opt[i] == 0 )
				if( cd->debug ) printf( "%-26s: %g\n", pd->var_id[i], ( *pd ).var[i] );
	}
	/*
		if( (*cd).problem_type != 0 )
		{
			// Read binary data if available
			sprintf( buf, "%s.bin", filename );
			if ( ( infileb = fopen( buf, "r" ) ) == NULL )
				printf( "Binary file %s cannot be opened to read problem information!\n", buf );
			else
			{
				printf( "Binary file %s is found and parameter values are read\n", buf );
				fread( (*pd).var, sizeof((*pd).var[i]), (*pd).nParam, infileb);
				for( i = 0; i < (*pd).nParam; i++ )
						printf( "%-26s: binary init %15.12g\n", pd->var_id[i], (*pd).var[i] );
				fclose( infileb );
			}
		}
	*/
	pd->var_current = ( double * ) malloc( ( *pd ).nOptParam * sizeof( double ) );
	pd->var_best = ( double * ) malloc( ( *pd ).nOptParam * sizeof( double ) );
	if( cd->solution_type == EXTERNAL )
	{
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &od->nObs );
		od->nTObs = od->nObs;
		od->obs_id = char_matrix( ( *od ).nObs, 50 );
		od->obs_target = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_weight = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_min = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_max = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_current = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_best = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->res = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_log = ( int * ) malloc( ( *od ).nObs * sizeof( int ) );
		for( i = 0; i < od->nObs; i++ )
		{
			od->obs_min[i] = -1e6; od->obs_max[i] = 1e6; od->obs_weight[i] = 1; od->obs_log[i] = 0;
			fscanf( infile, "%s %lf %lf %d %lf %lf\n", od->obs_id[i], &od->obs_target[i], &od->obs_weight[i], &od->obs_log[i], &od->obs_min[i], &od->obs_max[i] );
			if( od->obs_max[i] < od->obs_target[i] || od->obs_min[i] > od->obs_target[i] )
			{
				printf( "ERROR: Observation target is outside the specified min/max range! " );
				printf( "Observation %s: %g min %g max %g\n", od->obs_id[i], od->obs_target[i], od->obs_min[i], od->obs_max[i] );
				bad_data = 1;
			}
			if( od->obs_max[i] <= od->obs_min[i] )
			{
				printf( "ERROR: Calibration range is not correctly specified! " );
				printf( "Observation %s: min %g max %g\n", od->obs_id[i], od->obs_min[i], od->obs_max[i] );
				bad_data = 1;
			}
		}
		if( bad_data ) return( 0 );
		if( cd->ologtrans == 1 ) od->obs_log[i] = 1;
		else if( cd->ologtrans == 0 ) od->obs_log[i] = 0;
		if( cd->oweight == 1 )( *od ).obs_weight[i] = 1;
		else if( cd->oweight == 0 )( *od ).obs_weight[i] = 0;
		else if( cd->oweight == 2 ) { if( abs( ( *od ).obs_target[i] ) > DBL_EPSILON )( *od ).obs_weight[i] = ( double ) 1.0 / ( *od ).obs_target[i]; else( *od ).obs_weight[i] = HUGE_VAL; }
		if( cd->debug ) printf( "\n" );
		printf( "Number of calibration targets = %d\n", ( *od ).nObs );
		if( cd->debug )
		{
			for( i = 0; i < od->nObs; i++ )
			{
				if( od->nObs < 50 || ( i < 20 || i > od->nObs - 20 ) )
					printf( "%-20s: %15g weight %12g log %1d acceptable range: min %15g max %15g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
				if( od->nObs > 50 && i == 21 ) printf( "...\n" );
			}
		}
		ed->cmdline = ( char * ) malloc( 80 * sizeof( char ) );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": " ); fgets( ed->cmdline, 80, infile );
		ed->cmdline[strlen( ed->cmdline ) - 1] = 0;
		printf( "Execution command: %s\n", ed->cmdline );
		if( sscanf( ed->cmdline, "%i", &i ) == -1 )
		{
			printf( "ERROR: Execution command is not valid!\n" );
			sprintf( buf, "rm -f %s.running", op->root ); system( buf );
			exit( 1 );
		}
		strcpy( buf, ed->cmdline );
		file = &buf[0];
		file = strsep( &file, " \t" );
		k = 1;
		if( access( file, X_OK ) == -1 )
		{
			k = 0;
			if( file[0] == '~' )
			{
				sprintf( exec, "%s/%s", getenv( "HOME" ), &file[1] );
				if( access( exec, X_OK ) == 0 )
					k = 1;
			}
			if( k == 0 )
			{
				path = shellpath();
				for( i = 0; k == 0 && path[i]; i++ )
				{
					if( cd->debug > 2 ) printf( "%s\n", path[i] );
					sprintf( exec, "%s/%s", path[i], file );
					if( access( exec, X_OK ) == 0 )
						k = 1;
				}
			}
		}
		if( k == 0 )
		{
			printf( "WARNING: Program \'%s\' does not exist or cannot be executed!\n", file );
			exit( 1 );
		}
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %d\n", &ed->ntpl );
		ed->fn_tpl = char_matrix( ed->ntpl, 80 );
		ed->fn_out = char_matrix( ed->ntpl, 80 );
		for( i = 0; i < ed->ntpl; i++ )
		{
			fscanf( infile, "%s %s\n", ed->fn_tpl[i], ed->fn_out[i] );
			if( access( ed->fn_tpl[i], R_OK ) == -1 )
			{
				printf( "ERROR: File \'%s\' does not exist!\n", ed->fn_tpl[i] );
				exit( 1 );
			}
		}
		printf( "External files:\n" );
		printf( "- to provide current model parameters:\n" );
		for( i = 0; i < ed->ntpl; i++ )
			printf( "%s -> %s\n", ed->fn_tpl[i], ed->fn_out[i] );
		fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %d\n", &ed->nins );
		ed->fn_ins = char_matrix( ed->nins, 80 );
		ed->fn_obs = char_matrix( ed->nins, 80 );
		for( i = 0; i < ed->nins; i++ )
		{
			fscanf( infile, "%s %s\n", ed->fn_ins[i], ed->fn_obs[i] );
			if( access( ed->fn_ins[i], R_OK ) == -1 )
			{
				printf( "ERROR: File \'%s\' does not exist!\n", ed->fn_ins[i] );
				exit( 1 );
			}
		}
		printf( "- to read current model predictions:\n" );
		for( i = 0; i < ed->nins; i++ )
			printf( "%s <- %s\n", ed->fn_ins[i], ed->fn_obs[i] );
		fclose( infile );
		( *gd ).min_t = ( *gd ).time = 0;
		return( 1 ); // EXIT; Done with external problem
	}
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i\n", &( *wd ).nW );
	wd->id = char_matrix( ( *wd ).nW, 40 );
	wd->x = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->y = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->z1 = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->z2 = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->xa = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->ya = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->za1 = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->za2 = ( double * ) malloc( ( *wd ).nW * sizeof( double ) );
	wd->nWellObs = ( int * ) malloc( ( *wd ).nW * sizeof( int ) );
	wd->obs_target = ( double ** ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->obs_log = ( int ** ) malloc( ( *wd ).nW * sizeof( int * ) );
	wd->obs_time = ( double ** ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->obs_weight = ( double ** ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->obs_min = ( double ** ) malloc( ( *wd ).nW * sizeof( double * ) );
	wd->obs_max = ( double ** ) malloc( ( *wd ).nW * sizeof( double * ) );
	( *od ).nObs = ( *preds ).nObs = ( *od ).nTObs = 0;
	if( cd->debug ) printf( "\nObservation data:\n" );
	for( i = 0; i < ( *wd ).nW; i++ )
	{
		fscanf( infile, "%s %lf %lf %lf %lf %i ", wd->id[i], &( *wd ).x[i], &( *wd ).y[i], &( *wd ).z1[i], &( *wd ).z2[i], &( *wd ).nWellObs[i] );
		if( cd->debug ) printf( "Well %-6s x %8g y %8g z0 %6g z1 %6g nObs %2i ", wd->id[i], wd->x[i], ( *wd ).y[i], ( *wd ).z1[i], ( *wd ).z2[i], ( *wd ).nWellObs[i] );
		if( ( *wd ).nWellObs[i] <= 0 ) { if( cd->debug ) printf( "Warning: no observations\n" ); fscanf( infile, "%lf %lf %lf %i %lf %lf\n", &d, &d, &d, &j, &d, &d ); continue; }
		wd->obs_target[i] = ( double * ) malloc( ( *wd ).nWellObs[i] * sizeof( double ) );
		wd->obs_time[i] = ( double * ) malloc( ( *wd ).nWellObs[i] * sizeof( double ) );
		wd->obs_log[i] = ( int * ) malloc( ( *wd ).nWellObs[i] * sizeof( int ) );
		wd->obs_weight[i] = ( double * ) malloc( ( *wd ).nWellObs[i] * sizeof( double ) );
		wd->obs_min[i] = ( double * ) malloc( ( *wd ).nWellObs[i] * sizeof( double ) );
		wd->obs_max[i] = ( double * ) malloc( ( *wd ).nWellObs[i] * sizeof( double ) );
		for( j = 0; j < ( *wd ).nWellObs[i]; j++ )
		{
			wd->obs_min[i][j] = -1e6; wd->obs_max[i][j] = 1e6; wd->obs_weight[i][j] = 1; wd->obs_log[i][j] = 0;
			fscanf( infile, "%lf %lf %lf %i %lf %lf\n", &( *wd ).obs_time[i][j], &( *wd ).obs_target[i][j], &( *wd ).obs_weight[i][j], &( *wd ).obs_log[i][j], &( *wd ).obs_min[i][j], &( *wd ).obs_max[i][j] );
			if( cd->ologtrans == 1 )( *wd ).obs_log[i][j] = 1;
			else if( cd->ologtrans == 0 )( *wd ).obs_log[i][j] = 0;
			if( cd->oweight == 1 )( *wd ).obs_weight[i][j] = 1;
			else if( cd->oweight == 0 )( *wd ).obs_weight[i][j] = 0;
			else if( cd->oweight == 2 ) { if( abs( ( *wd ).obs_target[i][j] ) > DBL_EPSILON )( *wd ).obs_weight[i][j] = ( double ) 1.0 / ( *wd ).obs_target[i][j]; else( *wd ).obs_weight[i][j] = HUGE_VAL; }
			if( cd->debug ) printf( "t %5g c %5g weight %12g log %1d acceptable range: min %5g max %5g\n", ( *wd ).obs_time[i][j], ( *wd ).obs_target[i][j], ( *wd ).obs_weight[i][j], ( *wd ).obs_log[i][j], ( *wd ).obs_min[i][j], ( *wd ).obs_max[i][j] );
			if( wd->obs_max[i][j] < wd->obs_target[i][j] || wd->obs_min[i][j] > wd->obs_target[i][j] )
			{
				printf( "ERROR: Observation target is outside the specified min/max range! " );
				printf( "Observation %s(%g): %g min %g max %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
				bad_data = 1;
			}
			if( wd->obs_max[i][j] <= wd->obs_min[i][j] )
			{
				printf( "ERROR: Calibration range is not correctly specified! " );
				printf( "Observation %s(%g): min %g max %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
				bad_data = 1;
			}
			if( ( *wd ).obs_weight[i][j] > 0 )( *od ).nObs++;
			if( ( *wd ).obs_weight[i][j] < 0 )( *preds ).nObs++;
			( *od ).nTObs++;
			if( j + 1 < ( *wd ).nWellObs[i] ) { fscanf( infile, "\t\t" ); if( cd->debug ) printf( "\t\t\t\t\t\t\t      " ); }
		}
	}
	if( bad_data ) return( 0 );
	if( pd->nParam == 0 || ( pd->nOptParam == 0 && pd->nFlgParam == 0 ) ) { printf( "\nERROR: Number of model parameters is zero!\n\n" ); return( 0 ); }
	if( od->nObs == 0 )
	{
		if( cd->problem_type != FORWARD && cd->problem_type != MONTECARLO )
		{ printf( "\nERROR: Number of calibration targets is equal to zero!\n\n" ); return( 0 ); }
		else printf( "\nWARNING: Number of calibration targets is equal to zero!\n\n" );
	}
	if( cd->debug > 2 )
	{
		d = ( -pd->var[FLOW_ANGLE] * M_PI ) / 180;
		alpha = cos( d );
		beta = sin( d );
		printf( "\nCoordinate transformation of the observation points relative to the source:\n" );
		for( i = 0; i < ( *wd ).nW; i++ )
		{
			x0 = wd->x[i] - pd->var[SOURCE_X];
			y0 = wd->y[i] - pd->var[SOURCE_Y];
			x = x0 * alpha - y0 * beta;
			y = x0 * beta  + y0 * alpha;
			printf( "Well %10s %.15g %.15g : %.15g %.15g\n", wd->id[i], wd->x[i], ( *wd ).y[i], x, y );
		}
	}
	od->obs_target = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_current = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_best = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_id = char_matrix( ( *od ).nObs, 50 );
	od->preds_id = char_matrix( ( *od ).nObs, 50 );
	od->obs_weight = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_min = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_max = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->res = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_log = ( int * ) malloc( ( *od ).nObs * sizeof( int ) );
	od->well_index = ( int * ) malloc( ( *od ).nObs * sizeof( int ) );
	od->time_index = ( int * ) malloc( ( *od ).nObs * sizeof( int ) );
	if( cd->debug ) printf( "\n" );
	printf( "Number of calibration targets = %d\n", ( *od ).nObs );
	for( k = i = 0; i < ( *wd ).nW; i++ )
		for( j = 0; j < ( *wd ).nWellObs[i]; j++ )
			if( ( *wd ).obs_weight[i][j] > 0 )
			{
				od->obs_target[k] = ( *wd ).obs_target[i][j];
				od->obs_weight[k] = ( *wd ).obs_weight[i][j];
				od->obs_min[k] = ( *wd ).obs_min[i][j];
				od->obs_max[k] = ( *wd ).obs_max[i][j];
				od->obs_log[k] = ( *wd ).obs_log[i][j];
				od->well_index[k] = i;
				od->time_index[k] = j;
				sprintf( od->obs_id[k], "%s(%g)", wd->id[i], wd->obs_time[i][j] );
				if( cd->debug ) printf( "%s(%g): %g weight %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_target[i][j], ( *wd ).obs_weight[i][j] );
				k++;
			}
	if( ( *preds ).nObs > 0 )
	{
		printf( "Number of performance criterion predictions for info-gap analysis = %d\n", ( *preds ).nObs );
		preds->obs_target = ( double * ) malloc( ( *preds ).nObs * sizeof( double ) );
		preds->obs_current = ( double * ) malloc( ( *preds ).nObs * sizeof( double ) );
		preds->obs_best = ( double * ) malloc( ( *preds ).nObs * sizeof( double ) );
		preds->well_index = ( int * ) malloc( ( *preds ).nObs * sizeof( int ) );
		preds->time_index = ( int * ) malloc( ( *preds ).nObs * sizeof( int ) );
		preds->obs_id = char_matrix( ( *preds ).nObs, 50 );
		preds->obs_weight = ( double * ) malloc( ( *preds ).nObs * sizeof( double ) );
		preds->obs_min = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		preds->obs_max = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		preds->res = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		for( k = i = 0; i < ( *wd ).nW; i++ )
			for( j = 0; j < ( *wd ).nWellObs[i]; j++ )
				if( ( *wd ).obs_weight[i][j] < 0 )
				{
					preds->obs_target[k] = ( *wd ).obs_target[i][j];
					preds->obs_weight[k] = 1.0;
					preds->obs_min[k] = ( *wd ).obs_target[i][j];
					preds->obs_max[k] = ( *wd ).obs_target[i][j];
					preds->obs_log[k] = ( *wd ).obs_log[i][j];
					preds->well_index[k] = i;
					preds->time_index[k] = j;
					sprintf( preds->obs_id[k], "%s(%g)", wd->id[i], wd->obs_time[i][j] );
					if( cd->debug ) printf( "%s(%g): %g weight %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_target[i][j], ( *wd ).obs_weight[i][j] );
					k++;
				}
	}
	else if( ( ( *preds ).nObs <= 0 ) && ( ( *cd ).problem_type == INFOGAP ) )
	{ printf( "\nWeight of at least one observation must be set as performance criterion prediction\nby setting weight to -1 for infogap analysis\n\n" ); exit( 0 ); }
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf\n", &( *gd ).time );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %i %i %i\n", &( *gd ).nx, &( *gd ).ny, &( *gd ).nz );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &( *gd ).min_x, &( *gd ).min_y, &( *gd ).min_z );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &( *gd ).max_x, &( *gd ).max_y, &( *gd ).max_z );
	fscanf( infile, "%[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &( *gd ).min_t, &( *gd ).max_t, &( *gd ).dt );
	fclose( infile );
	if( cd->debug )
	{
		printf( "\nGrid Time: %g\n", ( *gd ).time );
		printf( "Grid lines: %i %i %i\n", ( *gd ).nx, ( *gd ).ny, ( *gd ).nz );
		printf( "Grid Minimums: %g %g %g\n", ( *gd ).min_x, ( *gd ).min_y, ( *gd ).min_z );
		printf( "Grid Maximums: %g %g %g\n", ( *gd ).max_x, ( *gd ).max_y, ( *gd ).max_z );
	}
	if( ( *gd ).nx == 1 )( *gd ).dx = 0;
	else( *gd ).dx = ( ( *gd ).max_x - ( *gd ).min_x ) / ( ( *gd ).nx - 1 );
	if( ( *gd ).ny == 1 )( *gd ).dy = 0;
	( *gd ).dy = ( ( *gd ).max_y - ( *gd ).min_y ) / ( ( *gd ).ny - 1 );
//	if(( *gd ).nz == 1 )( *gd ).dz = ( *gd ).max_z - ( *gd ).min_z ); // In this way compute_grid computed for min_z
	if( ( *gd ).nz == 1 )( *gd ).dz = 0;
	else( *gd ).dz = ( ( *gd ).max_z - ( *gd ).min_z ) / ( ( *gd ).nz - 1 );
	if( cd->debug ) printf( "Breakthrough-curve time window: %g %g %g\n", ( *gd ).min_t, ( *gd ).max_t, ( *gd ).dt );
	gd->nt = 1 + ( ( *gd ).max_t - ( *gd ).min_t ) / ( *gd ).dt;
	return( 1 );
}

int save_problem( char *filename, struct opt_data *op )
{
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od;
	struct well_data *wd;
	struct grid_data *gd;
	struct extrn_data *ed;
	cd = op->cd;
	pd = op->pd;
	od = op->od;
	wd = op->wd;
	gd = op->gd;
	ed = op->ed;
	FILE *outfile;
//	FILE *outfileb;
//	char buf[100];
	int  i, j;
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		printf( "File \'%s\' cannot be opened to save the problem information!\n", filename );
		return( 0 );
	}
	/*
		sprintf( buf, "%s.bin", filename );
		if ( ( outfileb = fopen( buf, "wb" ) ) == NULL )
			printf( "Binary file %s cannot be opened to save problem information!\n", buf );
		else
		{
			fwrite( (void *) (*pd).var, sizeof((*pd).var[i]), (*pd).nParam, outfileb );
			fclose( outfileb );
		}
	*/
	i = 0;
	fprintf( outfile, "Problem type: " );
	switch( cd->problem_type )
	{
		case CREATE: fprintf( outfile, "create" ); break;
		case FORWARD: fprintf( outfile, "forward" ); break;
		case CALIBRATE: fprintf( outfile, "calibration" ); break;
		case LOCALSENS: fprintf( outfile, "lsens" ); break;
		case GLOBALSENS: fprintf( outfile, "gsens" ); break;
		case EIGEN: fprintf( outfile, "eigen" ); break;
		case MONTECARLO: fprintf( outfile, "montecarlo real=%d", cd->nreal ); break;
		case ABAGUS: fprintf( outfile, " abagus energy=%d", cd->energy ); break;
		case POSTPUA: fprintf( outfile, " postpua" ); break;
	}
	if( cd->debug > 0 ) fprintf( outfile, " debug=%d", cd->debug );
	if( cd->fdebug > 0 ) fprintf( outfile, " fdebug=%d", cd->fdebug );
	if( cd->ldebug > 0 ) fprintf( outfile, " ldebug=%d", cd->ldebug );
	if( cd->pdebug > 0 ) fprintf( outfile, " pdebug=%d", cd->pdebug );
	if( cd->mdebug > 0 ) fprintf( outfile, " mdebug=%d", cd->mdebug );
	if( cd->odebug > 0 ) fprintf( outfile, " odebug=%d", cd->odebug );
	if( cd->insdebug > 0 ) fprintf( outfile, " insdebug=%d", cd->insdebug );
	if( cd->tpldebug > 0 ) fprintf( outfile, " tpldebug=%d", cd->tpldebug );
	if( cd->pardebug > 0 ) fprintf( outfile, " pardebug=%d", cd->pardebug );
	if( cd->test_func >= 0 ) fprintf( outfile, " test=%d", cd->test_func );
	if( cd->test_func_dim > 2 ) fprintf( outfile, " dim=%d", cd->test_func_dim );
	if( cd->phi_cutoff > 0 ) fprintf( outfile, " cutoff=%g", cd->phi_cutoff );
	fprintf( outfile, " sindx=%g", cd->sindx );
	fprintf( outfile, " lindx=%g", cd->lindx );
	fprintf( outfile, " pardx=%g", cd->lindx );
	if( cd->check_success ) fprintf( outfile, " success" );
	fprintf( outfile, " " );
	switch( cd->calib_type )
	{
		case SIMPLE: fprintf( outfile, "single" ); break;
		case PPSD: fprintf( outfile, "ppsd" ); break;
		case IGRND: fprintf( outfile, "igrnd real=%d", cd->nreal ); break;
		case IGPD: fprintf( outfile, "igpd" ); break;
	}
	if( cd->opt_method[0] != 0 ) fprintf( outfile, " opt=%s", cd->opt_method );
	if( cd->nretries > 0 ) fprintf( outfile, " retry=%d", cd->nretries );
	if( cd->init_particles > 1 ) fprintf( outfile, " particles=%d", cd->init_particles );
	else if( cd->init_particles < 0 ) fprintf( outfile, " particles" );
	if( cd->niter > 0 ) fprintf( outfile, " iter=%d", cd->niter );
	fprintf( outfile, " eval=%d", cd->maxeval );
	if( cd->smp_method[0] != 0 ) fprintf( outfile, " rnd=%s", cd->smp_method );
	if( cd->paran_method[0] != 0 ) fprintf( outfile, " paran=%s", cd->paran_method );
	fprintf( outfile, " " );
	switch( cd->objfunc_type )
	{
		case SSR: fprintf( outfile, "ssr" ); break;
		case SSDR: fprintf( outfile, "ssdr" ); break;
		case SSD0: fprintf( outfile, "ssd0" ); break;
		case SSDA: fprintf( outfile, "ssda" ); break;
	}
	fprintf( outfile, "\n" );
	fprintf( outfile, "Solution type: %s\n", ( *cd ).solution_id );
	fprintf( outfile, "Number of parameters: %i\n", ( *pd ).nParam );
	for( i = 0; i < ( *pd ).nParam; i++ )
	{
		if( ( *pd ).var_opt[i] >= 1 && ( *pd ).var_log[i] == 1 )
			fprintf( outfile, "%s: %.15g %d %d %g %g %g\n", pd->var_id[i], pow( 10, ( *pd ).var[i] ), ( *pd ).var_opt[i], ( *pd ).var_log[i], pow( 10, ( *pd ).var_dx[i] ), pow( 10, pd->var_min[i] ), pow( 10, pd->var_max[i] ) );
		else
			fprintf( outfile, "%s: %.15g %d %d %g %g %g\n", pd->var_id[i], ( *pd ).var[i], ( *pd ).var_opt[i], ( *pd ).var_log[i], ( *pd ).var_dx[i], ( *pd ).var_min[i], ( *pd ).var_max[i] );
	}
	if( cd->solution_type != EXTERNAL )
	{
		fprintf( outfile, "Number of wells: %i\n", ( *wd ).nW );
		for( i = 0; i < ( *wd ).nW; i++ )
		{
			fprintf( outfile, "%s %.15g %.15g %g %g %i ", wd->id[i], ( *wd ).x[i], ( *wd ).y[i], ( *wd ).z1[i], ( *wd ).z2[i], ( *wd ).nWellObs[i] );
			if( ( *wd ).nWellObs[i] > 1 ) { fprintf( outfile, "\n" ); }
			for( j = 0; j < ( *wd ).nWellObs[i]; j++ )
				fprintf( outfile, "%g %g %g %i %g %g\n", ( *wd ).obs_time[i][j], ( *wd ).obs_target[i][j], ( *wd ).obs_weight[i][j], ( *wd ).obs_log[i][j], ( *wd ).obs_min[i][j], ( *wd ).obs_max[i][j] );
		}
		fprintf( outfile, "Grid Time: %g\n", ( *gd ).time );
		fprintf( outfile, "Grid lines: %i %i %i\n", ( *gd ).nx, ( *gd ).ny, ( *gd ).nz );
		fprintf( outfile, "Grid Minimums: %g %g %g\n", ( *gd ).min_x, ( *gd ).min_y, ( *gd ).min_z );
		fprintf( outfile, "Grid Maximums: %g %g %g\n", ( *gd ).max_x, ( *gd ).max_y, ( *gd ).max_z );
		fprintf( outfile, "Time Window: %g %g %g\n", ( *gd ).min_t, ( *gd ).max_t, gd->dt );
	}
	else
	{
		fprintf( outfile, "Number of observations: %i\n", od->nObs );
		for( i = 0; i < od->nObs; i++ )
			fprintf( outfile, "%s %g %g %d %g %g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
		fprintf( outfile, "Execution command: %s\n", ed->cmdline );
		fprintf( outfile, "Number of execution templates: %d\n", ed->ntpl );
		for( i = 0; i < ed->ntpl; i++ )
			fprintf( outfile, "%s %s\n", ed->fn_tpl[i], ed->fn_out[i] );
		fprintf( outfile, "Number of execution instructions: %d\n", ed->ntpl );
		for( i = 0; i < ed->nins; i++ )
			fprintf( outfile, "%s %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	}
	fclose( outfile );
	return( 1 );
}

void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd )
{
	FILE *outfile;
	double x, y, z, t, c, dx, dy, max_conc, max_x, max_y, max_z;
	int  i, j, k, nx, ny;
	if( gd->time <= 0 ) return;
	t = gd->time;
	nx = ny = 15;
	dx = ( gd->max_x - gd->min_x ) / nx;
	dy = ( gd->max_y - gd->min_y ) / ny;
	z = gd->min_z;
	for( j = 0; j < ny; j++ )
	{
		y = gd->max_y - j * dy;
		for( i = 0; i < nx; i++ )
		{
			x = i * dx + gd->min_x;
			c = func_solver1( x, y, z, t, ( void * ) cd );
			printf( "%5.0g ", c );
		}
		printf( "\n" );
	}
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		printf( "Output file %s cannot be opened!\n", filename );
		return;
	}
	fprintf( outfile, "# vtk DataFile Version 2.0\n" );
	fprintf( outfile, "Anasol output\n" );
	fprintf( outfile, "ASCII\n" );
	fprintf( outfile, "DATASET STRUCTURED_POINTS\n" );
	fprintf( outfile, "DIMENSIONS %d %d %d\n", gd->nx, gd->ny, gd->nz );
	fprintf( outfile, "ORIGIN %g %g %g\n", gd->min_x, gd->min_y, gd->min_z ); /* negative z */
//	fprintf( outfile, "ORIGIN %g %g %g\n", gd->min_x, gd->min_y, -gd->max_z ); /* negative z */
	fprintf( outfile, "SPACING %g %g %g\n", gd->dx, gd->dy, gd->dz );
	fprintf( outfile, "POINT_DATA %d\n", gd->nx * gd->ny * gd->nz );
	fprintf( outfile, "SCALARS Cr_concentration_ppb float 1\n" );
	fprintf( outfile, "LOOKUP_TABLE default\n" );
	max_conc = max_x = max_y = max_z = 0;
	for( k = 0; k < gd->nz; k++ )
	{
//		z = gd->max_z - k * gd->dz; /* reverse z */
		z = k * gd->dz + gd->min_z;
		for( j = 0; j < gd->ny; j++ )
		{
			y = j * gd->dy + gd->min_y;
			for( i = 0; i < gd->nx; i++ )
			{
				x = i * gd->dx + gd->min_x;
				c = func_solver1( x, y, z, t, ( void * ) cd );
				fprintf( outfile, "%g ", c );
				if( max_conc < c ) { max_conc = c; max_x = x; max_y = y; max_z = z; }
			}
			fprintf( outfile, "\n" );
		}
		fprintf( outfile, "\n" );
	}
	printf( "Max Conc = %g @ (%g,%g,%g)\n", max_conc, max_x, max_y, max_z );
	fclose( outfile );
	printf( "Spatial concentration data saved in %s.\n", filename );
}

void compute_btc2( char *filename, char *filename2, struct opt_data *op )
{
	FILE *outfile;
	double time, time_expected, c, max_source_conc, max_source_time, *max_conc, *max_time, d, x0, y0, alpha, beta, xe, ye, v, v_apparent;
	int  i, k;
	struct grid_data *gd;
	gd = op->gd;
	if( gd->min_t <= 0 ) return;
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		printf( "Output file %s cannot be opened!\n", filename );
		return;
	}
	printf( "\n" );
	max_time = ( double * ) malloc( op->wd->nW * sizeof( double ) );
	max_conc = ( double * ) malloc( op->wd->nW * sizeof( double ) );
	fprintf( outfile, "variable = \"Time [a]\"" );
	for( i = 0; i < op->wd->nW; i++ )
	{
		fprintf( outfile, " \"%s\"", op->wd->id[i] );
		max_conc[i] = 0;
		max_time[i] = -1;
	}
	fprintf( outfile, "\n" );
	for( k = 0; k < gd->nt; k++ )
	{
		time = gd->min_t + gd->dt * k;
		fprintf( outfile, "%g", time );
		for( i = 0; i < op->wd->nW; i++ )
		{
			c = func_solver( op->wd->x[i], op->wd->y[i], op->wd->z1[i], op->wd->z2[i], time, ( void * ) op->cd );
			fprintf( outfile, " %g", c );
			if( max_conc[i] < c ) { max_conc[i] = c; max_time[i] = time; }
		}
		fprintf( outfile, "\n" );
	}
	fclose( outfile );
	printf( "Concentration breakthrough data saved in %s\n", filename );
	if( ( outfile = fopen( filename2, "w" ) ) == NULL )
	{
		printf( "Output file %s cannot be opened!\n", filename );
		return;
	}
	printf( "\nCompute source max concentration ..." );
	max_source_conc = 0; max_source_time = op->pd->var[TIME_INIT];
	for( k = 0; k < gd->nt; k++ )
	{
		time = gd->min_t + gd->dt * k;
		c = func_solver( op->pd->var[SOURCE_X], op->pd->var[SOURCE_Y], op->pd->var[SOURCE_Z], op->pd->var[SOURCE_Z], time, ( void * ) op->cd );
		if( max_source_conc < c ) { max_source_conc = c; max_source_time = time; }
	}
	printf( "\nPeak source concentration = %g @ %g\n", max_source_conc, max_source_time );
	for( i = 0; i < op->pd->nParam; i++ ) // IMPORTANT: Take care of log transformed variable
		if( op->pd->var_opt[i] >= 1 && op->pd->var_log[i] == 1 )
			op->pd->var[i] = pow( 10, op->pd->var[i] );
	v = sqrt( op->pd->var[VX] * op->pd->var[VX] + op->pd->var[VY] * op->pd->var[VY] + op->pd->var[VZ] * op->pd->var[VZ] ); // Flow velocity
	printf( "Flow velocity = %g (%g %g %g)\n", v, op->pd->var[VX], op->pd->var[VY], op->pd->var[VZ] );
	for( i = 0; i < op->wd->nW; i++ )
	{
		x0 = ( op->wd->x[i] - op->pd->var[SOURCE_X] );
		y0 = ( op->wd->y[i] - op->pd->var[SOURCE_Y] );
		d = ( -op->pd->var[FLOW_ANGLE] * M_PI ) / 180;
		alpha = cos( d );
		beta = sin( d );
		xe = x0 * alpha - y0 * beta;
		ye = x0 * beta  + y0 * alpha;
		d = sqrt( xe * xe + ye * ye );
		time = max_time[i] - max_source_time;
		if( time > DBL_EPSILON ) v_apparent = d / time; else { v_apparent = -1; if( time < 0 ) time = -1; };
		c = max_conc[i] / max_source_conc; // Normalized concentration
		if( v > DBL_EPSILON ) time_expected = d / v; else time_expected = -1;
		printf( "%s\tPeak Concentration = %12g (%12g) @ time %12g (%12g) velocity = %12g (%g) distance = %g\n", op->wd->id[i], c, max_conc[i], time, time_expected, v_apparent, v, d );
		fprintf( outfile, "%s\tPeak Conc = %12g @ time %12g (%g) velocity = %12g (%g) distance = %12g\n", op->wd->id[i], c, time, time_expected, v_apparent, v, d );
	}
	fclose( outfile );
	printf( "Concentration peak data saved in %s\n", filename2 );
	for( i = 0; i < op->pd->nParam; i++ ) // IMPORTANT: Take care of log transformed variable
		if( op->pd->var_opt[i] >= 1 && op->pd->var_log[i] == 1 )
			op->pd->var[i] = log10( op->pd->var[i] );
	free( max_conc );
	free( max_time );
}

void compute_btc( char *filename, struct opt_data *op )
{
	FILE *outfile;
	double time, c, max_conc, max_time;
	int  i, k;
	struct grid_data *gd;
	gd = op->gd;
	if( gd->min_t <= 0 ) return;
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		printf( " Output file %s cannot be opened!\n", filename );
		return;
	}
	printf( "\n" );
	for( i = 0; i < op->wd->nW; i++ )
	{
		max_conc = max_time = 0;
		for( k = 0; k < gd->nt; k++ )
		{
			time = gd->min_t + gd->dt * k;
			c = func_solver( op->wd->x[i], op->wd->y[i], op->wd->z1[i], op->wd->z2[i], time, ( void * ) op->cd );
			fprintf( outfile, "%s %g %g\n", op->wd->id[i], time, c );
			if( max_conc < c ) { max_conc = c; max_time = time; }
		}
		printf( "%s\tPeak Conc = %12.4g @ time %12g\tConc = %12.4f @ time %12g\n", op->wd->id[i], max_conc, max_time, c, time );
	}
	fclose( outfile );
	printf( "Concentration breakthrough data saved in %s\n", filename );
}

static char *strsave( const char *s, const char *lim )
{
	if( lim == NULL )
		lim = s + strlen( s );
	char *p = ( char * ) malloc_check( "save string", lim - s + 1 );
	strncpy( p, s, lim - s );
	p[lim - s] = '\0';
	return p;
}

char **shellpath( void )
{
	const char *path = getenv( "PATH" );
	if( !path )
		path = "/bin:/usr/bin:/usr/local/bin";
	char **vector = // size is overkill
		( char ** ) malloc_check( "hold path elements", strlen( path ) * sizeof( *vector ) );
	const char *p = path;
	int next = 0;
	while( p )
	{
		char *q = strchr( p, ':' );
		vector[next++] = strsave( p, q );
		p = q ? q + 1 : NULL;
	}
	vector[next] = NULL;
	return vector;
}

void freeshellpath( char *shellpath[] )
{
	int i;
	for( i = 0; shellpath[i]; i++ )
		free( shellpath[i] );
	free( shellpath );
}

unsigned maxpathlen( char *path[], const char *base )
{
	unsigned blen = strlen( base );
	unsigned n = 0;
	int i;
	for( i = 0; path[i]; i++ )
	{
		unsigned pn = strlen( path[i] );
		if( pn > n ) n = pn;
	}
	return blen + n + 1;
}

void execvepath( char *path[], const char *base, char *const argv[], char *const envp[] )
{
	int i;
	if( strchr( base, '/' ) )
		execve( base, argv, envp );
	else
	{
		size_t maxlen = maxpathlen( path, base ) + 1;
		char *buf = ( char * ) malloc_check( "hold path", maxlen );
		for( i = 0; path[i]; i++ )
		{
			snprintf( buf, maxlen, "%s/%s", path[i], base );
			execve( buf, argv, envp );
		}
	}
}

int count_lines( char *filename )
{
	int nol = 0;
	FILE *fl;
	char buf[1000];
	fl = fopen( filename, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", filename ); exit( 0 ); }
	while( ( fgets( buf, sizeof buf, fl ) ) != NULL ) nol++;
	fclose( fl );
	return nol;
}

int count_cols( char *filename, int row )
{
	int ncol = 0, i, n = 0;
	FILE *fl;
	char buf[1000], entry[16], *ln;
	fl = fopen( filename, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", filename ); exit( 0 ); }
	for( i = 1; i < row; i++ ) ln = fgets( buf, sizeof buf, fl );
	while( sscanf( ln, "%10s%n", entry, &n ) == 1 )
	{
		ncol++;
		ln += n;
	}
	fclose( fl );
	return ncol;
}
