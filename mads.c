// MADS: Model Analyses & Decision Support (v1.1) 2012
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
// Brianeisha Eure
// Leif Zinn-Bjorkman
//
// http://mads.lanl.gov
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_interp.h>
#include "misc/levmar-2.5/levmar.h"
#include "mads.h"
#ifdef MATHEVAL
#include <matheval.h>
#endif

#define FIT(i) gsl_vector_get(solver->x, i)
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
// Model analyses
int optimize_lm( struct opt_data *op ); // LM (Levenberg-Marquardt) optimization
int optimize_pso( struct opt_data *op ); // PSO optimization
int eigen( struct opt_data *op, double *f_x, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar ); // Eigen analysis
int check( struct opt_data *op );
int igrnd( struct opt_data *op );
int igpd( struct opt_data *op );
int ppsd( struct opt_data *op );
int montecarlo( struct opt_data *op );
//
void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op, int debug ); // Random sampling
void print_results( struct opt_data *op, int verbosity ); // Print final results
void save_final_results( char *filename, struct opt_data *op, struct grid_data *gd ); // Save final results
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var );
void ave_sorted( double data[], int n, double *ave, double *ep );
int sort_int( const void *x, const void *y );
double sort_double( const void *a, const void *b );

/* Functions elsewhere */
// Optimization strategies
int pso_tribes( struct opt_data *op );
int pso_std( struct opt_data *op );
int mopso( struct opt_data *op );
int lm_opt( int func( double x[], void *data, double f[] ), int func_dx( double *x, double *f_x, void *data, double *jacobian ), void *data,
			int nObs, int nParam, int nsig, double eps, double delta, int max_eval, int max_iter,
			int iopt, double parm[], double x[], double *phi, double f[],
			double jacobian[], int nian, double jacTjac[], int *infer );
int zxssqch( int func( double x[], void *, double f[] ), void *func_data,
			 int m, int n, int nsig, double eps, double delta, int maxfn,
			 int iopt, double parm[], double x[], double *phi, double f[],
			 double xjac[], int ixjac, double xjtj[], int *infer );
int lm_gsl( gsl_vector *x, struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *covar );
// Info
void mads_info();
// IO
char *timestamp(); // create time stamp
char *datestamp(); // create date stamp
int parse_cmd( char *buf, struct calc_data *cd );
int load_problem( char *filename, int argn, char *argv[], struct opt_data *op );
int load_pst( char *filename, struct opt_data *op );
int save_problem( char *filename, struct opt_data *op );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc( char *filename, struct opt_data *op );
void compute_btc2( char *filename, char *filename2, struct opt_data *op );
int Ftest( char *filename );
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );
FILE *Fappend( char *filename );
char *Fdatetime( char *filename, int debug );
time_t Fdatetime_t( char *filename, int debug );
// External IO
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
// Random sampling
double epsilon();
void lhs_imp_dist( int nvar, int npoint, int d, int *seed, double x[] );
void lhs_center( int nvar, int npoint, int *seed, double x[] );
void lhs_edge( int nvar, int npoint, int *seed, double x[] );
void lhs_random( int nvar, int npoint, int *seed, double x[] );
void smp_random( int nvar, int npoint, int *seed, double x[] );
int get_seed( );
// YAML
#ifdef YAML
int load_yaml_problem( char *filename, int argn, char *argv[], struct opt_data *op );
#endif
// Memory
double **double_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
void zero_double_matrix( double **matrix, int maxCols, int maxRows );
char *white_trim( char *x );
char **char_matrix( int maxCols, int maxRows );
// Func
int func_gsl_dx( const gsl_vector *x, void *data, gsl_matrix *J );
int func_gsl_xdx( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J );
double func_gsl_deriv( double x, void *data );
int func_gsl_deriv_dx( const gsl_vector *x, void *data, gsl_matrix *J );
int func_extrn( double *x, void *data, double *f );
int func_intrn( double *x, void *data, double *f );
void func_levmar( double *x, double *f, int m, int n, void *data );
void func_dx_levmar( double *x, double *f, double *jacobian, int m, int n, void *data );
int func_dx( double *x, double *f_x, void *data, double *jacobian );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
int func_extrn_write( int ieval, double *x, void *data );
int func_extrn_exec_serial( int ieval, void *data );
int func_extrn_read( int ieval, void *data, double *f );
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );
// Parallel
int mprun( int nJob, void *data );
char *dir_hosts( void *data, char *timedate_stamp );
// Methods
int sa_sobol( struct opt_data *op );
int sa_saltelli( struct opt_data *op );
int sa_moat( struct opt_data *op );
int abagus( struct opt_data *op );
int infogap_obs( struct opt_data *op );
int infogap( struct opt_data *op );
int postpua( struct opt_data *op );
int glue( struct opt_data *op );

int main( int argn, char *argv[] )
{
	// TODO return status of the function calls is not always checked; needs to be checked
	int i, j, k, ier, status, success, success_all, count,  predict = 0, compare, bad_data = 0;
	double c, err, min, max, dx, phi, *opt_params;
	struct calc_data cd;
	struct param_data pd;
	struct regul_data rd;
	struct obs_data od;
	struct obs_data preds;
	struct well_data wd;
	struct extrn_data ed;
	struct grid_data gd;
	struct opt_data op;
	struct anal_data ad;
	struct source_data sd;
	struct aquifer_data qd;
	char filename[255], filename2[255], root[255], extension[255], buf[255], *dot, *cwd;
	int ( *optimize_func )( struct opt_data * op ); // function pointer to optimization function (LM or PSO)
	char *host, *nodelist, *hostlist, *proclist, *lsblist, *beowlist; // parallel variables
	FILE *in, *out, *out2;
	time_t time_start, time_end, time_elapsed;
	pid_t pid;
	struct tm *ptr_ts;
	time_start = time( NULL );
	op.datetime_stamp = datestamp(); // create execution date stamp
	op.pd = &pd; // create opt_data structures ...
	op.rd = &rd;
	op.od = &od;
	op.preds = &preds;
	op.wd = &wd;
	op.cd = &cd;
	op.gd = &gd;
	op.ed = &ed;
	op.ad = &ad;
	op.sd = &sd;
	op.qd = &qd;
	cd.neval = 0;
	cd.njac = 0;
	cd.nlmo = 0;
	cd.lmstandalone = 1; // LM variable; LM is stand-alone if not part of tribes optimization
	cd.compute_phi = 1; // function calls compute OF (phi); turned off only when the jacobians are computed for LM
	cd.pderiv = cd.oderiv = -1; // internal flags; do not compute parameter and observation derivatives
	cd.c_background = 0;
	op.phi = HUGE_VAL;
	op.success = op.global_success = 0;
	op.f_ofe = NULL;
	printf( "MADS: Model Analyses & Decision Support (v1.1) 2012\n" );
	printf( "---------------------------------------------------\n" );
	if( argn < 2 )
	{
		mads_info(); // print short mads help manual
		exit( 1 );
	}
	else
		printf( "Velimir Vesselinov (monty) vvv@lanl.gov -:- velimir.vesselinov@gmail.com\nhttp://mads.lanl.gov -:- http://www.ees.lanl.gov/staff/monty/codes/mads\n\n" );
	op.label = ( char * ) malloc( 10 * sizeof( char ) ); op.label[0] = 0;
	if( cd.debug ) printf( "Argument[1]: %s\n", argv[1] );
	else
	{
		if( cd.debug > 1 && argn > 2 )
			for( i = 2; i < argn; i++ )
				printf( "Argument[%d]: %s\n", i, argv[i] );
	}
	strcpy( root, argv[1] ); // Defined problem name (root)
	dot = strrchr( root, '.' );
	if( dot != NULL && dot[1] != '/' )
	{
		strcpy( filename, argv[1] );
		strcpy( extension, &dot[1] );
		dot[0] = 0;
	}
	else
	{
		sprintf( filename, "%s.mads", argv[1] );
		extension[0] = 0;
	}
	printf( "Input file name: %s ", filename );
	if( strcasestr( filename, "yaml" ) || strcasestr( filename, "yml" ) )
	{
		printf( "(YAML format expected)\n" );
		op.yaml = 1; // YAML format
	}
	else if( strcasecmp( extension, "pst" ) == 0 )
	{
		printf( "(PEST format expected)\n" );
		op.yaml = 0; // PEST format
	}
	else
	{
		printf( "(Plain text MADS format expected)\n" );
		op.yaml = 0; // MADS plain text format
	}
	cd.time_infile = Fdatetime_t( filename, 0 );
	cd.datetime_infile = Fdatetime( filename, 0 );
	printf( "Problem root name: %s", root );
	if( extension[0] != 0 )	printf( " Extension: %s\n", extension );
	else printf( "\n" );
	op.root = root;
	op.filename = filename;
	op.counter = 0;
	sprintf( filename2, "%s.mads_output", op.root );
	if( Ftest( filename2 ) == 0 ) // If file already exists quit ...
	{
		sprintf( buf, "mv %s.mads_output %s.mads_output_%s >& /dev/null", op.root, op.root, Fdatetime( filename2, 0 ) );  // Move existing output file
		system( buf );
	}
	mads_output = Fwrite( filename2 );
	fprintf( mads_output, "MADS: Model Analyses & Decision Support (v1.1) 2012\n" );
	fprintf( mads_output, "---------------------------------------------------\n" );
	fprintf( mads_output, "Velimir Vesselinov (monty) vvv@lanl.gov -:- velimir.vesselinov@gmail.com\nhttp://mads.lanl.gov -:- http://www.ees.lanl.gov/staff/monty/codes/mads\n\n" );
	fprintf( mads_output, "Input file name: %s\n", filename );
	fprintf( mads_output, "Problem root name: %s\n", root );
	sprintf( buf, "%s.running", op.root ); // File named root.running is used to prevent simultaneous execution of multiple problems
	if( Ftest( buf ) == 0 ) // If file already exists quit ...
	{
		tprintf( "WARNING: Potentially another MADS run is currently performed for problem \'%s\' since file %s exists!\n\n", op.root, buf );
		// tprintf( "ERROR: Potentially another MADS run is currently performed for problem \'%s\' since file %s exists!\n", op.root, buf );
		// tprintf( "Delete %s to execute (sorry for the inconvenience)!\n", buf );
		// exit( 0 );
	}
	sprintf( buf, "touch %s.running", op.root ); system( buf ); // Create a file named root.running to prevent simultaneous execution of multiple problems
	/*
	 *  Read input data
	 */
	if( strcasecmp( extension, "pst" ) == 0 ) // PEST Problem
	{
		tprintf( "PEST problem:\n" );
		if( ( ier = load_pst( filename, &op ) ) <= 0 )
		{
			tprintf( "\nMADS quits! Data input problem!\nExecute \'mads\' without any arguments to check the acceptable command-line keywords and options.\n" );
			if( ier == 0 )
			{
				sprintf( filename, "%s-error.mads", op.root );
				save_problem( filename, &op );
				tprintf( "MADS problem file named %s-error.mads is created to debug.\n", op.root );
			}
			sprintf( buf, "rm -f %s.running", op.root ); system( buf ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			exit( 0 );
		}
		if( cd.opt_method[0] == 0 ) { strcpy( cd.opt_method, "lm" ); cd.calib_type = SIMPLE; cd.problem_type = CALIBRATE; }
		cd.solution_type[0] = EXTERNAL; func_global = func_extrn;
		buf[0] = 0;
		for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
		if( parse_cmd( buf, &cd ) == -1 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	else // MADS Problem
	{
		if( op.yaml ) // YAML format
		{
			// TODO parse command line after parsing the Problem and Solution YAML classes; the lines below should be removed
			buf[0] = 0;
			for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
			cd.solution_type = ( int * ) malloc( sizeof( int ) );
			if( parse_cmd( buf, &cd ) == -1 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
#ifdef YAML
			ier = load_yaml_problem( filename, argn, argv, &op );
#else
			tprintf( "\nMADS quits! YAML format is not supported in the compiled version of MADS. Recompile with YAML libraries.\n" );
			exit( 0 );
#endif
		}
		else			ier = load_problem( filename, argn, argv, &op ); // MADS plain text format
		if( ier <= 0 )
		{
			tprintf( "\nMADS quits! Data input problem!\nExecute \'mads\' without any arguments to check the acceptable command-line keywords and options.\n" );
			if( ier == 0 )
			{
				sprintf( filename, "%s-error.mads", op.root );
				save_problem( filename, &op );
				tprintf( "MADS problem file named %s-error.mads is created to debug.\n", op.root );
			}
			sprintf( buf, "rm -f %s.running", op.root ); system( buf ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			exit( 0 );
		}
		if( cd.solution_type[0] == EXTERNAL ) func_global = func_extrn;
		else func_global = func_intrn;
	}
	/*
	 *  Check for parallel environment
	 */
	cd.paral_hosts = NULL;
	hostlist = NULL;
	if( ( nodelist = getenv( "NODELIST" ) ) != NULL )
	{
		if( cd.debug ) tprintf( "\nParallel environment is detected (environmental variable NODELIST is defined)\n" );
		if( cd.debug ) tprintf( "Node list %s\n", nodelist );
		hostlist = nodelist;
	}
	if( ( beowlist = getenv( "BEOWULF_JOB_MAP" ) ) != NULL )
	{
		if( cd.debug ) tprintf( "\nParallel environment is detected (environmental variable BEOWULF_JOB_MAP is defined)\n" );
		if( cd.debug ) tprintf( "Node list %s\n", beowlist );
		hostlist = beowlist;
	}
	if( ( lsblist = getenv( "LSB_HOSTS" ) ) != NULL )
	{
		if( cd.debug ) tprintf( "\nParallel environment is detected (environmental variable LSB_HOSTS is defined)\n" );
		if( cd.debug ) tprintf( "Node list %s\n", lsblist );
		hostlist = lsblist;
		if( ( proclist = getenv( "LSB_MCPU_HOSTS" ) ) != NULL && cd.debug ) tprintf( "LSB_MCPU_HOSTS Processors list %s\n", proclist );
	}
	if( hostlist != NULL )
	{
		if( cd.debug == 0 ) tprintf( "\nParallel environment is detected.\n" );
		if( ( host = getenv( "HOSTNAME" ) ) == NULL ) host = getenv( "HOST" );
		tprintf( "Host: %s\n", host );
		k = strlen( hostlist );
		i = count = 0;
		tprintf( "Nodes:" );
		while( i <= k )
		{
			sscanf( &hostlist[i], "%s", buf );
			tprintf( " \'%s\'", buf );
			i += strlen( buf ) + 1;
			count++;
		}
		tprintf( "\n" );
		tprintf( "Number of available nodes for parallel execution: %i\n", count );
		if( count <= 0 )
			tprintf( "ERROR: There is problem with the description of execution nodes!\n" );
		else
		{
			cd.num_proc = count;
			cd.paral_hosts = char_matrix( cd.num_proc, 95 );
			k = strlen( hostlist );
			i = 0;
			j = 0;
			while( i <= k )
			{
				sscanf( &hostlist[i], "%s", cd.paral_hosts[j] );
				i += strlen( cd.paral_hosts[j++] ) + 1;
			}
		}
	}
	else if( cd.num_proc > 1 )
	{
		tprintf( "\nLocal parallel execution is requested using %d processors (np=%d)\n", cd.num_proc, cd.num_proc );
		cwd = getenv( "OSTYPE" ); tprintf( "OS type: %s\n", cwd );
		if( strncasecmp( cwd, "darwin", 6 ) == 0 )
			system( "\\rm -f num_proc >& /dev/null; ( sysctl hw.logicalcpu | cut -d : -f 2 ) > num_proc" ); // MAC OS
		else
			system( "\\rm -f num_proc >& /dev/null; ( cat /proc/cpuinfo | grep processor | wc -l ) > num_proc" ); // LINUX
		in = Fread( "num_proc" );
		fscanf( in, "%d", &k );
		fclose( in );
		system( "\\rm -f num_proc >& /dev/null" );
		tprintf( "Number of local processors available for parallel execution: %i\n", k );
		if( k < cd.num_proc ) tprintf( "WARNING: Number of requested processors exceeds the available resources!\n" );
	}
	if( cd.num_proc > 1 ) // Parallel job
	{
		pid = getpid();
		if( cd.debug ) tprintf( "Parent ID [%d]\n", pid );
		cwd = getenv( "PWD" ); dot = strrchr( cwd, '/' );
		cd.mydir = &dot[1];
		if( cd.debug ) tprintf( "Working directory: %s (%s)\n", cwd, cd.mydir );
		cd.mydir_hosts = dir_hosts( &op, op.datetime_stamp ); // Directories for parallel execution have unique name based on the execution time
	}
	tprintf( "\n" );
	/*
	 *  Problem based on external model
	 */
	if( cd.solution_type[0] == EXTERNAL ) // Check the files for external execution
	{
		tprintf( "Checking the template files for errors ...\n" );
		bad_data = 0;
		for( i = 0; i < pd.nParam; i++ ) cd.var[i] = ( double ) - 1;
		for( i = 0; i < ed.ntpl; i++ ) // Check template files ...
			if( check_par_tpl( pd.nParam, pd.var_id, cd.var, ed.fn_tpl[i], cd.tpldebug ) == -1 )
				bad_data = 1;
		for( i = 0; i < pd.nParam; i++ )
		{
			if( cd.var[i] < 0 )
			{
				tprintf( "ERROR: Model parameter \'%s\' is not represented in the template file(s)!\n", pd.var_id[i] );
				bad_data = 1;
			}
			else if( cd.var[i] > 1.5 )
				tprintf( "WARNING: Model parameter \'%s\' is represented more than once (%d times) in the template file(s)!\n", pd.var_id[i], ( int ) cd.var[i] );
		}
		if( !bad_data ) tprintf( "Template files are ok.\n\n" );
		tprintf( "Checking the instruction files for errors ...\n" );
		for( i = 0; i < od.nObs; i++ ) od.obs_current[i] = ( double ) - 1;
		for( i = 0; i < ed.nins; i++ )
			if( check_ins_obs( od.nObs, od.obs_id, od.obs_current, ed.fn_ins[i], cd.insdebug ) == -1 ) // Check instruction files.
				bad_data = 1;
		for( i = 0; i < od.nObs; i++ )
		{
			if( od.obs_current[i] < 0 )
			{
				tprintf( "ERROR: Observation \'%s\' is not defined in the instruction file(s)!\n", od.obs_id[i] );
				bad_data = 1;
			}
			else if( od.obs_current[i] > 1.5 )
				tprintf( "WARNING: Observation \'%s\' is defined more than once (%d times) in the instruction file(s)! Arithmetic average will be computed!\n", od.obs_id[i], ( int ) od.obs_current[i] );
		}
		if( !bad_data ) tprintf( "Instruction files are ok.\n" );
		if( bad_data )
		{
			sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			system( buf );
			exit( 0 );
		}
		if( cd.sintrans ) { if( cd.sindx < DBL_EPSILON ) cd.sindx = 0.1; else if( cd.sindx < 1e-3 ) tprintf( "WARNING: sindx (%g) is potentially too small for external problems; consider increasing sindx (add sindx=1e-2)\n", cd.sindx ); }
		else { if( cd.lindx < DBL_EPSILON ) cd.lindx = 0.01; }
	}
	else
	{
		if( cd.sintrans ) { if( cd.sindx < DBL_EPSILON ) cd.sindx = 1e-7; else if( cd.sindx > 1e-5 ) tprintf( "WARNING: sindx (%g) is potentially too large for internal problems; consider decreasing sindx (add sindx=1e-6)\n", cd.sindx ); }
		else { if( cd.lindx < DBL_EPSILON ) cd.lindx = 0.001; }
	}
	/*
	 *  Check for restart conditions
	 */
	tprintf( "\nExecution date & time stamp: %s\n", op.datetime_stamp ); // Stamp will be applied to name / rename various output files
	if( cd.solution_type[0] == EXTERNAL && cd.num_proc > 1 )
	{
		if( cd.restart == 1 ) // Restart by default
		{
			strcpy( buf, filename ); // Temporarily preserve the input file name
			sprintf( filename, "%s.restart_%s.zip", op.root, cd.datetime_infile );
			strcpy( cd.restart_zip_file, filename );
			if( Ftest( cd.restart_zip_file ) != 0 ) { if( cd.pardebug ) tprintf( "ZIP file (%s) with restart information is not available.\n", cd.restart_zip_file ); cd.restart = 0; }
			else
			{
				time_elapsed = cd.time_infile - Fdatetime_t( cd.restart_zip_file, 0 ); // time_infile - time_zipfile ...
				if( time_elapsed >= 0 ) { if( cd.pardebug ) tprintf( "No restart: the zip file (%s) with restart information is older than the MADS input file (%s)\n(restart can be enforced using \'restart=-1\' or \'rstfile=%s\')\n", cd.restart_zip_file, buf, cd.restart_zip_file ); cd.restart = 0; } // No restart
				else cd.restart = 1; // Attempt restart
			}
			if( cd.restart )
				tprintf( "DEFAULT Restart: zip file %s is consistent with date/time stamp of the MADS input file\n(IMPORTANT: to avoid restart either delete zip file %s, or use keyword \'restart=0\')\n", filename, filename );
		}
		else if( cd.restart == -1 ) // Forced restart
		{
			sprintf( filename, "%s.restart_%s.zip", op.root, cd.datetime_infile );
			if( cd.restart_zip_file[0] == 0 ) strcpy( cd.restart_zip_file, filename );
			if( Ftest( cd.restart_zip_file ) != 0 ) { tprintf( "Restart is requested but a zip file (%s) with restart information is not available.\n", cd.restart_zip_file ); cd.restart = 0; }
			else tprintf( "FORCED Restart: using zip file %s ...\n", cd.restart_zip_file );
		}
		if( cd.restart )
		{
			tprintf( "MADS  input  file \'%40s\' last modified on %s\n", buf, Fdatetime( buf, 0 ) );
			sprintf( filename, "%s.results", op.root ); if( Ftest( filename ) != 0 ) tprintf( "MADS results file \'%40s\' last modified on %s\n", filename, Fdatetime( filename, 0 ) );
			tprintf( "MADS restart file \'%40s\' last modified on %s\n", cd.restart_zip_file, Fdatetime( cd.restart_zip_file, 0 ) );
			tprintf( "ZIP file (%s) with restart information is unzipped ... \n", cd.restart_zip_file );
			sprintf( buf, "rm -fR ../%s* %s.restart_info; unzip -o -u -: %s >& /dev/null", cd.mydir_hosts, op.root, cd.restart_zip_file ); // the input file name was temporarily in buf; not any more ...
			system( buf );
			sprintf( filename, "%s.restart_info", op.root );
			in = Fread( filename );
			fgets( buf, 255, in );
			white_trim( buf );
			cd.mydir_hosts = dir_hosts( &op, buf ); // Directories for parallel execution have unique name based on the old execution time (when restart files were created)
			fclose( in );
			tprintf( "Date & time stamp of the previous run: %s\n", buf );
		}
		// Preserve the existing restart zip file
		if( Ftest( cd.restart_zip_file ) == 0 )
		{
			if( cd.pardebug ) tprintf( "Previous restart file (%s) exists!\n", cd.restart_zip_file );
			if( cd.restart ) sprintf( buf, "cp %s %s.restart_%s_%s.zip >& /dev/null", cd.restart_zip_file, op.root, cd.datetime_infile, Fdatetime( cd.restart_zip_file, 0 ) );  // Copy if restart
			else sprintf( buf, "mv %s %s.restart_%s_%s.zip >& /dev/null", cd.restart_zip_file, op.root, cd.datetime_infile, Fdatetime( cd.restart_zip_file, 0 ) );  // Move if no restart
			system( buf );
		}
		if( cd.restart == 0 )
		{
			sprintf( filename, "%s.restart_info", op.root );
			out = Fwrite( filename );
			fprintf( out, "%s\n", op.datetime_stamp );
			for( i = 0; i < argn; i++ )
				fprintf( out, "%s ", argv[i] );
			fprintf( out, "\n" );
			fclose( out );
			sprintf( buf, "zip %s %s.restart_info >& /dev/null", cd.restart_zip_file, op.root );
			system( buf );
		}
	}
	sprintf( filename, "%s.cmdline_hist", op.root );
	out = Fappend( filename );
	fprintf( out, "%s :", op.datetime_stamp );
	for( i = 0; i < argn; i++ )
		fprintf( out, " %s", argv[i] );
	fprintf( out, "\n" );
	fclose( out );
	//
	// DONE with file reading and problem setup
	//
	// Model analyses are performed below based on provided inputs
	//
	status = 1;
	// ------------------------------------------------------------------------------------------------ CHECK
	if( cd.problem_type == CHECK ) /* Check model input files */
		status = check( &op );
	// ------------------------------------------------------------------------------------------------ INFO-GAP
	// ------------------------------------------------------------------------------------------------ IGRND
	if( cd.problem_type == CALIBRATE && cd.calib_type == IGRND ) /* Calibration analysis using random initial guessed */
		status = igrnd( &op );
	// ------------------------------------------------------------------------------------------------ IGPD
	if( cd.problem_type == CALIBRATE && cd.calib_type == IGPD ) /* Calibration analysis using discretized initial guesses */
		status = igpd( &op );
	// ------------------------------------------------------------------------------------------------ PPSD
	if( cd.problem_type == CALIBRATE && cd.calib_type == PPSD ) /* Calibration analysis using discretized parameters */
		status = ppsd( &op );
	// ------------------------------------------------------------------------------------------------ MONTECARLO
	if( cd.problem_type == MONTECARLO ) /* Monte Carlo analysis */
		status = montecarlo( &op );
	// ------------------------------------------------------------------------------------------------ GLOBALSENS
	if( cd.problem_type == GLOBALSENS ) /* Global sensitivity analysis */
	{
		switch( cd.gsa_type )
		{
			case SOBOL: status = sa_sobol( &op ); break;
			case SALTELLI: status = sa_saltelli( &op ); break;
			case MOAT: status = sa_moat( &op ); break;
			default: tprintf( "WARNING: unknown global sensitivity analysis type; Sobol's method assumed" ); status = sa_sobol( &op ); break;
		}
	}
	// ------------------------------------------------------------------------------------------------ SIMPLE CALIBRATION
	if( cd.problem_type == CALIBRATE && cd.calib_type == SIMPLE ) /* Inverse analysis */
	{
		if( cd.nretries > 1 || cd.paranoid ) tprintf( "\nMULTI-START CALIBRATION using a series of random initial guesses:\n" );
		else tprintf( "\nSINGLE CALIBRATION: single optimization based on initial guesses provided in the input file:\n" );
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
		else optimize_func = optimize_pso; // Define optimization method: PSO
		for( i = 0; i < pd.nParam; i++ ) cd.var[i] = pd.var[i]; // Set all the initial values
		success = optimize_func( &op ); // Optimize
		if( success == 0 ) { tprintf( "ERROR: Optimization did not start!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		if( cd.debug == 0 ) tprintf( "\n" );
		print_results( &op, 1 );
		save_final_results( "", &op, &gd );
		predict = 0;
	}
	if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	strcpy( op.label, "" ); // No labels needed below
	// ------------------------------------------------------------------------------------------------ EIGEN || LOCALSENS
	if( cd.problem_type == EIGEN || cd.problem_type == LOCALSENS )
		status = eigen( &op, NULL, NULL, NULL ); // Eigen or sensitivity analysis run
	// ------------------------------------------------------------------------------------------------ ABAGUS
	if( cd.problem_type == ABAGUS ) // Particle swarm sensitivity analysis run
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = abagus( &op );
	}
	// ------------------------------------------------------------------------------------------------ POSTPUA
	if( cd.problem_type == POSTPUA ) // Predictive uncertainty analysis of sampling results
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = postpua( &op );
	}
	// ------------------------------------------------------------------------------------------------ GLUE
	if( cd.problem_type == GLUE ) // Generalized Likelihood Uncertainty Estimation
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = glue( &op );
	}
	// ------------------------------------------------------------------------------------------------ INFOGAP
	if( cd.problem_type == INFOGAP ) // Info-gap decision analysis
	{
		if( fabs( cd.obsstep ) > DBL_EPSILON )
			status = infogap_obs( &op );
		else
		{
			if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
			status = infogap( &op );
		}
	}
	if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	// ------------------------------------------------------------------------------------------------ FORWARD
	if( cd.problem_type == FORWARD ) // Forward run
	{
		if( cd.resultsfile[0] != 0 )
		{
			FILE *infile2;
			char bigbuffer[5000];
			infile2 = Fread( cd.resultsfile );
			int lines = 0, caseid;
			while( !feof( infile2 ) )
			{
				if( fgets( bigbuffer, sizeof bigbuffer, infile2 ) == NULL )
					break;
				if( sscanf( bigbuffer, "%d", &caseid ) ) lines++;
			}
			tprintf( "\nModel parameters will be initiated based on previously saved results in file %s (total number of cases %d)\n", cd.resultsfile, lines );
			if( cd.calib_type == PPSD ) tprintf( "PPSD format assumed\n" );
			else if( cd.calib_type == IGRND ) tprintf( "IGRND format assumed\n" );
			if( cd.resultscase < 0 ) tprintf( "Forward runs based on first %d cases.", -cd.resultscase );
			else
			{
				int implemented = 0;
				if( cd.phi_cutoff > DBL_EPSILON ) { tprintf( "Model analyses for cases with phi < %g.", cd.phi_cutoff ); implemented = 1; }
				if( cd.obsrange ) { tprintf( "Model analyses for successful cases." ); implemented = 1; }
				if( !implemented )
				{
					if( cd.resultscase == 0 ) cd.resultscase = 1;
					if( cd.resultscase > 0 ) tprintf( "Model analyses for case #%d", cd.resultscase );
					else tprintf( "Model analyses for first %d cases", -cd.resultscase );
				}
				else cd.resultscase = -lines;
			}
			tprintf( "\n\n" );
		}
		if( cd.resultscase < 0 )
		{
			FILE *infile2;
			char bigbuffer[5000], bigbufferorig[5000], dummy[50], *start, *word, *separator = " ";
			int cases = cd.resultscase * -1;
			int caseid;
			bad_data = 0;
			infile2 = Fread( cd.resultsfile );
			double *res;
			if( ( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( res = ( double * ) malloc( od.nTObs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			for( j = 0; j < cases; j++ )
			{
				if( feof( infile2 ) )
				{
					bad_data = 1;
					tprintf( "ERROR File end reached (%s; case %d)\n", cd.resultsfile, j + 1 );
					break;
				}
				if( fgets( bigbuffer, sizeof bigbuffer, infile2 ) == NULL )
				{
					bad_data = 1;
					tprintf( "ERROR reading model parameters initiated based on previously saved results in file %s (case %d)\n", cd.resultsfile, j + 1 );
					break;
				}
				strcpy( bigbufferorig, bigbuffer );
				sscanf( bigbuffer, "%d", &caseid );
				// if( cd.debug ) tprintf( "\n" );
				tprintf( "Case ID %d in %s (case %d)\n", caseid, cd.resultsfile, j + 1 );
				k = 0;
				if( cd.calib_type == PPSD )
				{
					start = strstr( bigbuffer, ":" );
					word = strtok( start, separator ); // skip :
					for( i = 0; i < pd.nParam; i++ )
					{
						while( pd.var_opt[i] != 2 ) { i++; if( i >= pd.nParam ) break; }
						if( !( i < pd.nParam ) ) break;
						word = strtok( NULL, separator );
						// tprintf( "Par #%d Word %s\n", i, word );
						sscanf( word, "%lf", &pd.var[i] );
						k++;
						if( pd.var_log[i] ) pd.var[i] = log10( pd.var[i] );
						cd.var[i] = pd.var[i];
						if( cd.debug ) tprintf( "%s %g\n", pd.var_id[i], pd.var[i] );
					}
					if( cd.debug ) tprintf( "Number of initialized fixed parameters in previous PPSD simulation = %d\n", k );
					if( pd.nFlgParam != k )
					{
						bad_data = 1;
						tprintf( "ERROR Number of flagged (%d) and initialized (%d) parameters in %s (case %d) do not match\n", pd.nFlgParam, k, cd.resultsfile, j + 1 );
						break;
					}
				}
				strcpy( bigbuffer, bigbufferorig );
				start = strstr( bigbuffer, ": OF " );
				if( start == NULL ) tprintf( "WARNING Objective function value is missing in %s (case %d)\n", cd.resultsfile, j + 1 );
				else
				{
					sscanf( start, ": OF %lg", &phi );
					tprintf( "Objective function = %g\n", phi );
				}
				if( cd.phi_cutoff > DBL_EPSILON && phi > cd.phi_cutoff )
				{
					tprintf( "Case skipped: phi %g > cutoff %g\n", phi, cd.phi_cutoff );
					continue;
				}
				// tprintf( "%s\n", bigbuffer );
				start = strcasestr( bigbuffer, "success" );
				if( start == NULL ) tprintf( "WARNING Success value is missing in %s (case %d)\n", cd.resultsfile, j + 1 );
				else
				{
					sscanf( start, "%s %d", dummy, &success );
					tprintf( "Success = %d\n", success );
				}
				if( cd.obsrange && !success )
				{
					tprintf( "Case skipped: no success\n" );
					continue;
				}
				start = strcasestr( bigbuffer, "final var" );
				if( start == NULL )
				{
					bad_data = 1;
					tprintf( "ERROR Final model parameters cannot be located in %s (case %d)\n", cd.resultsfile, j + 1 );
					break;
				}
				// tprintf( "%s\n", start );
				strcpy( bigbuffer, start );
				// c = -2 final var
				// c - word count
				// k - number of processed variables
				// i - variable index
				for( i = 0, k = 0, c = -2, word = strtok( bigbuffer, separator ); word; c++, word = strtok( NULL, separator ) )
				{
					if( c > -1 )
					{
						// tprintf( "Par #%d %s\n", c + 1, word );
						if( cd.calib_type == PPSD ) while( pd.var_opt[i] == 2 || pd.var_opt[i] == 0 ) i++;
						else while( pd.var_opt[i] == 0 ) i++;
						if( i >= pd.nParam ) break;
						sscanf( word, "%lf", &pd.var[i] );
						k++;
						if( pd.var_log[i] ) pd.var[i] = log10( pd.var[i] );
						cd.var[i] = pd.var[i];
						if( cd.debug ) tprintf( "%s %g\n", pd.var_id[i], pd.var[i] );
						i++;
					}
				}
				if( cd.debug ) tprintf( "Number of initialized parameters = %d\n\n", k );
				if( pd.nOptParam != k )
				{
					bad_data = 1;
					tprintf( "ERROR Number of optimized (%d) and initialized (%d) parameters in %s (case %d) do not match\n", pd.nOptParam, k, cd.resultsfile, j + 1 );
					break;
				}
				for( i = 0; i < pd.nOptParam; i++ )
					opt_params[i] = pd.var[pd.var_index[i]];
				// for( i = 0; i < pd.nOptParam; i++ )
				// tprintf( "%s %g\n", pd.var_id[pd.var_index[i]], opt_params[i] );
				Transform( opt_params, &op, opt_params );
				func_global( opt_params, &op, res );
				tprintf( "Rerun results:\n" );
				tprintf( "Objective function = %g\n", op.phi );
				tprintf( "Success = %d\n", op.success );
				if( cd.debug )
				{
					tprintf( "\n" );
					print_results( &op, 1 );
				}
				if( cd.phi_cutoff > DBL_EPSILON && phi > cd.phi_cutoff )
				{
					tprintf( "Case skipped: phi %g > cutoff %g\n", phi, cd.phi_cutoff );
					continue;
				}
				if( cd.obsrange && !op.success )
				{
					tprintf( "Case skipped: no success\n" );
					continue;
				}
				if( cd.save )
				{
					if( !cd.debug )
					{
						tprintf( "\n" );
						print_results( &op, 1 );
					}
					sprintf( filename, "%d", caseid );
					save_final_results( filename, &op, &gd );
				}
				tprintf( "\n" );
			}
			free( opt_params );
			free( res );
			fclose( infile2 );
			predict = 0;
		}
		else
		{
			if( cd.resultscase ) sprintf( filename, "%s.%d.results", op.root, cd.resultscase );
			else sprintf( filename, "%s.results", op.root );
			out = Fwrite( filename );
			tprintf( "\nModel parameter values:\n" );
			fprintf( out, "Model parameter values:\n" );
			for( i = 0; i < pd.nParam; i++ )
			{
				if( pd.var_opt[i] && pd.var_log[i] ) cd.var[i] = pow( 10, pd.var[i] );
				else cd.var[i] = pd.var[i];
				tprintf( "%s %g\n", pd.var_id[i], cd.var[i] );
				fprintf( out, "%s %g\n", pd.var_id[i], cd.var[i] );
			}
			fflush( out );
			if( od.nTObs > 0 || wd.nW > 0 )
			{
				tprintf( "\nModel predictions (forward run; no calibration):\n" );
				fprintf( out, "\nModel predictions (forward run; no calibration):\n" );
				fflush( out );
				if( cd.resultscase ) sprintf( filename, "%s.%d.forward", op.root, cd.resultscase );
				else sprintf( filename, "%s.forward", op.root );
				out2 = Fwrite( filename );
				predict = 1;
			}
			else
			{
				tprintf( "\nNo model predictions!\n" );
				exit( 1 );
			}
		}
	}
	//
	// ------------------------------------------------------------------------------------------------ CREATE
	//
	if( cd.problem_type == CREATE ) // Create a MADS file based on a forward run
	{
		tprintf( "\nModel predictions (forward run; no calibration):\n" );
		predict = 1;
	}
	//
	// ------------------------------------------------------------------------------------------------ PREDICT
	//
	if( predict )
	{
		success_all = 1;
		compare = 0;
		phi = 0;
		if( cd.solution_type[0] == EXTERNAL )
		{
			double *res;
			if( ( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( res = ( double * ) malloc( od.nTObs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			for( i = 0; i < pd.nOptParam; i++ )
				opt_params[i] = pd.var[pd.var_index[i]];
			Transform( opt_params, &op, opt_params );
			func_extrn( opt_params, &op, res );
			free( opt_params );
			free( res );
			for( i = 0; i < od.nTObs; i++ )
			{
				if( cd.problem_type == CALIBRATE && od.obs_weight[i] != 0 ) { if( od.nTObs > 50 && i == 21 ) tprintf( "...\n" ); continue; }
				compare = 1;
				c = od.obs_current[i];
				err = od.obs_target[i] - c;
				// err = od.res[i];
				min = od.obs_min[i];
				max = od.obs_max[i];
				if( cd.objfunc_type != SSR )
				{
					if( cd.objfunc_type == SSDA )
					{
						err = sqrt( fabs( err ) );
						if( c < od.obs_target[i] ) err *= -1;
					}
					else err = 0; // SSD0 & SSDX
					if( cd.objfunc_type == SSDX ) { dx = max - min; if( cd.obsdomain > DBL_EPSILON && cd.obsdomain < dx ) dx = cd.obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
					if( c < min ) err += min - c;
					else if( c > max ) err += c - max;
					if( cd.objfunc_type == SSDX ) { min = od.obs_min[i]; max = od.obs_max[i]; }
				}
				phi += ( err * err ) * od.obs_weight[i];
				if( ( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) && ( od.obs_weight[i] > 0.0 ) ) { success_all = 0; success = 0; }
				else success = 1;
				if( od.nTObs < 50 || ( i < 20 || i > od.nTObs - 20 ) ) tprintf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", od.obs_id[i], od.obs_target[i], c, err, err * od.obs_weight[i], success, min, max );
				if( od.nTObs > 50 && i == 21 ) tprintf( "...\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", od.obs_id[i], od.obs_target[i], c, err, err * od.obs_weight[i], success, min, max );
				else od.obs_target[i] = c; // Save computed values as calibration targets
			}
		}
		else if( cd.solution_type[0] != TEST )
			for( i = 0; i < wd.nW; i++ )
				for( j = 0; j < wd.nWellObs[i]; j++ )
				{
					if( cd.problem_type == CALIBRATE && wd.obs_weight[i][j] > DBL_EPSILON ) continue;
					compare = 1;
					c = func_solver( wd.x[i], wd.y[i], wd.z1[i], wd.z2[i], wd.obs_time[i][j], &cd );
					err = wd.obs_target[i][j] - c;
					min = wd.obs_min[i][j];
					max = wd.obs_max[i][j];
					if( cd.objfunc_type != SSR )
					{
						if( cd.objfunc_type == SSDA )
						{
							err = sqrt( fabs( err ) );
							if( c < wd.obs_target[i][j] ) err *= -1;
						}
						else err = 0; // SSD0 & SSDX
						if( cd.objfunc_type == SSDX ) { dx = max - min; if( cd.obsdomain > DBL_EPSILON && cd.obsdomain < dx ) dx = cd.obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
						if( c < min ) err += min - c;
						else if( c > max ) err += c - max;
						if( cd.objfunc_type == SSDX ) { min = od.obs_min[i]; max = od.obs_max[i]; }
					}
					if( cd.problem_type != CALIBRATE ) phi += ( err * err ) * wd.obs_weight[i][j];
					else phi += ( err * err );
					if( ( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) && ( wd.obs_weight[i][j] > 0.0 ) ) { success_all = 0; success = 0; }
					else success = 1;
					if( cd.problem_type != CALIBRATE )
						tprintf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err * wd.obs_weight[i][j], success, min, max );
					else
						tprintf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err, success, min, max );
					if( cd.problem_type != CREATE ) fprintf( out, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err * wd.obs_weight[i][j], success );
					else wd.obs_target[i][j] = c; // Save computed values as calibration targets
					if( cd.problem_type == FORWARD ) fprintf( out2, "%s(%g) %g\n", wd.id[i], wd.obs_time[i][j], c ); // Forward run
				}
		if( cd.problem_type == FORWARD && od.nTObs > 0 ) fclose( out2 );
		cd.neval++;
		if( compare )
		{
			op.phi = phi;
			tprintf( "Objective function: %g Success: %d\n", op.phi, success_all );
			if( cd.problem_type != CREATE ) fprintf( out, "Objective function = %g Success: %d\n", op.phi, success_all );
			if( success_all )
			{
				tprintf( "All the predictions are within acceptable ranges!\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "All the predictions are within acceptable ranges!\n" );
			}
			else
			{
				tprintf( "At least one of the predictions is outside acceptable ranges!\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "At least one of the predictions is outside acceptable ranges!\n" );
			}
		}
		else    tprintf( "No calibration targets!\n" );
		if( cd.problem_type != CREATE ) fclose( out );
	}
	if( predict && cd.problem_type == FORWARD ) save_final_results( "", &op, &gd );
	if( predict ) // Write phi in a separate file
	{
		if( od.nTObs > 0 )
		{
			if( cd.problem_type == FORWARD && cd.resultscase > 0 ) sprintf( filename, "%s.%d.phi", op.root, cd.resultscase );
			else sprintf( filename, "%s.phi", op.root );
			out2 = Fwrite( filename );
			fprintf( out2, "%g\n", op.phi );
			fclose( out2 );
		}
	}
	if( cd.problem_type == CREATE ) /* Create a file with calibration targets equal to the model predictions */
	{
		cd.problem_type = CALIBRATE;
		sprintf( filename, "%s-truth.mads", op.root );
		save_problem( filename, &op );
		tprintf( "\nMADS problem file named %s-truth.mads is created; modify the file if needed\n\n", op.root );
		cd.problem_type = CREATE;
	}
	//
	// Finalize the run
	//
	time_end = time( NULL );
	time_elapsed = time_end - time_start;
	if( time_elapsed > 86400 ) tprintf( "Simulation time = %g days\n", ( ( double ) time_elapsed / 86400 ) );
	else if( time_elapsed > 3600 ) tprintf( "Simulation time = %g hours\n", ( ( double ) time_elapsed / 3600 ) );
	else if( time_elapsed > 60 ) tprintf( "Simulation time = %g minutes\n", ( ( double ) time_elapsed / 60 ) );
	else tprintf( "Simulation time = %ld seconds\n", time_elapsed );
	tprintf( "Functional evaluations = %d\n", cd.neval );
	if( cd.njac > 0 ) tprintf( "Jacobian evaluations = %d\n", cd.njac );
	if( cd.nlmo > 0 ) tprintf( "Levenberg-Marquardt optimizations = %d\n", cd.nlmo );
	if( time_elapsed > 0 )
	{
		c = cd.neval / time_elapsed;
		if( c < ( ( double ) 1 / 86400 ) ) tprintf( "Functional evaluations per day = %g\n", c * 86400 );
		else if( c < ( ( double ) 1 / 3600 ) ) tprintf( "Functional evaluations per hour = %g\n", c * 3600 );
		else if( c < ( ( double ) 1 / 60 ) ) tprintf( "Functional evaluations per minute = %g\n", c * 60 );
		else tprintf( "Functional evaluations per second = %g\n", c );
	}
	if( op.cd->seed_init > 0 ) tprintf( "Seed = %d\n", op.cd->seed_init );
	ptr_ts = localtime( &time_start );
	tprintf( "Execution  started  on %s", asctime( ptr_ts ) );
	ptr_ts = localtime( &time_end );
	tprintf( "Execution completed on %s", asctime( ptr_ts ) );
	tprintf( "Execution date & time stamp: %s\n", op.datetime_stamp );
	sprintf( buf, "rm -f %s.running", op.root ); system( buf );
	if( op.f_ofe != NULL ) { fclose( op.f_ofe ); op.f_ofe = NULL; }
	free( op.cd->solution_id ); free( op.cd->solution_type );
	fclose( mads_output );
	exit( 0 ); // DONE
}

int optimize_pso( struct opt_data *op )
{
	if( op->cd->debug ) tprintf( "\nParticle-Swarm Optimization:" );
	if( strncasecmp( op->cd->opt_method, "pso", 3 ) == 0 || strncasecmp( op->cd->opt_method, "swarm", 5 ) == 0 )
	{
		if( op->cd->debug ) tprintf( " Standard (2006)\n" );
		pso_std( op );
	}
	else if( strncasecmp( op->cd->opt_method, "tribes", 5 ) == 0 && strcasestr( op->cd->opt_method, "std" ) != NULL )
	{
		if( op->cd->debug ) tprintf( " TRIBES (Clerc, 2006)\n" );
		mopso( op );
	}
	else
	{
		if( op->cd->debug )
		{
			if( strncasecmp( op->cd->opt_method, "squads", 5 ) == 0 || strcasestr( op->cd->opt_method, "lm" ) != NULL )
				tprintf( " SQUADS\n" );
			else
				tprintf( " TRIBES (modified based on Clerc, 2006)\n" );
		}
		pso_tribes( op );
	}
	if( op->cd->debug )
	{
		tprintf( "\n------------------------- Optimization Results:\n" );
		print_results( op, 1 );
	}
	if( op->cd->lm_eigen && op->cd->solution_type[0] != TEST )
		if( eigen( op, NULL, NULL, NULL ) == 0 ) // Execute eigen analysis of the final results
			return( 0 );
	return( 1 );
}

int optimize_lm( struct opt_data *op )
{
	double phi, phi_min;
	double *opt_params, *opt_params_best, *res, *x_c;
	int   nsig, maxfn, maxiter, maxiter_levmar, iopt, infer, ier, debug, standalone;
	int   i, j, k, debug_level, count, count_set, npar;
	double opt_parm[4], *jacobian, *jacTjac, *covar, *work, eps, delta, *var_lhs;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	char buf[80];
	if( op->cd->maxeval <= op->cd->neval )
	{
		tprintf( "WARNING: LM optimization cannot be performed! Number of the maximum evaluations has been exceeded (%d<%d)\n", op->cd->maxeval, op->cd->neval );
		return( 1 );
	}
	debug = MAX( op->cd->debug, op->cd->ldebug );
	standalone = op->cd->lmstandalone;
	if( debug == 0 && standalone && op->cd->calib_type != PPSD ) op->cd->lmstandalone = 2;
	if( op->cd->squads ) standalone = op->cd->lmstandalone = 0;
	if( op->od->nTObs == 0 ) { tprintf( "ERROR: Number of observations is equal to zero! Levenberg-Marquardt Optimization cannot be performed!\n" ); return( 0 ); }
	if( op->pd->nOptParam == 0 ) { tprintf( "ERROR: Number of optimized model parameters is equal to zero! Levenberg-Marquardt Optimization cannot be performed!\n" ); return( 0 ); }
	if( ( op->pd->nOptParam > op->od->nTObs ) && ( !op->cd->squads && op->cd->calib_type == SIMPLE ) ) { tprintf( "WARNING: Number of optimized model parameters is greater than number of observations (%d>%d)\n", op->pd->nOptParam, op->od->nTObs ); }
	gsl_matrix *gsl_jacobian = gsl_matrix_alloc( op->od->nTObs, op->pd->nOptParam );
	gsl_matrix *gsl_covar = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam );
	gsl_vector *gsl_opt_params = gsl_vector_alloc( op->pd->nOptParam );
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( opt_params_best = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( x_c = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( res = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( op->cd->niter <= 0 )
	{
		if( op->cd->squads ) maxiter = 8;
		else maxiter = 50;
	}
	else maxiter = op->cd->niter;
	if( op->cd->ldebug && standalone ) tprintf( "Number of Levenberg-Marquardt iterations = %d\n", maxiter );
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]];
	//	for( i = 0; i < op->pd->nOptParam; i++ )
	//		tprintf( "lmi %g\n", opt_params[i] );
	if( !op->cd->squads ) Transform( opt_params, op, opt_params ); // No need to transform if part of SQUADS runs
	//	for( i = 0; i < op->pd->nOptParam; i++ )
	//		tprintf( "lmi %g\n", opt_params[i] );
	if( op->cd->paranoid )
	{
		tprintf( "Multi-Start Levenberg-Marquardt (MSLM) Optimization ... " );
		npar = op->pd->nOptParam;
		if( op->cd->nretries <= 0 ) op->cd->nretries = ( int )( ( double )( op->cd->maxeval - op->cd->neval ) / ( maxiter * npar / 10 ) );
		if( debug ) tprintf( "\nRandom sampling for MSLM optimization (variables %d; realizations %d) using ", npar, op->cd->nretries );
		if( ( var_lhs = ( double * ) malloc( npar * op->cd->nretries * sizeof( double ) ) ) == NULL )
		{ tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( op->cd->seed < 0 ) { op->cd->seed *= -1; if( debug ) tprintf( "Imported seed: %d\n", op->cd->seed ); }
		else if( op->cd->seed == 0 ) { if( debug ) tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
		else if( debug ) tprintf( "Current seed: %d\n", op->cd->seed );
		if( op->cd->paran_method[0] != 0 ) { strcpy( buf, op->cd->smp_method ); strcpy( op->cd->smp_method, op->cd->paran_method ); }
		sampling( npar, op->cd->nretries, &op->cd->seed, var_lhs, op, debug );
		if( op->cd->paran_method[0] != 0 ) strcpy( op->cd->smp_method, buf );
		if( debug ) tprintf( "done.\n" );
		op->cd->retry_ind = count = count_set = 0;
	}
	else if( standalone ) tprintf( "Levenberg-Marquardt Optimization ... " );
	phi_min = HUGE_VAL;
	do // BEGIN Paranoid loop
	{
		if( op->cd->maxeval <= op->cd->neval ) { if( debug || op->cd->paranoid ) tprintf( "Maximum number of evaluations is exceeded (%d <= %d)!\n", op->cd->maxeval, op->cd->neval ); break; }
		op->cd->nlmo++;
		if( op->cd->paranoid )
		{
			count++;
			op->cd->retry_ind = count;
			if( count == 1 )
			{
				if( op->cd->calib_type == IGRND )
				{
					if( debug ) tprintf( "\nCALIBRATION %d: initial guesses from IGRND random set: ", count );
					for( i = 0; i < op->pd->nOptParam; i++ )
					{
						k = op->pd->var_index[i];
						opt_params[i] = op->pd->var[k];
						if( debug > 1 )
						{
							if( op->pd->var_log[k] ) tprintf( "%s %.15g\n", op->pd->var_id[k], pow( 10, opt_params[i] ) );
							else tprintf( "%s %.15g\n", op->pd->var_id[k], opt_params[i] );
						}
					}
				}
				else
				{
					if( debug ) tprintf( "\nCALIBRATION %d: initial guesses from MADS input file #: ", count );
					for( i = 0; i < op->pd->nOptParam; i++ )
					{
						k = op->pd->var_index[i];
						opt_params[i] = op->pd->var[k];
						if( debug > 1 )
						{
							if( op->pd->var_log[k] ) tprintf( "%s %.15g\n", op->pd->var_id[k], pow( 10, opt_params[i] ) );
							else tprintf( "%s %.15g\n", op->pd->var_id[k], opt_params[i] );
						}
					}
				}
			}
			else if( count > 1 )
			{
				if( op->cd->ldebug ) tprintf( "\n********************************************************************\n" );
				if( debug ) tprintf( "CALIBRATION %d: initial guesses from internal MSLM random set #%d: ", count, count_set + 1 );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					opt_params[i] = var_lhs[i + count_set * npar] * op->pd->var_range[k] + op->pd->var_min[k];
					if( debug > 1 )
					{
						if( op->pd->var_log[k] ) tprintf( "%s %.15g\n", op->pd->var_id[k], pow( 10, opt_params[i] ) );
						else tprintf( "%s %.15g\n", op->pd->var_id[k], opt_params[i] );
					}
				}
				count_set++;
			}
			if( debug > 1 ) tprintf( "\n" );
			Transform( opt_params, op, opt_params );
		}
		if( ( op->cd->debug > 1 || op->cd->ldebug > 5 ) && standalone )
		{
			tprintf( "\n-------------------- Initial state:\n" );
			op->cd->pderiv = op->cd->oderiv = -1;
			debug_level = op->cd->fdebug; op->cd->fdebug = 3;
			func_global( opt_params, op, res );
			op->cd->fdebug = debug_level;
		}
		// LM optimization ...
		if( strcasestr( op->cd->opt_method, "mon" ) != NULL || strcasestr( op->cd->opt_method, "chav" ) != NULL ) // Monty/Chavo versions
		{
			if( debug > 1 && standalone ) tprintf( "\nLevenberg-Marquardt Optimization:\n" );
			else if( op->cd->ldebug ) tprintf( "\n" );
			if( ( jacobian = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->od->nTObs ) ) == NULL ) { tprintf( "ERROR: Not enough memory!\n" ); return( 0 ); }
			if( ( jacTjac = ( double * ) malloc( sizeof( double ) * ( ( op->pd->nOptParam + 1 ) * op->pd->nOptParam / 2 ) ) ) == NULL ) { tprintf( "ERROR: Not enough memory!\n" ); return( 0 ); }
			iopt = 2; /*    iopt=0 Brown's algorithm without strict descent
		       		iopt=1 strict descent and default values for input vector parm
					iopt=2 strict descent with user parameter choices in input vector parm */
			opt_parm[0] = 10; /* initial value of the Marquardt parameter */
			opt_parm[1] = 2.0; /* scaling factor used to modify the Marquardt parameter */
			opt_parm[2] = 1e30; /* upper bound for increasing the Marquardt parameter */
			opt_parm[3] = 500; /* value for indicating when central differencing is to be used for calculating the jacobian */
			/* First convergence criterion */
			nsig = 8; /* parameter estimates agree to nsig digits */
			/* Second convergence criterion */
			eps = epsilon(); /* two successive iterations the residual sum of squares estimates have relative difference less than or equal to eps */
			/* Third convergence criterion */
			delta = 0; /* norm of the approximate gradient is less than or equal to delta */
			maxfn = op->cd->maxeval - op->cd->neval; /* maximum number of function evaluations; remove the number of evaluation already performed */
			if( strcasestr( op->cd->opt_method, "chav" ) != NULL ) ier = zxssqch( func_global, op, op->od->nTObs, op->pd->nOptParam, nsig, eps, delta, maxfn, iopt, opt_parm, opt_params, &phi, res, jacobian, op->od->nTObs, jacTjac, &infer ); // Chavo's version
			else ier = lm_opt( func_global, func_dx, op, op->od->nTObs, op->pd->nOptParam, nsig, eps, delta, maxfn, maxiter, iopt, opt_parm, opt_params, &phi, res, jacobian, op->od->nTObs, jacTjac, &infer ); // Monty's version
			for( k = i = 0; i < op->pd->nOptParam; i++ )
				for( j = 0; j < op->od->nTObs; j++, k++ )
					gsl_matrix_set( gsl_jacobian, j, i, jacobian[k] );
			gsl_multifit_covar( gsl_jacobian, 0.0, gsl_covar );
			free( jacTjac ); free( jacobian );
			op->phi = phi;
		}
		else if( strcasestr( op->cd->opt_method, "gsl" ) != NULL ) // GSL version of LM
		{
			if( debug > 1 && standalone ) tprintf( "\nLevenberg-Marquardt Optimization using GSL library:\n" );
			else if( op->cd->ldebug ) tprintf( "\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
				gsl_vector_set( gsl_opt_params, i, opt_params[i] );
			lm_gsl( gsl_opt_params, op, gsl_jacobian, gsl_covar );
			for( i = 0; i < op->pd->nOptParam; i++ )
				opt_params[i] = gsl_vector_get( gsl_opt_params, i );
			phi = op->phi;
		}
		else if( strcasestr( op->cd->opt_method, "tra" ) != NULL )// Transtrum version of LM
		{
			if( debug > 1 && standalone ) tprintf( "\nTranstrum version of Levenberg-Marquardt Optimization:\n" ); // TODO add explicit transtrum
			else if( op->cd->ldebug ) tprintf( "\n" );
		}
		else // DEFAULT LevMar version of LM
		{
			if( debug > 1 && standalone ) tprintf( "\nLevenberg-Marquardt Optimization using LevMar library:\n" );
			else if( op->cd->ldebug ) tprintf( "\n" );
			if( ( covar = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->pd->nOptParam ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			// LM_DIF_WORKSZ(m,n) = 4*n+4*m + n*m + m*m
			if( ( work = ( double * ) malloc( sizeof( double ) * LM_DIF_WORKSZ( op->pd->nOptParam, op->od->nTObs ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			for( i = 0; i < op->od->nTObs; i++ ) res[i] = 0;
			jacobian = work + op->pd->nOptParam + 2 * op->od->nTObs;
			opts[0] = op->cd->lm_error * 100; opts[1] = op->cd->lm_error; opts[2] = op->cd->lm_error;
			opts[3] = op->cd->phi_cutoff;
			if( op->cd->sintrans == 0 ) opts[4] = op->cd->lindx; // Forward difference; Central difference if negative; DO NOT USE CENTRAL DIFFERENCE
			else opts[4] = op->cd->sindx;
			while( op->cd->maxeval > op->cd->neval )
			{
				// Levmar has no termination criteria based on the number of functional evaluations or number of jacobian evaluations
				// ... it uses the number of LM iterations that count the jacobian and lambda evaluations
				// ... assuming about 10 lambda searches per jacobian iteration
				// Levmar has now adder terminaitonal based on op->cd->maxeval and op->cd->niter
				if( opts[4] > 0 ) maxiter_levmar = ( int )( ( double )( ( op->cd->maxeval - op->cd->neval ) / ( op->pd->nOptParam + 1 ) + 1 ) ); // Forward derivatives
				else              maxiter_levmar = ( int )( ( double )( ( op->cd->maxeval - op->cd->neval ) / ( 2 * op->pd->nOptParam + 1 ) + 1 ) ); // Central derivatives
				if( maxiter_levmar > maxiter ) maxiter_levmar = maxiter;
				// dlevmar_der called by DEFAULT
				if( strcasestr( op->cd->opt_method, "dif" ) != NULL ) ier = dlevmar_dif( func_levmar, opt_params, res, op->pd->nOptParam, op->od->nTObs, maxiter_levmar, opts, info, work, covar, op );
				else ier = dlevmar_der2( func_levmar, func_dx_levmar, opt_params, res, op->pd->nOptParam, op->od->nTObs, maxiter_levmar, opts, info, work, covar, op );
				if( info[6] == 4 || info[6] == 5 )
				{
					if( op->cd->maxeval > op->cd->neval )
					{
						opts[0] *= 10;
						tprintf( "\n\nIMPORTANT: LM optimization rerun with larger initial lambda (%g)\n\n", opts[0] );
					}
					else tprintf( "\n\nWARNING: LM optimization may benefit from rerun with larger initial lambda (%g)\n\n", opts[0] );
				}
				else break;
			}
			if( op->cd->ldebug > 5 )
			{
				tprintf( "Levenberg-Marquardt Optimization completed after %g iteration (reason %g) (returned value %d)\n", info[5], info[6], ier );
				tprintf( "initial phi %g final phi %g ||J^T e||_inf %g ||Dp||_2 %g mu/max[J^T J]_ii %g\n", info[0], info[1], info[2], info[3], info[4] );
				tprintf( "function evaluation %g jacobian evaluations %g linear systems solved %g\n", info[7], info[8], info[9] );
			}
			op->cd->njac += ( int ) info[8];
			for( k = j = 0; j < op->od->nTObs; j++ )
				for( i = 0; i < op->pd->nOptParam; i++ )
					gsl_matrix_set( gsl_jacobian, j, i, jacobian[k++] );
			for( i = 0; i < op->pd->nOptParam; i++ )
				for( j = 0; j < op->pd->nOptParam; j++ )
					gsl_matrix_set( gsl_covar, i, j, covar[i * op->pd->nOptParam + j] );
			op->phi = phi = info[1];
			free( work ); free( covar );
		}
		if( !op->cd->squads ) // if LM is not part of PSO (SQUADS) run
		{
			if( debug > 1 )
			{
				tprintf( "\n------------------------- Final state:\n" );
				op->cd->pderiv = op->cd->oderiv = -1;
				debug_level = op->cd->fdebug; op->cd->fdebug = 3;
				func_global( opt_params, op, op->od->res ); // opt_params are already transformed
				op->cd->fdebug = debug_level;
			}
			else
			{
				// Make a Forward run with the best results
				func_global( opt_params, op, op->od->res ); // opt_params are already transformed
			}
			DeTransform( opt_params, op, x_c );
			for( i = 0; i < op->pd->nOptParam; i++ )
				op->pd->var[op->pd->var_index[i]] = x_c[i]; // Save the obtained results
			if( debug > 1 )
			{
				tprintf( "\n------------------------- LM Optimization Results:\n" );
				print_results( op, 1 );
			}
		}
		else // if LM is part of PSO (SQUADS) run
			for( i = 0; i < op->pd->nOptParam; i++ )
				op->pd->var[op->pd->var_index[i]] = opt_params[i];
		if( op->cd->paranoid )
		{
			if( op->phi < phi_min ) { phi_min = op->phi; for( i = 0; i < op->pd->nOptParam; i++ ) opt_params_best[i] = op->pd->var[op->pd->var_index[i]]; }
			if( debug ) tprintf( "Objective function: %g Success: %d\n", op->phi, op->success );
			if( phi_min < op->cd->phi_cutoff )
			{
				if( debug ) tprintf( "MSLM optimization objective function is below the cutoff value after %d random initial guess attempts\n", count );
				break;
			}
			if( op->cd->check_success && op->success )
			{
				if( debug ) tprintf( "MSLM optimization within calibration ranges after %d random initial guess attempts\n", count );
				phi_min = op->phi;
				for( i = 0; i < op->pd->nOptParam; i++ ) opt_params_best[i] = op->pd->var[op->pd->var_index[i]];
				break;
			}
			if( op->cd->maxeval <= op->cd->neval )
			{ if( debug ) tprintf( "MSLM optimization terminated after evaluations %d (max evaluations %d)\n", op->cd->neval, op->cd->maxeval ); break; }
			if( count == op->cd->nretries )
			{ if( debug ) tprintf( "MSLM optimization terminated after %d attempts (evaluations %d; max evaluations %d)\n", count, op->cd->neval, op->cd->maxeval ); break; }
		}
		else break; // Quit if not Paranoid run
	}
	while( 1 ); // END Paranoid loop
	if( op->cd->paranoid ) // Recompute for the best results
	{
		if( !debug ) tprintf( "(retries=%d) ", count );
		for( i = 0; i < op->pd->nOptParam; i++ )
			op->pd->var[op->pd->var_index[i]] = opt_params_best[i];
		Transform( opt_params_best, op, opt_params );
		func_global( opt_params, op, op->od->res );
	}
	if( ( op->cd->lm_eigen || op->cd->ldebug || op->cd->debug ) && standalone && op->cd->calib_type == SIMPLE )
		if( eigen( op, op->od->res, gsl_jacobian, NULL ) == 0 ) // Eigen analysis
			return( 0 );
	if( !debug && standalone && op->cd->calib_type == SIMPLE ) tprintf( "\n" );
	if( op->cd->paranoid ) free( var_lhs );
	free( opt_params ); free( opt_params_best ); free( x_c ); free( res );
	gsl_matrix_free( gsl_jacobian ); gsl_matrix_free( gsl_covar ); gsl_vector_free( gsl_opt_params );
	return( 1 );
}

int eigen( struct opt_data *op, double *f_x, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar )
{
	FILE *out;
	double phi, stddev_scale, gf;
	double *opt_params, *x_u, *x_d, *stddev, *jacobian;
	double aopt, copt, eopt, dopt, aic, bic, cic, kic, ln_det_v, ln_det_weight, sml, tt;
	int   debug, compute_center, compute_jacobian, compute_covar;
	int   i, j, k, ier, debug_level, status, dof;
	double eps;
	char filename[200], buf[20];
	static double student_dist[34] = {12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042, 2.021, 2.000, 1.980, 1.960 };
	tprintf( "Eigen analysis ...\n" );
	gsl_vector *gsl_opt_params = gsl_vector_alloc( op->pd->nOptParam );
	gsl_matrix *eigenvec = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam );
	gsl_vector *eigenval = gsl_vector_alloc( op->pd->nOptParam );
	gsl_eigen_symmv_workspace *eigenwork = gsl_eigen_symmv_alloc( op->pd->nOptParam );
	compute_center = compute_jacobian = compute_covar = 0;
	if( f_x == NULL ) { compute_center = 1; if( ( f_x = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } }
	if( gsl_jacobian == NULL )
	{
		compute_jacobian = 1;
		gsl_jacobian = gsl_matrix_alloc( op->od->nTObs, op->pd->nOptParam );
		if( ( jacobian = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	}
	if( gsl_covar == NULL ) { gsl_covar = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam ); compute_covar = 1; }
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( x_u = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( x_d = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( stddev = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	debug = MAX( op->cd->lm_eigen, op->cd->debug );
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]];
	Transform( opt_params, op, opt_params );
	for( i = 0; i < op->pd->nOptParam; i++ )
		gsl_vector_set( gsl_opt_params, i, opt_params[i] );
	op->cd->pderiv = op->cd->oderiv = -1;
	if( compute_center )
	{
		if( debug )
		{
			tprintf( "Analyzed state:\n" );
			debug_level = op->cd->fdebug; op->cd->fdebug = 3;
		}
		func_global( opt_params, op, op->od->res );
		for( j = 0; j < op->od->nTObs; j++ )
			f_x[j] = op->od->res[j];
		phi = op->phi;
		if( debug ) op->cd->fdebug = debug_level;
	}
	else // Do not need to compute the center state
	{
		if( debug ) tprintf( "Analyzed state provided externally.\n" );
		phi = op->phi;
	}
	if( compute_jacobian )
	{
		/*
		op->pd->var_current_gsl = gsl_vector_alloc( op->pd->nOptParam );
		op->od->obs_current_gsl = gsl_vector_alloc(op->od->nObs);
		func_gsl_deriv_dx( gsl_opt_params, op, gsl_jacobian ); // Compute Jacobian using GSL function
		gsl_vector_free( op->od->obs_current_gsl );
		gsl_vector_free( op->pd->var_current_gsl ); */
		// func_gsl_dx( gsl_opt_params, op, gsl_jacobian ); // Compute Jacobian using forward difference
		func_dx( opt_params, f_x, op, jacobian ); // Compute Jacobian using forward difference
		for( k = i = 0; i < op->pd->nOptParam; i++ )
			for( j = 0; j < op->od->nTObs; j++, k++ )
				gsl_matrix_set( gsl_jacobian, j, i, jacobian[k] );
	}
	if( debug )
	{
		if( compute_jacobian ) tprintf( "\nJacobian matrix\n" ); // Print Jacobian
		else tprintf( "\nJacobian matrix (provided externally)\n" );
		tprintf( "%-25s :", "Observations" );
		for( k = 0; k < op->od->nTObs; k++ )
		{
			if( op->od->nTObs < 30 || ( k < 10 || k > op->od->nTObs - 10 ) ) tprintf( " %s", op->od->obs_id[k] );
			if( op->od->nTObs >= 30 && k == 11 ) tprintf( " ..." );
		}
		tprintf( "\n" );
		for( k = i = 0; i < op->pd->nOptParam; i++ )
		{
			tprintf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->od->nTObs; j++ )
			{
				eps = gsl_matrix_get( gsl_jacobian, j, i );
				if( fabs( eps ) > 1e3 ) sprintf( buf, "%6.0e", eps );
				else sprintf( buf, "%6.1f", eps );
				if( op->od->nTObs < 30 || ( j < 10 || j > op->od->nTObs - 10 ) ) tprintf( " %s", buf );
				if( op->od->nTObs >= 30 && j == 11 ) tprintf( " ..." );
			}
			tprintf( "\n" );
			/*
			tprintf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->od->nObs; j++, k++ )
			{
				eps = jacobian[k];
				if( fabs( eps ) > 1e3 ) sprintf( buf, " %6.0e", eps );
				else sprintf( buf, " %6.2f", eps );
				if( op->od->nObs < 30 || ( j < 10 || j > op->od->nObs - 10 ) ) tprintf( " %s", buf );
				if( op->od->nObs > 30 && j == 11 ) tprintf( " ..." );
			}
			tprintf( "\n" );
			 */
		}
	}
	sprintf( filename, "%s", op->root );
	if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
	if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
	strcat( filename, ".jacobian" );
	out = Fwrite( filename );
	fprintf( out, "%-25s :", "Parameters" );
	for( i = 0; i < op->pd->nOptParam; i++ )
		fprintf( out, " \"%s\"", op->pd->var_id[op->pd->var_index[i]] );
	fprintf( out, "\n" );
	for( j = 0; j < op->od->nTObs; j++ )
	{
		fprintf( out, "%-25s :", op->od->obs_id[j] );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fprintf( out, " %g", gsl_matrix_get( gsl_jacobian, j, i ) );
		fprintf( out, "\n" );
	}
	fclose( out );
	tprintf( "Jacobian matrix stored (%s)\n", filename );
	if( compute_covar ) // Standalone eigen analysis
	{
		// gsl_matrix_set_zero( gsl_covar );
		ier = gsl_multifit_covar( gsl_jacobian, 0.0, gsl_covar );
		if( ier != GSL_SUCCESS ) { tprintf( "Problem computing covariance matrix!\n" ); ier = 1; }
		else ier = 0;
	}
	if( debug )
	{
		if( compute_covar ) tprintf( "\nCovariance matrix\n" ); // Print Jacobian
		else tprintf( "\nCovariance matrix (provided externally)\n" );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			tprintf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->pd->nOptParam; j++ )
				tprintf( " %7.0e", gsl_matrix_get( gsl_covar, i, j ) );
			tprintf( "\n" );
		}
	}
	sprintf( filename, "%s", op->root );
	if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
	if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
	strcat( filename, ".covariance" );
	out = Fwrite( filename );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
		for( j = 0; j < op->pd->nOptParam; j++ )
			fprintf( out, " %g", gsl_matrix_get( gsl_covar, i, j ) );
		fprintf( out, "\n" );
	}
	fclose( out );
	tprintf( "Covariance matrix stored (%s)\n", filename );
	// Compute A optimality
	aopt = 0;
	for( i = 0; i < op->pd->nOptParam; i++ )
		aopt += gsl_matrix_get( gsl_covar, i, i );
	ier = 0;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		stddev[i] = sqrt( gsl_matrix_get( gsl_covar, i, i ) ); // compute standard deviations before the covariance matrix is destroyed by eigen functions
		if( stddev[i] < DBL_EPSILON ) ier = 1;
	}
	if( ier == 0 )
	{
		if( debug )
		{
			tprintf( "\nCorrelation matrix\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				tprintf( "%-25s : ", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					tprintf( " %6.3f", gsl_matrix_get( gsl_covar, i, j ) / ( stddev[i] * stddev[j] ) );
				tprintf( "\n" );
			}
		}
		sprintf( filename, "%s", op->root );
		if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
		if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
		strcat( filename, ".correlation" );
		out = Fwrite( filename );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->pd->nOptParam; j++ )
				fprintf( out, " %g", gsl_matrix_get( gsl_covar, i, j ) / ( stddev[i] * stddev[j] ) );
			fprintf( out, "\n" );
		}
		fclose( out );
		tprintf( "Correlation matrix stored (%s)\n", filename );
		// GSL_COVAR is destroyed during eigen computation
		gsl_eigen_symmv( gsl_covar, eigenval, eigenvec, eigenwork );
		if( debug )
		{
			tprintf( "\nEigenvectors (sorted by absolute values of eigenvalues)\n" );
			gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_ABS_ASC );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				tprintf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					tprintf( " %6.3f", gsl_matrix_get( eigenvec, i, j ) );
				tprintf( "\n" );
			}
			tprintf( "%-25s :", "Eigenvalues" );
			for( j = 0; j < op->pd->nOptParam; j++ )
				tprintf( " %6.0e", gsl_vector_get( eigenval, j ) );
			tprintf( "\n" );
			tprintf( "\nEigenvectors (sorted by eigenvalues)\n" );
			gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_VAL_ASC );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				tprintf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					tprintf( " %6.3f", gsl_matrix_get( eigenvec, i, j ) );
				tprintf( "\n" );
			}
			tprintf( "%-25s :", "Eigenvalues" );
			for( j = 0; j < op->pd->nOptParam; j++ )
				tprintf( " %6.0e", gsl_vector_get( eigenval, j ) );
			tprintf( "\n" );
		}
		else gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_VAL_ASC );
		sprintf( filename, "%s", op->root );
		if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
		if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
		strcat( filename, ".eigen" );
		out = Fwrite( filename );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = op->pd->nOptParam - 1; j >= 0; j-- )
				fprintf( out, " %g", gsl_matrix_get( eigenvec, i, j ) );
			fprintf( out, "\n" );
		}
		fprintf( out, "%-25s :", "Eigenvalues" );
		for( j = 0; j < op->pd->nOptParam; j++ )
			fprintf( out, " %g", gsl_vector_get( eigenval, j ) );
		fprintf( out, "\n" );
		fclose( out );
		tprintf( "Eigen vactors and eigen values stored (%s)\n", filename );
		// compute performance metrics
		copt = fabs( gsl_vector_get( eigenval, op->pd->nOptParam - 1 ) ) / fabs( gsl_vector_get( eigenval, 0 ) );
		eopt = fabs( gsl_vector_get( eigenval, op->pd->nOptParam - 1 ) );
		dopt = 1;
		for( j = op->pd->nOptParam - 1; j >= 0; j-- ) // reverse ordering to enhance accuracy small to big
			dopt *= fabs( gsl_vector_get( eigenval, j ) );
	}
	else
		tprintf( "Correlation matrix and eigen vectors cannot be computed!\n" );
	dof = op->od->nCObs - op->pd->nOptParam;
	if( dof > 0 )
	{
		stddev_scale = sqrt( phi / dof );
		for( i = 0; i < op->pd->nOptParam; i++ )
			stddev[i] *= stddev_scale;
		gf = phi / dof;
		ln_det_weight = 0;
		for( i = 0; i < op->od->nTObs; i++ )
			if( op->od->obs_weight[i] > 0 )
				ln_det_weight += log( op->od->obs_weight[i] );
		ln_det_v = ln_det_weight + op->pd->nOptParam * log( gf );
		sml = ( double ) dof + ln_det_v + op->od->nCObs * 1.837877;
		aic = sml + ( double ) 2 * op->pd->nOptParam;
		bic = sml + ( double ) op->pd->nOptParam * log( op->od->nCObs );
		cic = sml + ( double ) 2 * op->pd->nOptParam * log( log( op->od->nCObs ) );
		kic = sml + ( double ) op->pd->nOptParam * log( op->od->nCObs * 0.159154943 ) - log( dopt );
	}
	if( op->cd->problem_type == EIGEN || debug )
	{
		tprintf( "\nNumber of parameters           : %d\n", op->pd->nOptParam );
		tprintf( "Number of observations         : %d\n", op->od->nCObs );
		tprintf( "Number of degrees of freedom   : %d\n", dof );
		tprintf( "Objective function             : %g\n", phi );
		if( dof > 0 ) tprintf( "Posterior measurement variance : %g\n", gf );
		tprintf( "\nOptimality metrics based on covariance matrix of observation errors:\n" );
		tprintf( "A-optimality (matrix trace)               : %g\n", aopt );
		tprintf( "C-optimality (matrix conditioning number) : %g\n", copt );
		tprintf( "E-optimality (matrix maximum eigenvalue)  : %g\n", eopt );
		tprintf( "D-optimality (matrix determinant)         : %g\n", dopt );
		tprintf( "\nDeterminant of covariance matrix of observation errors : %-15g ( ln(det S) = %g )\n", dopt, log( dopt ) );
		if( dof > 0 )
		{
			tprintf( "Determinant of observation weight matrix               : %-15g ( ln(det W) = %g )\n", exp( ln_det_weight ) , ln_det_weight );
			tprintf( "Determinant of covariance matrix of measurement errors : %-15g ( ln(det V) = %g )\n", exp( ln_det_v ), ln_det_v );
			tprintf( "\nLog likelihood function             : %g\n", -sml / 2 );
			tprintf( "Maximum likelihood                  : %g\n", sml );
			tprintf( "AIC (Akaike information criterion)  : %g\n", aic );
			tprintf( "BIC                                 : %g\n", bic );
			tprintf( "CIC                                 : %g\n", cic );
			tprintf( "KIC (Kashyap Information Criterion) : %g\n", kic );
		}
	}
	if( dof <= 0 ) tt = 1;
	else if( dof <= 30 )  tt = student_dist[dof - 1]; // Critical values of Student's t distribution for 0.975
	else if( dof <= 40 )  tt = student_dist[30] + ( dof - 30 ) * ( student_dist[31] - student_dist[30] ) / 10;
	else if( dof <= 60 )  tt = student_dist[31] + ( dof - 40 ) * ( student_dist[32] - student_dist[31] ) / 20;
	else if( dof <= 120 ) tt = student_dist[32] + ( dof - 60 ) * ( student_dist[33] - student_dist[32] ) / 60;
	else tt = student_dist[33];
	if( dof > 0 )
	{
		tprintf( "\nObtained fit is " );
		if( gf > 200 ) tprintf( "not very good (chi^2/dof = %g > 200)\n", gf );
		else tprintf( "relatively good (chi^2/dof = %g < 200)\n", gf );
	}
	tprintf( "Optimized model parameters:\n" );
	if( debug && op->cd->sintrans == 1 ) tprintf( "Transformed space (applied during optimization):\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->cd->sintrans == 1 ) opt_params[i] = asin( sin( opt_params[i] ) );
		x_u[i] = opt_params[i] + ( double ) tt * stddev[i];
		x_d[i] = opt_params[i] - ( double ) tt * stddev[i];
		status = 0;
		if( op->cd->sintrans == 1 )
		{
			if( x_d[i] < -M_PI / 2 ) { status = 1; x_d[i] = -M_PI / 2; }
			if( x_u[i] > M_PI / 2 ) { status = 1; x_u[i] = M_PI / 2; }
		}
		else
		{
			if( x_d[i] < op->pd->var_min[k] ) { status = 1; x_d[i] = op->pd->var_min[k]; }
			if( x_u[i] > op->pd->var_max[k] ) { status = 1; x_u[i] = op->pd->var_max[k]; }
		}
		if( debug )
		{
			if( dof > 0 )
			{
				tprintf( "%-40s : %12g stddev %12g (%12g - %12g)", op->pd->var_id[k], opt_params[i], stddev[i], x_d[i], x_u[i] );
				if( status ) tprintf( " Uncertainty ranges constrained by prior bounds\n" );
				else tprintf( "\n" );
			}
			else
				tprintf( "%-40s : %12g -- Uncertainty ranges cannot be estimated\n", op->pd->var_id[k], opt_params[i] );
		}
	}
	DeTransform( x_u, op, x_u );
	DeTransform( x_d, op, x_d );
	if( op->cd->sintrans == 1 )
	{
		if( debug ) tprintf( "Untransformed space:\n" );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			tprintf( "%-40s : ", op->pd->var_id[k] );
			if( dof > 0 )
			{
				if( op->pd->var_log[k] == 0 ) tprintf( "%12g stddev %12g (%12g - %12g)", op->pd->var[k], stddev[i], x_d[i], x_u[i] );
				else tprintf( "%12g stddev %12g (%12g - %12g)", pow( 10, op->pd->var[k] ), stddev[i], pow( 10, x_d[i] ), pow( 10, x_u[i] ) );
				status = 0;
				if( x_d[i] <= op->pd->var_min[k] ) { status = 1; x_d[i] = op->pd->var_min[k]; }
				if( x_u[i] >= op->pd->var_max[k] ) { status = 1; x_u[i] = op->pd->var_max[k]; }
				if( status ) tprintf( " Uncertainty ranges constrained by prior bounds\n" );
				else tprintf( "\n" );
			}
			else
			{
				if( op->pd->var_log[k] == 0 ) tprintf( "%12g", op->pd->var[k] );
				else tprintf( "%12g", pow( 10, op->pd->var[k] ) );
				tprintf( " -- Uncertainty ranges cannot be estimated\n" );
			}
		}
	}
	free( opt_params ); free( stddev ); free( x_u ); free( x_d );
	gsl_vector_free( gsl_opt_params );
	gsl_matrix_free( eigenvec ); gsl_vector_free( eigenval ); gsl_eigen_symmv_free( eigenwork );
	if( compute_center ) free( f_x );
	if( compute_jacobian ) { gsl_matrix_free( gsl_jacobian ); free( jacobian ); }
	if( compute_covar ) gsl_matrix_free( gsl_covar );
	return( 1 );
}

// Check consistency of input and output files
int check( struct opt_data *op )
{
	struct opt_data *p = ( struct opt_data * )op;
	int i, bad_data = 0;
	for( i = 0; i < p->ed->ntpl; i++ )
		if( par_tpl( p->pd->nParam, p->pd->var_id, p->cd->var, p->ed->fn_tpl[i], p->ed->fn_out[i], p->cd->tpldebug + 1 ) == -1 )
			bad_data = 1;
	for( i = 0; i < p->od->nTObs; i++ ) p->od->obs_current[i] = p->od->res[i] = 0;
	for( i = 0; i < p->ed->nins; i++ )
		if( ins_obs( p->od->nTObs, p->od->obs_id, p->od->obs_current, p->od->res, p->ed->fn_ins[i], p->ed->fn_obs[i], p->cd->insdebug + 1 ) == -1 )
			bad_data = 1;
	for( i = 0; i < p->od->nTObs; i++ )
	{
		if( p->od->res[i] < 0 )
		{
			tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
	}
	if( bad_data ) return( 0 );
	return( 1 );
}

// Initial guesses -- random
int igrnd( struct opt_data *op )
{
	int i, k, m, q1, q2, npar, status, phi_global, success_global, count, debug_level, solution_found, no_memory = 0, neval_total, njac_total;
	int *eval_success, *eval_total;
	double c, phi_min, *orig_params, *opt_params,
		   *opt_params_min, *opt_params_max, *opt_params_avg,
		   *sel_params_min, *sel_params_max, *sel_params_avg,
		   *var_lhs, v;
	char filename[255], buf[255];
	int ( *optimize_func )( struct opt_data * op ); // function pointer to optimization function (LM or PSO)
	FILE *out, *out2;
	char ESC = 27; // Escape
	strcpy( op->label, "igrnd" );
	if( ( orig_params = ( double * ) malloc( op->pd->nParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( no_memory ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( op->cd->nreal > 1 )
	{
		if( ( opt_params_min = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
		if( ( opt_params_max = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
		if( ( opt_params_avg = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
		if( no_memory ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			opt_params_min[i] = HUGE_VAL;
			opt_params_max[i] = opt_params_avg[i] = 0;
		}
		if( op->cd->phi_cutoff > DBL_EPSILON || op->cd->check_success )
		{
			if( ( sel_params_min = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
			if( ( sel_params_max = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
			if( ( sel_params_avg = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
			if( no_memory ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				sel_params_min[i] = HUGE_VAL;
				sel_params_max[i] = opt_params_avg[i] = 0;
			}
		}
	}
	tprintf( "\nSEQUENTIAL RUNS using random initial guesses for model parameters (realizations = %d):\n", op->cd->nreal );
	if( op->pd->nFlgParam != 0 ) { tprintf( "Only flagged parameters are randomized\n" ); npar = op->pd->nFlgParam; }
	else if( op->pd->nOptParam != 0 ) { tprintf( "No flagged parameters; all optimizable parameters are randomized\n" ); npar = op->pd->nOptParam; }
	else { tprintf( "No flagged or optimizable parameters; all parameters are randomized\n" ); npar = op->pd->nParam; }
	if( ( var_lhs = ( double * ) malloc( npar * op->cd->nreal * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( ( eval_success = ( int * ) malloc( op->cd->nreal * sizeof( int ) ) ) == NULL ) no_memory = 1;
	if( ( eval_total = ( int * ) malloc( op->cd->nreal * sizeof( int ) ) ) == NULL ) no_memory = 1;
	if( no_memory )
	{
		tprintf( "Not enough memory!\n" );
		return( 0 );
	}
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	tprintf( "Random sampling (variables %d; realizations %d) using ", npar, op->cd->nreal );
	sampling( npar, op->cd->nreal, &op->cd->seed, var_lhs, op, 1 );
	tprintf( "done.\n" );
	if( op->cd->mdebug )
	{
		sprintf( filename, "%s.igrnd_set", op->root );
		out = Fwrite( filename );
		for( count = 0; count < op->cd->nreal; count ++ )
		{
			for( k = 0; k < npar; k++ )
				fprintf( out, "%.15g ", var_lhs[k + count * npar] );
			fprintf( out, "\n" );
		}
		fclose( out );
		tprintf( "Random sampling set saved in %s.igrnd_set\n", op->root );
		sprintf( filename, "%s.igrnd_param", op->root );
		out = Fwrite( filename );
		for( count = 0; count < op->cd->nreal; count ++ )
		{
			for( k = i = 0; i < op->pd->nParam; i++ )
				if( op->pd->var_opt[i] == 2 || ( op->pd->var_opt[i] == 1 && op->pd->nFlgParam == 0 ) )
				{
					v = var_lhs[k + count * npar] * op->pd->var_range[i] + op->pd->var_min[i];
					if( op->pd->var_log[i] == 0 ) fprintf( out, "%.15g ", v );
					else fprintf( out, "%.15g ", pow( 10, v ) );
					k++;
				}
			fprintf( out, "\n" );
		}
		fclose( out );
		tprintf( "Randomly sampled parameters saved in %s.mcrnd_param\n", op->root );
	}
	for( i = 0; i < op->pd->nParam; i++ )
		orig_params[i] = op->pd->var[i]; // Save original initial values for all parameters
	if( strncasecmp( op->cd->opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
	else optimize_func = optimize_pso; // Define optimization method: PSO
	// File management
	sprintf( filename, "%s.igrnd.zip", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.igrnd.zip %s.igrnd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	sprintf( buf, "zip -m %s.igrnd.zip %s.igrnd-[0-9]*.* >& /dev/null", op->root, op->root ); system( buf );
	sprintf( buf, "mv %s.igrnd.zip %s.igrnd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
	sprintf( filename, "%s.igrnd.results", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.igrnd_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	out = Fwrite( filename );
	sprintf( filename, "%s.igrnd-opt=%s_eval=%d_real=%d", op->root, op->cd->opt_method, op->cd->maxeval, op->cd->nreal );
	out2 = Fwrite( filename );
	if( op->pd->nOptParam == 0 )
		tprintf( "WARNING: No parameters to optimize! Forward runs performed instead (ie Monte Carlo analysis)\n" );
	phi_min = HUGE_VAL;
	solution_found = 0;
	phi_global = success_global = neval_total = njac_total = 0;
	if( op->cd->ireal != 0 ) k = op->cd->ireal - 1; // applied if execution of a specific realization is requested (ncase)
	else k = 0;
	if( op->cd->debug || op->cd->mdebug ) op->cd->lmstandalone = 1;
	else op->cd->lmstandalone = 0;
	for( count = k; count < op->cd->nreal; count++ )
	{
		op->cd->neval = op->cd->njac = 0;
		fprintf( out, "%d : init var", count + 1 );
		tprintf( "\nRandom set #%d: ", count + 1 );
		if( op->cd->mdebug || op->cd->nreal == 1 ) tprintf( "\n" );
		op->counter = count + 1;
		for( k = i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] == 2 || ( op->pd->var_opt[i] == 1 && op->pd->nFlgParam == 0 ) )
			{
				op->pd->var[i] = var_lhs[k + count * npar] * op->pd->var_range[i] + op->pd->var_min[i];
				if( op->pd->var_log[i] )
				{
					if( op->cd->mdebug || op->cd->nreal == 1 ) tprintf( "%s %.15g\n", op->pd->var_id[i], pow( 10, op->pd->var[i] ) );
					fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
				}
				else
				{
					if( op->cd->mdebug || op->cd->nreal == 1 ) tprintf( "%s %.15g\n", op->pd->var_id[i], op->pd->var[i] );
					fprintf( out, " %.15g", op->pd->var[i] );
				}
				k++;
			}
			else op->pd->var[i] = orig_params[i];
		if( op->pd->nOptParam > 0 )
		{
			if( op->cd->pargen )
			{
				if( op->cd->solution_type[0] != TEST && op->cd->solution_type[0] != EXTERNAL )
				{
					sprintf( filename, "%s-igrnd.%d.mads", op->root, count + 1 );
					op->cd->calib_type = SIMPLE;
					save_problem( filename, op );
					op->cd->calib_type = IGRND;
					continue;
				}
			}
			else
			{
				status = optimize_func( op ); // Optimize
				if( status == 0 ) return( 0 );
			}
		}
		else
		{
			for( i = 0; i < op->pd->nParam; i++ )
				op->pd->var[i] = var_lhs[i + count * npar] * op->pd->var_range[i] + op->pd->var_min[i];
			if( op->cd->mdebug ) { tprintf( "Forward run ... \n" ); debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			func_global( op->pd->var, op, op->od->res ); // op->pd->var is a dummy variable because op->pd->nOptParam == 0
			if( op->cd->mdebug ) op->cd->fdebug = debug_level;
		}
		if( op->cd->debug > 1 )
			tprintf( "\n" );
		else
		{
			tprintf( "Evaluations: %d ", op->cd->neval );
			if( op->cd->njac > 0 ) tprintf( "Jacobians: %d ", op->cd->njac );
			tprintf( "Objective function: %g Success: %d", op->phi, op->success );
		}
		if( op->phi < phi_min && ( ( op->cd->check_success && op->success ) || !op->cd->check_success ) )
		{
			solution_found = 1;
			phi_min = op->phi;
			for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]];
			for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
		}
		if( op->cd->mdebug || op->cd->nreal == 1 )
		{
			tprintf( "\n" );
			print_results( op, 0 );
		}
		neval_total += eval_total[count] = op->cd->neval;
		njac_total += op->cd->njac;
		if( op->cd->check_success && op->success )
		{
			eval_success[success_global] = op->cd->neval;
			success_global++;
		}
		else if( op->phi < op->cd->phi_cutoff )
		{
			eval_success[phi_global] = op->cd->neval;
			phi_global++;
		}
		if( op->cd->nreal > 1 )
		{
			for( i = 0; i < op->pd->nOptParam; i++ ) // Posterior parameter statistics for all simulations
			{
				c = op->pd->var[op->pd->var_index[i]];
				if( op->pd->var_log[op->pd->var_index[i]] ) c = pow( 10, c );
				if( c < opt_params_min[i] ) opt_params_min[i] = c;
				if( c > opt_params_max[i] ) opt_params_max[i] = c;
				opt_params_avg[i] += c;
				if( ( op->cd->check_success && op->success ) || op->phi < op->cd->phi_cutoff )
				{
					if( c < sel_params_min[i] ) sel_params_min[i] = c;
					if( c > sel_params_max[i] ) sel_params_max[i] = c;
					sel_params_avg[i] += c;
				}
			}
		}
		if( op->cd->pdebug || op->cd->ldebug ) tprintf( "\n" ); // extra new line if the optimization process is debugged
		fprintf( out2, "%g %d %d\n", op->phi, op->success, op->cd->neval );
		fflush( out2 );
		fprintf( out, " : OF %g success %d : final var", op->phi, op->success );
		for( i = 0; i < op->pd->nOptParam; i++ ) // Print only optimized parameters (including flagged); ignore fixed parameters
		{
			k = op->pd->var_index[i];
			if( op->pd->var_log[k] ) fprintf( out, " %.15g", pow( 10, op->pd->var[k] ) );
			else fprintf( out, " %.15g", op->pd->var[k] );
		}
		fprintf( out, "\n" );
		fflush( out );
		if( op->success && op->cd->nreal > 1 && op->cd->save )
			save_final_results( "igrnd", op, op->gd );
		if( op->f_ofe != NULL ) { fclose( op->f_ofe ); op->f_ofe = NULL; }
		if( op->cd->ireal != 0 ) break;
	}
	op->counter = 0;
	free( var_lhs );
	op->cd->neval = neval_total; // provide the correct number of total evaluations
	op->cd->njac = njac_total; // provide the correct number of total evaluations
	tprintf( "\nTotal number of evaluations = %d\n", neval_total );
	tprintf( "Total number of jacobians = %d\n", njac_total );
	if( !solution_found )
		tprintf( "WARNING: No IGRND solution has been identified matching required criteria!\n" );
	op->phi = phi_min; // get the best phi
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]] = op->pd->var_current[i] = op->pd->var_best[i]; // get the best estimate
	for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_current[i] = op->od->obs_best[i] ; // get the best observations
	fprintf( out, "Minimum objective function: %g\n", phi_min );
	tprintf( "Minimum objective function: %g\n", phi_min );
	tprintf( "Repeat the run producing the best results ...\n" );
	if( op->cd->debug ) { debug_level = op->cd->fdebug; op->cd->fdebug = 1; }
	Transform( opt_params, op, opt_params );
	func_global( opt_params, op, op->od->res );
	if( op->cd->debug ) op->cd->fdebug = debug_level;
	if( op->cd->nreal > 1 )
	{
		if( op->cd->phi_cutoff > DBL_EPSILON )
		{
			if( phi_global == 0 ) tprintf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", op->cd->nreal, op->cd->phi_cutoff );
			else tprintf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
		}
		if( op->cd->obsrange > DBL_EPSILON || op->cd->obserror > DBL_EPSILON )
		{
			if( success_global == 0 ) tprintf( "None of the %d sequential calibration runs produced successful calibration ranges!\n", op->cd->nreal );
			else tprintf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		}
		if( op->cd->parerror > DBL_EPSILON )
		{
			if( success_global == 0 ) tprintf( "None of the %d sequential calibration runs produced acceptable model parameters!\n", op->cd->nreal );
			else tprintf( "Number of the sequential calibration runs producing acceptable model parameters = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		}
		qsort( eval_success, success_global, sizeof( int ), sort_int );
		qsort( eval_total, op->cd->nreal, sizeof( int ), sort_int );
		q1 = ( int )( ( double ) success_global / 4 - 0.25 ) ;
		m = ( int )( ( double ) success_global / 2 - 0.5 );
		q2 = ( int )( ( double ) success_global * 3 / 4 - 0.25 );
		if( success_global > 0 )
		{
			tprintf( "Statistics of successful number of evaluations : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_success[0], eval_success[q1], ESC, eval_success[m], ESC, eval_success[q2], eval_success[success_global - 1], success_global );
			fprintf( out2, "Statistics of successful number of evaluations : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_success[0], eval_success[q1], ESC, eval_success[m], ESC, eval_success[q2], eval_success[success_global - 1], success_global );
		}
		q1 = ( int )( ( double ) op->cd->nreal / 4 - 0.25 );
		m = ( int )( ( double ) op->cd->nreal / 2 - 0.5 );
		q2 = ( int )( ( double ) op->cd->nreal * 3 / 4 - 0.25 );
		tprintf( "Statistics of total number of evaluations      : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_total[0], eval_total[q1], ESC, eval_total[m], ESC, eval_total[q2], eval_total[op->cd->nreal - 1], op->cd->nreal );
		fprintf( out2, "Statistics of total number of evaluations      : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_total[0], eval_total[q1], ESC, eval_total[m], ESC, eval_total[q2], eval_total[op->cd->nreal - 1], op->cd->nreal );
		tprintf( "Statistics of all the model parameter estimates:\n" );
		fprintf( out2, "Statistics of all the model parameter estimates:\n" );
		for( i = 0; i < op->pd->nOptParam; i++ ) // Posterior parameter statistics for all simulations
		{
			tprintf( "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], opt_params_avg[i] / op->cd->nreal, opt_params_min[i], opt_params_max[i] );
			fprintf( out2, "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], opt_params_avg[i] / op->cd->nreal, opt_params_min[i], opt_params_max[i] );
		}
		if( success_global > 0 || phi_global > 0 )
		{
			tprintf( "Statistics of all the successful model parameter estimates:\n" );
			fprintf( out2, "Statistics of all the successful model parameter estimates:\n" );
			k = success_global + phi_global;
			for( i = 0; i < op->pd->nOptParam; i++ ) // Posterior parameter statistics for all simulations
			{
				tprintf( "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], sel_params_avg[i] / k, sel_params_min[i], sel_params_max[i] );
				fprintf( out2, "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], sel_params_avg[i] / k, sel_params_min[i], sel_params_max[i] );
			}
		}
	}
	fprintf( out, "Total number of evaluations = %d\n", neval_total );
	if( op->cd->nreal > 1 )
	{
		if( op->cd->phi_cutoff > DBL_EPSILON ) fprintf( out, "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
		if( op->cd->obsrange > DBL_EPSILON || op->cd->obserror > DBL_EPSILON )
			fprintf( out, "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		if( op->cd->parerror > DBL_EPSILON )
			fprintf( out, "Number of the sequential calibration runs producing acceptable model parameters = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
	}
	fprintf( out2, "OF min: %g\n", phi_min );
	fprintf( out2, "Total number of evaluations: %d\n", neval_total );
	fprintf( out2, "Success rate %g\n", ( double ) success_global / op->cd->nreal );
	fclose( out ); fclose( out2 );
	tprintf( "Results are saved in %s.igrnd.results and %s.igrnd-opt=%s_eval=%d_real=%d\n", op->root, op->root, op->cd->opt_method, op->cd->maxeval, op->cd->nreal );
	tprintf( "\nFinal results:\n" );
	print_results( op, 1 );
	free( opt_params ); free( orig_params ); free( eval_success ); free( eval_total );
	if( op->cd->nreal > 1 )
	{
		free( opt_params_min ); free( opt_params_max ); free( opt_params_avg );
		if( op->cd->phi_cutoff > DBL_EPSILON || op->cd->check_success ) { free( sel_params_min ); free( sel_params_max ); free( sel_params_avg ); }
	}
	save_final_results( "", op, op->gd );
	if( solution_found ) return( 1 );
	else return( 0 );
}

// Initial guesses -- distributed to the parameter space
int igpd( struct opt_data *op )
{
	int i, j, k, status, phi_global, success_global, count, debug_level, no_memory = 0, neval_total, njac_total;
	double phi_min, *orig_params, *opt_params;
	char filename[255], buf[255];
	int ( *optimize_func )( struct opt_data * op ); // function pointer to optimization function (LM or PSO)
	FILE *out, *out2;
	strcpy( op->label, "igpd" );
	if( ( orig_params = ( double * ) malloc( op->pd->nParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( no_memory ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nSEQUENTIAL CALIBRATIONS using discretized initial guesses for model parameters:\n" );
	if( op->pd->nFlgParam == 0 )
		tprintf( "WARNING: No flagged parameters! Discretization of the initial guesses cannot be performed! Forward run will be performed instead.\n" );
	if( op->pd->nOptParam == 0 )
		tprintf( "WARNING: No parameters to optimize! Forward run will be performed instead.\n" );
	for( i = 0; i < op->pd->nParam; i++ ) orig_params[i] = op->pd->var[i]; // Save original initial values for all parameters
	if( strncasecmp( op->cd->opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
	else optimize_func = optimize_pso; // Define optimization method: PSO
	// File management
	sprintf( filename, "%s.igpd.zip", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.igpd.zip %s.igpd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	sprintf( buf, "zip -m %s.igpd.zip %s.igpd-[0-9]*.* >& /dev/null", op->root, op->root ); system( buf );
	sprintf( buf, "mv %s.igpd.zip %s.igpd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
	sprintf( filename, "%s.igpd.results", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.igpd_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	out = Fwrite( filename );
	k = 1;
	for( i = 0; i < op->pd->nParam; i++ )
		if( op->pd->var_opt[i] == 2 )
		{
			j = ( int )( ( double )( op->pd->var_max[i] - op->pd->var_min[i] ) / op->pd->var_dx[i] + 2 );
			if( op->pd->var_dx[i] > ( op->pd->var_max[i] - op->pd->var_min[i] ) ) j++;
			k *= j;
		}
	tprintf( "Total number of sequential calibrations will be %i\n", k );
	op->cd->nreal = k;
	sprintf( filename, "%s.igpd-opt=%s_eval=%d_real=%d", op->root, op->cd->opt_method, op->cd->maxeval, op->cd->nreal );
	out2 = Fwrite( filename );
	for( i = 0; i < op->pd->nParam; i++ )
		if( op->pd->var_opt[i] == 2 )
			orig_params[i] = op->pd->var_min[i];
	phi_global = success_global = 0;
	phi_min = HUGE_VAL;
	count = neval_total = njac_total = 0;
	do
	{
		op->cd->neval = op->cd->njac = 0;
		count++;
		if( op->cd->ireal == 0 || op->cd->ireal == count )
		{
			fprintf( out, "%d : init var", count ); // counter
			tprintf( "SEQUENTIAL CALIBRATIONS #%d: ", count );
			op->counter = count;
			if( op->cd->debug == 0 ) tprintf( "\n" );
			for( i = 0; i < op->pd->nParam; i++ )
			{
				op->pd->var[i] = orig_params[i];
				if( op->pd->var_opt[i] == 2 ) // Print flagged parameters
				{
					tprintf( "%s %g\n", op->pd->var_id[i], orig_params[i] );
					if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, orig_params[i] ) );
					else fprintf( out, " %.15g", orig_params[i] );
				}
			}
			if( op->pd->nOptParam > 0 && op->pd->nFlgParam > 0 )
			{
				status = optimize_func( op ); // Optimize
				if( status == 0 ) { return( 0 ); }
			}
			else
			{
				if( op->cd->debug ) { tprintf( "Forward run ... \n" ); debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
				func_global( op->pd->var, op, op->od->res ); // op->pd->var is dummy because op->pd->nOptParam == 0
				if( op->cd->debug ) op->cd->fdebug = debug_level;
			}
			neval_total += op->cd->neval;
			njac_total += op->cd->njac;
			if( op->phi < op->cd->phi_cutoff ) phi_global++;
			if( op->cd->debug )
			{
				tprintf( "\n" );
				print_results( op, 0 );
			}
			else tprintf( "Objective function: %g Success: %d\n", op->phi, op->success );
			success_global += op->success;
			if( op->phi < phi_min )
			{
				phi_min = op->phi;
				for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]];
				for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
			}
			fprintf( out, " : OF %g success %d : final var ", op->phi, op->success );
			for( i = 0; i < op->pd->nParam; i++ )
				if( op->pd->var_opt[i] >= 1 ) // Print only optimized parameters (including flagged); ignore fixed parameters
				{
					if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
					else fprintf( out, " %.15g", op->pd->var[i] );
				}
			fprintf( out, "\n" );
			fflush( out );
			if( op->f_ofe != NULL ) { fclose( op->f_ofe ); op->f_ofe = NULL; }
			if( op->success && op->cd->nreal > 1 && op->cd->save )
				save_final_results( "igpd", op, op->gd );
			if( op->cd->ireal != 0 ) break;
		}
		if( op->pd->nFlgParam == 0 || op->pd->nOptParam == 0 ) break;
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] == 2 )
			{
				if( orig_params[i] < op->pd->var_max[i] )
				{
					orig_params[i] += op->pd->var_dx[i];
					if( orig_params[i] > op->pd->var_max[i] ) orig_params[i] = op->pd->var_max[i];
					break;
				}
				else orig_params[i] = op->pd->var_min[i];
			}
		if( i == op->pd->nParam ) break;
	}
	while( 1 );
	op->counter = 0;
	op->cd->neval = neval_total; // provide the correct number of total evaluations
	op->cd->njac = njac_total; // provide the correct number of total evaluations
	tprintf( "\nTotal number of evaluations = %d\n", neval_total );
	tprintf( "Total number of jacobians = %d\n", njac_total );
	op->phi = phi_min; // get the best phi
	for( i = 0; i < op->pd->nOptParam; i++ ) opt_params[i] = op->pd->var[op->pd->var_index[i]] = op->pd->var_current[i] = op->pd->var_best[i]; // get the best estimate
	for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_current[i] = op->od->obs_best[i] ; // get the best observations
	tprintf( "Minimum objective function: %g\n", phi_min );
	tprintf( "Repeat the run producing the best results ...\n" );
	if( op->cd->debug ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
	Transform( opt_params, op, opt_params );
	func_global( opt_params, op, op->od->res );
	if( op->cd->debug ) op->cd->fdebug = debug_level;
	print_results( op, 1 );
	fprintf( out, "Minimum objective function: %g\n", phi_min );
	if( op->cd->nreal > 1 )
	{
		if( success_global == 0 ) tprintf( "None of the %d sequential calibration runs produced predictions within calibration ranges!\n", op->cd->nreal );
		else tprintf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		if( op->cd->phi_cutoff > DBL_EPSILON )
		{
			if( phi_global == 0 ) tprintf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", op->cd->nreal, op->cd->phi_cutoff );
			else tprintf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
		}
	}
	fprintf( out, "Number of evaluations = %d\n", neval_total );
	if( op->cd->nreal > 1 )
	{
		fprintf( out, "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		if( op->cd->phi_cutoff > DBL_EPSILON ) fprintf( out, "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
	}
	fprintf( out2, "OF min %g\n", phi_min );
	fprintf( out2, "eval %d\n", neval_total );
	fprintf( out2, "success %d\n", success_global );
	fclose( out ); fclose( out2 );
	tprintf( "Results are saved in %s.igpd.results\n", op->root );
	free( opt_params ); free( orig_params );
	save_final_results( "", op, op->gd );
	return( 1 );
}

// Partial parameter space discretizations
int ppsd( struct opt_data *op )
{
	int i, j, k, status, phi_global, success_global, count, debug_level, no_memory = 0, neval_total, njac_total;
	double phi_min, *orig_params;
	int *orig_opt;
	char filename[255], buf[255];
	int ( *optimize_func )( struct opt_data * op ); // function pointer to optimization function (LM or PSO)
	FILE *out;
	strcpy( op->label, "ppsd" );
	op->cd->lmstandalone = 1;
	if( ( orig_params = ( double * ) malloc( op->pd->nParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( no_memory ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nSEQUENTIAL RUNS using partial parameter-space discretization (PPSD):\n" );
	if( op->pd->nFlgParam == 0 )
		tprintf( "WARNING: No flagged parameters! Discretization of the initial guesses cannot be performed!\n" );
	if( op->pd->nOptParam == 0 )
		tprintf( "WARNING: No parameters to optimize! Forward runs performed instead\n" );
	for( i = 0; i < op->pd->nParam; i++ ) orig_params[i] = op->pd->var[i]; // Save original initial values for all parameters
	if( strncasecmp( op->cd->opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
	else optimize_func = optimize_pso; // Define optimization method: PSO
	sprintf( filename, "%s.ppsd.zip", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.ppsd.zip %s.ppsd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	sprintf( buf, "zip -m %s.ppsd.zip %s.ppsd-[0-9]*.* >& /dev/null", op->root, op->root ); system( buf );
	sprintf( buf, "mv %s.ppsd.zip %s.ppsd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
	sprintf( filename, "%s.ppsd.results", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.ppsd_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	out = Fwrite( filename );
	k = 1;
	for( i = 0; i < op->pd->nParam; i++ )
		if( op->pd->var_opt[i] == 2 )
		{
			j = ( int )( ( double )( op->pd->var_max[i] - op->pd->var_min[i] ) / op->pd->var_dx[i] + 2 );
			if( op->pd->var_dx[i] > ( op->pd->var_max[i] - op->pd->var_min[i] ) ) j++;
			k *= j;
		}
	tprintf( "Total number of sequential runs will be %i\n", k );
	op->cd->nreal = k;
	for( i = 0; i < op->pd->nParam; i++ )
		if( op->pd->var_opt[i] == 2 ) orig_params[i] = op->cd->var[i] = op->pd->var_min[i];
	phi_min = HUGE_VAL;
	count = neval_total = njac_total = phi_global = success_global = 0;
	if( op->cd->pargen )
	{
		orig_opt = ( int * ) malloc( op->pd->nParam * sizeof( int ) );
		for( i = 0; i < op->pd->nParam; i++ )
		{
			orig_opt[i] = op->pd->var_opt[i];
			if( op->pd->var_opt[i] == 2 ) op->pd->var_opt[i] = 0;
		}
	}
	do
	{
		op->cd->neval = op->cd->njac = 0;
		count++;
		if( op->cd->ireal == 0 || op->cd->ireal == count )
		{
			fprintf( out, "%d : ", count );
			tprintf( "\nSEQUENTIAL RUN #%d:\n", count );
			op->counter = count;
			tprintf( "Discretized parameters:\n" );
			for( i = 0; i < op->pd->nParam; i++ )
			{
				op->cd->var[i] = op->pd->var[i] = orig_params[i]; // these are the true original parameters
				if( op->pd->var_opt[i] == 2 ) // Print only flagged parameters
				{
					tprintf( "%s %g\n", op->pd->var_id[i], op->cd->var[i] );
					fprintf( out, "%g ", op->cd->var[i] );
				}
			}
			if( op->pd->nOptParam > 0 )
			{
				if( op->cd->pargen )
				{
					if( op->cd->solution_type[0] != TEST && op->cd->solution_type[0] != EXTERNAL )
					{
						sprintf( filename, "%s-ppsd.%d.mads", op->root, count + 1 );
						op->cd->calib_type = SIMPLE;
						save_problem( filename, op );
						op->cd->calib_type = PPSD;
						for( i = 0; i < op->pd->nParam; i++ )
							if( orig_opt[i] == 2 )
							{
								tprintf( "%s %g\n", op->pd->var_id[i], op->cd->var[i] );
								fprintf( out, "%g ", op->cd->var[i] );
								if( orig_params[i] < op->pd->var_max[i] )
								{
									orig_params[i] += op->pd->var_dx[i];
									if( orig_params[i] > op->pd->var_max[i] ) orig_params[i] = op->pd->var_max[i];
									break;
								}
								else orig_params[i] = op->pd->var_min[i];
							}
						if( i == op->pd->nParam ) break;
						else continue;
					}
				}
				else
				{
					status = optimize_func( op ); // Optimize
					if( status == 0 ) return( 0 );
				}
			}
			else
			{
				if( op->cd->debug ) { tprintf( "Forward run ... \n" ); debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
				func_global( op->pd->var, op, op->od->res ); // op->pd->var is dummy because op->pd->nOptParam == 0
				if( op->cd->debug ) op->cd->fdebug = debug_level;
			}
			if( op->phi < op->cd->phi_cutoff ) phi_global++;
			neval_total += op->cd->neval;
			njac_total += op->cd->njac;
			if( op->cd->debug > 2 )
			{
				tprintf( "\n" );
				print_results( op, 0 );
			}
			else tprintf( "Objective function: %g Success: %d", op->phi, op->success );
			success_global += op->success;
			if( op->phi < phi_min )
			{
				phi_min = op->phi;
				for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]];
				for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
			}
			if( op->pd->nOptParam > 0 )
			{
				fprintf( out, " : OF %g Success %d : final var", op->phi, op->success );
				for( i = 0; i < op->pd->nParam; i++ )
					if( op->pd->var_opt[i] == 1 ) // Print only optimized parameters; ignore fixed and flagged parameters
					{
						if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
						else fprintf( out, " %.15g", op->pd->var[i] );
					}
				fprintf( out, "\n" );
			}
			else
				fprintf( out, " : OF %g Success %d\n", op->phi, op->success );
			fflush( out );
			if( op->f_ofe != NULL ) { fclose( op->f_ofe ); op->f_ofe = NULL; }
			if( op->success && op->cd->save )
			{
				op->cd->calib_type = SIMPLE;
				sprintf( filename, "%s-ppsd.%d.mads", op->root, count + 1 );
				save_problem( filename, op );
				op->cd->calib_type = PPSD;
				save_final_results( "ppsd", op, op->gd );
			}
			if( op->cd->ireal != 0 ) break;
		}
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] == 2 )
			{
				if( orig_params[i] < op->pd->var_max[i] )
				{
					orig_params[i] += op->pd->var_dx[i];
					if( orig_params[i] > op->pd->var_max[i] ) orig_params[i] = op->pd->var_max[i];
					break;
				}
				else orig_params[i] = op->pd->var_min[i];
			}
		if( i == op->pd->nParam ) break;
	}
	while( 1 );
	fclose( out );
	if( op->cd->pargen )
	{
		for( i = 0; i < op->pd->nParam; i++ )
			op->pd->var_opt[i] = orig_opt[i];
		free( orig_opt );
	}
	op->cd->neval = neval_total; // provide the correct number of total evaluations
	op->cd->njac = njac_total; // provide the correct number of total evaluations
	tprintf( "\nTotal number of evaluations = %d\n", neval_total );
	tprintf( "Total number of jacobians = %d\n", njac_total );
	if( success_global == 0 ) tprintf( "None of the %d sequential calibration runs produced predictions within calibration ranges!\n", op->cd->nreal );
	else tprintf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
	if( op->cd->phi_cutoff > DBL_EPSILON )
	{
		if( phi_global == 0 ) tprintf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", op->cd->nreal, op->cd->phi_cutoff );
		else tprintf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
	}
	op->counter = 0;
	tprintf( "Results are saved in %s.ppsd.results\n", op->root );
	free( orig_params );
	return( 1 );
}

int montecarlo( struct opt_data *op )
{
	int i, j, k, npar, phi_global, success_global, success_all, count, debug_level, bad_data = 0;
	double phi_min, *opt_params, *var_lhs, v;
	char filename[255], buf[255];
	FILE *out;
	strcpy( op->label, "mcrnd" );
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	npar = op->pd->nOptParam;
	if( ( var_lhs = ( double * ) malloc( npar * op->cd->nreal * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nMonte Carlo analysis using latin-hyper cube sampling:\n" );
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	tprintf( "Random sampling (variables %d; realizations %d) using ", npar, op->cd->nreal );
	sampling( npar, op->cd->nreal, &op->cd->seed, var_lhs, op, 1 );
	tprintf( "done.\n" );
	if( op->cd->mdebug )
	{
		sprintf( filename, "%s.mcrnd_set", op->root );
		out = Fwrite( filename );
		for( count = 0; count < op->cd->nreal; count ++ )
		{
			for( k = 0; k < npar; k++ )
				fprintf( out, "%.15g ", var_lhs[k + count * npar] );
			fprintf( out, "\n" );
		}
		fclose( out );
		tprintf( "Random sampling set saved in %s.mcrnd_set\n", op->root );
		sprintf( filename, "%s.mcrnd_param", op->root );
		out = Fwrite( filename );
		for( count = 0; count < op->cd->nreal; count ++ )
		{
			for( i = 0; i < npar; i++ )
			{
				k = op->pd->var_index[i];
				v = var_lhs[i + count * npar] * op->pd->var_range[k] + op->pd->var_min[k];
				if( op->pd->var_log[k] == 0 ) fprintf( out, "%.15g ", v );
				else fprintf( out, "%.15g ", pow( 10, v ) );
			}
			fprintf( out, "\n" );
		}
		fclose( out );
		tprintf( "Randomly sampled parameters saved in %s.mcrnd_param\n", op->root );
	}
	sprintf( filename, "%s.mcrnd.zip", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.mcrnd.zip %s.mcrnd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	sprintf( buf, "zip -m %s.mcrnd.zip %s.mcrnd-[0-9]*.* >& /dev/null", op->root, op->root ); system( buf );
	sprintf( buf, "mv %s.mcrnd.zip %s.mcrnd_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
	sprintf( filename, "%s.mcrnd.results", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.mcrnd_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	out = Fwrite( filename );
	phi_global = success_global = success_all = 0;
	phi_min = HUGE_VAL;
	if( op->cd->ireal != 0 ) k = op->cd->ireal - 1;
	else k = 0;
	if( op->cd->solution_type[0] == EXTERNAL && op->cd->num_proc > 1 && k == 0 ) // Parallel job
	{
		if( op->cd->debug || op->cd->mdebug ) tprintf( "Parallel execution of external jobs ...\n" );
		if( op->cd->debug || op->cd->mdebug ) tprintf( "Generation of all the model input files ...\n" );
		for( count = 0; count < op->cd->nreal; count ++ ) // Write all the files
		{
			fprintf( out, "%d : ", count + 1 ); // counter
			if( op->cd->mdebug ) tprintf( "\n" );
			tprintf( "Random set #%d: ", count + 1 );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				k = op->pd->var_index[i];
				opt_params[i] = op->pd->var[k] = var_lhs[i + count * npar] * op->pd->var_range[k] + op->pd->var_min[k];
			}
			if( op->cd->mdebug ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			Transform( opt_params, op, opt_params );
			func_extrn_write( count + 1, opt_params, op );
			tprintf( "external model input file(s) generated ...\n" );
			if( op->cd->mdebug ) op->cd->fdebug = debug_level;
			if( op->cd->mdebug )
			{
				tprintf( "\nRandom parameter values:\n" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					if( op->pd->var_log[op->pd->var_index[i]] == 0 ) tprintf( "%s %g\n", op->pd->var_id[op->pd->var_index[i]], op->pd->var[op->pd->var_index[i]] );
					else tprintf( "%s %g\n", op->pd->var_id[op->pd->var_index[i]], pow( 10, op->pd->var[op->pd->var_index[i]] ) );
			}
			for( i = 0; i < op->pd->nParam; i++ )
				if( op->pd->var_opt[i] >= 1 )
				{
					if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
					else fprintf( out, " %.15g", op->pd->var[i] );
				}
			fprintf( out, "\n" );
		}
		fclose( out );
		if( op->cd->pardebug > 4 )
		{
			for( count = 0; count < op->cd->nreal; count ++ ) // Perform all the runs in serial model (for testing)
			{
				tprintf( "Execute model #%d ... ", count + 1 );
				func_extrn_exec_serial( count + 1, op );
				tprintf( "done!\n" );
			}
		}
		else if( mprun( op->cd->nreal, op ) < 0 ) // Perform all the runs in parallel
		{
			tprintf( "ERROR: there is a problem with the parallel execution!\n" );
			return( 0 );
		}
		out = Fwrite( filename ); // rewrite results file including the results
		for( count = 0; count < op->cd->nreal; count ++ ) // Read all the files
		{
			if( op->cd->debug || op->cd->mdebug ) tprintf( "Reading all the model output files ...\n" );
			op->counter = count + 1;
			fprintf( out, "%d : ", op->counter ); // counter
			for( i = 0; i < op->pd->nOptParam; i++ ) // re
			{
				k = op->pd->var_index[i];
				op->pd->var[k] = var_lhs[i + count * npar] * op->pd->var_range[k] + op->pd->var_min[k];
			}
			if( op->cd->mdebug ) tprintf( "\n" );
			tprintf( "Model results #%d: ", op->counter );
			bad_data = 0;
			bad_data = func_extrn_read( count + 1, op, op->od->res );
			if( bad_data ) return( 0 );
			if( ( op->cd->check_success && op->success ) || op->phi < op->cd->phi_cutoff ) success_all = 1;
			else success_all = 0;
			if( op->cd->mdebug > 1 ) { tprintf( "\n" ); print_results( op, 0 ); }
			else if( op->cd->mdebug )
			{
				tprintf( "Objective function: %g Success: %d\n", op->phi, success_all );
				if( success_all ) tprintf( "All the predictions are within calibration ranges!\n" );
				else tprintf( "At least one of the predictions is outside calibration ranges!\n" );
			}
			else
				tprintf( "Objective function: %g Success = %d\n", op->phi, success_all );
			if( op->phi < phi_min )
			{
				phi_min = op->phi;
				for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]];
				for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
			}
			if( success_all ) success_global++;
			if( op->phi < op->cd->phi_cutoff ) phi_global++;
			for( i = 0; i < op->pd->nParam; i++ )
				if( op->pd->var_opt[i] >= 1 )
				{
					if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
					else fprintf( out, " %.15g", op->pd->var[i] );
				}
			if( op->od->nTObs > 0 ) fprintf( out, " OF %g success %d\n", op->phi, success_all );
			else fprintf( out, "\n" );
			fflush( out );
			if( ( success_all || op->od->nTObs == 0 ) && op->cd->save )
				save_final_results( "mcrnd", op, op->gd );
		}
	}
	else // Serial job
	{
		for( count = k; count < op->cd->nreal; count ++ )
		{
			op->counter = count + 1;
			fprintf( out, "%d : ", count + 1 ); // counter
			if( op->cd->mdebug ) tprintf( "\n" );
			tprintf( "Random set #%d: ", count + 1 );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				j = op->pd->var_index[i];
				opt_params[i] = op->pd->var[j] = var_lhs[i + count * npar] * op->pd->var_range[j] + op->pd->var_min[j];
			}
			if( op->cd->mdebug ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			if( ( op->cd->check_success && op->success ) || op->phi < op->cd->phi_cutoff ) success_all = 1;
			else success_all = 0;
			if( op->cd->mdebug ) op->cd->fdebug = debug_level;
			if( op->cd->mdebug )
			{
				tprintf( "\nRandom parameter values:\n" );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					j = op->pd->var_index[i];
					if( op->pd->var_log[j] == 0 ) tprintf( "%s %g\n", op->pd->var_id[j], op->pd->var[j] );
					else tprintf( "%s %g\n", op->pd->var_id[j], pow( 10, op->pd->var[j] ) );
				}
			}
			if( op->cd->mdebug > 1 ) { tprintf( "\nPredicted calibration targets:\n" ); print_results( op, 1 ); }
			else if( op->cd->mdebug )
			{
				tprintf( "Objective function: %g Success: %d\n", op->phi, success_all );
				if( success_all ) tprintf( "All the predictions are within calibration ranges!\n" );
				else tprintf( "At least one of the predictions is outside calibration ranges!\n" );
			}
			else
				tprintf( "Objective function: %g Success = %d\n", op->phi, success_all );
			if( op->phi < phi_min )
			{
				phi_min = op->phi;
				for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]];
				for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
			}
			if( success_all ) success_global++;
			if( op->phi < op->cd->phi_cutoff ) phi_global++;
			for( i = 0; i < op->pd->nParam; i++ )
				if( op->pd->var_opt[i] >= 1 )
				{
					if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
					else fprintf( out, " %.15g", op->pd->var[i] );
				}
			if( op->od->nTObs > 0 ) fprintf( out, " OF %g success %d\n", op->phi, success_all );
			else fprintf( out, "\n" );
			fflush( out );
			if( ( success_all || op->od->nTObs == 0 ) && op->cd->save )
				save_final_results( "mcrnd", op, op->gd );
			if( op->f_ofe != NULL ) { fclose( op->f_ofe ); op->f_ofe = NULL; }
			if( op->cd->ireal != 0 ) break;
		}
	}
	op->counter = 0;
	free( var_lhs );
	fclose( out );
	op->phi = phi_min; // get the best phi
	for( i = 0; i < op->pd->nOptParam; i++ ) opt_params[i] = op->pd->var[op->pd->var_index[i]] = op->pd->var_current[i] = op->pd->var_best[i]; // get the best estimate
	for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_current[i] = op->od->obs_best[i] ; // get the best observations
	tprintf( "\nMinimum objective function: %g\n", phi_min );
	tprintf( "Repeat the run producing the best results ...\n" );
	if( op->cd->debug ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
	Transform( opt_params, op, opt_params );
	func_global( opt_params, op, op->od->res );
	if( op->cd->debug ) op->cd->fdebug = debug_level;
	print_results( op, 1 );
	tprintf( "Results are saved in %s.mcrnd.results\n", op->root );
	tprintf( "\nMinimum objective function: %g\n", phi_min );
	if( success_global == 0 ) tprintf( "None of the Monte-Carlo runs produced predictions within calibration ranges!\n" );
	else tprintf( "Number of Monte-Carlo runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
	if( op->cd->phi_cutoff > DBL_EPSILON )
	{
		if( phi_global == 0 ) tprintf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", op->cd->nreal, op->cd->phi_cutoff );
		else tprintf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
	}
	free( opt_params );
	save_final_results( "", op, op->gd );
	return( 1 );
}

void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op, int debug )
{
	if( debug ) tprintf( "%s\n", op->cd->smp_method );
	if( nreal == 1 || strncasecmp( op->cd->smp_method, "random", 6 ) == 0 )
	{
		if( debug )
			tprintf( "Pure random sampling method ... " );
		smp_random( npar, nreal, seed, var_lhs );
	}
	else if( ( nreal <= 500 && op->cd->smp_method[0] == 0 ) || strncasecmp( op->cd->smp_method, "idlhs", 5 ) == 0 )
	{
		if( debug )
		{
			tprintf( "Improved Distributed LHS method " );
			if( strncasecmp( op->cd->smp_method, "idlhs", 5 ) != 0 ) tprintf( "(number of realizations < 500) " );
			tprintf( "... " );
		}
		lhs_imp_dist( npar, nreal, 5, seed, var_lhs );
	}
	else if( ( nreal > 500 && op->cd->smp_method[0] == 0 ) || strncasecmp( op->cd->smp_method, "lhs", 3 ) == 0 )
	{
		if( debug )
		{
			tprintf( "Standard LHS method " );
			if( strncasecmp( op->cd->smp_method, "lhs", 3 ) != 0 ) tprintf( "(number of realizations > 500) " );
			tprintf( "... " );
		}
		lhs_random( npar, nreal, seed, var_lhs );
	}
}

void print_results( struct opt_data *op, int verbosity )
{
	int i, j, k, success, success_all, predict = 0;
	double c, err, min, max, dx;
	success_all = 1;
	if( verbosity > 0 ) tprintf( "Optimized model parameters:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->pd->var_log[k] == 0 ) op->cd->var[k] = op->pd->var[k];
		else op->cd->var[k] = pow( 10, op->pd->var[k] );
		tprintf( "%s %g\n", op->pd->var_id[k], op->cd->var[k] );
	}
	if( verbosity > 0 && op->pd->nExpParam > 0 ) tprintf( "Tied model parameters:\n" );
#ifdef MATHEVAL
	for( i = 0; i < op->pd->nExpParam; i++ )
	{
		k = op->pd->param_expressions_index[i];
		tprintf( "%s = ", op->pd->var_id[k] );
		tprintf( "%s", evaluator_get_string( op->pd->param_expressions[i] ) );
		if( op->cd->solution_type[0] == EXTERNAL ) op->pd->var[k] = evaluator_evaluate( op->pd->param_expressions[i], op->pd->nParam, op->pd->var_id, op->cd->var );
		else op->pd->var[k] = evaluator_evaluate( op->pd->param_expressions[i], op->pd->nParam, op->pd->var_id_short, op->cd->var );
		tprintf( " = %g\n", op->pd->var[k] );
	}
#else
	if( op->pd->nExpParam > 0 )
	{
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		exit( 0 );
	}
#endif
	if( verbosity == 0 ) return;
	if( op->cd->solution_type[0] != TEST && op->od->nTObs > 0 )
	{
		tprintf( "\nModel calibration targets:\n" );
		if( op->cd->solution_type[0] == EXTERNAL )
		{
			for( i = 0; i < op->od->nObs; i++ )
			{
				if( op->od->obs_weight[i] == 0 ) { predict = 1; if( op->od->nCObs > 50 && i == 21 ) tprintf( "...\n" ); continue; }
				c = op->od->obs_current[i];
				err = op->od->obs_target[i] - c;
				// err = op->od->res[i];
				min = op->od->obs_min[i];
				max = op->od->obs_max[i];
				if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success_all = 0; success = 0; }
				else success = 1;
				if( op->od->nCObs < 50 || ( i < 20 || i > op->od->nCObs - 20 ) )
					tprintf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, min, max );
				if( op->od->nCObs > 50 && i == 21 ) tprintf( "...\n" );
			}
			if( op->rd->nRegul > 0 ) tprintf( "Model regularization terms:\n" );
			for( i = op->od->nObs; i < op->od->nTObs; i++ )
			{
				c = op->od->obs_current[i];
				err = op->od->obs_target[i] - c;
				// err = op->od->res[i];
				min = op->od->obs_min[i];
				max = op->od->obs_max[i];
				if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success_all = 0; success = 0; }
				else success = 1;
				tprintf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, min, max );
			}
		}
		else
		{
			for( k = 0, i = 0; i < op->wd->nW; i++ )
				for( j = 0; j < op->wd->nWellObs[i]; j++ )
				{
					if( op->wd->obs_weight[i][j] == 0 ) { predict = 1; continue; }
					c = op->od->obs_current[k];
					// err = op->wd->obs_target[i][j] - c;
					err = op->od->res[k];
					k++;
					min = op->wd->obs_min[i][j];
					max = op->wd->obs_max[i][j];
					if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success_all = 0; success = 0; }
					else success = 1;
					tprintf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], success, min, max );
				}
			if( op->rd->nRegul > 0 )  tprintf( "Model regularization terms:\n" );
			for( i = op->od->nObs; i < op->od->nTObs; i++ )
			{
				c = op->od->obs_current[i];
				// err = op->od->obs_target[i] - c;
				err = op->od->res[i];
				min = op->od->obs_min[i];
				max = op->od->obs_max[i];
				if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success_all = 0; success = 0; }
				else success = 1;
				tprintf( "%-17s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, min, max );
			}
		}
	}
	else
		for( i = 0; i < op->pd->nOptParam; i++ )
			if( fabs( op->pd->var[i] - op->pd->var_truth[i] ) > op->cd->parerror ) success_all = 0;
	if( op->phi < op->cd->phi_cutoff ) success_all = 1;
	tprintf( "Objective function: %g Success: %d (%d)\n", op->phi, op->success, success_all );
	if( op->cd->check_success > 0 && op->cd->obserror < 0 && op->cd->parerror < 0 )
	{
		if( op->success ) tprintf( "SUCCESS: All the model predictions are within calibration ranges!\n" );
		else tprintf( "At least one of the model predictions is outside calibration ranges!\n" );
	}
	if( op->cd->check_success > 0 && op->cd->obserror > 0 )
	{
		if( op->success ) tprintf( "SUCCESS: All the model predictions are within a predefined absolute error %g!\n", op->cd->obserror );
		else tprintf( "At least one of the model predictions has an absolute error greater than %g!\n", op->cd->obserror );
	}
	if( op->cd->check_success > 0 && op->cd->parerror > 0 )
	{
		if( op->success ) tprintf( "SUCCESS: All the estimated model parameters have an absolute error from the true parameters less than %g!\n", op->cd->parerror );
		else tprintf( "At least one of the estimated model parameters has an absolute error from the true parameters greater than %g!\n", op->cd->parerror );
	}
	if( op->cd->phi_cutoff > DBL_EPSILON && op->phi < op->cd->phi_cutoff )
		tprintf( "SUCCESS: Objective function is below the predefined cutoff value (%g < %g)!\n", op->phi, op->cd->phi_cutoff );
	if( op->cd->solution_type[0] != TEST && op->od->nTObs > 0 && predict )
	{
		tprintf( "\nModel predictions for not calibration targets:\n" );
		if( op->cd->solution_type[0] == EXTERNAL )
		{
			predict = op->od->nTObs - op->od->nCObs;
			j = 0;
			for( i = 0; i < op->od->nTObs; i++ )
			{
				if( op->od->obs_weight[i] != 0 ) { if( predict > 50 && j == 21 ) tprintf( "...\n" ); continue; }
				c = op->od->obs_current[i];
				err = op->od->obs_target[i] - c;
				// err = op->od->res[i];
				min = op->od->obs_min[i];
				max = op->od->obs_max[i];
				if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) success = 0;
				else success = 1;
				if( predict < 50 || ( j < 20 || j > predict - 20 ) )
					tprintf( "%-20s:%12g - %12g = %12g success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, success, min, max );
				if( predict > 50 && j == 21 ) tprintf( "...\n" );
				j++;
			}
		}
		else
		{
			for( i = 0; i < op->wd->nW; i++ )
				for( j = 0; j < op->wd->nWellObs[i]; j++ )
				{
					if( op->wd->obs_weight[i][j] != 0 ) continue;
					c = func_solver( op->wd->x[i], op->wd->y[i], op->wd->z1[i], op->wd->z2[i], op->wd->obs_time[i][j], op->cd );
					err = op->wd->obs_target[i][j] - c;
					min = op->wd->obs_min[i][j];
					max = op->wd->obs_max[i][j];
					if( op->cd->objfunc_type != SSR )
					{
						if( op->cd->objfunc_type == SSDA )
						{
							err = sqrt( fabs( err ) );
							if( c < op->wd->obs_target[i][j] ) err *= -1;
						}
						else err = 0; // SSD0 & SSDX
						if( op->cd->objfunc_type == SSDX ) { dx = max - min; if( op->cd->obsdomain > DBL_EPSILON && op->cd->obsdomain < dx ) dx = op->cd->obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
						if( c < min ) err += min - c;
						else if( c > max ) err += c - max;
						// tprintf( "%g %g %g %g\n", err, c, min - c, c - max );
						if( op->cd->objfunc_type == SSDX ) { min = op->wd->obs_min[i][j]; max = op->wd->obs_max[i][j]; }
					}
					if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) success = 0;
					else success = 1;
					tprintf( "%-10s(%5g):%12g - %12g = %12g success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, success, min, max );
				}
		}
	}
}

void save_final_results( char *label, struct opt_data *op, struct grid_data *gd )
{
	FILE *out, *out2;
	int i, j, k, success, success_all;
	double c, err, min, max, dx;
	char filename[255], filename2[255], fileroot[255];
	success_all = 1;
	// Generate general filename
	i = strlen( op->root ); // check for previous version number in the root name
	k = -1;
	if( label[0] == 0 && op->root[i - 4] == '-' && op->root[i - 3] == 'v' )
	{
		sscanf( &op->root[i - 2], "%d", &k );
		tprintf( "Current MADS analysis version: %d\n", k );
		strncpy( filename2, op->root, i - 4 );
		filename2[i - 4] = 0;
	}
	strcpy( filename, op->root );
	if( label[0] != 0 ) sprintf( filename, "%s.%s", filename, label );
	if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
	else if( op->cd->resultscase > 1 ) sprintf( filename, "%s.%d", filename, op->cd->resultscase );
	strcpy( fileroot, filename ); // Save filename root
	// Save MADS rerun file
	if( k == -1 ) sprintf( filename, "%s-rerun.mads", filename );
	else sprintf( filename, "%s-v%02d.mads", filename2, k + 1 ); // create new version mads file
	if( op->cd->solution_type[0] != TEST ) save_problem( filename, op );
	// Save results file
	strcpy( filename, fileroot );
	strcat( filename, ".results" );
	out = Fwrite( filename );
	fprintf( out, "Model parameters:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->pd->var_log[k] == 0 ) op->cd->var[k] = op->pd->var[k];
		else op->cd->var[k] = pow( 10, op->pd->var[k] );
		fprintf( out, "%s %g\n", op->pd->var_id[k], op->cd->var[k] );
	}
	if( op->pd->nExpParam > 0 ) fprintf( out, "Tied model parameters:\n" );
#ifdef MATHEVAL
	for( i = 0; i < op->pd->nExpParam; i++ )
	{
		k = op->pd->param_expressions_index[i];
		fprintf( out, "%s = ", op->pd->var_id[k] );
		fprintf( out, "%s", evaluator_get_string( op->pd->param_expressions[i] ) );
		if( op->cd->solution_type[0] == EXTERNAL ) op->pd->var[k] = evaluator_evaluate( op->pd->param_expressions[i], op->pd->nParam, op->pd->var_id, op->cd->var );
		else op->pd->var[k] = evaluator_evaluate( op->pd->param_expressions[i], op->pd->nParam, op->pd->var_id_short, op->cd->var );
		fprintf( out, " = %g\n", op->pd->var[k] );
	}
#else
	if( op->pd->nExpParam > 0 )
	{
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		fprintf( out, "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		exit( 0 );
	}
#endif
	if( op->cd->solution_type[0] != TEST && op->od->nTObs > 0 )
	{
		fprintf( out, "\nModel predictions:\n" );
		// Save residuals file
		strcpy( filename, fileroot );
		strcat( filename, ".residuals" );
		out2 = Fwrite( filename );
		if( op->cd->solution_type[0] == EXTERNAL )
			for( i = 0; i < op->od->nTObs; i++ )
			{
				c = op->od->obs_current[i];
				err = op->od->obs_target[i] - c;
				// err = op->od->res[i];
				min = op->od->obs_min[i];
				max = op->od->obs_max[i];
				if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { if( op->od->obs_weight[i] != 0 ) success_all = 0; success = 0; }
				else success = 1;
				fprintf( out, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, min, max );
				fprintf( out2, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, min, max );
			}
		else
		{
			for( k = 0, i = 0; i < op->wd->nW; i++ )
				for( j = 0; j < op->wd->nWellObs[i]; j++ )
				{
					min = op->wd->obs_min[i][j];
					max = op->wd->obs_max[i][j];
					if( op->wd->obs_weight[i][j] != 0 )
					{
						c = op->od->obs_current[k++];
						if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success_all = 0; success = 0; }
						else success = 1;
					}
					else
					{
						c = func_solver( op->wd->x[i], op->wd->y[i], op->wd->z1[i], op->wd->z2[i], op->wd->obs_time[i][j], op->cd );
						if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) success = 0;
						else success = 1;
					}
					err = op->wd->obs_target[i][j] - c;
					if( op->cd->objfunc_type != SSR )
					{
						if( op->cd->objfunc_type == SSDA )
						{
							err = sqrt( fabs( err ) );
							if( c < op->wd->obs_target[i][j] ) err *= -1;
						}
						else err = 0; // SSD0 & SSDX
						if( op->cd->objfunc_type == SSDX ) { dx = max - min; if( op->cd->obsdomain > DBL_EPSILON && op->cd->obsdomain < dx ) dx = op->cd->obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
						if( c < min ) err += min - c;
						else if( c > max ) err += c - max;
						if( op->cd->objfunc_type == SSDX ) { min = op->wd->obs_min[i][j]; max = op->wd->obs_max[i][j]; }
					}
					fprintf( out, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], success, min, max );
					fprintf( out2, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], success, min, max );
				}
		}
		fclose( out2 );
	}
	else
	{
		if( op->cd->check_success > 0 && op->cd->parerror > 0 )
		{
			for( i = 0; i < op->pd->nOptParam; i++ )
				if( fabs( op->pd->var[i] - op->pd->var_truth[i] ) > op->cd->parerror ) success_all = 0;
		}
	}
	if( op->phi < op->cd->phi_cutoff ) success_all = 1;
	fprintf( out, "Objective function: %g Success: %d (%d)\n", op->phi, op->success, success_all );
	if( op->cd->check_success > 0 && op->cd->obserror < 0 && op->cd->parerror < 0 )
	{
		if( op->success ) fprintf( out, "SUCCESS: All the model predictions are within calibration ranges!\n" );
		else fprintf( out, "At least one of the model predictions is outside calibration ranges!\n" );
	}
	if( op->cd->check_success > 0 && op->cd->obserror > 0 )
	{
		if( op->success ) fprintf( out, "SUCCESS: All the model predictions are within a predefined absolute error %g!\n", op->cd->obserror );
		else fprintf( out, "At least one of the model predictions has an absolute error greater than %g!\n", op->cd->obserror );
	}
	if( op->cd->check_success > 0 && op->cd->parerror > 0 )
	{
		if( op->success ) fprintf( out, "SUCCESS: All the estimated model parameters have an absolute error from the true parameters less than %g!\n", op->cd->parerror );
		else fprintf( out, "At least one of the estimated model parameters has an absolute error from the true parameters greater than %g!\n", op->cd->parerror );
	}
	if( op->cd->phi_cutoff > DBL_EPSILON && op->phi < op->cd->phi_cutoff )
		fprintf( out, "SUCCESS: Objective function is below the predefined cutoff value (%g < %g)!\n", op->phi, op->cd->phi_cutoff );
	fprintf( out, "Number of function evaluations = %d\n", op->cd->neval );
	if( op->cd->seed > 0 ) fprintf( out, "Seed = %d\n", op->cd->seed_init );
	fclose( out );
	// Save breakthrough files
	if( gd->min_t > 0 && op->cd->solution_type[0] != TEST )
	{
		tprintf( "\nCompute breakthrough curves at all the wells ..." );
		sprintf( filename, "%s.btc", fileroot );
		sprintf( filename2, "%s.btc-peak", fileroot );
		compute_btc2( filename, filename2, op );
		//			compute_btc( filename, &op, &gd );
	}
	// Save grid files (VTK)
	if( gd->time > 0 && op->cd->solution_type[0] != TEST )
	{
		tprintf( "\nCompute spatial distribution of predictions at t = %g ...\n", gd->time );
		sprintf( filename, "%s.vtk", fileroot );
		compute_grid( filename, op->cd, gd );
	}
}

// Modified from Numerical Recipes in C: The Art of Scientific Computing (ISBN 0-521-43108-5)
// corrected three-pass algorithm to minimize roundoff error in variance
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var )
{
	int j;
	double dev2[n];
	// First pass to calculate mean
	//	s = 0.0;
	//	for( j = 0; j < n; j++ ) s += data[j];
	//	*ave = s/n;
	// Second pass to calculate absolute deviations
	for( j = 0; j < n; j++ )
		dev2[j] = ( data[j] - ave ) * ( datb[j] - ave );
	// Sort devs
	gsl_sort( dev2, 1, n );
	// Third pass to calculate first (absolute) and second moments
	*var = 0.0;
	for( j = 0; j < n; j++ )
		*var += dev2[j];
	*var = ( *var - ep * ep / n ) / ( n - 1 );
}

void ave_sorted( double data[], int n, double *ave, double *ep )
{
	int j;
	double s, dev[n];
	// First pass to calculate mean
	s = 0.0;
	for( j = 0; j < n; j++ ) s += data[j];
	*ave = s / n;
	// Second pass to calculate absolute deviations
	for( j = 0; j < n; j++ )
		dev[j] = data[j] - *ave;
	// Sort devs
	gsl_sort( dev, 1, n );
	// Third pass to calculate first (absolute) moment
	*ep = 0.0;
	for( j = 0; j < n; j++ )
		*ep += dev[j];
}

int sort_int( const void *x, const void *y )
{
	return ( *( int * ) x - * ( int * ) y );
}

double sort_double( const void *x, const void *y )
{
	return ( *( double * ) x - * ( double * ) y );
}
