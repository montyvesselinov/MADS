// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
// Brianeisha Eure
// Leif
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

#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
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
#include "levmar-2.5/levmar.h"
#include "mads.h"

#define FIT(i) gsl_vector_get(solver->x, i)
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
int optimize_lm( struct opt_data *op ); // LM (Levenberg-Marquardt) optimization
int optimize_pso( struct opt_data *op ); // PSO optimization
int eigen( struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar ); // Eigen analysis
void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op, int debug ); // Random sampling
void print_results( struct opt_data *op, int verbosity ); // Print final results
void save_results( char *filename, struct opt_data *op, struct grid_data *gd ); // Save final results
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var );
void ave_sorted( double data[], int n, double *ave, double *ep );
char *timestamp(); // create time stamp
char *datestamp(); // create date stamp
int sort_int( const void *x, const void *y );
int sort_double( const void *a, const void *b );

/* Functions elsewhere */
// Model analyses
int pssa( struct opt_data *op );
int postpua( struct opt_data *op );
int glue( struct opt_data *op );
int infogap( struct opt_data *op );
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
// IO
void mads_info();
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
int count_lines( char *filename );
int count_cols( char *filename, int row );
// Pesting
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
// Random sampling
double epsilon();
void lhs_imp_dist( int nvar, int npoint, int d, int *seed, double x[] );
void lhs_center( int nvar, int npoint, int *seed, double x[] );
void lhs_edge( int nvar, int npoint, int *seed, double x[] );
void lhs_random( int nvar, int npoint, int *seed, double x[] );
void smp_random( int nvar, int npoint, int *seed, double x[] );
int get_seed( );
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

int main( int argn, char *argv[] )
{
	// TODO return status of the function calls is not always checked; needs to be checked
	int i, j, k, ier, npar, status, success, phi_global, success_global, success_all, count, debug_level, predict = 0, compare, bad_data = 0, no_memory = 0;
	int *eval_success, *eval_total;
	unsigned long neval_total, njac_total;
	double c, err, phi, phi_min, *orig_params, *opt_params,
		   *var_lhs, *var_a_lhs, *var_b_lhs;
	struct calc_data cd;
	struct param_data pd;
	struct obs_data od;
	struct obs_data preds;
	struct well_data wd;
	struct extrn_data ed;
	struct grid_data gd;
	struct opt_data op;
	struct gsens_data gs;
	char filename[255], root[255], extension[255], buf[255], *dot, *cwd;
	int ( *optimize_func )( struct opt_data * op ); // function pointer to optimization function (LM or PSO)
	char *host, *nodelist, *hostlist, *proclist, *lsblist, *beowlist; // parallel variables
	FILE *in, *out, *out2;
	time_t time_start, time_end, time_elapsed;
	pid_t pid;
	struct tm *ptr_ts;
	time_start = time( NULL );
	op.datetime_stamp = datestamp();
	op.pd = &pd;
	op.od = &od;
	op.preds = &preds;
	op.wd = &wd;
	op.cd = &cd;
	op.gd = &gd;
	op.ed = &ed;
	cd.neval = 0;
	cd.njac = 0;
	cd.nlmo = 0;
	cd.standalone = 1; // LM variable; LM is stand-alone if not part of tribes optimization
	cd.compute_phi = 1; // function calls compute OF (phi); turned off only when the jacobians are computed for LM
	cd.pderiv = cd.oderiv = -1; // internal flags; do not compute parameter and observation derivatives
	op.phi = HUGE_VAL;
	op.success = 0;
	op.f_ofe = NULL;
	printf( "MADS: Model Analyses & Decision Support (v1.1) 2011\n" );
	printf( "---------------------------------------------------\n" );
	if( argn < 2 )
	{
		mads_info(); // print short mads help manual
		exit( 1 );
	}
	else
		printf( "Velimir Vesselinov (monty) vvv@lanl.gov\nhttp://mads.lanl.gov -:- http://www.ees.lanl.gov/staff/monty/codes/mads\n\n" );
	if( cd.debug ) printf( "Argument[1]: %s\n", argv[1] );
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
	if( cd.debug ) printf( "Input file name: %s\n", filename );
	cd.time_infile = Fdatetime_t( filename, 0 );
	cd.datetime_infile = Fdatetime( filename, 0 );
	printf( "Problem root name: %s", root );
	op.root = root;
	op.counter = 0;
	if( cd.debug && extension[0] != 0 ) printf( " Extension: %s", extension );
	printf( "\n" );
	sprintf( buf, "%s.running", op.root ); // File named root.running is used to prevent simultaneous execution of multiple problems
	if( Ftest( buf ) == 0 ) // If file already exists quit ...
	{
		printf( "WARNING: Potentially another MADS run is currently performed for problem \'%s\' since file %s exists!\n\n", op.root, buf );
		//		printf( "ERROR: Potentially another MADS run is currently performed for problem \'%s\' since file %s exists!\n", op.root, buf );
		//		printf( "Delete %s to execute (sorry for the inconvenience)!\n", buf );
		//		exit( 0 );
	}
	sprintf( buf, "touch %s.running", op.root ); system( buf ); // Create a file named root.running to prevent simultaneous execution of multiple problems
	/*
	 *  Read input data
	 */
	if( strcasecmp( extension, "pst" ) == 0 ) // PEST Problem
	{
		printf( "PEST problem:\n" );
		if( ( ier = load_pst( filename, &op ) ) <= 0 )
		{
			printf( "MADS quits! Data input problem!\nExecute \'mads\' without any arguments to check the acceptable command-line keywords and options.\n" );
			if( ier == 0 )
			{
				sprintf( filename, "%s-error.mads", op.root );
				save_problem( filename, &op );
				printf( "\nMADS problem file named %s-error.mads is created to debug.\n", op.root );
			}
			sprintf( buf, "rm -f %s.running", op.root ); system( buf ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			exit( 0 );
		}
		if( cd.opt_method[0] == 0 ) { strcpy( cd.opt_method, "lm" ); cd.calib_type = SIMPLE; cd.problem_type = CALIBRATE; }
		cd.solution_type = EXTERNAL; func_global = func_extrn;
		buf[0] = 0;
		for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
		if( parse_cmd( buf, &cd ) == -1 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	else // MADS Problem
	{
		if( ( ier = load_problem( filename, argn, argv, &op ) ) <= 0 )
		{
			printf( "MADS quits! Data input problem!\nExecute \'mads\' without any arguments to check the acceptable command-line keywords and options.\n" );
			if( ier == 0 )
			{
				sprintf( filename, "%s-error.mads", op.root );
				save_problem( filename, &op );
				printf( "\nMADS problem file named %s-error.mads is created to debug.\n", op.root );
			}
			sprintf( buf, "rm -f %s.running", op.root ); system( buf ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			exit( 0 );
		}
		if( cd.solution_type == EXTERNAL ) func_global = func_extrn;
		else func_global = func_intrn;
	}
	/*
	 *  Check for parallel environment
	 */
	cd.paral_hosts = NULL;
	hostlist = NULL;
	if( ( nodelist = getenv( "NODELIST" ) ) != NULL )
	{
		if( cd.debug ) printf( "\nParallel environment is detected (environmental variable NODELIST is defined)\n" );
		if( cd.debug ) printf( "Node list %s\n", nodelist );
		hostlist = nodelist;
	}
	if( ( beowlist = getenv( "BEOWULF_JOB_MAP" ) ) != NULL )
	{
		if( cd.debug ) printf( "\nParallel environment is detected (environmental variable BEOWULF_JOB_MAP is defined)\n" );
		if( cd.debug ) printf( "Node list %s\n", beowlist );
		hostlist = beowlist;
	}
	if( ( lsblist = getenv( "LSB_HOSTS" ) ) != NULL )
	{
		if( cd.debug ) printf( "\nParallel environment is detected (environmental variable LSB_HOSTS is defined)\n" );
		if( cd.debug ) printf( "Node list %s\n", lsblist );
		hostlist = lsblist;
		if( ( proclist = getenv( "LSB_MCPU_HOSTS" ) ) != NULL && cd.debug ) printf( "LSB_MCPU_HOSTS Processors list %s\n", proclist );
	}
	if( hostlist != NULL )
	{
		if( cd.debug == 0 ) printf( "\nParallel environment is detected.\n" );
		if( ( host = getenv( "HOSTNAME" ) ) == NULL ) host = getenv( "HOST" );
		printf( "Host: %s\n", host );
		k = strlen( hostlist );
		i = count = 0;
		printf( "Nodes:" );
		while( i <= k )
		{
			sscanf( &hostlist[i], "%s", buf );
			printf( " \'%s\'", buf );
			i += strlen( buf ) + 1;
			count++;
		}
		printf( "\n" );
		printf( "Number of available nodes for parallel execution: %i\n", count );
		if( count <= 0 )
			printf( "ERROR: There is problem with the description of execution nodes!\n" );
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
		printf( "\nLocal parallel execution is requested using %d processors (np=%d)\n", cd.num_proc, cd.num_proc );
		cwd = getenv( "OSTYPE" ); printf( "OS type: %s\n", cwd );
		if( strncasecmp( cwd, "darwin", 6 ) == 0 )
			system( "\\rm -f num_proc >& /dev/null; ( sysctl hw.logicalcpu | cut -d : -f 2 ) > num_proc" ); // MAC OS
		else
			system( "\\rm -f num_proc >& /dev/null; ( cat /proc/cpuinfo | grep processor | wc -l ) > num_proc" ); // LINUX
		in = Fread( "num_proc" );
		fscanf( in, "%d", &k );
		fclose( in );
		system( "\\rm -f num_proc >& /dev/null" );
		printf( "Number of local processors available for parallel execution: %i\n", k );
		if( k < cd.num_proc ) printf( "WARNING: Number of requested processors exceeds the available resources!\n" );
	}
	if( cd.num_proc > 1 ) // Parallel job
	{
		pid = getpid();
		if( cd.debug ) printf( "Parent ID [%d]\n", pid );
		cwd = getenv( "PWD" ); dot = strrchr( cwd, '/' );
		cd.mydir = &dot[1];
		if( cd.debug ) printf( "Working directory: %s (%s)\n", cwd, cd.mydir );
		cd.mydir_hosts = dir_hosts( &op, op.datetime_stamp ); // Directories for parallel execution have unique name based on the execution time
	}
	op.label = ( char * ) malloc( 10 * sizeof( char ) ); op.label[0] = 0;
	if( ( orig_params = ( double * ) malloc( pd.nParam * sizeof( double ) ) ) == NULL ) { printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	/*
	 *  Problem based on external model
	 */
	if( cd.solution_type == EXTERNAL ) // Check the files for external execution
	{
		if( cd.debug || cd.tpldebug || cd.insdebug ) printf( "Checking the template files for errors ...\n" );
		bad_data = 0;
		for( i = 0; i < pd.nParam; i++ ) orig_params[i] = ( double ) - 1;
		for( i = 0; i < ed.ntpl; i++ ) // Check template files ...
			if( check_par_tpl( pd.nParam, pd.var_id, orig_params, ed.fn_tpl[i], cd.tpldebug ) == -1 )
				bad_data = 1;
		for( i = 0; i < pd.nParam; i++ )
		{
			if( orig_params[i] < 0 )
			{
				printf( "ERROR: Model parameter \'%s\' is not represented in the template file(s)!\n", pd.var_id[i] );
				bad_data = 1;
			}
			else if( orig_params[i] > 1.5 )
				printf( "WARNING: Model parameter \'%s\' is represented more than once (%d) in the template file(s)!\n", pd.var_id[i], ( int ) orig_params[i] );
		}
		if( cd.debug || cd.tpldebug || cd.insdebug ) printf( "Checking the instruction files for errors ...\n" );
		for( i = 0; i < od.nObs; i++ ) od.obs_current[i] = ( double ) - 1;
		for( i = 0; i < ed.nins; i++ )
			if( check_ins_obs( od.nObs, od.obs_id, od.obs_current, ed.fn_ins[i], cd.insdebug ) == -1 ) // Check instruction files.
				bad_data = 1;
		for( i = 0; i < od.nObs; i++ )
		{
			if( od.obs_current[i] < 0 )
			{
				printf( "ERROR: Observation \'%s\' is not defined in the instruction file(s)!\n", od.obs_id[i] );
				bad_data = 1;
			}
			else if( od.obs_current[i] > 1.5 )
				printf( "WARNING: Observation \'%s\' is defined more than once (%d) in the instruction file(s)! Arithmetic average will be computed!\n", od.obs_id[i], ( int ) od.obs_current[i] );
		}
		if( bad_data )
		{
			sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			system( buf );
			exit( 0 );
		}
	}
	/*
	 *  Check for restart conditions
	 */
	printf( "\nExecution date & time stamp: %s\n", op.datetime_stamp ); // Stamp will be applied to name / rename various output files
	if( cd.solution_type == EXTERNAL && cd.num_proc > 1 )
	{
		if( cd.restart == 1 ) // Restart by default
		{
			strcpy( buf, filename ); // Temporarily preserve the input file name
			sprintf( filename, "%s.restart_%s.zip", op.root, cd.datetime_infile );
			strcpy( cd.restart_zip_file, filename );
			if( Ftest( cd.restart_zip_file ) != 0 ) { if( cd.pardebug ) printf( "ZIP file (%s) with restart information is not available.\n", cd.restart_zip_file ); cd.restart = 0; }
			else
			{
				time_elapsed = cd.time_infile - Fdatetime_t( cd.restart_zip_file, 0 ); // time_infile - time_zipfile ...
				if( time_elapsed >= 0 ) { if( cd.pardebug ) printf( "No restart: the zip file (%s) with restart information is older than the MADS input file (%s)\n(restart can be enforced using \'restart=-1\' or \'rstfile=%s\')\n", cd.restart_zip_file, buf, cd.restart_zip_file ); cd.restart = 0; } // No restart
				else cd.restart = 1; // Attempt restart
			}
			if( cd.restart )
				printf( "DEFAULT Restart: zip file %s is consistent with date/time stamp of the MADS input file\n(to ignore the zip file, either delete the zip file, or use keyword \'restart=0\' ... \n", filename );
		}
		else if( cd.restart == -1 ) // Forced restart
		{
			sprintf( filename, "%s.restart_%s.zip", op.root, cd.datetime_infile );
			if( cd.restart_zip_file[0] == 0 ) strcpy( cd.restart_zip_file, filename );
			if( Ftest( cd.restart_zip_file ) != 0 ) { printf( "Restart is requested but a zip file (%s) with restart information is not available.\n", cd.restart_zip_file ); cd.restart = 0; }
			else printf( "FORCED Restart: using zip file %s ...\n", cd.restart_zip_file );
		}
		if( cd.restart )
		{
			printf( "MADS  input  file \'%40s\' last modified on %s\n", buf, Fdatetime( buf, 0 ) );
			sprintf( filename, "%s.results", op.root ); if( Ftest( filename ) != 0 ) printf( "MADS results file \'%40s\' last modified on %s\n", filename, Fdatetime( filename, 0 ) );
			printf( "MADS restart file \'%40s\' last modified on %s\n", cd.restart_zip_file, Fdatetime( cd.restart_zip_file, 0 ) );
			printf( "ZIP file (%s) with restart information is unzipped ... \n", cd.restart_zip_file );
			sprintf( buf, "rm -fR ../%s* %s.restart_info; unzip -o -u -: %s >& /dev/null", cd.mydir_hosts, op.root, cd.restart_zip_file ); // the input file name was temporarily in buf; not any more ...
			system( buf );
			sprintf( filename, "%s.restart_info", op.root );
			in = Fread( filename );
			fgets( buf, 255, in );
			white_trim( buf );
			cd.mydir_hosts = dir_hosts( &op, buf ); // Directories for parallel execution have unique name based on the old execution time (when restart files were created)
			fclose( in );
			printf( "Date & time stamp of the previous run: %s\n", buf );
		}
		// Preserve the existing restart zip file
		if( Ftest( cd.restart_zip_file ) == 0 )
		{
			if( cd.pardebug ) printf( "Previous restart file (%s) exists!\n", cd.restart_zip_file );
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
	//
	// DONE with file reading and problem setup
	//
	// Model analyses are performed below based on provided inputs
	//
	// ------------------------ IGRND
	//
	if( cd.problem_type == CALIBRATE && cd.calib_type == IGRND ) /* Calibration analysis using random initial guessed */
	{
		status = igrnd( &op );
		if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	//
	// ------------------------ IGPD
	//
	if( cd.problem_type == CALIBRATE && cd.calib_type == IGPD ) /* Calibration analysis using discretized initial guesses */
	{
		strcpy( op.label, "igpd" );
		if( ( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		printf( "\nSEQUENTIAL CALIBRATIONS using discretized initial guesses for model parameters:\n" );
		if( pd.nFlgParam == 0 )
			printf( "WARNING: No flagged parameters! Discretization of the initial guesses cannot be performed! Forward run will be performed instead.\n" );
		if( pd.nOptParam == 0 )
			printf( "WARNING: No parameters to optimize! Forward run will be performed instead.\n" );
		for( i = 0; i < pd.nParam; i++ ) orig_params[i] = pd.var[i]; // Save original initial values for all parameters
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
		else optimize_func = optimize_pso; // Define optimization method: PSO
		// File management
		sprintf( filename, "%s.igpd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.igpd.zip %s.igpd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.igpd.zip %s.igpd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.igpd.zip %s.igpd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.igpd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.igpd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		k = 1;
		for( i = 0; i < pd.nParam; i++ )
			if( pd.var_opt[i] == 2 )
			{
				j = ( double )( pd.var_max[i] - pd.var_min[i] ) / pd.var_dx[i] + 1;
				if( pd.var_dx[i] > ( pd.var_max[i] - pd.var_min[i] ) ) j++;
				k *= j;
			}
		printf( "Total number of sequential calbrations will be %i\n", k );
		cd.nreal = k;
		sprintf( filename, "%s.igpd-opt=%s_eval=%d_real=%d", op.root, cd.opt_method, cd.maxeval, cd.nreal );
		out2 = Fwrite( filename );
		for( i = 0; i < pd.nParam; i++ )
			if( pd.var_opt[i] == 2 )
				orig_params[i] = pd.var_min[i];
		phi_global = success_global = 0;
		phi_min = HUGE_VAL;
		count = neval_total = 0;
		do
		{
			cd.neval = 0;
			count++;
			if( cd.ireal == 0 || cd.ireal == count )
			{
				fprintf( out, "%d : init var", count ); // counter
				printf( "SEQUENTIAL CALIBRATIONS #%d: ", count );
				op.counter = count;
				if( cd.debug == 0 ) printf( "\n" );
				for( i = 0; i < pd.nParam; i++ )
				{
					pd.var[i] = orig_params[i];
					if( pd.var_opt[i] == 2 ) // Print flagged parameters
					{
						printf( "%s %g\n", pd.var_id[i], orig_params[i] );
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, orig_params[i] ) );
						else fprintf( out, " %.15g", orig_params[i] );
					}
				}
				if( pd.nOptParam > 0 )
				{
					status = optimize_func( &op ); // Optimize
					if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
				}
				else
				{
					if( cd.debug ) { printf( "Forward run ... \n" ); debug_level = cd.fdebug; cd.fdebug = 3; }
					func_global( pd.var, &op, od.res ); // pd.var is dummy because pd.nOptParam == 0
					if( cd.debug ) cd.fdebug = debug_level;
				}
				neval_total += cd.neval;
				if( op.phi < op.cd->phi_cutoff ) phi_global++;
				if( cd.debug )
				{
					printf( "\n" );
					print_results( &op, 0 );
				}
				else printf( "Objective function: %g Success: %d\n", op.phi, op.success );
				success_global += op.success;
				if( op.phi < phi_min )
				{
					phi_min = op.phi;
					for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]];
					for( i = 0; i < od.nObs; i++ ) od.obs_best[i] = od.obs_current[i];
				}
				fprintf( out, " : OF %g success %d : final var ", op.phi, op.success );
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 ) // Print only optimized parameters (including flagged); ignore fixed parameters
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				fprintf( out, "\n" );
				fflush( out );
				if( op.f_ofe != NULL ) { fclose( op.f_ofe ); op.f_ofe = NULL; }
				if( op.success && cd.nreal > 1 && cd.odebug > 1 ) save_results( "igpd", &op, &gd );
				if( cd.ireal != 0 ) break;
			}
			if( pd.nFlgParam == 0 || pd.nOptParam == 0 ) break;
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] == 2 )
				{
					if( orig_params[i] < pd.var_max[i] )
					{
						orig_params[i] += pd.var_dx[i];
						if( orig_params[i] > pd.var_max[i] ) orig_params[i] = pd.var_max[i];
						break;
					}
					else orig_params[i] = pd.var_min[i];
				}
			if( i == pd.nParam ) break;
		}
		while( 1 );
		op.counter = 0;
		cd.neval = neval_total; // provide the correct number of total evaluations
		printf( "\nTotal number of evaluations = %lu\n", neval_total );
		op.phi = phi_min; // get the best phi
		for( i = 0; i < pd.nOptParam; i++ ) opt_params[i] = pd.var[pd.var_index[i]] = pd.var_current[i] = pd.var_best[i]; // get the best estimate
		for( i = 0; i < od.nObs; i++ ) od.obs_current[i] = od.obs_best[i] ; // get the best observations
		printf( "Minimum objective function: %g\n", phi_min );
		print_results( &op, 1 );
		if( cd.debug )
		{
			printf( "Repeat the run producing the best results ...\n" );
			debug_level = cd.fdebug; cd.fdebug = 3;
			Transform( opt_params, &op, opt_params );
			func_global( opt_params, &op, od.res );
			cd.fdebug = debug_level;
		}
		fprintf( out, "Minimum objective function: %g\n", phi_min );
		if( cd.nreal > 1 )
		{
			if( success_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions within calibration ranges!\n", cd.nreal );
			else printf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, cd.nreal, ( double ) success_global / cd.nreal );
			if( op.cd->phi_cutoff > DBL_EPSILON )
			{
				if( phi_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", cd.nreal, op.cd->phi_cutoff );
				else printf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op.cd->phi_cutoff, phi_global, cd.nreal, ( double ) phi_global / cd.nreal );
			}
		}
		fprintf( out, "Number of evaluations = %lu\n", neval_total );
		if( cd.nreal > 1 )
		{
			fprintf( out, "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, cd.nreal, ( double ) success_global / cd.nreal );
			if( op.cd->phi_cutoff > DBL_EPSILON ) fprintf( out, "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op.cd->phi_cutoff, phi_global, cd.nreal, ( double ) phi_global / cd.nreal );
		}
		fprintf( out2, "OF min %g\n", phi_min );
		fprintf( out2, "eval %lu\n", neval_total );
		fprintf( out2, "success %d\n", success_global );
		fclose( out ); fclose( out2 );
		printf( "Results are saved in %s.igpd.results\n", op.root );
		free( opt_params );
		save_results( "", &op, &gd );
	}
	//
	// ------------------------ PPSD
	//
	if( cd.problem_type == CALIBRATE && cd.calib_type == PPSD ) /* Calibration analysis using discretized parameters */
	{
		strcpy( op.label, "ppsd" );
		printf( "\nSEQUENTIAL RUNS using partial parameter-space discretization (PPSD):\n" );
		if( pd.nFlgParam == 0 )
			printf( "WARNING: No flagged parameters! Discretization of the initial guesses cannot be performed!\n" );
		if( pd.nOptParam == 0 )
			printf( "WARNING: No parameters to optimize! Forward runs performed instead\n" );
		for( i = 0; i < pd.nParam; i++ ) orig_params[i] = pd.var[i]; // Save original initial values for all parameters
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
		else optimize_func = optimize_pso; // Define optimization method: PSO
		sprintf( filename, "%s.igpd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.ppsd.zip %s.ppsd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.ppsd.zip %s.ppsd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.ppsd.zip %s.ppsd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.ppsd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.ppsd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		k = 1;
		for( i = 0; i < pd.nParam; i++ )
			if( pd.var_opt[i] == 2 )
			{
				j = ( double )( pd.var_max[i] - pd.var_min[i] ) / pd.var_dx[i] + 1;
				if( pd.var_dx[i] > ( pd.var_max[i] - pd.var_min[i] ) ) j++;
				k *= j;
			}
		printf( "Total number of sequential runs will be %i\n", k );
		cd.nreal = k;
		for( i = 0; i < pd.nParam; i++ )
			if( pd.var_opt[i] == 2 ) orig_params[i] = cd.var[i] = pd.var_min[i];
		phi_min = HUGE_VAL;
		count = neval_total = phi_global = success_global = 0;
		do
		{
			cd.neval = 0;
			count++;
			if( cd.ireal == 0 || cd.ireal == count )
			{
				fprintf( out, "%d : ", count );
				printf( "\nSEQUENTIAL RUN #%d:\n", count );
				op.counter = count;
				printf( "Discretized parameters:\n" );
				for( i = 0; i < pd.nParam; i++ )
				{
					cd.var[i] = pd.var[i] = orig_params[i]; // these are the true original parameters
					if( pd.var_opt[i] == 2 ) // Print only flagged parameters
					{
						printf( "%s %g\n", pd.var_id[i], cd.var[i] );
						fprintf( out, "%g ", cd.var[i] );
					}
				}
				if( pd.nOptParam > 0 )
				{
					status = optimize_func( &op ); // Optimize
					if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
				}
				else
				{
					if( cd.debug ) { printf( "Forward run ... \n" ); debug_level = cd.fdebug; cd.fdebug = 3; }
					func_global( pd.var, &op, od.res ); // pd.var is dummy because pd.nOptParam == 0
					if( cd.debug ) cd.fdebug = debug_level;
				}
				if( op.phi < op.cd->phi_cutoff ) phi_global++;
				neval_total += cd.neval;
				if( cd.debug > 2 )
				{
					printf( "\n" );
					print_results( &op, 0 );
				}
				else printf( "Objective function: %g Success: %d", op.phi, op.success );
				success_global += op.success;
				if( op.phi < phi_min )
				{
					phi_min = op.phi;
					for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]];
					for( i = 0; i < od.nObs; i++ ) od.obs_best[i] = od.obs_current[i];
				}
				if( pd.nOptParam > 0 )
				{
					fprintf( out, " : OF %g Success %d : final var", op.phi, op.success );
					for( i = 0; i < pd.nParam; i++ )
						if( pd.var_opt[i] == 1 ) // Print only optimized parameters; ignore fixed and flagged parameters
						{
							if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
							else fprintf( out, " %.15g", pd.var[i] );
						}
					fprintf( out, "\n" );
				}
				else
					fprintf( out, " : OF %g Success %d\n", op.phi, op.success );
				fflush( out );
				if( op.f_ofe != NULL ) { fclose( op.f_ofe ); op.f_ofe = NULL; }
				if( op.success && cd.odebug > 1 ) save_results( "ppsd", &op, &gd );
				if( cd.ireal != 0 ) break;
			}
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] == 2 )
				{
					if( orig_params[i] < pd.var_max[i] )
					{
						orig_params[i] += pd.var_dx[i];
						if( orig_params[i] > pd.var_max[i] ) orig_params[i] = pd.var_max[i];
						break;
					}
					else orig_params[i] = pd.var_min[i];
				}
			if( i == pd.nParam ) break;
		}
		while( 1 );
		fclose( out );
		cd.neval = neval_total; // provide the correct number of total evaluations
		printf( "\nTotal number of evaluations = %lu\n", neval_total );
		if( success_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions within calibration ranges!\n", cd.nreal );
		else printf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, cd.nreal, ( double ) success_global / cd.nreal );
		if( op.cd->phi_cutoff > DBL_EPSILON )
		{
			if( phi_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", cd.nreal, op.cd->phi_cutoff );
			else printf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op.cd->phi_cutoff, phi_global, cd.nreal, ( double ) phi_global / cd.nreal );
		}
		op.counter = 0;
		printf( "Results are saved in %s.ppsd.results\n", op.root );
	}
	//
	// ------------------------ MONTECARLO
	//
	if( cd.problem_type == MONTECARLO ) /* Monte Carlo analysis */
	{
		strcpy( op.label, "mcrnd" );
		if( ( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		npar = pd.nOptParam;
		if( ( var_lhs = ( double * ) malloc( npar * cd.nreal * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		printf( "\nMonte Carlo analysis using latin-hyper cube sampling:\n" );
		if( cd.seed < 0 ) { cd.seed *= -1; printf( "Imported seed: %d\n", cd.seed ); }
		else if( cd.seed == 0 ) { printf( "New " ); cd.seed_init = cd.seed = get_seed(); }
		else printf( "Current seed: %d\n", cd.seed );
		printf( "Random sampling (variables %d; realizations %d) using ", npar, cd.nreal );
		sampling( npar, cd.nreal, &cd.seed, var_lhs, &op, 1 );
		printf( "done.\n" );
		if( cd.mdebug )
		{
			sprintf( filename, "%s.mcrnd_set", op.root );
			out = Fwrite( filename );
			for( count = 0; count < cd.nreal; count ++ )
			{
				for( k = 0; k < npar; k++ )
					fprintf( out, "%.15g ", var_lhs[k + count * npar] );
				fprintf( out, "\n" );
			}
			fclose( out );
			printf( "Random sampling set saved in %s.mcrnd_set\n", op.root );
		}
		sprintf( filename, "%s.mcrnd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.mcrnd.zip %s.mcrnd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.mcrnd.zip %s.mcrnd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.mcrnd.zip %s.mcrnd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.mcrnd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.mcrnd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		phi_global = success_global = 0;
		phi_min = HUGE_VAL;
		if( cd.ireal != 0 ) k = cd.ireal - 1;
		else k = 0;
		if( cd.solution_type == EXTERNAL && cd.num_proc > 1 && k == 0 ) // Parallel job
		{
			if( cd.debug || cd.mdebug ) printf( "Parallel execution of external jobs ...\n" );
			if( cd.debug || cd.mdebug ) printf( "Generation of all the model input files ...\n" );
			for( count = 0; count < cd.nreal; count ++ ) // Write all the files
			{
				fprintf( out, "%d : ", count + 1 ); // counter
				if( cd.mdebug ) printf( "\n" );
				printf( "Random set #%d: ", count + 1 );
				fflush( stdout );
				for( i = 0; i < pd.nOptParam; i++ )
				{
					k = pd.var_index[i];
					opt_params[i] = pd.var[k] = var_lhs[i + count * npar] * pd.var_range[k] + pd.var_min[k];
				}
				if( cd.mdebug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
				Transform( opt_params, &op, opt_params );
				func_extrn_write( count + 1, opt_params, &op );
				printf( "external model input file(s) generated ...\n" );
				if( cd.mdebug ) cd.fdebug = debug_level;
				if( cd.mdebug )
				{
					printf( "\nRandom parameter values:\n" );
					for( i = 0; i < pd.nOptParam; i++ )
						if( pd.var_log[pd.var_index[i]] == 0 ) printf( "%s %g\n", pd.var_id[pd.var_index[i]], pd.var[pd.var_index[i]] );
						else printf( "%s %g\n", pd.var_id[pd.var_index[i]], pow( 10, pd.var[pd.var_index[i]] ) );
				}
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 )
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				fprintf( out, "\n" );
			}
			fclose( out );
			if( cd.pardebug > 4 )
			{
				for( count = 0; count < cd.nreal; count ++ ) // Perform all the runs in serial model (for testing)
				{
					printf( "Execute model #%d ... ", count + 1 );
					fflush( stdout );
					func_extrn_exec_serial( count + 1, &op );
					printf( "done!\n" );
				}
			}
			else if( mprun( cd.nreal, &op ) < 0 ) // Perform all the runs in parallel
			{
				printf( "ERROR: there is a problem with the parallel execution!\n" );
				sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
				system( buf );
			}
			out = Fwrite( filename ); // rewrite results file including the results
			for( count = 0; count < cd.nreal; count ++ ) // Read all the files
			{
				if( cd.debug || cd.mdebug ) printf( "Reading all the model output files ...\n" );
				op.counter = count + 1;
				fprintf( out, "%d : ", op.counter ); // counter
				for( i = 0; i < pd.nOptParam; i++ ) // re
				{
					k = pd.var_index[i];
					pd.var[k] = var_lhs[i + count * npar] * pd.var_range[k] + pd.var_min[k];
				}
				if( cd.mdebug ) printf( "\n" );
				printf( "Model results #%d: ", op.counter );
				bad_data = 0;
				bad_data = func_extrn_read( count + 1, &op, od.res );
				if( bad_data )
				{
					sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
					system( buf );
					exit( -1 );
				}
				if( cd.mdebug > 1 ) { printf( "\n" ); print_results( &op, 0 ); }
				else if( cd.mdebug )
				{
					printf( "Objective function: %g Success: %d\n", op.phi, success_all );
					if( success_all ) printf( "All the predictions are within calibration ranges!\n" );
					else printf( "At least one of the predictions is outside calibration ranges!\n" );
				}
				else
					printf( "Objective function: %g Success = %d\n", op.phi, success_all );
				if( op.phi < phi_min )
				{
					phi_min = op.phi;
					for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]];
					for( i = 0; i < od.nObs; i++ ) od.obs_best[i] = od.obs_current[i];
				}
				if( success_all ) success_global++;
				if( op.phi < op.cd->phi_cutoff ) phi_global++;
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 )
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				if( od.nObs > 0 ) fprintf( out, " OF %g success %d\n", op.phi, success_all );
				fflush( out );
				if( ( success_all || od.nObs == 0 ) && cd.odebug > 1 ) save_results( "mcrnd", &op, &gd );
			}
		}
		else // Serial job
		{
			for( count = k; count < cd.nreal; count ++ )
			{
				op.counter = count + 1;
				fprintf( out, "%d : ", count + 1 ); // counter
				if( cd.mdebug ) printf( "\n" );
				printf( "Random set #%d: ", count + 1 );
				fflush( stdout );
				for( i = 0; i < pd.nOptParam; i++ )
				{
					j = pd.var_index[i];
					opt_params[i] = pd.var[j] = var_lhs[i + count * npar] * pd.var_range[j] + pd.var_min[j];
				}
				if( cd.mdebug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
				Transform( opt_params, &op, opt_params );
				func_global( opt_params, &op, od.res );
				if( cd.mdebug ) cd.fdebug = debug_level;
				if( cd.mdebug )
				{
					printf( "\nRandom parameter values:\n" );
					for( i = 0; i < pd.nOptParam; i++ )
					{
						j = pd.var_index[i];
						if( pd.var_log[j] == 0 ) printf( "%s %g\n", pd.var_id[j], pd.var[j] );
						else printf( "%s %g\n", pd.var_id[j], pow( 10, pd.var[j] ) );
					}
				}
				if( cd.mdebug > 1 ) { printf( "\nPredicted calibration targets:\n" ); print_results( &op, 1 ); }
				else if( cd.mdebug )
				{
					printf( "Objective function: %g Success: %d\n", op.phi, success_all );
					if( success_all ) printf( "All the predictions are within calibration ranges!\n" );
					else printf( "At least one of the predictions is outside calibration ranges!\n" );
				}
				else
					printf( "Objective function: %g Success = %d\n", op.phi, success_all );
				if( op.phi < phi_min )
				{
					phi_min = op.phi;
					for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]];
					for( i = 0; i < od.nObs; i++ ) od.obs_best[i] = od.obs_current[i];
				}
				if( success_all ) success_global++;
				if( op.phi < op.cd->phi_cutoff ) phi_global++;
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 )
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				if( od.nObs > 0 ) fprintf( out, " OF %g success %d\n", op.phi, success_all );
				fflush( out );
				if( ( success_all || od.nObs == 0 ) && cd.odebug > 1 ) save_results( "mcrnd", &op, &gd );
				if( op.f_ofe != NULL ) { fclose( op.f_ofe ); op.f_ofe = NULL; }
				if( cd.ireal != 0 ) break;
			}
		}
		op.counter = 0;
		free( var_lhs );
		fclose( out );
		op.phi = phi_min; // get the best phi
		for( i = 0; i < pd.nOptParam; i++ ) opt_params[i] = pd.var[pd.var_index[i]] = pd.var_current[i] = pd.var_best[i]; // get the best estimate
		for( i = 0; i < od.nObs; i++ ) od.obs_current[i] = od.obs_best[i] ; // get the best observations
		printf( "\nMinimum objective function: %g\n", phi_min );
		print_results( &op, 1 );
		if( cd.debug )
		{
			printf( "Repeat the run producing the best results ...\n" );
			debug_level = cd.fdebug; cd.fdebug = 3;
			Transform( opt_params, &op, opt_params );
			func_global( opt_params, &op, od.res );
			cd.fdebug = debug_level;
		}
		printf( "Results are saved in %s.mcrnd.results\n", op.root );
		printf( "\nMinimum objective function: %g\n", phi_min );
		if( success_global == 0 ) printf( "None of the Monte-Carlo runs produced predictions within calibration ranges!\n" );
		else printf( "Number of Monte-Carlo runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, cd.nreal, ( double ) success_global / cd.nreal );
		if( op.cd->phi_cutoff > DBL_EPSILON )
		{
			if( phi_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", cd.nreal, op.cd->phi_cutoff );
			else printf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op.cd->phi_cutoff, phi_global, cd.nreal, ( double ) phi_global / cd.nreal );
		}
		free( opt_params );
		save_results( "", &op, &gd );
	}
	//
	// ------------------------ GLOBALSENS
	//
	// TODO gsens needs to be a separate function
	if( cd.problem_type == GLOBALSENS ) // Global sensitivity analysis run
	{
		strcpy( op.label, "gsens" );
		double fhat, fhat2, *phis_full, *phis_half;
		int n_sub; //! number of samples for subsets a and b
		//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_sobol, pd.nOptParam );
		n_sub = cd.nreal / 2;	// set to half of user specified reals
		if( ( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Temporary variable to store cd.nreal phis
		if( ( phis_full = ( double * ) malloc( cd.nreal * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Temporary variable to store m_sub phis
		if( ( phis_half = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Temporary variable to store random sample a
		if( ( var_a_lhs = ( double * ) malloc( pd.nOptParam * n_sub * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Sample a phis
		if( ( gs.f_a = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Sample b phis
		if( ( gs.f_b = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Temporary variable to store random sample b
		if( ( var_b_lhs = ( double * ) malloc( pd.nOptParam * n_sub * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// matrices to store lhs samples
		gs.var_a_lhs = double_matrix( n_sub, pd.nOptParam );
		gs.var_b_lhs = double_matrix( n_sub, pd.nOptParam );
		// Matrices to store phis with different combinations of parameters from samples a and b
		if( ( gs.fmat_a = double_matrix( pd.nOptParam, n_sub ) ) == NULL )
		{ printf( "Error creating 3D matrix\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		if( ( gs.fmat_b = double_matrix( pd.nOptParam, n_sub ) ) == NULL )
		{ printf( "Error creating 3D matrix\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Vector of variances for individual component contribution
		if( ( gs.D_hat = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Vector of variances for total component contribution
		if( ( gs.D_hat_n = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		printf( "\nGlobal sensitivity analysis using random sampling:\n" );
		// Create samples
		if( cd.seed < 0 ) { cd.seed *= -1; printf( "Imported seed: %d\n", cd.seed ); }
		else if( cd.seed == 0 ) { printf( "New " ); cd.seed_init = cd.seed = get_seed(); }
		else printf( "Current seed: %d\n", cd.seed );
		printf( "Random sampling set 1 (variables %d; realizations %d) using ", pd.nOptParam, cd.nreal );
		sampling( pd.nOptParam, n_sub, &cd.seed, var_a_lhs, &op, 1 );
		printf( "done.\n" );
		printf( "Random sampling set 2 (variables %d; realizations %d) using ", pd.nOptParam, cd.nreal );
		sampling( pd.nOptParam, n_sub, &cd.seed, var_b_lhs, &op, 1 );
		printf( "done.\n" );
		// Create samples using Sobol's quasi-random sequence
		/*		for( count = 0; count < n_sub; count++ )
				{
					double v[ pd.nOptParam ];
					gsl_qrng_get( q, v);
					for( i = 0; i < pd.nOptParam; i++ )
					{
						k = pd.var_index[i];
						gs.var_a_lhs[count][i] = v[i] * pd.var_range[k] + pd.var_min[k];
					}
				}

				for( count = 0; count < n_sub; count++ )
				{
					double v[ pd.nOptParam ];
					gsl_qrng_get( q, v);
					for( i = 0; i < pd.nOptParam; i++ )
					{
						k = pd.var_index[i];
						gs.var_b_lhs[count][i] = v[i] * pd.var_range[k] + pd.var_min[k];
					}
				}*/
		// Copy temp lhs vectors to matrices
		for( count = 0; count < n_sub; count++ )
			for( i = 0; i < pd.nOptParam; i++ )
			{
				k = pd.var_index[i];
				gs.var_a_lhs[count][i] = var_a_lhs[i + count * pd.nOptParam] * pd.var_range[k] + pd.var_min[k];
				gs.var_b_lhs[count][i] = var_b_lhs[i + count * pd.nOptParam] * pd.var_range[k] + pd.var_min[k];
			}
		free( var_a_lhs );
		free( var_b_lhs );
		// Output samples to files
		if( cd.mdebug )
		{
			sprintf( filename, "%s.gsens.zip", op.root );
			if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.gsens.zip %s.gsens_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
			sprintf( buf, "zip -m %s.gsens.zip %s.gsens_set_* >& /dev/null", op.root, op.root ); system( buf );
			sprintf( buf, "mv %s.gsens.zip %s.gsens_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
			sprintf( filename, "%s.gsens_set_a", op.root ); out = Fwrite( filename );
			sprintf( filename, "%s.gsens_set_b", op.root ); out2 = Fwrite( filename );
			for( count = 0; count < n_sub; count ++ )
			{
				for( k = 0; k < pd.nOptParam; k++ )
				{
					fprintf( out, "%.15g ", gs.var_a_lhs[count][k] );
					fprintf( out2, "%.15g ", gs.var_b_lhs[count][k] );
				}
				fprintf( out, "\n" );
				fprintf( out2, "\n" );
			}
			fclose( out );
			fclose( out2 );
			printf( "Random sampling sets a and b saved in %s.mcrnd_set_a and %s.mcrnd_set_b\n", op.root, op.root );
		}
		sprintf( filename, "%s.gsens.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.gsens_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		// Accumulate phis into fhat and fhat2 for total output mean and variance
		fhat = fhat2 = 0;
		printf( "Computing phis to calculate total output mean and variance...\n" );
		// Compute sample a phis
		for( count = 0; count < n_sub; count ++ )
		{
			for( i = 0; i < pd.nOptParam; i++ )
			{
				k = pd.var_index[i];
				opt_params[i] = pd.var[k] = gs.var_a_lhs[count][i];
			}
			Transform( opt_params, &op, opt_params );
			func_global( opt_params, &op, od.res );
			// Sum phi and phi^2
			fhat += op.phi;
			fhat2 += pow( op.phi, 2 );
			// Save sample a phis
			gs.f_a[count] = op.phi;
			phis_full[count] = op.phi;
			// save to results file
			fprintf( out, "%d : ", count + 1 ); // counter
			fprintf( out, "%g :", op.phi );
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] >= 1 )
					fprintf( out, " %.15g", pd.var[i] );
			fprintf( out, "\n" );
			fflush( out );
		}
		// Compute sample b phis
		for( count = 0; count < n_sub; count ++ )
		{
			for( i = 0; i < pd.nOptParam; i++ )
			{
				k = pd.var_index[i];
				opt_params[i] = pd.var[k] = gs.var_b_lhs[count][i];
			}
			Transform( opt_params, &op, opt_params );
			func_global( opt_params, &op, od.res );
			// Sum phi and phi^2
			fhat += op.phi;
			fhat2 += pow( op.phi, 2 );
			// Save sample b phis
			gs.f_b[count] = op.phi;
			phis_full[ n_sub + count ] = op.phi;
			// save to results file
			fprintf( out, "%d : ", n_sub + count ); // counter
			fprintf( out, "%g :", op.phi );
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] >= 1 )
					fprintf( out, " %.15g", pd.var[i] );
			fprintf( out, "\n" );
			fflush( out );
		}
		fclose( out );
		printf( "Global Sensitivity MC results are saved in %s.gsens.results\n", op.root );
		// Calculate total output mean and variance based on sample a
		gs.f_hat_0 = fhat / ( 2 * n_sub );
		gs.D_hat_t = fhat2 / ( 2 * n_sub ) - gs.f_hat_0;
		printf( "Total output mean = %g\n", gs.f_hat_0 );
		printf( "Total output variance = %g\n", gs.D_hat_t );
		gs.f_hat_0 = gsl_stats_mean( phis_full, 1, cd.nreal );
		gs.D_hat_t = gsl_stats_variance( phis_full, 1, cd.nreal );
		printf( "Total output mean = %g\n", gs.f_hat_0 );
		printf( "Total output variance = %g\n", gs.D_hat_t );
		gs.f_hat_0 = gs.D_hat_t = 0.0;
		ave_sorted( phis_full, cd.nreal, &gs.f_hat_0, &gs.ep );
		printf( "Total output mean = %g abs 1st moment = %g\n", gs.f_hat_0, gs.ep );
		var_sorted( phis_full, phis_full, cd.nreal, gs.f_hat_0, gs.ep, &gs.D_hat_t );
		printf( "Total output variance = %g\n", gs.D_hat_t );
		/*		// Subtract f_hat_0 from phis and recalculate total output variance
				fhat2 = 0;
				for( count = 0; count < n_sub; count++ )
				{
					gs.f_a[count] -= gs.f_hat_0;
					gs.f_b[count] -= gs.f_hat_0;
					phis_full[ count ] = gs.f_a[count];
					phis_full[ n_sub + count ] = gs.f_b[count];
					fhat2 += pow( gs.f_a[count], 2 );
					fhat2 += pow( gs.f_b[count], 2 );
				}
				gs.D_hat_t = fhat2 / (2 * n_sub);
				printf( "Total output variance = %g\n", gs.D_hat_t );
				gs.D_hat_t = gsl_stats_variance( phis_full, 1, cd.nreal );
				printf( "Total output variance = %g\n", gs.D_hat_t );
		 */		free( phis_full );
		// Collect matrix of phis for fmat_a
		printf( "Computing phis for calculation of individual output variances:\n" );
		fflush( stdout );
		for( i = 0; i < pd.nOptParam; i++ )
		{
			printf( "Parameter %d...\n", i + 1 );
			for( count = 0; count < n_sub; count ++ )
			{
				for( j = 0; j < pd.nOptParam; j++ )
				{
					k = pd.var_index[j];
					if( i == j ) // then select from sample a
						opt_params[j] = pd.var[k] = gs.var_a_lhs[count][j];
					else // else select from sample b
						opt_params[j] = pd.var[k] = gs.var_b_lhs[count][j];
				}
				Transform( opt_params, &op, opt_params );
				func_global( opt_params, &op, od.res );
				// Save phi to fmat_a
				gs.fmat_a[i][count] = op.phi;
			}
		}
		// Collect matrix of phis for fmat_b
		printf( "Computing phis for calculation of individual plus interaction output variances:\n" );
		for( i = 0; i < pd.nOptParam; i++ )
		{
			printf( "Parameter %d...\n", i + 1 );
			for( count = 0; count < n_sub; count ++ )
			{
				for( j = 0; j < pd.nOptParam; j++ )
				{
					k = pd.var_index[j];
					if( i == j ) // then select from sample b
						opt_params[j] = pd.var[k] = gs.var_b_lhs[count][j];
					else // else select from sample a
						opt_params[j] = pd.var[k] = gs.var_a_lhs[count][j];
				}
				Transform( opt_params, &op, opt_params );
				func_global( opt_params, &op, od.res );
				// Save phi to fmat_b
				gs.fmat_b[i][count] = op.phi;
			}
		}
		printf( "done.\n" );
		// Calculate individual and interaction output variances
		for( i = 0; i < pd.nOptParam; i++ )
		{
			fhat2 = 0;
			for( j = 0; j < n_sub; j++ )
			{
				fhat2 += ( gs.f_a[j] * gs.fmat_a[i][j] );
				phis_half[ j ] = ( gs.f_a[j] * gs.fmat_a[i][j] );
			}
			gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			gs.D_hat[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			gs.D_hat[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_a[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			var_sorted( gs.f_a, gs.fmat_a[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat[i] );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			//gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
			fhat2 = 0;
			for( j = 0; j < n_sub; j++ )
			{
				fhat2 += ( gs.f_a[j] * gs.fmat_b[i][j] );
				phis_half[ j ] = ( gs.f_a[j] * gs.fmat_b[i][j] );
			}
			gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			gs.D_hat_n[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			gs.D_hat_n[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_b[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			var_sorted( gs.f_a, gs.fmat_b[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_n[i] );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			//gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		}
		// Print sensitivity indices
		printf( "\nParameter sensitivity indices:\n" );
		printf( "parameter individual interaction\n" );
		for( i = 0; i < pd.nOptParam; i++ ) printf( "%d %g %g\n", i + 1, gs.D_hat[i] / gs.D_hat_t, 1 - ( gs.D_hat_n[i] / gs.D_hat_t ) );
		printf( "\n" );
		free( opt_params ); free( phis_half ); free( gs.f_a ); free( gs.f_b ); free( gs.D_hat ); free( gs.D_hat_n );
		free_matrix( ( void ** ) gs.var_a_lhs, n_sub );
		free_matrix( ( void ** ) gs.var_b_lhs, n_sub );
		free_matrix( ( void ** ) gs.fmat_a, pd.nOptParam );
		free_matrix( ( void ** ) gs.fmat_b, pd.nOptParam );
	}
	//
	// ------------------------ SIMPLE CALIBRATION
	//
	if( cd.problem_type == CALIBRATE && cd.calib_type == SIMPLE ) /* Inverse analysis */
	{
		if( cd.nretries > 1 ) printf( "\nMULTI-START CALIBRATION using a series of random initial guesses:\n" );
		else printf( "\nSINGLE CALIBRATION: single optimization based on initial guesses provided in the input file:\n" );
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
		else optimize_func = optimize_pso; // Define optimization method: PSO
		for( i = 0; i < pd.nParam; i++ ) cd.var[i] = pd.var[i]; // Set all the initial values
		success = optimize_func( &op ); // Optimize
		if( success == 0 ) { printf( "ERROR: Optimization did not start!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		sprintf( filename, "%s-rerun.mads", op.root );
		if( cd.solution_type != TEST ) save_problem( filename, &op );
		if( cd.debug == 0 ) printf( "\n" );
		print_results( &op, 1 );
		save_results( "", &op, &gd );
		if( od.nObs < od.nTObs )
		{
			predict = 1; // Produce outputs for predictions that are not calibration targets (see below PREDICT)
			sprintf( filename, "%s.results", op.root );
			out = Fappend( filename ); // Reopen results file
			printf( "\nModel predictions that are not calibration targets:\n" );
			fprintf( out, "\nModel predictions that are not calibration targets:\n" );
		}
		else predict = 0; // There are no observations that are not calibration targets
	}
	strcpy( op.label, "" ); // No labels needed below
	//
	// ------------------------ EIGEN || LOCALSENS
	//
	if( cd.problem_type == EIGEN || cd.problem_type == LOCALSENS )
	{
		status = eigen( &op, NULL, NULL ); // Eigen or sensitivity analysis run
		if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	//
	//------------------------- ABAGUS
	//
	if( cd.problem_type == ABAGUS ) // Particle swarm sensitivity analysis run
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = pssa( &op );
		if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	//
	//------------------------ POSTPUA
	//
	if( cd.problem_type == POSTPUA ) // Predictive uncertainty analysis of sampling results
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = postpua( &op );
		if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	//
	//------------------------ GLUE
	//
	if( cd.problem_type == GLUE ) // Generalized Likelihood Uncertainty Estimation
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = glue( &op );
		if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	//
	//------------------------ INFOGAP
	//
	if( cd.problem_type == INFOGAP ) // Info-gap decision analysis
	{
		if( cd.pardx < DBL_EPSILON ) cd.pardx = 0.1;
		status = infogap( &op );
		if( status == 0 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
	}
	//
	// ------------------------ FORWARD
	//
	if( cd.problem_type == FORWARD ) // Forward run
	{
		sprintf( filename, "%s.results", op.root );
		out = Fwrite( filename );
		fprintf( out, "Model parameter values:\n" );
		for( i = 0; i < pd.nParam; i++ )
			fprintf( out, "%s %g\n", pd.var_id[i], pd.var[i] );
		if( od.nObs > 0 )
		{
			printf( "\nModel predictions (forward run; no calibration):\n" );
			fprintf( out, "\nModel predictions (forward run; no calibration):\n" );
			sprintf( filename, "%s.forward", op.root );
			out2 = Fwrite( filename );
		}
	}
	//
	// ------------------------ CREATE
	//
	if( cd.problem_type == CREATE ) // Create a MADS file based on a forward run
		printf( "\nModel predictions (forward run; no calibration):\n" );
	//
	// ------------------------ PREDICT
	//
	if( cd.problem_type == FORWARD || cd.problem_type == CREATE || predict )
	{
		success_all = 1;
		compare = 0;
		phi = 0;
		if( cd.solution_type == EXTERNAL )
			for( i = 0; i < od.nObs; i++ )
			{
				if( cd.problem_type == CALIBRATE && od.obs_weight[i] != 0 ) continue;
				compare = 1;
				c = od.obs_current[i];
				err = od.obs_target[i] - c;
				phi += ( err * err ) * od.obs_weight[i];
				if( ( c < od.obs_min[i] || c > od.obs_max[i] ) && ( wd.obs_weight[i][j] > 0.0 ) ) { success_all = 0; success = 0; }
				else success = 1;
				if( od.nObs < 50 || ( i < 20 || i > od.nObs - 20 ) ) printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", od.obs_id[i], od.obs_target[i], c, err, err * od.obs_weight[i], success, od.obs_min[i], od.obs_max[i] );
				if( od.nObs > 50 && i == 21 ) printf( "...\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", od.obs_id[i], od.obs_target[i], c, err, err * od.obs_weight[i], success, od.obs_min[i], od.obs_max[i] );
				else od.obs_target[i] = c; // Save computed values as calibration targets
			}
		else if( cd.solution_type != TEST )
			for( i = 0; i < wd.nW; i++ )
				for( j = 0; j < wd.nWellObs[i]; j++ )
				{
					if( cd.problem_type == CALIBRATE && wd.obs_weight[i][j] > DBL_EPSILON ) continue;
					compare = 1;
					c = func_solver( wd.x[i], wd.y[i], wd.z1[i], wd.z2[i], wd.obs_time[i][j], &cd );
					err = wd.obs_target[i][j] - c;
					if( cd.problem_type != CALIBRATE ) phi += ( err * err ) * wd.obs_weight[i][j];
					else phi += ( err * err );
					if( ( c < wd.obs_min[i][j] || c > wd.obs_max[i][j] ) && ( wd.obs_weight[i][j] > 0.0 ) ) { success_all = 0; success = 0; }
					else success = 1;
					if( cd.problem_type != CALIBRATE )
						printf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err * wd.obs_weight[i][j], success, wd.obs_min[i][j], wd.obs_max[i][j] );
					else
						printf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err, success, wd.obs_min[i][j], wd.obs_max[i][j] );
					if( cd.problem_type != CREATE ) fprintf( out, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err * wd.obs_weight[i][j], success );
					else wd.obs_target[i][j] = c; // Save computed values as calibration targets
					if( cd.problem_type == FORWARD ) fprintf( out2, "%s(%g) %g\n", wd.id[i], wd.obs_time[i][j], c ); // Forward run
				}
		cd.neval++;
		if( compare )
		{
			op.phi = phi;
			printf( "Objective function: %g Success: %d\n", op.phi, success_all );
			if( cd.problem_type != CREATE ) fprintf( out, "Objective function = %g Success: %d\n", op.phi, success_all );
			if( success_all )
			{
				printf( "All the predictions are within acceptable ranges!\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "All the predictions are within acceptable ranges!\n" );
			}
			else
			{
				printf( "At least one of the predictions is outside acceptable ranges!\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "At least one of the predictions is outside acceptable ranges!\n" );
			}
		}
		else    printf( "No calibration targets!\n" );
		fclose( out );
	}
	if( cd.problem_type == FORWARD && od.nObs > 0 ) fclose( out2 );
	if( cd.problem_type == FORWARD || cd.problem_type == CALIBRATE )
	{
		if( od.nObs > 0 )
		{
			sprintf( filename, "%s.phi", op.root );
			out2 = Fwrite( filename );
			fprintf( out2, "%g\n", op.phi ); // Write phi in a separate file
			fclose( out2 );
		}
	}
	if( cd.problem_type == CREATE ) /* Create a file with calibration targets equal to the model predictions */
	{
		cd.problem_type = CALIBRATE;
		sprintf( filename, "%s-truth.mads", op.root );
		save_problem( filename, &op );
		printf( "\nMADS problem file named %s-truth.mads is created; modify the file if needed\n\n", op.root );
	}
	free( orig_params );
	// Finalize the run
	time_end = time( NULL );
	time_elapsed = time_end - time_start;
	if( time_elapsed > 86400 ) printf( "Simulation time = %g days\n", ( ( double ) time_elapsed / 86400 ) );
	else if( time_elapsed > 3600 ) printf( "Simulation time = %g hours\n", ( ( double ) time_elapsed / 3600 ) );
	else if( time_elapsed > 60 ) printf( "Simulation time = %g minutes\n", ( ( double ) time_elapsed / 60 ) );
	else printf( "Simulation time = %ld seconds\n", time_elapsed );
	printf( "Functional evaluations = %d\n", cd.neval );
	if( cd.njac > 0 ) printf( "Jacobian evaluations = %d\n", cd.njac );
	if( cd.problem_type == CALIBRATE ) printf( "Levenberg-Marquardt optimizations = %d\n", cd.nlmo );
	if( time_elapsed > 0 )
	{
		c = cd.neval / time_elapsed;
		if( c < ( ( double ) 1 / 86400 ) ) printf( "Functional evaluations per day = %g\n", c * 86400 );
		else if( c < ( ( double ) 1 / 3600 ) ) printf( "Functional evaluations per hour = %g\n", c * 3600 );
		else if( c < ( ( double ) 1 / 60 ) ) printf( "Functional evaluations per minute = %g\n", c * 60 );
		else printf( "Functional evaluations per second = %g\n", c );
	}
	if( op.cd->seed_init > 0 ) printf( "Seed = %d\n", op.cd->seed_init );
	ptr_ts = localtime( &time_start );
	printf( "Execution  started  on %s", asctime( ptr_ts ) );
	ptr_ts = localtime( &time_end );
	printf( "Execution completed on %s", asctime( ptr_ts ) );
	printf( "Execution date & time stamp: %s\n", op.datetime_stamp );
	sprintf( buf, "rm -f %s.running", op.root ); system( buf );
	if( op.f_ofe != NULL ) { fclose( op.f_ofe ); op.f_ofe = NULL; }
	exit( 0 ); // DONE
}

int optimize_pso( struct opt_data *op )
{
	if( op->cd->debug ) printf( "\nParticle-Swarm Optimization:" );
	if( strncasecmp( op->cd->opt_method, "pso", 3 ) == 0 || strncasecmp( op->cd->opt_method, "swarm", 5 ) == 0 )
	{
		if( op->cd->debug ) printf( " Standard (2006)\n" );
		pso_std( op );
	}
	else if( strncasecmp( op->cd->opt_method, "tribes", 5 ) == 0 && strcasestr( op->cd->opt_method, "std" ) != NULL )
	{
		if( op->cd->debug ) printf( " TRIBES (Clerc, 2006)\n" );
		mopso( op );
	}
	else
	{
		if( op->cd->debug ) printf( " TRIBES\n" );
		pso_tribes( op );
	}
	if( op->cd->debug )
	{
		printf( "\n------------------------- Optimization Results:\n" );
		print_results( op, 1 );
	}
	if( op->cd->leigen && op->cd->solution_type != TEST )
		if( eigen( op, NULL, NULL ) == 0 ) // Execute eigen analysis of the final results
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
	debug = op->cd->debug;
	standalone = op->cd->standalone;
	if( !standalone ) op->cd->paranoid = 0;
	if( op->od->nObs == 0 ) { printf( "ERROR: Number of observations is equal to zero! Levenberg-Marquardt Optimization cannot be performed!\n" ); return( 0 ); }
	if( op->pd->nOptParam == 0 ) { printf( "ERROR: Number of optimized model parameters is equal to zero! Levenberg-Marquardt Optimization cannot be performed!\n" ); return( 0 ); }
	if( ( op->pd->nOptParam > op->od->nObs ) && ( standalone && op->cd->calib_type == SIMPLE ) ) { printf( "WARNING: Number of optimized model parameters is greater than number of observations (%d>%d)\n", op->pd->nOptParam, op->od->nObs ); }
	gsl_matrix *gsl_jacobian = gsl_matrix_alloc( op->od->nObs, op->pd->nOptParam );
	gsl_matrix *gsl_covar = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam );
	gsl_vector *gsl_opt_params = gsl_vector_alloc( op->pd->nOptParam );
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( ( opt_params_best = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( ( x_c = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( ( res = ( double * ) malloc( op->od->nObs * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( op->cd->niter <= 0 )
	{
		if( strcasestr( op->cd->opt_method, "lm" ) ) maxiter = 50;
		if( strcasestr( op->cd->opt_method, "squad" ) ) maxiter = 8;
	}
	else maxiter = op->cd->niter;
	if( op->cd->ldebug && standalone ) printf( "Number of Levenberg-Marquardt iterations = %d\n", maxiter );
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]];
	if( !standalone )
	{
		// for( i = 0; i < op->pd->nOptParam; i++ )
		// printf( "%g\n", opt_params[i] );
	}
	else
		Transform( opt_params, op, opt_params );
	if( op->cd->paranoid )
	{
		if( standalone ) printf( "Multi-Start Levenberg-Marquardt (MSLM) Optimization ... " ); fflush( stdout );
		npar = op->pd->nOptParam;
		if( op->cd->nretries <= 0 ) op->cd->nretries = ( double )( op->cd->maxeval - op->cd->neval ) / ( maxiter * npar / 10 );
		if( debug ) printf( "\nRandom sampling for MSLM optimization (variables %d; realizations %d) using ", npar, op->cd->nretries );
		if( ( var_lhs = ( double * ) malloc( npar * op->cd->nretries * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 0 ); }
		if( op->cd->seed < 0 ) { op->cd->seed *= -1; if( debug ) printf( "Imported seed: %d\n", op->cd->seed ); }
		else if( op->cd->seed == 0 ) { if( debug ) printf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
		else if( debug ) printf( "Current seed: %d\n", op->cd->seed );
		if( op->cd->paran_method[0] != 0 ) { strcpy( buf, op->cd->smp_method ); strcpy( op->cd->smp_method, op->cd->paran_method ); }
		sampling( npar, op->cd->nretries, &op->cd->seed, var_lhs, op, debug );
		if( op->cd->paran_method[0] != 0 ) strcpy( op->cd->smp_method, buf );
		if( debug ) printf( "done.\n" );
		op->cd->retry_ind = count = count_set = 0;
	}
	else
	{ if( standalone ) printf( "Levenberg-Marquardt Optimization ... " ); fflush( stdout ); }
	phi_min = HUGE_VAL;
	do // BEGIN Paranoid loop
	{
		if( op->cd->maxeval <= op->cd->neval ) { if( debug || op->cd->paranoid ) printf( "Maximum number of evaluations is exceeded (%d <= %d)!\n", op->cd->maxeval, op->cd->neval ); break; }
		op->cd->nlmo++;
		if( op->cd->paranoid )
		{
			count++;
			op->cd->retry_ind = count;
			if( op->cd->calib_type == IGRND && count == 1 && debug )
				printf( "CALIBRATION %d: initial guesses from IGRND random set: ", count );
			else if( count > 1 )
			{
				if( op->cd->ldebug ) printf( "\n********************************************************************\n" );
				if( debug ) printf( "CALIBRATION %d: initial guesses from internal MSLM random set #%d: ", count, count_set + 1 );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					opt_params[i] = var_lhs[i + count_set * npar] * op->pd->var_range[k] + op->pd->var_min[k];
					if( debug > 1 )
					{
						if( op->pd->var_log[k] ) printf( "%s %.15g\n", op->pd->var_id[k], pow( 10, opt_params[i] ) );
						else printf( "%s %.15g\n", op->pd->var_id[k], opt_params[i] );
					}
				}
				count_set++;
			}
			else
			{
				if( debug ) printf( "CALIBRATION %d: initial guesses from MADS input file #: ", count );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					opt_params[i] = op->pd->var[op->pd->var_index[i]];
					if( debug > 1 )
					{
						if( op->pd->var_log[k] ) printf( "%s %.15g\n", op->pd->var_id[k], pow( 10, opt_params[i] ) );
						else printf( "%s %.15g\n", op->pd->var_id[k], opt_params[i] );
					}
				}
			}
			if( debug > 1 ) printf( "\n" );
			Transform( opt_params, op, opt_params );
			fflush( stdout );
		}
		if( debug > 1 && standalone )
		{
			printf( "\n-------------------- Initial state:\n" );
			op->cd->pderiv = op->cd->oderiv = -1;
			debug_level = op->cd->fdebug; op->cd->fdebug = 3;
			func_global( opt_params, op, res );
			op->cd->fdebug = debug_level;
		}
		// LM optimization ...
		if( strcasestr( op->cd->opt_method, "mon" ) != NULL || strcasestr( op->cd->opt_method, "chav" ) != NULL ) // Monty/Chavo versions
		{
			if( debug > 1 && standalone ) printf( "\nLevenberg-Marquardt Optimization:\n" );
			else if( op->cd->ldebug ) printf( "\n" );
			if( ( jacobian = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->od->nObs ) ) == NULL )
			{ printf( "Not enough memory!\n" ); return( 0 ); }
			if( ( jacTjac = ( double * ) malloc( sizeof( double ) * ( ( op->pd->nOptParam + 1 ) * op->pd->nOptParam / 2 ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); return( 0 ); }
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
			if( strcasestr( op->cd->opt_method, "chav" ) != NULL ) ier = zxssqch( func_global, op, op->od->nObs, op->pd->nOptParam, nsig, eps, delta, maxfn, iopt, opt_parm, opt_params, &phi, res, jacobian, op->od->nObs, jacTjac, &infer ); // Chavo's version
			else ier = lm_opt( func_global, func_dx, op, op->od->nObs, op->pd->nOptParam, nsig, eps, delta, maxfn, maxiter, iopt, opt_parm, opt_params, &phi, res, jacobian, op->od->nObs, jacTjac, &infer ); // Monty's version
			for( k = i = 0; i < op->pd->nOptParam; i++ )
				for( j = 0; j < op->od->nObs; j++, k++ )
					gsl_matrix_set( gsl_jacobian, j, i, jacobian[k] );
			gsl_multifit_covar( gsl_jacobian, 0.0, gsl_covar );
			free( jacTjac ); free( jacobian );
			op->phi = phi;
		}
		else if( strcasestr( op->cd->opt_method, "gsl" ) != NULL ) // GSL version of LM
		{
			if( debug > 1 && standalone ) printf( "\nLevenberg-Marquardt Optimization using GSL library:\n" );
			else if( op->cd->ldebug ) printf( "\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
				gsl_vector_set( gsl_opt_params, i, opt_params[i] );
			lm_gsl( gsl_opt_params, op, gsl_jacobian, gsl_covar );
			for( i = 0; i < op->pd->nOptParam; i++ )
				opt_params[i] = gsl_vector_get( gsl_opt_params, i );
			phi = op->phi;
		}
		else if( strcasestr( op->cd->opt_method, "tra" ) != NULL )// Transtrum version of LM
		{
			if( debug > 1 && standalone ) printf( "\nTranstrum version of Levenberg-Marquardt Optimization:\n" );
			else if( op->cd->ldebug ) printf( "\n" );
		}
		else // DEFAULT LevMar version of LM
		{
			if( debug > 1 && standalone ) printf( "\nLevenberg-Marquardt Optimization using LevMar library:\n" );
			else if( op->cd->ldebug ) printf( "\n" );
			if( ( covar = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->pd->nOptParam ) ) == NULL )
			{ printf( "Not enough memory!\n" ); return( 0 ); }
			// LM_DIF_WORKSZ(m,n) = 4*n+4*m + n*m + m*m
			if( ( work = ( double * ) malloc( sizeof( double ) * LM_DIF_WORKSZ( op->pd->nOptParam, op->od->nObs ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); return( 0 ); }
			for( i = 0; i < op->od->nObs; i++ ) res[i] = 0;
			jacobian = work + op->pd->nOptParam + 2 * op->od->nObs;
			opts[0] = 1e-3; opts[1] = 1E-5; opts[2] = 1E-5;
			opts[3] = op->cd->phi_cutoff;
			if( op->cd->sintrans == 0 ) opts[4] = op->cd->lindx; // Forward difference; Central difference if negative; DO NOT USE CENTRAL DIFFERENCE
			else opts[4] = op->cd->sindx;
			while( op->cd->maxeval > op->cd->neval )
			{
				// Levmar has no termination criteria based on the number of functional evaluations or number of jacobian evaluations
				if( opts[4] > 0 ) maxiter_levmar = ( double )( ( op->cd->maxeval - op->cd->neval ) / ( op->pd->nOptParam + 10 ) + 1 ); // Forward derivatives
				else              maxiter_levmar = ( double )( ( op->cd->maxeval - op->cd->neval ) / ( 2 * op->pd->nOptParam + 10 ) + 1 ); // Central derivatives
				if( maxiter_levmar > maxiter ) maxiter_levmar = maxiter;
				maxiter_levmar *= 10; // Assuming about 10 lambda searches per iteration
				if( strcasestr( op->cd->opt_method, "dif" ) != NULL ) ier = dlevmar_dif( func_levmar, opt_params, res, op->pd->nOptParam, op->od->nObs, maxiter_levmar, opts, info, work, covar, op );
				else ier = dlevmar_der( func_levmar, func_dx_levmar, opt_params, res, op->pd->nOptParam, op->od->nObs, maxiter_levmar, opts, info, work, covar, op );
				if( info[6] == 4 || info[6] == 5 ) { opts[0] *= 10; if( op->cd->ldebug ) printf( "Rerun with larger initial lambda (%g)\n", opts[0] ); }
				else break;
			}
			if( op->cd->ldebug > 1 )
			{
				printf( "Levenberg-Marquardt Optimization completed after %g iteration (reason %g) (returned value %d)\n", info[5], info[6], ier );
				printf( "initial phi %g final phi %g ||J^T e||_inf %g ||Dp||_2 %g mu/max[J^T J]_ii %g\n", info[0], info[1], info[2], info[3], info[4] );
				printf( "function evaluation %g jacobian evaluations %g linear systems solved %g\n", info[7], info[8], info[9] );
			}
			op->cd->njac += info[8];
			for( k = j = 0; j < op->od->nObs; j++ )
				for( i = 0; i < op->pd->nOptParam; i++ )
					gsl_matrix_set( gsl_jacobian, j, i, jacobian[k++] );
			for( i = 0; i < op->pd->nOptParam; i++ )
				for( j = 0; j < op->pd->nOptParam; j++ )
					gsl_matrix_set( gsl_covar, i, j, covar[i * op->pd->nOptParam + j] );
			op->phi = phi = info[1];
			free( work ); free( covar );
		}
		if( standalone ) // if LM is stand alone (not part of PSO run)
		{
			if( debug > 1 )
			{
				printf( "\n------------------------- Final state:\n" );
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
				printf( "\n------------------------- LM Optimization Results:\n" );
				print_results( op, 1 );
			}
		}
		else // if LM is part of PSO run
			for( i = 0; i < op->pd->nOptParam; i++ )
				op->pd->var[op->pd->var_index[i]] = opt_params[i];
		if( op->cd->paranoid )
		{
			if( op->phi < phi_min ) { phi_min = op->phi; for( i = 0; i < op->pd->nOptParam; i++ ) opt_params_best[i] = op->pd->var[op->pd->var_index[i]]; }
			if( debug ) printf( "Objective function: %g Success %d\n", op->phi, op->success );
			if( phi_min < op->cd->phi_cutoff )
			{
				if( debug ) printf( "MSLM optimization objective function is below the cutoff value after %d random initial guess attempts\n", count );
				break;
			}
			if( op->cd->check_success && op->success )
			{
				if( debug ) printf( "MSLM optimization within calibration ranges after %d random initial guess attempts\n", count );
				phi_min = op->phi;
				for( i = 0; i < op->pd->nOptParam; i++ ) opt_params_best[i] = op->pd->var[op->pd->var_index[i]];
				break;
			}
			if( op->cd->maxeval <= op->cd->neval )
			{ if( debug ) printf( "MSLM optimization terminated after evaluations %d (max evaluations %d)\n", op->cd->neval, op->cd->maxeval ); break; }
			if( count == op->cd->nretries )
			{ if( debug ) printf( "MSLM optimization terminated after %d attempts (evaluations %d; max evaluations %d)\n", count, op->cd->neval, op->cd->maxeval ); break; }
		}
		else break; // Quit if not Paranoid run
	}
	while( 1 ); // END Paranoid loop
	if( op->cd->paranoid ) // Recompute for the best results
	{
		if( !debug ) printf( "(retries=%d) ", count );
		op->phi = phi_min;
		for( i = 0; i < op->pd->nOptParam; i++ )
			op->pd->var[op->pd->var_index[i]] = opt_params_best[i];
		Transform( opt_params_best, op, opt_params );
		func_global( opt_params, op, op->od->res );
	}
	if( ( op->cd->leigen || op->cd->ldebug || op->cd->debug ) && standalone && op->cd->calib_type == SIMPLE )
		if( eigen( op, gsl_jacobian, gsl_covar ) == 0 ) // Eigen analysis
			return( 0 );
	if( op->cd->paranoid ) free( var_lhs );
	free( opt_params ); free( opt_params_best ); free( x_c ); free( res );
	gsl_matrix_free( gsl_jacobian ); gsl_matrix_free( gsl_covar ); gsl_vector_free( gsl_opt_params );
	if( !debug && standalone && op->cd->calib_type == SIMPLE ) printf( "\n" );
	return( 1 );
}

int eigen( struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar )
{
	FILE *out;
	double phi, stddev_scale, gf;
	double *opt_params, *x_u, *x_d, *stddev, *jacobian;
	double aopt, copt, eopt, dopt, aic, bic, cic, kic, ln_det_v, ln_det_weight, sml, tt;
	int   debug, compute_covar, compute_jacobian;
	int   i, j, k, ier, debug_level, status, dof;
	double eps;
	char filename[200], buf[20];
	static double student_dist[34] = {12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042, 2.021, 2.000, 1.980, 1.960 };
	gsl_vector *gsl_opt_params = gsl_vector_alloc( op->pd->nOptParam );
	gsl_matrix *eigenvec = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam );
	gsl_vector *eigenval = gsl_vector_alloc( op->pd->nOptParam );
	gsl_eigen_symmv_workspace *eigenwork = gsl_eigen_symmv_alloc( op->pd->nOptParam );
	if( ( jacobian = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->od->nObs ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	compute_jacobian = compute_covar = 0;
	if( gsl_jacobian == NULL ) { gsl_jacobian = gsl_matrix_alloc( op->od->nObs, op->pd->nOptParam ); compute_jacobian = 1; }
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( ( x_u = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( ( x_d = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	if( ( stddev = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	debug = op->cd->debug;
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]];
	Transform( opt_params, op, opt_params );
	for( i = 0; i < op->pd->nOptParam; i++ )
		gsl_vector_set( gsl_opt_params, i, opt_params[i] );
	printf( "\nEigen analysis ...\n" );
	op->cd->pderiv = op->cd->oderiv = -1;
	if( compute_jacobian )
	{
		/*		op->pd->var_current_gsl = gsl_vector_alloc( op->pd->nOptParam );
				op->od->obs_current_gsl = gsl_vector_alloc(op->od->nObs);
				func_gsl_deriv_dx( gsl_opt_params, op, gsl_jacobian ); // Using GSL function
				gsl_vector_free( op->od->obs_current_gsl );
				gsl_vector_free( op->pd->var_current_gsl ); */
		func_gsl_dx( gsl_opt_params, op, gsl_jacobian ); // Compute Jacobian using forward difference
		func_dx( opt_params, NULL, op, jacobian ); // Compute Jacobian using forward difference
	}
	if( debug )
	{
		printf( "Analyzed state:\n" );
		debug_level = op->cd->fdebug; op->cd->fdebug = 3;
	}
	func_global( opt_params, op, op->od->res );
	if( debug ) op->cd->fdebug = debug_level;
	phi = op->phi;
	if( debug )
	{
		printf( "\nJacobian matrix\n" ); // Print Jacobian
		printf( "%-25s :", "Observations" );
		for( k = 0; k < op->od->nObs; k++ )
		{
			if( op->od->nObs < 30 || ( k < 10 || k > op->od->nObs - 10 ) ) printf( " %s", op->od->obs_id[k] );
			if( op->od->nObs >= 30 && k == 11 ) printf( " ..." );
		}
		printf( "\n" );
		for( k = i = 0; i < op->pd->nOptParam; i++ )
		{
			printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->od->nObs; j++ )
			{
				eps = gsl_matrix_get( gsl_jacobian, j, i );
				if( fabs( eps ) > 1e3 ) sprintf( buf, " %6.0e", eps );
				else sprintf( buf, " %6.2f", eps );
				if( op->od->nObs < 30 || ( j < 10 || j > op->od->nObs - 10 ) ) printf( " %s", buf );
				if( op->od->nObs >= 30 && j == 11 ) printf( " ..." );
			}
			printf( "\n" );
			/*
			printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->od->nObs; j++, k++ )
			{
				eps = jacobian[k];
				if( fabs( eps ) > 1e3 ) sprintf( buf, " %6.0e", eps );
				else sprintf( buf, " %6.2f", eps );
				if( op->od->nObs < 30 || ( j < 10 || j > op->od->nObs - 10 ) ) printf( " %s", buf );
				if( op->od->nObs > 30 && j == 11 ) printf( " ..." );
			}
			printf( "\n" );
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
	for( j = 0; j < op->od->nObs; j++ )
	{
		fprintf( out, "%-25s :", op->od->obs_id[j] );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fprintf( out, " %g", gsl_matrix_get( gsl_jacobian, j, i ) );
		fprintf( out, "\n" );
	}
	fclose( out );
	printf( "Jacobian matrix stored (%s)\n", filename );
	for( i = 0; i < op->pd->nOptParam; i++ )
		aopt = gsl_matrix_get( gsl_covar, i, i );
	if( op->cd->problem_type == EIGEN || op->cd->leigen )
	{
		if( gsl_covar == NULL ) { gsl_covar = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam ); compute_covar = 1; }
		if( compute_covar ) // Standalone eigen analysis
		{
			ier = gsl_multifit_covar( gsl_jacobian, 0.0, gsl_covar );
			if( ier != GSL_SUCCESS ) { printf( "Problem computing covariance matrix!\n" ); ier = 1; }
			else ier = 0;
		}
		if( debug )
		{
			printf( "\nCovariance matrix\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					printf( " %7.0e", gsl_matrix_get( gsl_covar, i, j ) );
				printf( "\n" );
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
		printf( "Covariance matrix stored (%s)\n", filename );
		for( ier = 0, i = 0; i < op->pd->nOptParam; i++ )
		{
			stddev[i] = sqrt( gsl_matrix_get( gsl_covar, i, i ) ); // compute standard deviations before the covariance matrix is destroyed by eigen functions
			if( stddev[i] < DBL_EPSILON ) ier = 1;
		}
		if( ier == 0 )
		{
			if( debug )
			{
				printf( "\nCorrelation matrix\n" );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					printf( "%-25s : ", op->pd->var_id[op->pd->var_index[i]] );
					for( j = 0; j < op->pd->nOptParam; j++ )
						printf( " %6.3f", gsl_matrix_get( gsl_covar, i, j ) / ( stddev[i] * stddev[j] ) );
					printf( "\n" );
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
			printf( "Correlation matrix stored (%s)\n", filename );
			// GSL_COVAR is destroyed during eigen computation
			gsl_eigen_symmv( gsl_covar, eigenval, eigenvec, eigenwork );
			if( debug )
			{
				printf( "\nEigenvectors (sorted by absolute values of eigenvalues)\n" );
				gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_ABS_ASC );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
					for( j = 0; j < op->pd->nOptParam; j++ )
						printf( " %6.3f", gsl_matrix_get( eigenvec, i, j ) );
					printf( "\n" );
				}
				printf( "%-25s :", "Eigenvalues" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					printf( " %6.0e", gsl_vector_get( eigenval, i ) );
				printf( "\n" );
				copt = fabs( gsl_vector_get( eigenval, op->pd->nOptParam - 1 ) ) / fabs( gsl_vector_get( eigenval, 0 ) );
				eopt = fabs( gsl_vector_get( eigenval, op->pd->nOptParam - 1 ) );
				dopt = 1;
				for( i = op->pd->nOptParam - 1; i >= 0; i-- )
					dopt *= fabs( gsl_vector_get( eigenval, i ) );
				printf( "\nEigenvectors (sorted by eigenvalues)\n" );
				gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_VAL_ASC );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
					for( j = 0; j < op->pd->nOptParam; j++ )
						printf( " %6.3f", gsl_matrix_get( eigenvec, i, j ) );
					printf( "\n" );
				}
				printf( "%-25s :", "Eigenvalues" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					printf( " %6.0e", gsl_vector_get( eigenval, i ) );
				printf( "\n" );
			}
			sprintf( filename, "%s", op->root );
			if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
			if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
			strcat( filename, ".eigen" );
			out = Fwrite( filename );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					fprintf( out, " %g", gsl_matrix_get( eigenvec, i, j ) );
				fprintf( out, "\n" );
			}
			fprintf( out, "%-25s :", "Eigenvalues" );
			for( i = 0; i < op->pd->nOptParam; i++ )
				fprintf( out, " %g", gsl_vector_get( eigenval, i ) );
			fprintf( out, "\n" );
			fclose( out );
			printf( "Eigen vactors and eigen values stored (%s)\n", filename );
		}
		else
			printf( "Correlation matrix and eigen vectors cannot be computed!\n" );
		dof = op->od->nObs - op->pd->nOptParam;
		stddev_scale = sqrt( phi / dof );
		gf = phi / dof;
		for( i = 0; i < op->od->nObs; i++ )
			if( op->od->obs_weight[i] > 0 )
				ln_det_weight += log( op->od->obs_weight[i] );
		ln_det_v = ln_det_weight + op->pd->nOptParam * log( gf );
		sml = ( double ) dof + ln_det_v + op->od->nObs * 1.837877;
		aic = sml + ( double ) 2 * op->pd->nOptParam;
		bic = sml + ( double ) op->pd->nOptParam * log( op->od->nObs );
		cic = sml + ( double ) 2 * op->pd->nOptParam * log( log( op->od->nObs ) );
		kic = sml + ( double ) op->pd->nOptParam * log( op->od->nObs * 0.159154943 ) - log( dopt );
		if( op->cd->problem_type == EIGEN || debug )
		{
			printf( "\nNumber of parameters           : %d\n", op->pd->nOptParam );
			printf( "Number of observations         : %d\n", op->od->nObs );
			printf( "Number of degrees of freedom   : %d\n", dof );
			printf( "Objective function             : %g\n", phi );
			printf( "Posterior measurement variance : %g\n", gf );
			printf( "\nOptimality metrics based on covariance matrix of observation errors:\n" );
			printf( "A-optimality (matrix trace)               : %g\n", aopt );
			printf( "C-optimality (matrix conditioning number) : %g\n", copt );
			printf( "E-optimality (matrix maximum eigenvalue)  : %g\n", eopt );
			printf( "D-optimality (matrix determinant)         : %g\n", dopt );
			printf( "\nDeterminant of covariance matrix of observation errors : %-15g ( ln(det S) = %g )\n", dopt, log( dopt ) );
			printf( "Determinant of observation weight matrix               : %-15g ( ln(det W) = %g )\n", exp( ln_det_weight ) , ln_det_weight );
			printf( "Determinant of covariance matrix of measurement errors : %-15g ( ln(det V) = %g )\n", exp( ln_det_v ), ln_det_v );
			printf( "\nLog likelihood function             : %g\n", -sml / 2 );
			printf( "Maximum likelihood                  : %g\n", sml );
			printf( "AIC (Akaike information criterion)  : %g\n", aic );
			printf( "BIC                                 : %g\n", bic );
			printf( "CIC                                 : %g\n", cic );
			printf( "KIC (Kashyap Information Criterion) : %g\n", kic );
		}
		if( dof < 0 ) tt = 1;
		else if( dof < 30 )  tt = student_dist[dof];
		else if( dof < 40 )  tt = student_dist[30] + ( dof - 30 ) * ( student_dist[31] - student_dist[30] ) / 10;
		else if( dof < 60 )  tt = student_dist[31] + ( dof - 40 ) * ( student_dist[32] - student_dist[31] ) / 20;
		else if( dof < 120 ) tt = student_dist[32] + ( dof - 60 ) * ( student_dist[33] - student_dist[32] ) / 60;
		else tt = student_dist[34];
		printf( "\nObtained fit is " );
		if( gf > 200 ) printf( "not very good (chi^2/dof = %g > 200)\n", gf );
		else printf( "relatively good (chi^2/dof = %g < 200)\n", gf );
		printf( "\nOptimized parameters:\n" );
		if( debug && op->cd->sintrans == 1 ) printf( "Transformed space (applied during optimization):\n" );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			if( op->cd->sintrans == 1 ) opt_params[i] = asin( sin( opt_params[i] ) );
			stddev[i] *= stddev_scale;
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
				printf( "%-40s : %12g stddev %12g (%12g - %12g)", op->pd->var_id[k], opt_params[i], stddev[i], x_d[i], x_u[i] );
				if( status ) printf( " Estimated ranges are constrained by prior uncertainty bounds\n" );
				else printf( "\n" );
			}
		}
		DeTransform( x_u, op, x_u );
		DeTransform( x_d, op, x_d );
		if( op->cd->sintrans == 1 )
		{
			if( debug ) printf( "Untransformed space:\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				k = op->pd->var_index[i];
				printf( "%-40s : ", op->pd->var_id[k] );
				if( op->pd->var_log[k] == 0 ) printf( "%12g stddev %12g (%12g - %12g)", op->pd->var[k], stddev[i], x_d[i], x_u[i] );
				else printf( "%12g stddev %12g (%12g - %12g)", pow( 10, op->pd->var[k] ), stddev[k], pow( 10, x_d[i] ), pow( 10, x_u[i] ) );
				status = 0;
				if( x_d[i] <= op->pd->var_min[k] ) { status = 1; x_d[i] = op->pd->var_min[k]; }
				if( x_u[i] >= op->pd->var_max[k] ) { status = 1; x_u[i] = op->pd->var_max[k]; }
				if( status ) printf( " Estimated ranges are constrained by prior uncertainty bounds\n" );
				else printf( "\n" );
			}
		}
	}
	free( opt_params ); free( stddev ); free( x_u ); free( x_d ); free( jacobian );
	gsl_vector_free( gsl_opt_params );
	gsl_matrix_free( eigenvec ); gsl_vector_free( eigenval ); gsl_eigen_symmv_free( eigenwork );
	if( compute_jacobian ) gsl_matrix_free( gsl_jacobian );
	if( compute_covar ) gsl_matrix_free( gsl_covar );
	return( 1 );
}

int infogap( struct opt_data *op )
{
	FILE *fl, *outfl;
	double *opt_params, of, maxof;
	char buf[80], filename[80];
	int i, j, k, n, npar, nrow, ncol, *nPreds, col;
	gsl_matrix *ig_mat; //! info gap matrix for sorting
	gsl_permutation *p;
	nPreds = &op->preds->nObs; // Set pointer to nObs for convenience
	if( op->cd->infile[0] == 0 ) { printf( "\nInfile must be specified for infogap run\n" ); exit( 0 );}
	nrow = count_lines( op->cd->infile ); nrow--; // Determine number of parameter sets in file
	npar = count_cols( op->cd->infile, 2 ); npar = npar - 2; // Determine number of parameter sets in file
	if( npar != op->pd->nOptParam ) { printf( "Number of optimization parameters in %s does not match input file\n", op->cd->infile ); exit( 0 ); } // Make sure MADS input file and PSSA file agree
	printf( "\n%s contains %d parameters and %d parameter sets\n", op->cd->infile, npar, nrow );
	ncol = npar + *nPreds + 1; // Number of columns for ig_mat = #pars + #preds + #ofs
	ig_mat = gsl_matrix_alloc( nrow, ncol );
	p = gsl_permutation_alloc( nrow );
	fl = fopen( op->cd->infile, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", op->cd->infile ); exit( 0 ); }
	printf( "Computing predictions for %s...", op->cd->infile );
	fflush( stdout );
	if( ( opt_params = ( double * ) malloc( npar * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); exit( 0 ); }
	fgets( buf, sizeof buf, fl ); // Skip header
	// Fill in ig_mat
	for( i = 0; i < nrow; i++ )
	{
		fscanf( fl, "%d %lf", &n, &of );
		gsl_matrix_set( ig_mat, i, *nPreds, of ); // Place of after predictions
		for( j = 0; j < npar; j++ )
		{
			fscanf( fl, "%lf", &opt_params[j] );
			col = *nPreds + 1 + j;
			gsl_matrix_set( ig_mat, i, col, opt_params[j] ); // Place after of
		}
		fscanf( fl, " \n" );
		func_global( opt_params, op, op->preds->res );
		for( j = 0; j < *nPreds; j++ )
		{
			gsl_matrix_set( ig_mat, i, j, op->preds->obs_current[j] ); // Place in first columns
		}
	}
	fclose( fl );
	for( k = 0; k < *nPreds; k++ )
	{
		gsl_vector_view column = gsl_matrix_column( ig_mat, k );
		gsl_sort_vector_index( p, &column.vector );
		// Print out ig_mat with headers
		sprintf( filename, "%s-pred%d.igap", op->root, k );
		outfl = fopen( filename , "w" );
		if( outfl == NULL ) { printf( "\nError opening %s\n", filename ); exit( 0 ); }
		fprintf( outfl, " %-12s", op->preds->obs_id[k] );
		fprintf( outfl, " OFmax OF" );
		for( i = 0; i < npar; i++ )
			fprintf( outfl, " (%-12s)", op->pd->var_id[i] );
		fprintf( outfl, "\n" );
		maxof = gsl_matrix_get( ig_mat, gsl_permutation_get( p, 0 ), *nPreds );
		for( i = 0; i < nrow; i++ )
		{
			if( maxof < gsl_matrix_get( ig_mat, gsl_permutation_get( p, i ), *nPreds ) )
				maxof = gsl_matrix_get( ig_mat, gsl_permutation_get( p, i ), *nPreds );
			fprintf( outfl, "%-12g", gsl_matrix_get( ig_mat, gsl_permutation_get( p, i ), k ) );
			fprintf( outfl, "%-12g", maxof );
			fprintf( outfl, "%-12g", gsl_matrix_get( ig_mat, gsl_permutation_get( p, i ), *nPreds ) );
			for( j = *nPreds + 1; j < ncol; j++ )
				fprintf( outfl, "%-12g", gsl_matrix_get( ig_mat, gsl_permutation_get( p, i ), j ) );
			fprintf( outfl, "\n" );
		}
		fclose( outfl );
		printf( "Done\n" );
		printf( "Results written to %s\n\n", filename );
	}
	gsl_matrix_free( ig_mat );
	return 1;
}

int postpua( struct opt_data *op )
{
	FILE *in, *out;
	double *opt_params, of;
	char buf[80], filename[80];
	int i, n;
	op->od = op->preds;
	if( op->cd->infile[0] == 0 ) { printf( "\nInfile (results file from abagus run) must be specified for postpua run\n" ); return( 0 );}
	in = Fread( op->cd->infile );
	// Create postpua output file
	sprintf( filename, "%s.pua", op->root );
	out = Fwrite( filename );
	printf( "\nComputing predictions for %s...", op->cd->infile );
	fflush( stdout );
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	fgets( buf, sizeof buf, in ); // Skip header
	fprintf( out, "Number       OF           " );
	for( i = 0; i < op->od->nObs; i++ )
		fprintf( out, " %-12s", op->od->obs_id[i] );
	fprintf( out, "\n" );
	while( fscanf( in, "%d %lf", &n, &of ) > 0 )
	{
		fprintf( out, "%-12d %-12lf ", n, of );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fscanf( in, "%lf", &opt_params[i] );
		fscanf( in, " \n" );
		func_global( opt_params, op, op->od->res );
		for( i = 0; i < op->od->nObs; i++ )
			fprintf( out, " %-12g", op->od->obs_current[i] );
		fprintf( out, "\n" );
	}
	fclose( in );
	fclose( out );
	printf( "Done.\n" );
	printf( "Results written to %s\n\n", filename );
	return 1;
}

int glue( struct opt_data *op )
{
	FILE *in, *out;
	double *phi, **preds, phi_temp, *percentile, *pred_temp, *sum;
	char buf[200], filename[80], pred_id[100][30];
	int num_lines = 0, j;
	gsl_matrix *glue_mat; // matrix for sorting predictions
	gsl_permutation *p1;
	// Open postpua output file
	if( op->cd->infile[0] == 0 ) { printf( "\nInfile (results file from postpua run) must be specified for glue run\n" ); return( 0 );}
	in = Fread( op->cd->infile );
	// Create glue output file
	sprintf( filename, "%s.glue", op->root );
	out = Fwrite( filename );
	fgets( buf, sizeof buf, in ); // Skip header
	// Count number of acceptable solutions
	while( fgets( buf, sizeof buf, in ) != NULL )
	{
		sscanf( buf, "%*d %lf ", &phi_temp );
		if( phi_temp <= op->cd->phi_cutoff ) num_lines++; // Count lines
	}
	printf( "\nNumber of solutions with phi <= %g: %d\n", op->cd->phi_cutoff, num_lines );
	printf( "\nPerforming GLUE analysis for %s...", op->cd->infile );
	fflush( stdout );
	// Read in data
	rewind( in );
	fscanf( in, "%*s %*s %[^\n]s", buf ); // Skip first part of header ("Number OF")
	int i = 0;
	// Read in names of predictions (e.g. combination of well names and times)
	while( sscanf( buf, " %s %[^\n]s", pred_id[i], buf ) > 1 ) { i++; }
	int num_preds = i + 1;
	// Allocate memory for phis and predictions
	if( ( phi = ( double * ) malloc( num_lines * sizeof( double ) ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 0 ); }
	//phi = ( double * ) malloc( num_lines * sizeof( double ) );
	preds = double_matrix( num_lines, num_preds );
	// Collect acceptable solutions
	glue_mat = gsl_matrix_alloc( num_lines, num_preds + 1 );
	int phi_index = num_preds; // phi_index indicates column of phis in glue_mat
	i = 0;
	//	printf( "\n\nAcceptable lines from %s:\n", op->cd->infile );
	while( fgets( buf, sizeof buf, in ) != NULL )
	{
		sscanf( buf, "%*d %lf %[^\n]s", &phi_temp, buf );
		if( phi_temp <= op->cd->phi_cutoff )
		{
			//			printf( "%lf %s\n", phi_temp, buf );
			phi[i] = phi_temp;
			gsl_matrix_set( glue_mat, i, phi_index, phi_temp ); // Place phi after predictions
			for( j = 0; j < num_preds; j++ )
			{
				sscanf( buf, " %lf %[^\n]s", &preds[i][j], buf );
				gsl_matrix_set( glue_mat, i, j, preds[i][j] ); // Place phi after predictions
			}
			i++;
		}
	}
	fclose( in );
	// printf( "\nglue_mat:\n" );
	// gsl_matrix_fprintf( stdout, glue_mat, "%g" );
	// Calculate weighted percentile of each phi; note: low phis imply high percentile
	p1 = gsl_permutation_alloc( num_lines );
	percentile = ( double * ) malloc( num_lines * sizeof( double ) );
	pred_temp = ( double * ) malloc( num_lines * sizeof( double ) );
	sum = ( double * ) malloc( num_lines * sizeof( double ) );
	double p05, p95; // 5th and 95th percentiles
	printf( "\n\nprediction p05 p95\n" );
	int count;
	for( i = 0; i < num_preds; i++ )
	{
		gsl_vector_view column = gsl_matrix_column( glue_mat, i );
		gsl_sort_vector_index( p1, &column.vector );
		sum[0] = gsl_matrix_get( glue_mat, gsl_permutation_get( p1, 0 ), num_preds );
		pred_temp[0] = gsl_matrix_get( glue_mat, gsl_permutation_get( p1, 0 ), i );
		// printf( "\nSample sum prediction:\n" );
		// printf( "0 %g %g\n", sum[0], pred_temp[0] );
		// Collect summation of weights and ordered predictions
		for( j = 1; j < num_lines; j++ )
		{
			sum[j] = sum[j - 1] + gsl_matrix_get( glue_mat, gsl_permutation_get( p1, j ), num_preds );
			pred_temp[j] = gsl_matrix_get( glue_mat, gsl_permutation_get( p1, j ), i );
			// printf( "%d %g %g\n", j, sum[j], pred_temp[j] );
		}
		// printf( "\nno prediction percentile:\n" );
		for( j = 0; j < num_lines; j++ ) { percentile[j] = ( 1.0 / sum[num_lines - 1] ) * ( sum[j] - pred_temp[j] / 2.0 ); /*printf( "%d %g %g\n", j+1, pred_temp[j], percentile[j]);*/ }
		if( percentile[0] > 0.05 ) p05 = pred_temp[0];
		else if( percentile[num_lines - 1] < 0.05 ) p05 = pred_temp[num_lines - 1];
		else
		{
			count = 0;
			for( j = 1; j < num_lines; j++ ) { if( percentile[j] < 0.05 ) count++; else {break;} }
			p05 = pred_temp[j - 1] + ( ( 0.05 - percentile[j - 1] ) / ( percentile[j] - percentile[j - 1] ) ) * ( pred_temp[j] - pred_temp[j - 1] );
		}
		//		printf( "\n%d\n", j );
		//		printf( "\n%g %g %g %g %g\n", pred_temp[j], pred_temp[j-1], percentile[j], percentile[j-1], p05 );
		if( percentile[0] > 0.95 ) p95 = pred_temp[0];
		else if( percentile[num_lines - 1] < 0.95 ) p95 = pred_temp[num_lines - 1];
		else
		{
			count = 0;
			for( j = 1; j < num_lines; j++ ) { if( percentile[j] < 0.95 ) count++; else {break;}  }
			p95 = pred_temp[j - 1] + ( ( 0.95 - percentile[j - 1] ) / ( percentile[j] - percentile[j - 1] ) ) * ( pred_temp[j] - pred_temp[j - 1] );
		}
		//		printf( "\n%d\n", j );
		//		printf( "\n%g %g %g %g %g\n", pred_temp[j], pred_temp[j-1], percentile[j], percentile[j-1], p95 );
		printf( "%d %g %g\n", i + 1, p05, p95 );
		//printf( "%g ", gsl_interp_eval( pred_interp, pred_temp, percentile, 0.95, accelerator ) );
		//printf( " %g\n", gsl_interp_eval( pred_interp, pred_temp, percentile, 0.05, accelerator ) );
	}
	printf( "\n" );
	fclose( out );
	printf( "Done.\n" );
	printf( "Results written to %s\n\n", filename );
	gsl_matrix_free( glue_mat );
	return 1;
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

void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op, int debug )
{
	if( debug ) printf( "%s\n", op->cd->smp_method );
	if( nreal == 1 || strncasecmp( op->cd->smp_method, "random", 6 ) == 0 )
	{
		if( debug )
		{
			printf( "Pure random sampling method ... " );
			fflush( stdout );
		}
		smp_random( npar, nreal, seed, var_lhs );
	}
	else if( ( nreal <= 500 && op->cd->smp_method[0] == 0 ) || strncasecmp( op->cd->smp_method, "idlhs", 5 ) == 0 )
	{
		if( debug )
		{
			printf( "Improved Distributed LHS method " );
			if( strncasecmp( op->cd->smp_method, "idlhs", 5 ) != 0 ) printf( "(number of realizations < 500) " );
			printf( "... " );
			fflush( stdout );
		}
		lhs_imp_dist( npar, nreal, 5, seed, var_lhs );
	}
	else if( ( nreal > 500 && op->cd->smp_method[0] == 0 ) || strncasecmp( op->cd->smp_method, "lhs", 3 ) == 0 )
	{
		if( debug )
		{
			printf( "Standard LHS method " );
			if( strncasecmp( op->cd->smp_method, "lhs", 3 ) != 0 ) printf( "(number of realizations > 500) " );
			printf( "... " );
			fflush( stdout );
		}
		lhs_random( npar, nreal, seed, var_lhs );
	}
}

void print_results( struct opt_data *op, int verbosity )
{
	int i, j, k, success, success_all;
	double c, err;
	success_all = 1;
	if( verbosity > 0 ) printf( "Model parameters:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->pd->var_log[k] == 0 ) printf( "%s %g\n", op->pd->var_id[k], op->pd->var[k] );
		else printf( "%s %g\n", op->pd->var_id[k], pow( 10, op->pd->var[k] ) );
	}
	if( verbosity == 0 ) return;
	if( op->cd->solution_type != TEST && op->od->nObs > 0 )
	{
		printf( "\nModel predictions:\n" );
		if( op->cd->solution_type == EXTERNAL )
			for( i = 0; i < op->od->nObs; i++ )
			{
				if( op->od->obs_weight[i] == 0 ) continue;
				c = op->od->obs_current[i];
				err = op->od->obs_target[i] - c;
				if( c < op->od->obs_min[i] || c > op->od->obs_max[i] ) { success_all = 0; success = 0; }
				else success = 1;
				if( op->od->nObs < 50 || ( i < 20 || i > op->od->nObs - 20 ) )
					printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, op->od->obs_min[i], op->od->obs_max[i] );
				if( op->od->nObs > 50 && i == 21 ) printf( "...\n" );
			}
		else
		{
			for( k = 0, i = 0; i < op->wd->nW; i++ )
				for( j = 0; j < op->wd->nWellObs[i]; j++ )
				{
					if( op->wd->obs_weight[i][j] == 0 ) continue;
					c = op->od->obs_current[k++];
					err = op->wd->obs_target[i][j] - c;
					if( c < op->wd->obs_min[i][j] || c > op->wd->obs_max[i][j] ) { success_all = 0; success = 0; }
					else success = 1;
					printf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], success, op->wd->obs_min[i][j], op->wd->obs_max[i][j] );
				}
		}
	}
	else
		for( i = 0; i < op->pd->nOptParam; i++ )
			if( fabs( op->pd->var[i] - op->pd->var_truth[i] ) > op->cd->parerror ) success_all = 0;
	op->success = success_all;
	printf( "Objective function: %g Success: %d \n", op->phi, op->success );
	if( op->cd->check_success > 0 && op->cd->obserror < 0 && op->cd->parerror < 0 )
	{
		if( success_all ) printf( "SUCCESS: All the model predictions are within calibration ranges!\n" );
		else printf( "At least one of the model predictions is outside calibration ranges!\n" );
	}
	if( op->cd->check_success > 0 && op->cd->obserror > 0 )
	{
		if( success_all ) printf( "SUCCESS: All the model predictions are within a predefined absolute error %g!\n", op->cd->obserror );
		else printf( "At least one of the model predictions has an absolute error greater than %g!\n", op->cd->obserror );
	}
	if( op->cd->check_success > 0 && op->cd->parerror > 0 )
	{
		if( success_all ) printf( "SUCCESS: All the estimated model parameters have an absolute error from the true parameters less than %g!\n", op->cd->parerror );
		else printf( "At least one of the estimated model parameters has an absolute error from the true parameters greater than %g!\n", op->cd->parerror );
	}
}

void save_results( char *label, struct opt_data *op, struct grid_data *gd )
{
	FILE *out, *out2;
	int i, j, k, success, success_all;
	double c, err;
	char filename[255], filename2[255], f[255];
	success_all = 1;
	sprintf( filename, "%s", op->root );
	if( label[0] != 0 ) sprintf( filename, "%s.%s", filename, label );
	if( op->counter > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->counter );
	strcpy( f, filename );
	strcat( filename, ".results" );
	out = Fwrite( filename );
	fprintf( out, "Model parameters:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->pd->var_log[k] == 0 ) fprintf( out, "%s %g\n", op->pd->var_id[k], op->pd->var[k] );
		else fprintf( out, "%s %g\n", op->pd->var_id[k], pow( 10, op->pd->var[k] ) );
	}
	if( op->cd->solution_type != TEST && op->od->nObs > 0 )
	{
		fprintf( out, "\nModel predictions:\n" );
		strcpy( filename, f );
		strcat( filename, ".residuals" );
		out2 = Fwrite( filename );
		if( op->cd->solution_type == EXTERNAL )
			for( i = 0; i < op->od->nObs; i++ )
			{
				if( op->od->obs_weight[i] != 0 )
				{
					c = op->od->obs_current[i];
					err = op->od->obs_target[i] - c;
					if( c < op->od->obs_min[i] || c > op->od->obs_max[i] ) { success_all = 0; success = 0; }
					else success = 1;
				}
				else success = 0;
				fprintf( out, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, op->od->obs_min[i], op->od->obs_max[i] );
				fprintf( out2, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], success, op->od->obs_min[i], op->od->obs_max[i] );
			}
		else
		{
			for( k = 0, i = 0; i < op->wd->nW; i++ )
				for( j = 0; j < op->wd->nWellObs[i]; j++ )
				{
					if( op->wd->obs_weight[i][j] != 0 )
					{
						c = op->od->obs_current[k++];
						err = op->wd->obs_target[i][j] - c;
						if( c < op->wd->obs_min[i][j] || c > op->wd->obs_max[i][j] ) { success_all = 0; success = 0; }
						else success = 1;
					}
					else success = 0;
					fprintf( out, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], success, op->wd->obs_min[i][j], op->wd->obs_max[i][j] );
					fprintf( out2, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], success, op->wd->obs_min[i][j], op->wd->obs_max[i][j] );
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
	op->success = success_all;
	fprintf( out, "Objective function: %g Success: %d \n", op->phi, op->success );
	if( op->cd->check_success > 0 && op->cd->obserror < 0 && op->cd->parerror < 0 )
	{
		if( success_all ) fprintf( out, "SUCCESS: All the model predictions are within calibration ranges!\n" );
		else fprintf( out, "At least one of the model predictions is outside calibration ranges!\n" );
	}
	if( op->cd->check_success > 0 && op->cd->obserror > 0 )
	{
		if( success_all ) fprintf( out, "SUCCESS: All the model predictions are within a predefined absolute error %g!\n", op->cd->obserror );
		else fprintf( out, "At least one of the model predictions has an absolute error greater than %g!\n", op->cd->obserror );
	}
	if( op->cd->check_success > 0 && op->cd->parerror > 0 )
	{
		if( success_all ) fprintf( out, "SUCCESS: All the estimated model parameters have an absolute error from the true parameters less than %g!\n", op->cd->parerror );
		else fprintf( out, "At least one of the estimated model parameters has an absolute error from the true parameters greater than %g!\n", op->cd->parerror );
	}
	fprintf( out, "Number of function evaluations = %d\n", op->cd->neval );
	if( op->cd->seed > 0 ) fprintf( out, "Seed = %d\n", op->cd->seed_init );
	fclose( out );
	if( gd->min_t > 0 && op->cd->solution_type != TEST )
	{
		printf( "\nCompute breakthrough curves at all the wells ..." );
		fflush( stdout );
		sprintf( filename, "%s.btc", f );
		sprintf( filename2, "%s.btc-peak", f );
		compute_btc2( filename, filename2, op );
		//			compute_btc( filename, &op, &gd );
	}
	if( gd->time > 0 && op->cd->solution_type != TEST )
	{
		printf( "\nCompute spatial distribution of predictions at t = %g ...\n", gd->time );
		fflush( stdout );
		sprintf( filename, "%s.vtk", f );
		compute_grid( filename, op->cd, gd );
	}
}

char *timestamp()
{
	time_t raw_time;
	struct tm *ptr_ts;
	char *datetime;
	datetime = ( char * ) malloc( 10 * sizeof( char ) );
	time( &raw_time );
	ptr_ts = localtime( &raw_time );
	//	printf( "%s\n", asctime( ptr_ts ) );
	sprintf( datetime, "%02d:%02d:%02d", ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
	return( datetime );
}

char *datestamp()
{
	time_t raw_time;
	struct tm *ptr_ts;
	char *datetime;
	datetime = ( char * ) malloc( 16 * sizeof( char ) );
	time( &raw_time );
	ptr_ts = localtime( &raw_time );
	//	printf( "%s\n", asctime( ptr_ts ) );
	sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
	return( datetime );
}

int sort_int( const void *x, const void *y )
{
	return ( *( int * )x - * ( int * )y );
}

int sort_double( const void *x, const void *y )
{
	return ( *( double * )x - * ( double * )y );
}

int igrnd( struct opt_data *op )
{
	int i, k, m, q1, q2, npar, status, success, phi_global, success_global, success_all, count, debug_level, predict = 0, compare, bad_data = 0, no_memory = 0;
	int *eval_success, *eval_total;
	unsigned long neval_total, njac_total;
	double c, err, phi, phi_min, *orig_params, *opt_params,
		   *opt_params_min, *opt_params_max, *opt_params_avg,
		   *sel_params_min, *sel_params_max, *sel_params_avg,
		   *var_lhs;
	char filename[255], buf[255];
	int ( *optimize_func )( struct opt_data * op ); // function pointer to optimization function (LM or PSO)
	FILE *out, *out2;
	char ESC = 27; // Escape

	strcpy( op->label, "igrnd" );
	if( ( orig_params = ( double * ) malloc( op->pd->nParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( no_memory ) { printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); return( 0 ); }
	if( op->cd->nreal > 1 )
	{
		if( ( opt_params_min = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
		if( ( opt_params_max = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
		if( ( opt_params_avg = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
		if( no_memory ) { printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); return( 0 ); }
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
			if( no_memory ) { printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); return( 0 ); }
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				sel_params_min[i] = HUGE_VAL;
				sel_params_max[i] = opt_params_avg[i] = 0;
			}
		}
	}
	printf( "\nSEQUENTIAL RUNS using random initial guesses for model parameters (realizations = %d):\n", op->cd->nreal );
	if( op->pd->nFlgParam != 0 ) { printf( "Only flagged parameters are randomized\n" ); npar = op->pd->nFlgParam; }
	else if( op->pd->nOptParam != 0 ) { printf( "No flagged parameters; all optimizable parameters are randomized\n" ); npar = op->pd->nOptParam; }
	else { printf( "No flagged or optimizable parameters; all parameters are randomized\n" ); npar = op->pd->nParam; }
	if( ( var_lhs = ( double * ) malloc( npar * op->cd->nreal * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( ( eval_success = ( int * ) malloc( op->cd->nreal * sizeof( int ) ) ) == NULL ) no_memory = 1;
	if( ( eval_total = ( int * ) malloc( op->cd->nreal * sizeof( int ) ) ) == NULL ) no_memory = 1;
	if( no_memory )
	{
		printf( "Not enough memory!\n" );
		sprintf( buf, "rm -f %s.running", op->root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
		system( buf );
		return( 0 );
	}
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; printf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { printf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else printf( "Current seed: %d\n", op->cd->seed );
	printf( "Random sampling (variables %d; realizations %d) using ", npar, op->cd->nreal );
	sampling( npar, op->cd->nreal, &op->cd->seed, var_lhs, op, 1 );
	printf( "done.\n" );
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
		printf( "Random sampling set saved in %s.igrnd_set\n", op->root );
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
		printf( "WARNING: No parameters to optimize! Forward runs performed instead (ie Monte Carlo analysis)\n" );
	phi_min = HUGE_VAL;
	phi_global = success_global = neval_total = njac_total = 0;
	if( op->cd->ireal != 0 ) k = op->cd->ireal - 1; // applied if execution of a specific realization is requested (ncase)
	else k = 0;
	for( count = k; count < op->cd->nreal; count++ )
	{
		op->cd->neval = op->cd->njac = 0;
		fprintf( out, "%d : init var", count + 1 );
		printf( "\nRandom set #%d: ", count + 1 );
		if( op->cd->debug || op->cd->nreal == 1 ) printf( "\n" );
		op->counter = count + 1;
		for( k = i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] == 2 || ( op->pd->var_opt[i] == 1 && op->pd->nFlgParam == 0 ) )
			{
				op->pd->var[i] = var_lhs[k + count * npar] * op->pd->var_range[i] + op->pd->var_min[i];
				if( op->pd->var_log[i] )
				{
					if( op->cd->debug || op->cd->nreal == 1 ) printf( "%s %.15g\n", op->pd->var_id[i], pow( 10, op->pd->var[i] ) );
					fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
				}
				else
				{
					if( op->cd->debug || op->cd->nreal == 1 ) printf( "%s %.15g\n", op->pd->var_id[i], op->pd->var[i] );
					fprintf( out, " %.15g", op->pd->var[i] );
				}
				k++;
			}
			else op->pd->var[i] = orig_params[i];
		if( op->pd->nOptParam > 0 )
		{
			status = optimize_func( op ); // Optimize
			if( status == 0 ) { sprintf( buf, "rm -f %s.running", op->root ); system( buf ); return( 0 ); }
		}
		else
		{
			for( i = 0; i < op->pd->nParam; i++ )
				op->pd->var[i] = var_lhs[i + count * npar] * op->pd->var_range[i] + op->pd->var_min[i];
			if( op->cd->mdebug ) { printf( "Forward run ... \n" ); debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			func_global( op->pd->var, op, op->od->res ); // op->pd->var is a dummy variable because op->pd->nOptParam == 0
			if( op->cd->mdebug ) op->cd->fdebug = debug_level;
		}
		if( op->cd->debug > 1 )
			printf( "\n" );
		else
		{
			printf( "Evaluations: %d ", op->cd->neval );
			if( op->cd->njac > 0 ) printf( "Jacobians: %d ", op->cd->njac );
			printf( "Objective function: %g Success: %d", op->phi, op->success );
		}
		if( op->cd->debug || op->cd->nreal == 1 )
		{
			printf( "\n" );
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
		if( op->phi < phi_min )
		{
			phi_min = op->phi;
			for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]];
			for( i = 0; i < op->od->nObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
		}
		if( op->cd->pdebug || op->cd->ldebug ) printf( "\n" ); // extra new line if the optimization process is debugged
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
		if( op->success && op->cd->nreal > 1 && op->cd->odebug > 1 ) save_results( "igrnd", op, op->gd );
		if( op->f_ofe != NULL ) { fclose( op->f_ofe ); op->f_ofe = NULL; }
		if( op->cd->ireal != 0 ) break;
	}
	op->counter = 0;
	free( var_lhs );
	op->cd->neval = neval_total; // provide the correct number of total evaluations
	op->cd->njac = njac_total; // provide the correct number of total evaluations
	printf( "\nTotal number of evaluations = %lu\n", neval_total );
	printf( "Total number of jacobians = %lu\n", njac_total );
	op->phi = phi_min; // get the best phi
	for( i = 0; i < op->pd->nOptParam; i++ ) opt_params[i] = op->pd->var[op->pd->var_index[i]] = op->pd->var_current[i] = op->pd->var_best[i]; // get the best estimate
	for( i = 0; i < op->od->nObs; i++ ) op->od->obs_current[i] = op->od->obs_best[i] ; // get the best observations
	fprintf( out, "Minimum objective function: %g\n", phi_min );
	printf( "Minimum objective function: %g\n", phi_min );
	if( op->cd->debug )
	{
		printf( "Repeat the run producing the best results ...\n" );
		debug_level = op->cd->fdebug; op->cd->fdebug = 1;
		Transform( opt_params, &op, opt_params );
		func_global( opt_params, &op, op->od->res );
		op->cd->fdebug = debug_level;
	}
	if( op->cd->nreal > 1 )
	{
		if( op->cd->phi_cutoff > DBL_EPSILON )
		{
			if( phi_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions below predefined OF cutoff %g!\n", op->cd->nreal, op->cd->phi_cutoff );
			else printf( "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
		}
		if( op->cd->obsrange > DBL_EPSILON || op->cd->obserror > DBL_EPSILON )
		{
			if( success_global == 0 ) printf( "None of the %d sequential calibration runs produced successful calibration ranges!\n", op->cd->nreal );
			else printf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		}
		if( op->cd->parerror > DBL_EPSILON )
		{
			if( success_global == 0 ) printf( "None of the %d sequential calibration runs produced acceptable model parameters!\n", op->cd->nreal );
			else printf( "Number of the sequential calibration runs producing acceptable model parameters = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		}
		qsort( eval_success, success_global, sizeof( int ), sort_int );
		qsort( eval_total, op->cd->nreal, sizeof( int ), sort_int );
		q1 = ( double ) success_global / 4 - 0.25;
		m = ( double ) success_global / 2 - 0.5;
		q2 = ( double ) success_global * 3 / 4 - 0.25;
		if( success_global > 0 )
		{
			printf( "Statistics of successful number of evaluations : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_success[0], eval_success[q1], ESC, eval_success[m], ESC, eval_success[q2], eval_success[success_global - 1], success_global );
			fprintf( out2, "Statistics of successful number of evaluations : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_success[0], eval_success[q1], ESC, eval_success[m], ESC, eval_success[q2], eval_success[success_global - 1], success_global );
		}
		q1 = ( double ) op->cd->nreal / 4 - 0.25;
		m = ( double ) op->cd->nreal / 2 - 0.5;
		q2 = ( double ) op->cd->nreal * 3 / 4 - 0.25;
		printf( "Statistics of total number of evaluations      : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_total[0], eval_total[q1], ESC, eval_total[m], ESC, eval_total[q2], eval_total[op->cd->nreal - 1], op->cd->nreal );
		fprintf( out2, "Statistics of total number of evaluations      : %d - %d %c[1m%d%c[0m %d - %d : %d\n", eval_total[0], eval_total[q1], ESC, eval_total[m], ESC, eval_total[q2], eval_total[op->cd->nreal - 1], op->cd->nreal );
		printf( "Statistics of all the model parameter estimates:\n" );
		fprintf( out2, "Statistics of all the model parameter estimates:\n" );
		for( i = 0; i < op->pd->nOptParam; i++ ) // Posterior parameter statistics for all simulations
		{
			printf( "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], opt_params_avg[i] / op->cd->nreal, opt_params_min[i], opt_params_max[i] );
			fprintf( out2, "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], opt_params_avg[i] / op->cd->nreal, opt_params_min[i], opt_params_max[i] );
		}
		if( success_global > 0 || phi_global > 0 )
		{
			printf( "Statistics of all the successful model parameter estimates:\n" );
			fprintf( out2, "Statistics of all the successful model parameter estimates:\n" );
			k = success_global + phi_global;
			for( i = 0; i < op->pd->nOptParam; i++ ) // Posterior parameter statistics for all simulations
			{
				printf( "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], sel_params_avg[i] / k, sel_params_min[i], sel_params_max[i] );
				fprintf( out2, "%-35s : average %12g min %12g max %12g\n", op->pd->var_id[op->pd->var_index[i]], sel_params_avg[i] / k, sel_params_min[i], sel_params_max[i] );
			}
		}
	}
	fprintf( out, "Total number of evaluations = %lu\n", neval_total );
	if( op->cd->nreal > 1 )
	{
		if( op->cd->phi_cutoff > DBL_EPSILON ) fprintf( out, "Number of the sequential calibration runs producing predictions below predefined OF cutoff (%g) = %d (out of %d; success ratio %g)\n", op->cd->phi_cutoff, phi_global, op->cd->nreal, ( double ) phi_global / op->cd->nreal );
		if( op->cd->obsrange > DBL_EPSILON || op->cd->obserror > DBL_EPSILON )
			fprintf( out, "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
		if( op->cd->parerror > DBL_EPSILON )
			fprintf( out, "Number of the sequential calibration runs producing acceptable model parameters = %d (out of %d; success ratio %g)\n", success_global, op->cd->nreal, ( double ) success_global / op->cd->nreal );
	}
	fprintf( out2, "OF min: %g\n", phi_min );
	fprintf( out2, "Total number of evaluations: %lu\n", neval_total );
	fprintf( out2, "Success rate %g\n", ( double ) success_global / op->cd->nreal );
	fclose( out ); fclose( out2 );
	printf( "Results are saved in %s.igrnd.results and %s.igrnd-opt=%s_eval=%d_real=%d\n", op->root, op->root, op->cd->opt_method, op->cd->maxeval, op->cd->nreal );
	printf( "\nFinal results:\n" );
	print_results( op, 1 );
	free( opt_params ); free( eval_success ); free( eval_total );
	free( opt_params_min ); free( opt_params_max ); free( opt_params_avg );
	if( op->cd->phi_cutoff > DBL_EPSILON || op->cd->check_success ) { free( sel_params_min ); free( sel_params_max ); free( sel_params_avg ); }
	save_results( "", op, op->gd );
}

