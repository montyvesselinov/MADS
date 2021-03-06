// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
// http://www.ees.lanl.gov/staff/monty/codes/mads
// http://gitlab.com/monty/mads
//
// Licensing: GPLv3: http://www.gnu.org/licenses/gpl-3.0.html
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

#include <math.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <stdbool.h>
#include "mads.h"
#ifdef MATHEVAL
#include <matheval.h>
#endif
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
int func_extrn( double *x, void *data, double *f, double *o );
int func_extrn_write( int ieval, double *x, void *data );
int func_extrn_read( int ieval, void *data, double *f, double *o );
int func_extrn_check_read( int ieval, void *data );
int func_intrn( double *x, void *data, double *f, double *o );
void func_levmar( double *x, double *f, double *o, int m, int n, void *data );
void func_dx_levmar( double *x, double *f, double *o, double *jacobian, int m, int n, void *data ); // Jacobian order: obs / param
int func_dx( double *x, double *f, double *o, void *data, double *jacobian ); // Jacobian order: param / obs
int func_dx_set( double *x, double *f, double *o, void *data, double *jacobian );
int func_set( int n_sub, double *var_mat[], double *phi, double *f[], double *o[], FILE *out, struct opt_data *op );
int func_set_omp( int n_sub, double *var_mat[], double *phi, double *f[], double *o[], FILE *out, struct opt_data *op );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );

/* Functions elsewhere */
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );
int ins_obs( int nobs, char **obs_id, double *obs, int *check, char *fn_in_t, char *fn_in_d, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
double test_problems( int D, int function, double *x, int nObs, double *o );
double point_source( double x, double y, double z, double t, void *params );
double point_source_triangle_time( double x, double y, double z, double t, void *params );
double rectangle_source( double x, double y, double z, double t, void *params );
double gaussian_source_2d( double x, double y, double z, double t, void *params );
double gaussian_source_3d( double x, double y, double z, double t, void *params );
double rectangle_source_vz( double x, double y, double z, double t, void *params );
double box_source( double x, double y, double z, double t, void *params );
double box_source_levy_dispersion( double x, double y, double z, double t, void *params );
double box_source_sym_levy_dispersion( double x, double y, double z, double t, void *params );
int create_mprun_dir( char *dir );
int delete_mprun_dir( char *dir );
int mprun( int nJob, void *data );
int mprunall( int nJob, void *data, double *var_mat[], double *phi, double *f[], double *o[] ); // Write and execute all the files in parallel; Reading does not work
int Ftestread( char *filename );
time_t Fdatetime_t( char *filename, int debug );

int func_extrn( double *x, void *data, double *f, double *o )
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000];
	double c, t, w, min, max, dx, err, phi = 0.0;
	int i, k, success, success_all = 1, bad_data = 0;
	if( p->cd->parallel_type ) // Parallel execution of a serial job to archive all the intermediate results
	{
		func_extrn_write( p->cd->neval + 1, x, data );
		if( mprun( 1, data ) < 0 ) // Perform one (1) run in parallel
		{
			tprintf( "ERROR: there is a problem with the parallel execution!\n" );
			mads_quits( p->root );
		}
		bad_data = func_extrn_read( p->cd->neval, data, f, o ); // p->cd->eval was already incremented in mprun
		if( bad_data ) mads_quits( p->root );
		return GSL_SUCCESS; // DONE
	}
	p->cd->neval++;
	DeTransform( x, p, p->pd->var_current );
	if( p->cd->fdebug >= 3 ) tprintf( "Optimized model parameters (%d):\n", p->pd->nOptParam );
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) p->cd->var[k] = pow( 10, p->pd->var_current[i] );
		else p->cd->var[k] = p->pd->var_current[i];
		if( p->cd->fdebug >= 3 )
			tprintf( "%s %.12g log %d\n", p->pd->var_name[k], p->cd->var[k], p->pd->var_log[k] );
	}
	if( p->pd->nExpParam > 0 )
	{
		if( p->cd->fdebug >= 3 ) tprintf( "Tied model parameters (%d):\n", p->pd->nExpParam );
#ifdef MATHEVAL
		for( i = 0; i < p->pd->nExpParam; i++ )
		{
			k = p->pd->param_expressions_index[i];
			p->cd->var[k] = evaluator_evaluate( p->pd->param_expression[i], p->pd->nParam, p->pd->var_name, p->cd->var );
			if( p->cd->fdebug >= 3 ) tprintf( "%s = %s = %.12g\n", p->pd->var_name[k], evaluator_get_string( p->pd->param_expression[i] ), p->cd->var[k] );
		}
#else
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		mads_quits( p->root );
#endif
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nFixParam == 0 ) tprintf( "NO fixed parameters.\n" );
		else
		{
			tprintf( "Fixed model parameters (%d):\n", p->pd->nFixParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->analysis_type == PPSD ) )
					tprintf( "%s %.12g\n", p->pd->var_name[i], p->cd->var[i] );
		}
	}
	if( p->cd->fdebug >= 4 )
	{
		tprintf( "Objective function: " );
		switch( p->cd->objfunc_type )
		{
			case SSR: tprintf( "sum of squared residuals" ); break;
			case SSDR: tprintf( "sum of squared discrepancies and squared residuals" ); break;
			case SSDA: tprintf( "sum of squared discrepancies and absolute residuals" ); break;
			case SSDX: tprintf( "sum of squared discrepancies with increased to get within the bounds" ); break;
			case SSD0: tprintf( "sum of squared discrepancies" ); break;
			default: tprintf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
		}
		tprintf( "\n" );
	}
	for( i = 0; i < p->ed->ntpl; i++ )
		if( par_tpl( p->pd->nParam, p->pd->var_id, p->cd->var, p->ed->fn_tpl[i], p->ed->fn_out[i], p->cd->tpldebug ) == -1 )
			mads_quits( p->root );
	if( p->cd->tpldebug || p->cd->insdebug ) tprintf( "\nDelete the expected output files before execution ... \n" );
	if( p->ed->nins > 0 )
	{
		for( i = 0; i < p->ed->nins; i++ )
			remove( p->ed->fn_obs[i] );
	}
	if( p->cd->tpldebug || p->cd->insdebug ) tprintf( "Execute external model \'%s\' ... ", p->ed->cmdline );
	sprintf( buf, "%s \"%s %s\"", BASH, p->ed->cmdline, quiet_string );
	system( buf );
	if( p->cd->tpldebug || p->cd->insdebug ) tprintf( "done!\n" );
	int *obs_count;
	obs_count = ( int * ) malloc( p->od->nObs * sizeof( int ) );
	for( i = 0; i < p->od->nObs; i++ ) obs_count[i] = 0;
	for( i = 0; i < p->ed->nins; i++ )
		if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, obs_count, p->ed->fn_ins[i], p->ed->fn_obs[i], p->cd->insdebug ) == -1 )
			bad_data = 1;
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( obs_count[i] == 0 )
		{
			tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( obs_count[i] > 1 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug )
				tprintf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], obs_count[i] );
			p->od->obs_current[i] /= obs_count[i];
		}
	}
	free( obs_count );
	if( bad_data ) mads_quits( p->root );
#ifdef MATHEVAL
	// for( k = 0; k < p->rd->regul_nMap; k++ ) { tprintf( "%s %g\n", p->rd->regul_map_id[k], p->rd->regul_map_val[k] ); }
	for( i = p->od->nObs; i < p->od->nTObs; i++ )
		p->od->obs_current[i] = evaluator_evaluate( p->rd->regul_expression[i - p->od->nObs], p->rd->regul_nMap, p->rd->regul_map_id, p->rd->regul_map_val );
#else
	if( p->od->nObs < p->od->nTObs )
	{
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		mads_quits( p->root );
	}
#endif
	if( p->cd->fdebug >= 2 ) tprintf( "\nModel predictions:\n" );
	for( i = 0; i < p->od->nTObs; i++ )
	{
		c  = p->od->obs_current[i];
		if( o != NULL ) o[i] = c;
		t = p->od->obs_target[i];
		w = p->od->obs_weight[i];
		min = p->od->obs_min[i];
		max = p->od->obs_max[i];
		if( p->od->obs_log[i] == 0 )
		{
			err = t - c;
		}
		else
		{
			if( c < DBL_EPSILON ) c = DBL_EPSILON;
			if( t < DBL_EPSILON ) t = DBL_EPSILON;
			err = log10( t ) - log10( c );
		}
		if( p->cd->objfunc_type != SSR )
		{
			if( p->cd->objfunc_type == SSDA )
			{
				err = sqrt( fabs( err ) );
				if( t < c ) err *= -1;
			}
			else if( p->cd->objfunc_type != SSDR ) err = 0; // SSD0 & SSDX
			if( p->cd->objfunc_type == SSDX ) { dx = max - min; if( p->cd->obsdomain > DBL_EPSILON && p->cd->obsdomain < dx ) dx = p->cd->obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
			if( c < min ) err += min - c;
			else if( c > max ) err += max - c;
			// tprintf( "%g %g %g %g\n", err, c, min - c, max - c );
			if( p->cd->objfunc_type == SSDX ) { min = p->od->obs_min[i]; max = p->od->obs_max[i]; }
		}
		f[i] = err * w;
		if( p->cd->compute_phi ) phi += f[i] * f[i];
		if( p->cd->obserror > DBL_EPSILON )
		{
			if( fabs( c - t ) > p->cd->obserror ) { success = 0; if( w > DBL_EPSILON ) success_all = 0; }
			else success = 1;
		}
		else // if( p->cd->obsrange )
		{
			if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success = 0; if( w > DBL_EPSILON ) success_all = 0; }
			else success = 1;
		}
		if( p->cd->fdebug >= 2 )
		{
			if( p->od->nTObs < 50 || ( i < 20 || i > p->od->nTObs - 20 ) )
				tprintf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", p->od->obs_id[i], t, c, err, err * w, success, min, max );
			if( p->od->nTObs > 50 && i == 21 ) tprintf( "...\n" );
			if( !p->cd->compute_phi ) phi += f[i] * f[i];
		}
		if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
	}
	if( isnan( phi ) ) success_all = 0;
	p->success = success_all; // Just in case
	if( p->cd->fdebug >= 2 ) tprintf( "Objective function %g\n", phi );
	if( p->cd->compute_phi ) { p->phi = phi; p->success = success_all; }
	if( p->cd->phi_cutoff > DBL_EPSILON )
	{
		if( phi < p->cd->phi_cutoff )
		{
			p->success = 1;
			if( p->cd->fdebug ) tprintf( "SUCCESS: OF is below predefined cutoff value (%g < %g; func_intrn)!\n", phi, p->cd->phi_cutoff );
		}
		else p->success = 0;
	}
	if( p->cd->check_success ) { p->success = success_all; p->phi = phi; if( p->cd->fdebug && success_all ) tprintf( "SUCCESS: Model results are within the predefined ranges (func_extrn)!\n" ); }
	return GSL_SUCCESS;
}

int func_extrn_write( int ieval, double *x, void *data ) // Create a series of input files for parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	int i, k, bad_data = 0;
	double *opt_params, *cur_params;
	if( p->cd->sintrans )
	{
		if( ( opt_params = ( double * ) malloc( p->pd->nOptParam * sizeof( double ) ) ) == NULL ) printf( "Not enough memory!\n" ); // needed for parallel runs
		DeTransform( x, p, opt_params );
	}
	else
		opt_params = x;
	if( ( cur_params = ( double * ) malloc( p->pd->nParam * sizeof( double ) ) ) == NULL ) printf( "Not enough memory!\n" ); // needed for parallel runs
	if( p->cd->fdebug >= 3 ) tprintf( "Optimized model parameters (%d; model run = %d):\n", p->pd->nOptParam, ieval );
	for( i = 0; i < p->pd->nParam; i++ )
		cur_params[i] = p->cd->var[i];
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) cur_params[k] = pow( 10, opt_params[i] );
		else cur_params[k] = opt_params[i];
		if( p->cd->fdebug >= 3 )
			tprintf( "%s %.12g log %d\n", p->pd->var_name[k], cur_params[k], p->pd->var_log[k] );
	}
	if( p->pd->nExpParam > 0 )
	{
		if( p->cd->fdebug >= 3 ) tprintf( "Tied model parameters (%d):\n", p->pd->nExpParam );
#ifdef MATHEVAL
		for( i = 0; i < p->pd->nExpParam; i++ )
		{
			k = p->pd->param_expressions_index[i];
			cur_params[k] = evaluator_evaluate( p->pd->param_expression[i], p->pd->nParam, p->pd->var_name, cur_params );
			if( p->cd->fdebug >= 3 ) tprintf( "%s = %s = %.12g\n", p->pd->var_name[k], evaluator_get_string( p->pd->param_expression[i] ), cur_params[k] );
		}
#else
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		mads_quits( p->root );
#endif
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nFixParam == 0 ) tprintf( "NO fixed parameters.\n" );
		else
		{
			tprintf( "Fixed model parameters (%d):\n", p->pd->nFixParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->analysis_type == PPSD ) )
					tprintf( "%s %.12g\n", p->pd->var_name[i], p->cd->var[i] );
		}
	}
	// Delete expected output files in the root directory to prevent the creation of links to these files in the "child" directories
	if( p->ed->nins > 0 )
	{
		for( i = 0; i < p->ed->nins; i++ )
			remove( p->ed->fn_obs[i] );
		if( p->cd->pardebug > 3 ) tprintf( "Delete the expected output files before execution ...\n" );
	}
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval ); // Name of directory for parallel runs
	create_mprun_dir( dir ); // Create the child directory for parallel runs with link to the files in the working root directory
	for( i = 0; i < p->ed->ntpl; i++ ) // Create all the model input files
	{
		sprintf( buf, "../%s/%s", dir, p->ed->fn_out[i] );
		if( par_tpl( p->pd->nParam, p->pd->var_id, cur_params, p->ed->fn_tpl[i], buf, p->cd->tpldebug ) == -1 )
			bad_data = 1;
	}
	if( p->cd->sintrans )
		free( opt_params );
	free( cur_params );
	if( bad_data ) return( 0 );
	if( p->cd->restart && !p->cd->bin_restart ) // Update model input files in zip restart files
	{
		// sprintf( buf, "%s \'WAIT_TIME=0; until [ $WAIT_TIME -eq 15 ]; do zip -q -u %s ", BASH, p->cd->restart_container ); // Archive input files
		sprintf( buf, "%s \'zip -q -u %s ", BASH, p->cd->restart_container ); // Archive input files
		for( i = 0; i < p->ed->ntpl; i++ )
			sprintf( &buf[( int ) strlen( buf )], "../%s/%s ", dir, p->ed->fn_out[i] );
		// strcat( buf, "; E=$?; echo ERROR=$E; if [[ $E -eq 0 || $E -eq 12 ]]; then break; fi; sleep $(( WAIT_TIME++ )); done" );
		if( p->cd->pardebug <= 3 || quiet ) strcat( buf, QUIET );
		strcat( buf, "\'" );
		if( p->cd->pardebug > 4 ) tprintf( "Execute: %s", buf );
		system( buf );
		if( p->cd->pardebug > 3 ) tprintf( "Input files for parallel run #%d are archived!\n", ieval );
	}
	if( p->cd->restart == 0 ) // Do not delete if restart is attempted
	{
		if( p->ed->nins > 0 )
		{
			if( p->cd->pardebug > 3 ) tprintf( "Delete the expected output files before execution (\'%s\')\n", buf );
			for( i = 0; i < p->ed->nins; i++ )
			{
				sprintf( buf, "../%s/%s", dir, p->ed->fn_obs[i] ); // Delete expected output files in the hosts directories
				remove( buf );
			}
		}
	}
	else if( !p->cd->bin_restart ) // Just in case; the restart file should have been already extracted
	{
		sprintf( buf, "%s \"unzip -u -: %s ", BASH, p->cd->restart_container ); // Archive input files
		for( i = 0; i < p->ed->nins; i++ )
			sprintf( &buf[( int ) strlen( buf )], "../%s/%s ", dir, p->ed->fn_obs[i] );
		if( p->cd->pardebug <= 3 || quiet ) strcat( buf, QUIET );
		strcat( buf, "\"" );
		if( p->cd->pardebug > 2 ) tprintf( "Extract the expected output files before execution (\'%s\')\n", buf );
		system( buf );
	}
	return GSL_SUCCESS;
}

int func_extrn_exec_serial( int ieval, void *data ) // Execute a series of external runs in serial (for testing only)
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	p->cd->neval++;
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval ); // Name of directory for parallel runs
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) tprintf( "\nWorking directory: ../%s\n", dir );
	if( p->cd->pardebug > 9 )
	{
		sprintf( buf, "%s \"cd ../%s; ls -altr\"", BASH, dir ); // Check directory content
		system( buf );
	}
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) tprintf( "Execute external model \'%s\' ... ", p->ed->cmdline );
	sprintf( buf, "%s \" ( cd ../%s; %s ) %s \"", BASH, dir, p->ed->cmdline, quiet_string );
	system( buf );
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) tprintf( "done!\n" );
	return GSL_SUCCESS;
}

int func_extrn_check_read( int ieval, void *data ) // Check a series of output files after parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	int i, bad_data;
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval );
	if( p->cd->pardebug > 9 )
	{
		sprintf( buf, "%s \"cd ../%s; ls -altr\"", BASH, dir ); // Check directory content
		system( buf );
	}
	int *obs_count;
	obs_count = ( int * ) malloc( p->od->nObs * sizeof( int ) );
	for( i = 0; i < p->od->nObs; i++ ) obs_count[i] = 0;
	for( i = 0; i < p->ed->nins; i++ )
	{
		sprintf( buf, "../%s/%s", dir, p->ed->fn_obs[i] );
		if( Ftestread( buf ) != 0 ) { if( p->cd->pardebug > 2 ) tprintf( "File %s cannot be opened to read.\n", buf ); return( 0 ); }
		else if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, obs_count, p->ed->fn_ins[i], buf, 0 ) == -1 )
			return( 0 );
	}
	bad_data = 0;
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( obs_count[i] == 0 )
		{
			if( p->cd->pardebug ) tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( obs_count[i] > 1 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug || p->cd->pardebug )
				tprintf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], obs_count[i] );
			p->od->obs_current[i] /= obs_count[i];
		}
	}
	free( obs_count );
	if( bad_data ) return( 0 );
	if( ( p->cd->time_infile - Fdatetime_t( buf, 0 ) ) > 0 && p->cd->restart != -1 )
	{
		tprintf( "RESTART ERROR: File %s is older than the MADS input file.\n", buf );
		return( 0 );
	}
	return( 1 );
}

int func_extrn_read( int ieval, void *data, double *f, double *o ) // Read a series of output files after parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	double c, t, w, min, max, dx, err, phi = 0.0;
	int i, success, success_all = 1, bad_data;
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval );
	int *obs_count;
	bool bin_read = false;
	if( p->cd->restart && p->cd->bin_restart )
	{
		FILE *infileb;
		sprintf( buf, "%s/%020d.obs", p->cd->restart_container, ieval ); // Archive model output
		if( ( infileb = fopen( buf, "rb" ) ) != NULL )
		{
			int obj_read = fread( ( void * ) f, sizeof( f[0] ), p->od->nTObs, infileb );
			fclose( infileb );
			if( obj_read != p->od->nTObs ) tprintf( "RESTART ERROR: Binary file %s cannot be applied to read model predictions; data mismatch!\n", buf );
			else
			{
				if( p->cd->pardebug > 1 )
					tprintf( "RESTART: Results (model residuals) from parallel run #%d are read from a file in directory %s!\n", ieval, p->cd->restart_container );
				if( p->cd->pardebug > 4 )
					for( i = 0; i < p->od->nTObs; i++ )
						tprintf( "%-27s: binary observations %15.12g\n", p->od->obs_id[i], f[i] );
				delete_mprun_dir( dir );
				bin_read = true;
			}
		}
	}
	if( !bin_read )
	{
		obs_count = ( int * ) malloc( p->od->nObs * sizeof( int ) );
		for( i = 0; i < p->od->nObs; i++ ) obs_count[i] = 0;
		bad_data = 0;
		for( i = 0; i < p->ed->nins; i++ )
		{
			sprintf( buf, "../%s/%s", dir, p->ed->fn_obs[i] );
			if( ins_obs( p->od->nObs, p->od->obs_id, f, obs_count, p->ed->fn_ins[i], buf, p->cd->insdebug ) == -1 )
				bad_data = 1;
		}
		for( i = 0; i < p->od->nObs; i++ )
		{
			if( obs_count[i] == 0 )
			{
				tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
				bad_data = 1;
			}
			else if( obs_count[i] > 1 )
			{
				if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug )
					tprintf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], obs_count[i] );
				f[i] /= obs_count[i];
			}
		}
		free( obs_count );
		if( bad_data ) return( bad_data );
	}
	if( p->cd->restart )
	{
		if( p->cd->bin_restart && !bin_read )
		{
			FILE *outfileb;
			sprintf( buf, "%s/%020d.par", p->cd->restart_container, ieval ); // Archive model inputs
			if( ( outfileb = fopen( buf, "wb" ) ) == NULL ) tprintf( "RESTART ERROR: Binary file %s cannot be opened to save problem information!\n", buf );
			else
			{
				fwrite( ( void * ) p->pd->var, sizeof( p->pd->var[0] ), p->pd->nParam, outfileb );
				fclose( outfileb );
			}
			sprintf( buf, "%s/%020d.obs", p->cd->restart_container, ieval ); // Archive model output
			if( ( outfileb = fopen( buf, "wb" ) ) == NULL ) tprintf( "RESTART ERROR: Binary file %s cannot be opened to save problem information!\n", buf );
			else
			{
				fwrite( ( void * ) f, sizeof( f[0] ), p->od->nTObs, outfileb );
				fclose( outfileb );
			}
			if( p->cd->pardebug > 3 ) tprintf( "RESTART: Results (model predictions) from parallel run #%d are archived in directory %s!\n", ieval, p->cd->restart_container );
		}
		else if( !p->cd->bin_restart )
		{
			// sprintf( buf, "%s \'WAIT_TIME=0; until [ $WAIT_TIME -eq 15 ]; do zip -q -u %s ", BASH, p->cd->restart_container ); // Archive input files
			sprintf( buf, "%s \'zip -q -u %s ", BASH, p->cd->restart_container ); // Archive input files
			for( i = 0; i < p->ed->nins; i++ )
				sprintf( &buf[strlen( buf )], "../%s/%s ", dir, p->ed->fn_obs[i] );
			// strcat( buf, "; E=$?; echo ERROR=$E; if [[ $E -eq 0 || $E -eq 12 ]]; then break; fi; sleep $(( WAIT_TIME++ )); done" );
			if( p->cd->pardebug <= 3 || quiet ) strcat( buf, QUIET );
			strcat( buf, "\'" );
			if( p->cd->pardebug > 4 ) tprintf( "Execute: %s", buf );
			system( buf );
			if( p->cd->pardebug > 3 ) tprintf( "RESTART: Results from parallel run #%d are archived in zip file %s!\n", ieval, p->cd->restart_container );
		}
	}
	delete_mprun_dir( dir ); // Delete directory for parallel runs
#ifdef MATHEVAL
	for( i = p->od->nObs; i < p->od->nTObs; i++ )
		p->od->obs_current[i] = evaluator_evaluate( p->rd->regul_expression[i - p->od->nObs], p->rd->regul_nMap, p->rd->regul_map_id, p->rd->regul_map_val );
#else
	if( p->od->nObs < p->od->nTObs )
	{
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		mads_quits( p->root );
	}
#endif
	if( p->cd->fdebug >= 2 ) tprintf( "\nModel predictions (model run = %d):\n", ieval );
	for( i = 0; i < p->od->nTObs; i++ )
	{
		c = f[i];
		if( o != NULL ) o[i] = c;
		p->od->obs_current[i] = c;
		t = p->od->obs_target[i];
		w = p->od->obs_weight[i];
		min = p->od->obs_min[i];
		max = p->od->obs_max[i];
		if( p->od->obs_log[i] == 0 )
		{
			err = t - c;
		}
		else
		{
			if( c < DBL_EPSILON ) c = DBL_EPSILON;
			if( t < DBL_EPSILON ) t = DBL_EPSILON;
			err = log10( t ) - log10( c );
		}
		if( p->cd->objfunc_type != SSR )
		{
			if( p->cd->objfunc_type == SSDA )
			{
				err = sqrt( fabs( err ) );
				if( t < c ) err *= -1;
			}
			else if( p->cd->objfunc_type != SSDR ) err = 0; // SSD0 & SSDX
			if( p->cd->objfunc_type == SSDX ) { dx = max - min; if( p->cd->obsdomain > DBL_EPSILON && p->cd->obsdomain < dx ) dx = p->cd->obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
			if( c < min ) err += min - c;
			else if( c > max ) err += max - c;
			// tprintf( "%g %g %g %g\n", err, c, min - c, max - c );
			if( p->cd->objfunc_type == SSDX ) { min = p->od->obs_min[i]; max = p->od->obs_max[i]; }
		}
		f[i] = err * w;
		if( p->cd->compute_phi ) phi += f[i] * f[i];
		if( p->cd->obserror > DBL_EPSILON )
		{
			if( fabs( c - t ) > p->cd->obserror ) { success = 0; if( w > DBL_EPSILON ) success_all = 0; }
			else success = 1;
		}
		else // if( p->cd->obsrange )
		{
			if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success = 0; if( w > DBL_EPSILON ) success_all = 0; }
			else success = 1;
		}
		if( p->cd->fdebug >= 2 )
		{
			if( p->od->nTObs < 50 || ( i < 20 || i > p->od->nTObs - 20 ) )
				tprintf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", p->od->obs_id[i], t, c, err, err * w, success, min, max );
			if( p->od->nTObs > 50 && i == 21 ) tprintf( "...\n" );
			if( !p->cd->compute_phi ) phi += f[i] * f[i];
		}
		if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
	}
	if( p->cd->restart && p->cd->bin_restart && !bin_read )
	{
		FILE *outfileb;
		sprintf( buf, "%s/%020d.res", p->cd->restart_container, ieval ); // Archive model output
		if( ( outfileb = fopen( buf, "wb" ) ) == NULL ) tprintf( "RESTART: Binary file %s cannot be opened to save problem information!\n", buf );
		else
		{
			fwrite( ( void * ) f, sizeof( f[0] ), p->od->nTObs, outfileb );
			fclose( outfileb );
		}
		if( p->cd->pardebug > 3 ) tprintf( "RESTART: Results (model residuals) from parallel run #%d are archived in directory %s!\n", ieval, p->cd->restart_container );
	}
	if( isnan( phi ) ) success_all = 0;
	p->success = success_all; // Just in case
	if( p->cd->fdebug >= 2 ) tprintf( "Objective function %g\n", phi );
	if( p->cd->compute_phi ) { p->phi = phi; p->success = success_all; }
	if( p->cd->phi_cutoff > DBL_EPSILON )
	{
		if( phi < p->cd->phi_cutoff )
		{
			p->success = 1;
			if( p->cd->fdebug ) tprintf( "SUCCESS: OF is below predefined cutoff value (%g < %g; func_intrn)!\n", phi, p->cd->phi_cutoff );
		}
		else p->success = 0;
	}
	if( p->cd->check_success ) { p->success = success_all; p->phi = phi; if( p->cd->fdebug && success_all ) tprintf( "SUCCESS: Model results are within the predefined ranges (func_extrn)!\n" ); }
	return GSL_SUCCESS;
}

int func_intrn( double *x, void *data, double *f, double *o ) /* forward run for LM */
{
	int i, j, k, p1, p2, l, s, success, success_all = 1;
	double c = 0, t, w, min, max, c1, c2, err, phi, dx, dy, dz, x1, y1, z1, dist;
	struct opt_data *p = ( struct opt_data * )data;
	char filename[255];
	dx = dy = dz = phi = 0.0;
	p->cd->neval++;
	if( p->cd->odebug && p->f_ofe == NULL )
	{
		if( p->counter > 0 ) sprintf( filename, "%s.%08d.ofe", p->root, p->counter );
		else sprintf( filename, "%s.ofe", p->root );
		if( p->cd->nretries > 1 && p->cd->retry_ind > 1 ) p->f_ofe = fopen( filename, "a" );
		else p->f_ofe = fopen( filename, "w" );
	}
	DeTransform( x, p, p->pd->var_current );
	if( p->cd->fdebug >= 3 ) tprintf( "Optimized model parameters (%d):\n", p->pd->nOptParam );
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) p->cd->var[k] = pow( 10, p->pd->var_current[i] );
		else p->cd->var[k] = p->pd->var_current[i];
		if( p->cd->fdebug >= 3 ) tprintf( "%s %.12g log %d\n", p->pd->var_name[k], p->cd->var[k], p->pd->var_log[k] );
	}
	if( p->pd->nExpParam > 0 )
	{
		if( p->cd->fdebug >= 3 ) tprintf( "Tied model parameters (%d):\n", p->pd->nExpParam );
#ifdef MATHEVAL
		for( i = 0; i < p->pd->nExpParam; i++ )
		{
			k = p->pd->param_expressions_index[i];
			p->cd->var[k] = evaluator_evaluate( p->pd->param_expression[i], p->pd->nParam, p->pd->var_id, p->cd->var );
			if( p->cd->fdebug >= 3 ) tprintf( "%s = %s = %.12g\n", p->pd->var_name[k], evaluator_get_string( p->pd->param_expression[i] ), p->cd->var[k] );
		}
#else
		tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
		mads_quits( p->root );
#endif
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nFixParam == 0 ) tprintf( "NO fixed parameters.\n" );
		else
		{
			tprintf( "Fixed model parameters (%d):\n", p->pd->nFixParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->analysis_type == PPSD ) )
					tprintf( "%s %.12g\n", p->pd->var_name[i], p->cd->var[i] );
		}
	}
	if( p->cd->fdebug >= 4 )
	{
		tprintf( "Objective function: " );
		if( p->cd->solution_type[0] != TEST )
		{
			switch( p->cd->objfunc_type )
			{
				case SSR: tprintf( "sum of squared residuals" ); break;
				case SSDR: tprintf( "sum of squared discrepancies and squared residuals" ); break;
				case SSDA: tprintf( "sum of squared discrepancies and absolute residuals" ); break;
				case SSDX: tprintf( "sum of squared discrepancies with increased to get within the bounds" ); break;
				case SSD0: tprintf( "sum of squared discrepancies" ); break;
				default: tprintf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
			}
			if( p->cd->compute_phi ) tprintf( " --- computed\n" );
		}
		else
			tprintf( "test function\n" );
	}
	// p->cd->compute_phi = 1;
	// if( p->cd->compute_phi ) tprintf( " --- computed!!!!\n" );
	if( p->cd->solution_type[0] == TEST )
	{
		p->phi = phi = test_problems( p->pd->nOptParam, p->cd->test_func, p->cd->var, p->od->nTObs, f );
		if( p->cd->check_success && p->cd->obsrange ) success_all = 0;
		if( p->cd->check_success && p->cd->parerror > DBL_EPSILON )
		{
			success_all = 1;
			for( k = 0; k < p->pd->nOptParam; k++ )
				if( fabs( p->cd->var[k] - p->pd->var_truth[k] ) > p->cd->parerror ) success_all = 0;
			if( p->cd->fdebug >= 4 ) tprintf( "Test OF %g Success %d\n", phi, success_all );
		}
		if( p->cd->test_func >= 40 && p->cd->test_func < 100 )
		{
			phi = 0;
			if( p->cd->check_success && p->cd->obserror > DBL_EPSILON ) success_all = 1;
			else success_all = 0;
			for( k = 0; k < p->od->nTObs; k++ )
			{
				c = f[k];
				f[k] = err = p->od->obs_target[k] - c;
				if( p->cd->check_success && p->cd->obserror > DBL_EPSILON )
				{
					if( c < p->od->obs_min[k] || c > p->od->obs_max[k] ) { success_all = success = 0; }
					else success = 1;
				}
				// tprintf( "e %g %g\n", err, p->od->obs_target[k] );
				phi += err * err;
			}
			p->phi = phi;
		}
		else
		{
			if( p->cd->check_success && p->cd->obserror > DBL_EPSILON ) success_all = 0;
			if( p->cd->fdebug >= 5 )
			{
				tprintf( "Test OF %d %g\n", success_all, phi );
				c = 0;
				for( k = 0; k < p->od->nTObs; k++ )
				{
					tprintf( "%s %g\n", p->od->obs_id[k], f[k] );
					c += f[k] * f[k];
				}
				tprintf( "Test OF ? %g\n", c );
			}
		}
	}
	else
	{
		for( p1 = p->cd->num_source_params * p->cd->num_sources, p2 = p->cd->num_source_params; p1 < p->pd->nAnalParam; p1++, p2++ )
			p->ad->var[p2] = p->cd->var[p1];
		if( p->cd->disp_tied && p->cd->disp_scaled == 0 ) // Tied dispersivities
		{
			if( p->cd->fdebug >= 5 )
			{
				tprintf( "Tied AY %.12g = %.12g / %.12g\n", p->ad->var[AX] / p->ad->var[AY], p->ad->var[AX], p->ad->var[AY] );
				tprintf( "Tied AZ %.12g = %.12g / %.12g\n", ( p->ad->var[AX] / p->ad->var[AY] ) / p->ad->var[AZ], p->ad->var[AX] / p->ad->var[AY], p->ad->var[AZ] );
			}
			p->ad->var[AY] = p->ad->var[AX] / p->ad->var[AY];
			p->ad->var[AZ] = p->ad->var[AY] / p->ad->var[AZ];
		}
		if( p->cd->fdebug > 10 )
		{
			for( i = 0; i < p->pd->nAnalParam; i++ )
			{
				tprintf( "%s %g\n", p->pd->var_name[i], p->cd->var[i] );
			}
			for( k = 0; k < p->rd->regul_nMap; k++ ) { tprintf( "%s %g\n", p->rd->regul_map_id[k], p->rd->regul_map_val[k] ); }
		}
		if( p->cd->fdebug >= 2 ) tprintf( "\nModel predictions:\n" );
		for( k = 0; k < p->od->nTObs; k++ ) // for computational efficiency performed only for observations with weight > DBL_EPSILON
		{
			if( p->cd->oderiv != -1 ) { k = p->cd->oderiv; }
			if( k >= p->od->nObs )
			{
				// regularization term
#ifdef MATHEVAL
				c = evaluator_evaluate( p->rd->regul_expression[k - p->od->nObs], p->rd->regul_nMap, p->rd->regul_map_id, p->rd->regul_map_val );
#else
				tprintf( "ERROR: MathEval is not installed; expressions cannot be evaluated. MADS Quits!\n" );
				mads_quits( p->root );
#endif
			}
			else
			{
				// model predicted calibration target
				i = p->od->obs_well_index[k];
				j = p->od->obs_time_index[k];
				c1 = c2 = p->cd->c_background;
				for( s = 0; s < p->cd->num_sources; s++ )
				{
					l = p->cd->num_source_params * ( s + 1 );
					p2 = 0;
					for( p1 = p->cd->num_source_params * s; p1 < l; p1++, p2++ )
						p->ad->var[p2] = p->cd->var[p1];
					if( p->cd->disp_scaled ) // Scaled dispersivities
					{
						dx = p->ad->var[AX]; dy = p->ad->var[AY]; dz = p->ad->var[AZ];
						x1 = p->wd->x[i] - p->ad->var[SOURCE_X];
						y1 = p->wd->y[i] - p->ad->var[SOURCE_Y];
						z1 = ( p->wd->z1[i] + p->wd->z2[i] ) / 2 - p->ad->var[SOURCE_Z];
						dist = sqrt( x1 * x1 + y1 * y1 + z1 * z1 );
						if( p->cd->fdebug >= 5 ) tprintf( "Scaled AX %.12g = %.12g * %.12g\n", p->ad->var[AX] * dist, p->ad->var[AX], dist );
						p->ad->var[AX] *= dist;
						if( p->cd->disp_scaled > 1 && !p->cd->disp_tied ) { p->ad->var[AY] *= dist; p->ad->var[AZ] *= dist; }
						else if( p->cd->disp_tied ) { p->ad->var[AY] = p->ad->var[AX] / p->ad->var[AY]; p->ad->var[AZ] = p->ad->var[AY] / p->ad->var[AZ]; };
						if( p->cd->fdebug >= 5 )
						{
							if( p->cd->disp_scaled > 1 && !p->cd->disp_tied ) tprintf( "Transverse dispersivities are scaled!\n" );
							else if( p->cd->disp_tied ) tprintf( "Transverse dispersivities are tied!\n" );
							else tprintf( "Transverse dispersivities are neither tied nor scaled!\n" );
							tprintf( "AY %.12g\n", p->ad->var[AY] );
							tprintf( "AZ %.12g\n", p->ad->var[AZ] );
						}
					}
					if( fabs( p->ad->var[TSCALE_DISP] - 1 ) < COMPARE_EPSILON || fabs( p->ad->var[TSCALE_DISP] ) < COMPARE_EPSILON ) p->ad->scaling_dispersion = 0;
					else p->ad->scaling_dispersion = 1;
					p->ad->time_step = p->cd->time_step;
					// TODO merge func_intrn and func_solver; func_intrn is called by methods (PE, UQ, ..); func_solver is called by forward and grid solvers
					if( p->cd->obs_int == 1 ) // TODO add other integration models ...
					{
						switch( p->cd->solution_type[s] )
						{
							case POINT:
								c1 += point_source( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case POINT_TRIANGLE_TIME:
								c1 += point_source_triangle_time( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case PLANE:
								c1 += rectangle_source( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case PLANE3D:
								c1 += rectangle_source_vz( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case GAUSSIAN2D:
								c1 += gaussian_source_2d( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case GAUSSIAN3D:
								c1 += gaussian_source_3d( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							default:
							case BOX:
								if( p->cd->levy == FULL_LEVY ) c1 += box_source_levy_dispersion( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								else if( p->cd->levy == SYM_LEVY ) c1 += box_source_sym_levy_dispersion( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								else c1 += box_source( p->wd->x[i], p->wd->y[i], ( p->wd->z1[i] + p->wd->z2[i] ) / 2, p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
						}
					}
					else
					{
						switch( p->cd->solution_type[s] )
						{
							case POINT:
								c1 += point_source( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								c2 += point_source( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case POINT_TRIANGLE_TIME:
								c1 += point_source_triangle_time( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								c2 += point_source_triangle_time( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case PLANE:
								c1 += rectangle_source( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								c2 += rectangle_source( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case PLANE3D:
								c1 += rectangle_source_vz( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								c2 += rectangle_source_vz( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case GAUSSIAN2D:
								c1 += gaussian_source_2d( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								c2 += gaussian_source_2d( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							case GAUSSIAN3D:
								c1 += gaussian_source_3d( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								c2 += gaussian_source_3d( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								break;
							default:
							case BOX:
								if( p->cd->levy == FULL_LEVY )
								{
									c1 += box_source_levy_dispersion( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
									c2 += box_source_levy_dispersion( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								}
								else if( p->cd->levy == SYM_LEVY )
								{
									c1 += box_source_sym_levy_dispersion( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
									c2 += box_source_sym_levy_dispersion( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								}
								else
								{
									c1 += box_source( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
									c2 += box_source( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
								}
								break;
						}
					}
					if( p->cd->disp_scaled ) // Scaled dispersivities
					{
						p->ad->var[AX] = dx; p->ad->var[AY] = dy; p->ad->var[AZ] = dz;
					}
				}
				if( p->cd->obs_int == 1 ) c = c1;
				else c = ( c1 + c2 ) / 2;
			}
			p->od->obs_current[k] = c;
			if( o != NULL ) o[k] = c;
			t = p->od->obs_target[k];
			w = p->od->obs_weight[k];
			min = p->od->obs_min[k];
			max = p->od->obs_max[k];
			if( p->od->obs_log[k] == 0 )
			{
				err = t - c;
			}
			else
			{
				if( c < DBL_EPSILON ) c = DBL_EPSILON;
				if( t < DBL_EPSILON ) t = DBL_EPSILON;
				err = log10( t ) - log10( c );
			}
			if( p->cd->objfunc_type != SSR )
			{
				if( p->cd->objfunc_type == SSDA )
				{
					err = sqrt( fabs( err ) );
					if( c < t ) err *= -1;
				}
				else if( p->cd->objfunc_type != SSDR ) err = 0; // SSD0 & SSDX
				if( p->cd->objfunc_type == SSDX ) { dx = max - min; if( p->cd->obsdomain > DBL_EPSILON && p->cd->obsdomain < dx ) dx = p->cd->obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
				if( c < min ) err += min - c;
				else if( c > max ) err += max - c;
				// tprintf( "%g %g %g %g\n", err, c, min - c, max - c );
				if( p->cd->objfunc_type == SSDX ) { min = p->od->obs_min[k]; max = p->od->obs_max[k]; }
			}
			f[k] = err * w;
			if( p->cd->compute_phi ) phi += f[k] * f[k];
			if( p->cd->obserror > DBL_EPSILON )
			{
				if( fabs( c - t ) > p->cd->obserror ) { success = success_all = 0; } // weight should be > DBL_EPSILON by default; if( w > DBL_EPSILON ) is not needed
				else success = 1;
			}
			else // if( p->cd->obsrange )
			{
				if( min - c > COMPARE_EPSILON || c - max > COMPARE_EPSILON ) { success = success_all = 0; } // weight should be > DBL_EPSILON by default; if( w > DBL_EPSILON ) is not needed
				else success = 1;
			}
			if( p->cd->fdebug >= 2 )
			{
				if( p->od->nTObs < 50 || ( i < 20 || i > p->od->nTObs - 20 ) )
					tprintf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", p->od->obs_id[k], t, c, err, err * w, success, min , max );
				if( p->od->nTObs > 50 && i == 21 ) tprintf( "...\n" );
				if( !p->cd->compute_phi ) phi += f[k] * f[k];
			}
			if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
		}
	}
	if( isnan( phi ) ) success_all = 0;
	p->success = success_all; // Just in case
	if( fabs( p->cd->obsstep ) > DBL_EPSILON && p->success ) // INFOGAP analysis
	{
		for( i = 0; i < p->preds->nTObs; i++ )
		{
			k = p->preds->obs_index[i];
			if( p->cd->obsstep >  DBL_EPSILON && p->preds->obs_best[i] < p->od->obs_current[k] ) p->preds->obs_best[i] = p->od->obs_current[k];
			if( p->cd->obsstep < -DBL_EPSILON && p->preds->obs_best[i] > p->od->obs_current[k] ) p->preds->obs_best[i] = p->od->obs_current[k];
		}
	}
	if( p->cd->odebug )
	{
		fprintf( p->f_ofe, "%d %g", p->cd->neval, phi ); // Print current best
		for( i = 0; i < p->pd->nOptParam; i++ )
			fprintf( p->f_ofe, " %g", p->cd->var[p->pd->var_index[i]] );
		fprintf( p->f_ofe, "\n" );
		fflush( p->f_ofe );
	}
	if( p->cd->fdebug >= 2 ) tprintf( "Objective function %g\n", phi );
	if( p->cd->compute_phi ) { p->phi = phi; p->success = success_all; }
	if( p->cd->phi_cutoff > DBL_EPSILON )
	{
		if( phi < p->cd->phi_cutoff )
		{
			p->success = 1;
			if( p->cd->fdebug ) tprintf( "SUCCESS: OF is below predefined cutoff value (%g < %g; func_intrn)!\n", phi, p->cd->phi_cutoff );
		}
		else p->success = 0;
	}
	if( p->cd->check_success ) { p->success = success_all; p->phi = phi; if( p->cd->fdebug && success_all ) tprintf( "SUCCESS: Model results are within the predefined ranges (func_intrn)!\n" ); }
	return GSL_SUCCESS;
}

void func_levmar( double *x, double *f, double *o, int m, int n, void *data ) /* forward run for LevMar */
{
	func_global( x, data, f, o );
}

void func_dx_levmar( double *x, double *f, double *o, double *jac, int m, int n, void *data ) /* forward run for LevMar */
{
	struct opt_data *p = ( struct opt_data * ) data;
	double *jacobian;
	int i, j, k;
	if( ( jacobian = ( double * ) malloc( sizeof( double ) * p->pd->nOptParam * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); mads_quits( p->root ); }
	if( p->cd->omp || p->cd->posix )
		func_dx_set( x, f, o, data, jacobian );
	else
		func_dx( x, f, o, data, jacobian );
	for( k = j = 0; j < p->pd->nOptParam; j++ ) // LEVMAR is using different jacobian order
		for( i = 0; i < p->od->nTObs; i++, k++ )
			jac[i * p->pd->nOptParam + j] = jacobian[k]; // order: obs / param
	free( jacobian );
}

int func_set( int n_sub, double *var_mat[], double *phi, double *f[], double *o[], FILE *out, struct opt_data *op ) // TODO use this function for executions in general
{
	int count;
	time_t time_start, time_end, time_elapsed;
	int ieval = op->cd->neval;
	if( op->cd->solution_type[0] == EXTERNAL && op->cd->posix ) // POSIX/MPRUN Parallel job
	{
		tprintf( "POSIX/MPRUN Parallel writing and execution of all external jobs ...\n" );
		time_start = time( NULL );
		if( mprunall( n_sub, op, var_mat, phi, f, o ) < 0 ) // Write and execute all the files in parallel; Reading does not work
		{
			tprintf( "ERROR: there is a problem with the MPRUNALL parallel execution!\n" );
			return( 0 );
		}
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( op->cd->tdebug ) tprintf( "Parallel set writing & execution PT = %ld seconds\n", time_elapsed );
		time_start = time_end;
		int read_error = 0;
		#pragma omp parallel for private(count)
		for( count = 0; count < n_sub; count++ ) // Read all the files in serial
		{
			double *opt_res, *opt_o;
			if( ( opt_res = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) tprintf( "Not enough memory!\n" );
			if( ( opt_o = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) tprintf( "Not enough memory!\n" );
			if( op->cd->pardebug ) tprintf( "OpenMP reading the model output files for case %d ...\n", ieval + count + 1 );
			if( func_extrn_read( ieval + count + 1, op, opt_res, opt_o ) )
				read_error = 1;
			else
			{
				if( f != NULL )
				{
					double lphi = 0;
					int j;
					for( j = 0; j < op->od->nTObs; j++ )
					{
						f[count][j] = opt_res[j];
						lphi += opt_res[j] * opt_res[j];
						if( o != NULL ) o[count][j] = opt_o[j];
					}
					if( phi != NULL ) phi[count] = lphi;
				}
			}
			free( opt_res );
			free( opt_o );
		}
		if( read_error == 1 ) return( 0 );
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( op->cd->tdebug ) tprintf( "Parallel set reading PT = %ld seconds\n", time_elapsed );
		return( 1 );
	}
	if( op->cd->solution_type[0] == EXTERNAL && !op->cd->posix && op->cd->omp ) // Pure OpenMP Parallel job
		return( func_set_omp( n_sub, var_mat, phi, f, o, out, op ) );
	if( op->cd->solution_type[0] == EXTERNAL && op->cd->parallel_type ) // Parallel job; potentially mix of OpenMP and POSIX threads
	{
		tprintf( "MPRUN Parallel execution of external jobs ...\n" );
		if( op->cd->omp ) tprintf( "OpenMP parallel" );
		else tprintf( "Serial" );
		tprintf( " generation of all the model input files ...\n" );
		time_start = time( NULL );
		#pragma omp parallel for private(count)
		for( count = 0; count < n_sub; count++ ) // Write all the files
		{
			if( out != NULL ) fprintf( out, "%d : ", count + 1 ); // counter
			if( op->cd->mdebug ) tprintf( "\n" );
			if( op->cd->debug > 1 || op->cd->mdebug )  tprintf( "Set #%d (%d): ", count + 1, omp_get_thread_num() );
			double *opt_params;
			if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) printf( "Not enough memory!\n" );
			int i, k;
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				k = op->pd->var_index[i];
				opt_params[i] = op->pd->var[k] = var_mat[count][i];
			}
			int debug_level = 0;
			if( op->cd->mdebug > 1 ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			func_extrn_write( ieval + count + 1, opt_params, op );
			free( opt_params );
			if( op->cd->debug > 1 || op->cd->mdebug ) tprintf( "external model input file(s) generated ...\n" );
			if( op->cd->mdebug > 1 ) op->cd->fdebug = debug_level;
			if( op->cd->mdebug )
			{
				tprintf( "Parameter values:\n" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					if( op->pd->var_log[op->pd->var_index[i]] == 0 ) tprintf( "%s %g\n", op->pd->var_name[op->pd->var_index[i]], op->pd->var[op->pd->var_index[i]] );
					else tprintf( "%s %g\n", op->pd->var_name[op->pd->var_index[i]], pow( 10, op->pd->var[op->pd->var_index[i]] ) );
			}
			if( out != NULL )
			{
				for( i = 0; i < op->pd->nParam; i++ )
					if( op->pd->var_opt[i] >= 1 )
					{
						if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
						else fprintf( out, " %.15g", op->pd->var[i] );
					}
				fprintf( out, "\n" );
			}
		}
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( op->cd->tdebug ) tprintf( "Parallel set writing PT = %ld seconds\n", time_elapsed );
		time_start = time_end;
		if( op->cd->pardebug > 9 )
		{
			tprintf( "Serial execution of external jobs ...\n" );
			for( count = 0; count < n_sub; count++ ) // Perform all the runs in serial model (for testing)
			{
				tprintf( "Execute model #%d ... ", ieval + count + 1 );
				func_extrn_exec_serial( ieval + count + 1, op );
				tprintf( "done!\n" );
			}
		}
		else
		{
			tprintf( "MPRUN parallel execution of external jobs ...\n" );
			if( mprun( n_sub, op ) < 0 )
			{
				tprintf( "ERROR: there is a problem with the parallel execution!\n" );
				return( 0 );
			}
		}
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( op->cd->tdebug ) tprintf( "Parallel set execution PT = %ld seconds\n", time_elapsed );
		if( op->cd->pardebug )
		{
			if( op->cd->omp ) tprintf( "OpenMP parallel" );
			else tprintf( "Serial" );
			tprintf( " reading of all the model input files ...\n" );
		}
		time_start = time_end;
		int read_error = 0;
		#pragma omp parallel for private(count)
		for( count = 0; count < n_sub; count++ ) // Read all the files in serial
		{
			double *opt_res, *opt_o;
			if( ( opt_res = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) tprintf( "Not enough memory!\n" );
			if( ( opt_o = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) tprintf( "Not enough memory!\n" );
			if( op->cd->debug > 1 || op->cd->mdebug ) tprintf( "Reading all the model output files ...\n" );
			if( out != NULL )
			{
				fprintf( out, "%d :\n", ieval + count + 1 ); // counter
				int i, k;
				for( i = 0; i < op->pd->nOptParam; i++ ) // re
				{
					k = op->pd->var_index[i];
					op->pd->var[k] = var_mat[count][i];
					fprintf( out, "%s %.12g\n", op->pd->var_name[k], op->pd->var[k] );
				}
				fflush( out );
			}
			if( func_extrn_read( ieval + count + 1, op, opt_res, opt_o ) )
				read_error = 1;
			else
			{
				if( f != NULL )
				{
					double lphi = 0;
					int j;
					for( j = 0; j < op->od->nTObs; j++ )
					{
						f[count][j] = opt_res[j];
						lphi += opt_res[j] * opt_res[j];
						if( o != NULL ) o[count][j] = opt_o[j];
					}
					if( phi != NULL ) phi[count] = lphi;
				}
			}
			free( opt_res );
			free( opt_o );
		}
		if( read_error == 1 ) return( 0 );
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( op->cd->tdebug ) tprintf( "Parallel set reading PT = %ld seconds\n", time_elapsed );
	}
	else // Serial job
	{
		int i, j, k;
		if( op->cd->pardebug ) tprintf( "Serial jobs ...\n" );
		double *opt_params;
		if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		for( count = 0; count < n_sub; count++ )
		{
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				k = op->pd->var_index[i];
				opt_params[i] = op->pd->var[k] = var_mat[count][i];
			}
			int debug_level = 0;
			if( op->cd->mdebug > 1 ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			func_global( opt_params, op, op->od->res, op->od->obs_current );
			if( op->cd->mdebug > 1 ) op->cd->fdebug = debug_level;
			if( phi != NULL ) phi[count] = op->phi;
			if( f != NULL )
				for( j = 0; j < op->od->nTObs; j++ )
				{
					f[count][j] = op->od->res[j];
					if( o != NULL ) o[count][j] = op->od->obs_current[j];
				}
			if( op->cd->mdebug > 1 && out != NULL )
			{
				// save to results file
				fprintf( out, "%d : ", count + 1 ); // counter
				fprintf( out, "%g :", op->phi );
				for( i = 0; i < op->pd->nParam; i++ )
					if( op->pd->var_opt[i] >= 1 )
						fprintf( out, " %g", op->pd->var[i] );
				fprintf( out, "\n" );
				fflush( out );
			}
		}
		free( opt_params );
	}
	return( 1 );
}

int func_set_omp( int n_sub, double *var_mat[], double *phi, double *f[], double *o[], FILE *out, struct opt_data *op ) // TODO use this function for executions in general
{
	int count, bad_data = 0;
	time_t time_start, time_end, time_elapsed;
	int ieval = op->cd->neval;
	tprintf( "OpenMP Parallel execution of external jobs ...\n" );
	time_start = time( NULL );
	#pragma omp parallel for private(count)
	for( count = 0; count < n_sub; count++ ) // Write all the files
	{
		int i, j, k;
		int rank = omp_get_thread_num();
		int done = 0;
		if( out != NULL ) fprintf( out, "%d : ", count + 1 ); // counter
		if( op->cd->mdebug ) tprintf( "\nSet #%d: ", count + 1 );
		double *opt_res, *opt_o;
		if( ( opt_res = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) tprintf( "Not enough memory!\n" );
		if( ( opt_o = ( double * ) malloc( op->od->nTObs * sizeof( double ) ) ) == NULL ) tprintf( "Not enough memory!\n" );
		if( op->cd->restart ) // Check for already computed jobs (smart restart)
		{
			done = func_extrn_check_read( ieval + count + 1, op );
			if( !done )
			{
				FILE *infileb;
				char buf[1000];
				sprintf( buf, "%s/%020d.res", op->cd->restart_container, ieval + count + 1 );
				if( ( infileb = fopen( buf, "rb" ) ) != NULL )
				{
					int obj_read = fread( ( void * ) opt_res, sizeof( opt_res[0] ), op->od->nTObs, infileb );
					fclose( infileb );
					if( obj_read != op->od->nTObs ) { tprintf( "RESTART ERROR: Binary file %s cannot be applied to read model predictions; data mismatch!\n", buf ); done = 0; }
					else { if( op->cd->pardebug ) tprintf( "RESTART: Binary file %s is applied to read model predictions\n", buf ); done = 2; }
					if( op->cd->pardebug > 4 )
						for( i = 0; i < op->od->nTObs; i++ )
							tprintf( "%-27s: binary observations %15.12g\n", op->od->obs_id[i], opt_res[i] );
				}
			}
			if( op->cd->pardebug > 1 )
			{
				if( done ) tprintf( "OPENMP #%d RESTART: Job %d is already completed; it will be skipped!\n", rank, ieval + count + 1 );
				else tprintf( "OPENMP #%d: Job %d will be executed!\n", rank, ieval + count + 1 );
			}
		}
		if( !done )
		{
			double *opt_params;
			if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) printf( "Not enough memory!\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				k = op->pd->var_index[i];
				opt_params[i] = op->pd->var[k] = var_mat[count][i];
			}
			int debug_level = 0;
			if( op->cd->mdebug > 1 ) { debug_level = op->cd->fdebug; op->cd->fdebug = 3; }
			func_extrn_write( ieval + count + 1, opt_params, op );
			free( opt_params );
			if( op->cd->mdebug > 1 ) op->cd->fdebug = debug_level;
			if( op->cd->mdebug )
			{
				tprintf( "Parameter values:\n" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					if( op->pd->var_log[op->pd->var_index[i]] == 0 ) tprintf( "%s %g\n", op->pd->var_name[op->pd->var_index[i]], op->pd->var[op->pd->var_index[i]] );
					else tprintf( "%s %g\n", op->pd->var_name[op->pd->var_index[i]], pow( 10, op->pd->var[op->pd->var_index[i]] ) );
			}
			if( out != NULL )
			{
				for( i = 0; i < op->pd->nParam; i++ )
					if( op->pd->var_opt[i] >= 1 )
					{
						if( op->pd->var_log[i] ) fprintf( out, " %.15g", pow( 10, op->pd->var[i] ) );
						else fprintf( out, " %.15g", op->pd->var[i] );
					}
				fprintf( out, "\n" );
			}
			func_extrn_exec_serial( ieval + count + 1, op );
			if( out != NULL )
			{
				fprintf( out, "%d :\n", ieval + count + 1 ); // counter
				for( i = 0; i < op->pd->nOptParam; i++ ) // re
				{
					k = op->pd->var_index[i];
					op->pd->var[k] = var_mat[count][i];
					fprintf( out, "%s %.12g\n", op->pd->var_name[k], op->pd->var[k] );
				}
				fflush( out );
			}
		}
		if( done != 2 )
		{
			if( func_extrn_read( ieval + count + 1, op, opt_res, opt_o ) )
				bad_data = 1;
		}
		if( bad_data == 0 && f != NULL )
		{
			double lphi = 0;
			for( j = 0; j < op->od->nTObs; j++ )
			{
				f[count][j] = opt_res[j];
				lphi += opt_res[j] * opt_res[j];
				if( o != NULL ) o[count][j] = opt_o[j];
			}
			if( phi != NULL ) phi[count] = lphi;
		}
		free( opt_res );
		free( opt_o );
	}
	time_end = time( NULL );
	time_elapsed = time_end - time_start;
	op->cd->neval = ieval + n_sub;
	if( op->cd->tdebug ) tprintf( "OpenMP Parallel execution PT = %ld seconds\n", time_elapsed );
	if( bad_data ) return( 0 ); // return without deleting the directories
	return( 1 );
}


int func_dx( double *x, double *f, double *o, void *data, double *jacobian ) /* Compute Jacobian using forward numerical derivatives */
{
	struct opt_data *p = ( struct opt_data * )data;
	double *f_xpdx;
	double x_old, dx;
	time_t time_start, time_end, time_elapsed;
	int i, j, k, old_phi, compute_center = 0, bad_data = 0, ieval;
	old_phi = p->cd->compute_phi;
	p->cd->compute_phi = 0;
	if( p->cd->parallel_type && p->cd->solution_type[0] == EXTERNAL ) // Parallel execution of external runs
	{
		time_start = time( NULL );
		ieval = p->cd->neval;
		if( f == NULL ) // Model predictions for x are not provided; need to compute
		{
			compute_center = 1;
			if( ( f = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
			if( ( o = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
			func_extrn_write( ++ieval, x, data );
		}
		else if( p->cd->compute_center ) // Model predictions have to be compute
		{
			compute_center = 1;
			func_extrn_write( ++ieval, x, data );
		}
		for( j = 0; j < p->pd->nOptParam; j++ )
		{
			x_old = x[j];
			if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			x[j] += dx;
			func_extrn_write( ++ieval, x, data );
			x[j] = x_old;
		}
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( p->cd->tdebug ) tprintf( "Parallel jacobian writing PT = %ld seconds\n", time_elapsed );
		time_start = time_end;
		if( mprun( p->pd->nOptParam + compute_center, data ) < 0 ) // Perform all the runs in parallel
		{
			tprintf( "ERROR: there is a problem with the parallel execution!\n" );
			mads_quits( p->root );
		}
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( p->cd->tdebug ) tprintf( "Parallel jacobian execution PT = %ld seconds\n", time_elapsed );
		// sleep( 1 );
		system( "sleep 0" ); // TODO investigate how much sleep is needed
		ieval -= ( p->pd->nOptParam + compute_center );
		time_start = time_end;
		if( compute_center ) func_extrn_read( ++ieval, data, f, o );
		if( ( f_xpdx = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			bad_data = func_extrn_read( ++ieval, data, f_xpdx, NULL );
			if( bad_data ) mads_quits( p->root );
			if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			for( i = 0; i < p->od->nTObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f[i] ) / dx;
		}
		free( f_xpdx );
		time_end = time( NULL );
		time_elapsed = time_end - time_start;
		if( p->cd->tdebug ) tprintf( "Parallel jacobian reading PT = %ld seconds\n", time_elapsed );
	}
	else
	{
		if( f == NULL ) // Model predictions for x are not provided; need to compute
		{
			compute_center = 1;
			if( ( f = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
			if( ( o = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
			func_global( x, data, f, o );
		}
		else if( p->cd->compute_center )
		{
			compute_center = 2;
			func_global( x, data, f, o );
		}
		if( ( f_xpdx = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			x_old = x[j];
			if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			x[j] += dx;
			func_global( x, data, f_xpdx, o );
			x[j] = x_old;
			for( i = 0; i < p->od->nTObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f[i] ) / dx;
		}
		free( f_xpdx );
		if( compute_center && !p->cd->compute_center ) { free( f ); free( o ); }
	}
	p->cd->compute_phi = old_phi;
	return GSL_SUCCESS;
}

int func_dx_set( double *x, double *f, double *o, void *data, double *jacobian ) /* Compute Jacobian using forward numerical derivatives */
{
	struct opt_data *p = ( struct opt_data * )data;
	double **par_mat, **obs_mat, **o_mat;
	double x_old, dx;
	time_t time_start, time_end, time_elapsed;
	int i, j, k, old_phi, compute_center = 0, ieval;
	if( p->cd->pardebug ) tprintf( "Parallel jacobian execution ...\n" );
	old_phi = p->cd->compute_phi;
	p->cd->compute_phi = 0;
	ieval = p->cd->neval;
	if( f == NULL || p->cd->compute_center ) compute_center = 1; // Model predictions for x are not provided; need to compute
	if( ( par_mat = double_matrix( p->pd->nOptParam + compute_center, p->pd->nOptParam ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( obs_mat = double_matrix( p->pd->nOptParam + compute_center, p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( o_mat = double_matrix( p->pd->nOptParam + compute_center, p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( compute_center )
	{
		for( k = 0; k < p->pd->nOptParam; k++ )
			par_mat[0][k] = x[k];
	}
	for( j = 0; j < p->pd->nOptParam; j++ )
	{
		x_old = x[j];
		if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
		else dx = p->cd->sindx;
		x[j] += dx;
		for( k = 0; k < p->pd->nOptParam; k++ )
			par_mat[compute_center + j][k] = x[k];
		x[j] = x_old;
	}
	/*
	for( j = 0; j < ( compute_center + p->pd->nOptParam ) ; j++ )
	{
		for( k = 0; k < p->pd->nOptParam; k++ )
			tprintf( "%g ", par_mat[compute_center + j][k] );
		tprintf( "\n" );
	}
	*/
	time_start = time( NULL );
	func_set( p->pd->nOptParam + compute_center, par_mat, NULL, obs_mat, o_mat, NULL, p );
	time_end = time( NULL );
	time_elapsed = time_end - time_start;
	if( p->cd->tdebug ) tprintf( "Parallel jacobian execution PT = %ld seconds\n", time_elapsed );
	if( p->cd->compute_center )
		for( i = 0; i < p->od->nTObs; i++ )
		{
			f[i] = obs_mat[0][i];
			if( o != NULL ) o[i] = o_mat[0][i];
		}
	else if( f == NULL ) f = obs_mat[0];
	ieval -= ( p->pd->nOptParam + compute_center );
	time_start = time_end;
	for( k = 0, j = compute_center; j < ( compute_center + p->pd->nOptParam ); j++ )
	{
		if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
		else dx = p->cd->sindx;
		for( i = 0; i < p->od->nTObs; i++, k++ )
			jacobian[k] = ( obs_mat[j][i] - f[i] ) / dx;
	}
	time_end = time( NULL );
	time_elapsed = time_end - time_start;
	if( p->cd->tdebug ) tprintf( "Parallel jacobian reading PT = %ld seconds\n", time_elapsed );
	free_matrix( ( void ** ) par_mat, compute_center + p->pd->nOptParam );
	free_matrix( ( void ** ) obs_mat, compute_center + p->pd->nOptParam );
	free_matrix( ( void ** ) o_mat, compute_center + p->pd->nOptParam );
	p->cd->compute_phi = old_phi;
	return GSL_SUCCESS;
}

double func_solver1( double x, double y, double z, double t, void *data ) // Compute for given (x, y, z, t)
{
	int i, j, k, s, num_params;
	double c, dx, dy, dz, x1, y1, z1, dist;
	struct calc_data *cd = ( struct calc_data * )data;
	struct anal_data ad;
	dx = dy = dz = 0.0;
	c = cd->c_background;
	num_params = cd->num_source_params + cd->num_aquifer_params;
	if( ( ad.var = ( double * ) malloc( ( num_params ) * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	k = cd->num_source_params * cd->num_sources + cd->num_aquifer_params;
	j = cd->num_source_params;
	for( i = cd->num_source_params * cd->num_sources; i < k; i++, j++ )
		ad.var[j] = cd->var[i];
	if( cd->disp_tied && cd->disp_scaled == 0 ) // Tied dispersivities
	{
		if( cd->fdebug >= 5 )
		{
			tprintf( "Tied AY %.12g = %.12g / %.12g\n", ad.var[AX] / ad.var[AY], ad.var[AX], ad.var[AY] );
			tprintf( "Tied AZ %.12g = %.12g / %.12g\n", ( ad.var[AX] / ad.var[AY] ) / ad.var[AZ],  ad.var[AX] / ad.var[AY], ad.var[AZ] );
		}
		ad.var[AY] = ad.var[AX] / ad.var[AY];
		ad.var[AZ] = ad.var[AY] / ad.var[AZ];
	}
	for( s = 0; s < cd->num_sources; s++ )
	{
		k = cd->num_source_params * ( s + 1 );
		j = 0;
		for( i = cd->num_source_params * s; i < k; i++, j++ )
			ad.var[j] = cd->var[i];
		if( cd->disp_scaled ) // Scaled dispersivities
		{
			dx = ad.var[AX]; dy = ad.var[AY]; dz = ad.var[AZ];
			x1 = x - ad.var[SOURCE_X];
			y1 = y - ad.var[SOURCE_Y];
			z1 = z - ad.var[SOURCE_Z];
			dist = sqrt( x1 * x1 + y1 * y1 + z1 * z1 );
			// tprintf( "func_solver1\n" );
			if( cd->fdebug >= 5 ) tprintf( "Scaled AX %.12g = %.12g * %.12g\n", ad.var[AX] * dist, ad.var[AX], dist );
			ad.var[AX] *= dist;
			if( cd->disp_scaled > 1 && !cd->disp_tied ) { ad.var[AY] *= dist; ad.var[AZ] *= dist; }
			else if( cd->disp_tied ) { ad.var[AY] = ad.var[AX] / ad.var[AY]; ad.var[AZ] = ad.var[AY] / ad.var[AZ]; };
			if( cd->fdebug >= 5 )
			{
				if( cd->disp_scaled > 1 && !cd->disp_tied ) tprintf( "Transverse dispersivities are scaled!\n" );
				else if( cd->disp_tied ) tprintf( "Transverse dispersivities are tied!\n" );
				else tprintf( "Transverse dispersivities are neither tied nor scaled!\n" );
				tprintf( "AY %.12g\n", ad.var[AY] );
				tprintf( "AZ %.12g\n", ad.var[AZ] );
			}
		}
		if( fabs( ad.var[TSCALE_DISP] - 1 ) < COMPARE_EPSILON || fabs( ad.var[TSCALE_DISP] ) < COMPARE_EPSILON ) ad.scaling_dispersion = 0;
		else ad.scaling_dispersion = 1;
		ad.time_step = cd->time_step;
		if( cd->fdebug > 6 )
			for( i = 0; i < num_params; i++ )
				tprintf( "func_solver1 source #%d parameter #%d %g\n", s + 1, i + 1, ad.var[i] );
		switch( cd->solution_type[s] )
		{
			case POINT:
				c += point_source( x, y, z, t, ( void * ) &ad );
				break;
			case POINT_TRIANGLE_TIME:
				c += point_source_triangle_time( x, y, z, t, ( void * ) &ad );
				break;
			case PLANE:
				c += rectangle_source( x, y, z, t, ( void * ) &ad );
				break;
			case PLANE3D:
				c += rectangle_source_vz( x, y, z, t, ( void * ) &ad );
				break;
			case GAUSSIAN2D:
				c += gaussian_source_2d( x, y, z, t, ( void * ) &ad );
				break;
			case GAUSSIAN3D:
				c += gaussian_source_3d( x, y, z, t, ( void * ) &ad );
				break;
			default:
			case BOX:
				if( cd->levy == FULL_LEVY ) c += box_source_levy_dispersion( x, y, z, t, ( void * ) &ad );
				else if( cd->levy == SYM_LEVY ) c += box_source_sym_levy_dispersion( x, y, z, t, ( void * ) &ad );
				else c += box_source( x, y, z, t, ( void * ) &ad );
				break;
		}
		if( cd->disp_scaled ) // Scaled dispersivities
		{
			ad.var[AX] = dx; ad.var[AY] = dy; ad.var[AZ] = dz;
		}
	}
	free( ad.var );
	return( c );
}

double func_solver( double x, double y, double z1, double z2, double t, void *data )
{
	struct calc_data *cd = ( struct calc_data * )data;
	if( cd->obs_int == 1 ) // TODO add other integration models ...
		return func_solver1( x, y, ( z1 + z2 ) / 2, t, data );
	else
		return( ( double )( func_solver1( x, y, z1, t, data ) + func_solver1( x, y, z2, t, data ) ) / 2 ); // Compute for (x, y, z1, t) and (x, y, z2, t) and average
}

void Transform( double *v, void *data, double *vt )
{
	struct opt_data *p = ( struct opt_data * )data;
	int i, k;
	if( p->cd->sintrans == 0 )
		for( i = 0; i < p->pd->nOptParam; i++ )
			vt[i] = v[i];
	else
		for( i = 0; i < p->pd->nOptParam; i++ )
		{
			k = p->pd->var_index[i];
			// tprintf( "trans %s %g %g %g -> ", p->pd->var_id[p->pd->var_index[i]], v[i], p->pd->var_range[k], p->pd->var_min[k] );
			vt[i] = ( v[i] - p->pd->var_min[k] ) / p->pd->var_range[k];
			vt[i] = asin( ( double ) vt[i] * 2.0 - 1.0 );
			// tprintf( "%g\n", vt[i] );
		}
}

void DeTransform( double *v, void *data, double *vt )
{
	struct opt_data *p = ( struct opt_data * )data;
	int i, k;
	if( p->cd->sintrans == 0 )
		for( i = 0; i < p->pd->nOptParam; i++ )
			vt[i] = v[i];
	else
		for( i = 0; i < p->pd->nOptParam; i++ )
		{
			k = p->pd->var_index[i];
			vt[i] = ( ( double ) sin( v[i] ) + 1.0 ) / 2.0;
			vt[i] = p->pd->var_min[k] + vt[i] * p->pd->var_range[k];
			// tprintf( "detrans %s %g -> %g\n", p->pd->var_id[p->pd->var_index[i]], v[i], vt[i] );
		}
}
