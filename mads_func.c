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

#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <matheval.h>
#include "mads.h"
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
int func_extrn( double *x, void *data, double *f );
int func_extrn_write( int ieval, double *x, void *data );
int func_extrn_read( int ieval, void *data, double *f );
int func_extrn_check_read( int ieval, void *data );
int func_intrn( double *x, void *data, double *f );
void func_levmar( double *x, double *f, int m, int n, void *data );
void func_dx_levmar( double *x, double *f, double *jacobian, int m, int n, void *data ); // Jacobian order: obs / param
int func_dx( double *x, double *f_x, void *data, double *jacobian ); // Jacobian order: param / obs
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );
/* Functions elsewhere */
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
double test_problems( int D, int function, double *x, int nObs, double *o );
double point_source( double x, double y, double z, double t, void *params );
double rectangle_source( double x, double y, double z, double t, void *params );
double rectangle_source_vz( double x, double y, double z, double t, void *params );
double box_source( double x, double y, double z, double t, void *params );
int create_mprun_dir( char *dir );
int delete_mprun_dir( char *dir );
int mprun( int nJob, void *data );
int Ftest( char *filename );
int Ftestread( char *filename );
time_t Fdatetime_t( char *filename, int debug );

int func_extrn( double *x, void *data, double *f )
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000];
	double c, t, w, min, max, err, phi = 0.0;
	int i, k, success, success_all = 1, bad_data = 0;
	if( p->cd->num_proc > 1 ) // Parallel execution of a serial job to archive all the intermediate results
	{
		func_extrn_write( p->cd->neval + 1, x, data );
		if( mprun( 1, data ) < 0 ) // Perform one (1) run in parallel
		{
			tprintf( "ERROR: there is a problem with the parallel execution!\n" );
			exit( 1 );
		}
		bad_data = func_extrn_read( p->cd->neval, data, f ); // p->cd->eval was already incremented in mprun
		if( bad_data ) exit( -1 );
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
			tprintf( "%s %.12g\n", p->pd->var_id[k], p->cd->var[k] );
	}
	if( p->pd->nExpParam > 0 )
	{
		if( p->cd->fdebug >= 3 ) tprintf( "Tied model parameters (%d):\n", p->pd->nExpParam );
		for( i = 0; i < p->pd->nExpParam; i++ )
		{
			k = p->pd->param_expressions_index[i];
			p->cd->var[k] = evaluator_evaluate( p->pd->param_expressions[i], p->pd->nParam, p->pd->var_id, p->cd->var );
			if( p->cd->fdebug >= 3 ) tprintf( "%s = %s = %.12g\n", p->pd->var_id[k], evaluator_get_string( p->pd->param_expressions[i] ), p->cd->var[k] );
		}
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nFixParam == 0 ) tprintf( "NO fixed parameters.\n" );
		else
		{
			tprintf( "Fixed model parameters (%d):\n", p->pd->nFixParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->calib_type == PPSD ) )
					tprintf( "%s %.12g\n", p->pd->var_id[i], p->cd->var[i] );
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
			case SSD0: tprintf( "sum of squared discrepancies" ); break;
			default: tprintf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
		}
	}
	for( i = 0; i < p->ed->ntpl; i++ )
		if( par_tpl( p->pd->nParam, p->pd->var_id, p->cd->var, p->ed->fn_tpl[i], p->ed->fn_out[i], p->cd->tpldebug ) == -1 )
			exit( -1 );
	strcpy( buf, "rm -f " );
	for( i = 0; i < p->ed->nins; i++ )
	{
		strcat( buf, p->ed->fn_obs[i] );
		strcat( buf, " " );
	}
	strcat( buf, " >& /dev/null" );
	if( p->cd->tpldebug || p->cd->insdebug ) tprintf( "\nDelete the expected output files before execution (\'%s\')\n", buf );
	system( buf );
	if( p->cd->tpldebug || p->cd->insdebug ) tprintf( "Execute external model \'%s\' ... ", p->ed->cmdline );
	system( p->ed->cmdline );
	if( p->cd->tpldebug || p->cd->insdebug ) tprintf( "done!\n" );
	for( i = 0; i < p->od->nObs; i++ ) p->od->res[i] = -1;
	for( i = 0; i < p->ed->nins; i++ )
		if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, p->od->res, p->ed->fn_ins[i], p->ed->fn_obs[i], p->cd->insdebug ) == -1 )
			exit( -1 );
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( p->od->res[i] < 0 )
		{
			tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( p->od->res[i] > 1.5 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug )
				tprintf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], ( int ) p->od->res[i] );
			p->od->obs_current[i] /= p->od->res[i];
		}
	}
	if( bad_data ) exit( -1 );
	for( i = p->od->nObs; i < p->od->nTObs; i++ )
		p->od->obs_current[i] = evaluator_evaluate( p->rd->regul_expressions[i - p->od->nObs], p->pd->nParam, p->pd->var_id, p->cd->var );
	if( p->cd->fdebug >= 2 ) tprintf( "\nModel predictions:\n" );
	for( i = 0; i < p->od->nTObs; i++ )
	{
		c = p->od->obs_current[i];
		t = p->od->obs_target[i];
		w = p->od->obs_weight[i];
		min = p->od->obs_min[i];
		max = p->od->obs_max[i];
		if( p->od->obs_log[i] == 0 )
		{
			err = c - t;
			if( p->cd->objfunc_type != SSR )
			{
				if( p->cd->objfunc_type == SSD0 ) err = 0;
				else if( p->cd->objfunc_type == SSDA )
				{
					err = sqrt( fabs( err ) );
					if( c < t ) err *= -1;
				}
				if( c < p->od->obs_min[i] ) err += p->od->obs_min[i] - c;
				else if( c > p->od->obs_max[i] ) err += c - p->od->obs_max[i];
			}
		}
		else
		{
			if( c < DBL_EPSILON ) c = DBL_EPSILON;
			if( t < DBL_EPSILON ) t = DBL_EPSILON;
			err = log10( c ) - log10( t );
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
	int i, k;
	DeTransform( x, p, p->pd->var_current );
	if( p->cd->fdebug >= 3 ) tprintf( "Optimized model parameters (%d; model run = %d):\n", p->pd->nOptParam, ieval );
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) p->cd->var[k] = pow( 10, p->pd->var_current[i] );
		else p->cd->var[k] = p->pd->var_current[i];
		if( p->cd->fdebug >= 3 )
			tprintf( "%s %.12g\n", p->pd->var_id[k], p->cd->var[k] );
	}
	if( p->pd->nExpParam > 0 )
	{
		if( p->cd->fdebug >= 3 ) tprintf( "Tied model parameters (%d):\n", p->pd->nExpParam );
		for( i = 0; i < p->pd->nExpParam; i++ )
		{
			k = p->pd->param_expressions_index[i];
			p->cd->var[k] = evaluator_evaluate( p->pd->param_expressions[i], p->pd->nParam, p->pd->var_id, p->cd->var );
			if( p->cd->fdebug >= 3 ) tprintf( "%s = %s = %.12g\n", p->pd->var_id[k], evaluator_get_string( p->pd->param_expressions[i] ), p->cd->var[k] );
		}
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nFixParam == 0 ) tprintf( "NO fixed parameters.\n" );
		else
		{
			tprintf( "Fixed model parameters (%d):\n", p->pd->nFixParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->calib_type == PPSD ) )
					tprintf( "%s %.12g\n", p->pd->var_id[i], p->cd->var[i] );
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
			case SSD0: tprintf( "sum of squared discrepancies" ); break;
			default: tprintf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
		}
	}
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval ); // Name of directory for parallel runs
	// Delete expected output files in the root directory to prevent the creation of links to these files in the "child" directories
	strcpy( buf, "rm -f " );
	for( i = 0; i < p->ed->nins; i++ )
	{
		strcat( buf, p->ed->fn_obs[i] );
		strcat( buf, " " );
	}
	if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
	if( p->cd->pardebug > 2 ) tprintf( "Delete the expected output files before execution (\'%s\')\n", buf );
	system( buf );
	create_mprun_dir( dir ); // Create the child directory for parallel runs with link to the files in the working root directory
	for( i = 0; i < p->ed->ntpl; i++ ) // Create all the model input files
	{
		sprintf( buf, "../%s/%s", dir, p->ed->fn_out[i] );
		if( par_tpl( p->pd->nParam, p->pd->var_id, p->cd->var, p->ed->fn_tpl[i], buf, p->cd->tpldebug ) == -1 )
			exit( -1 );
	}
	// Update model input files in zip restart files
	if( p->cd->restart )
		sprintf( buf, "zip -u %s ", p->cd->restart_zip_file ); // Archive input files
	for( i = 0; i < p->ed->ntpl; i++ )
		sprintf( &buf[( int ) strlen( buf )], "../%s/%s ", dir, p->ed->fn_out[i] );
	if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
	system( buf );
	if( p->cd->pardebug > 3 ) tprintf( "Input files for parallel run #%d are archived!\n", ieval );
	if( p->cd->restart == 0 ) // Do not delete if restart is attempted
	{
		sprintf( buf, "cd ../%s; rm -f ", dir ); // Delete expected output files in the hosts directories
		for( i = 0; i < p->ed->nins; i++ )
		{
			strcat( buf, p->ed->fn_obs[i] );
			strcat( buf, " " );
		}
		if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
		if( p->cd->pardebug > 2 ) tprintf( "Delete the expected output files before execution (\'%s\')\n", buf );
		system( buf );
	}
	else // Just in case; the restart file should have been already extracted
	{
		sprintf( buf, "unzip -u -: %s ", p->cd->restart_zip_file ); // Archive input files
		for( i = 0; i < p->ed->nins; i++ )
			sprintf( &buf[( int ) strlen( buf )], "../%s/%s ", dir, p->ed->fn_obs[i] );
		if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
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
	if( p->cd->pardebug > 2 )
	{
		sprintf( buf, "cd ../%s; ls -altr ", dir ); // Check directory content
		system( buf );
	}
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) tprintf( "Execute external model \'%s\' ... ", p->ed->cmdline );
	sprintf( buf, "cd ../%s; %s", dir, p->ed->cmdline );
	system( buf );
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) tprintf( "done!\n" );
	return GSL_SUCCESS;
}

int func_extrn_check_read( int ieval, void *data ) // Check a series of output files after parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	int i, bad_data;
	for( i = 0; i < p->od->nObs; i++ ) p->od->res[i] = -1;
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval );
	if( p->cd->pardebug > 2 )
	{
		sprintf( buf, "cd ../%s; ls -altr ", dir ); // Check directory content
		system( buf );
	}
	for( i = 0; i < p->ed->nins; i++ )
	{
		sprintf( buf, "../%s/%s", dir, p->ed->fn_obs[i] );
		if( Ftestread( buf ) == 1 ) { if( p->cd->pardebug ) tprintf( "File %s cannot be opened to read.", buf ); return( 0 ); }
		else if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, p->od->res, p->ed->fn_ins[i], buf, 0 ) == -1 )
			return( 0 );
	}
	bad_data = 0;
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( p->od->res[i] < 0 )
		{
			if( p->cd->pardebug ) tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( p->od->res[i] > 1.5 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug || p->cd->pardebug )
				tprintf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], ( int ) p->od->res[i] );
			p->od->obs_current[i] /= p->od->res[i];
		}
	}
	if( bad_data ) return( 0 );
	if( ( p->cd->time_infile - Fdatetime_t( buf, 0 ) ) > 0 )
	{
		if( p->cd->pardebug ) tprintf( "File %s is older than the MADS input file.\n", buf );
		if( p->cd->restart == -1 ) return( 1 ); else return( 0 );
	}
	return( 1 );
}

int func_extrn_read( int ieval, void *data, double *f ) // Read a series of output files after parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	double c, t, w, min, max, err, phi = 0.0;
	int i, success, success_all = 1, bad_data;
	for( i = 0; i < p->od->nObs; i++ ) p->od->res[i] = -1;
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval );
	for( i = 0; i < p->ed->nins; i++ )
	{
		sprintf( buf, "../%s/%s", dir, p->ed->fn_obs[i] );
		if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, p->od->res, p->ed->fn_ins[i], buf, p->cd->insdebug ) == -1 )
			exit( -1 );
	}
	bad_data = 0;
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( p->od->res[i] < 0 )
		{
			tprintf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( p->od->res[i] > 1.5 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug )
				tprintf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], ( int ) p->od->res[i] );
			p->od->obs_current[i] /= p->od->res[i];
		}
	}
	if( bad_data ) return( bad_data );
	sprintf( buf, "zip -u %s ", p->cd->restart_zip_file ); // Archive output files
	for( i = 0; i < p->ed->nins; i++ )
		sprintf( &buf[strlen( buf )], "../%s/%s ", dir, p->ed->fn_obs[i] );
	if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
	system( buf );
	if( p->cd->pardebug > 3 ) tprintf( "Results from parallel run #%d are archived!\n", ieval );
	delete_mprun_dir( dir ); // Delete directory for parallel runs
	for( i = p->od->nObs; i < p->od->nTObs; i++ )
		p->od->obs_current[i] = evaluator_evaluate( p->rd->regul_expressions[i - p->od->nObs], p->pd->nParam, p->pd->var_id, p->cd->var );
	if( p->cd->fdebug >= 2 ) tprintf( "\nModel predictions (model run = %d):\n", ieval );
	for( i = 0; i < p->od->nTObs; i++ )
	{
		c = p->od->obs_current[i];
		t = p->od->obs_target[i];
		w = p->od->obs_weight[i];
		min = p->od->obs_min[i];
		max = p->od->obs_max[i];
		if( p->od->obs_log[i] == 0 )
		{
			err = c - t;
			if( p->cd->objfunc_type != SSR )
			{
				if( p->cd->objfunc_type == SSD0 ) err = 0;
				else if( p->cd->objfunc_type == SSDA )
				{
					err = sqrt( fabs( err ) );
					if( c < t ) err *= -1;
				}
				if( c < p->od->obs_min[i] ) err += p->od->obs_min[i] - c;
				else if( c > p->od->obs_max[i] ) err += c - p->od->obs_max[i];
			}
		}
		else
		{
			if( c < DBL_EPSILON ) c = DBL_EPSILON;
			if( t < DBL_EPSILON ) t = DBL_EPSILON;
			err = log10( c ) - log10( t );
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

int func_intrn( double *x, void *data, double *f ) /* forward run for LM */
{
	int i, j, k, p1, p2, l, s, success, success_all = 1;
	double c, t, w, min, max, c1, c2, err, phi = 0.0, dx, dy, dz, x1, y1, z1, dist;
	struct opt_data *p = ( struct opt_data * )data;
	char filename[255];
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
		if( p->cd->fdebug >= 3 ) tprintf( "%s %.12g\n", p->pd->var_id[k], p->cd->var[k] );
	}
	if( p->pd->nExpParam > 0 )
	{
		if( p->cd->fdebug >= 3 ) tprintf( "Tied model parameters (%d):\n", p->pd->nExpParam );
		for( i = 0; i < p->pd->nExpParam; i++ )
		{
			k = p->pd->param_expressions_index[i];
			p->cd->var[k] = evaluator_evaluate( p->pd->param_expressions[i], p->pd->nParam, p->pd->var_id_short, p->cd->var );
			if( p->cd->fdebug >= 3 ) tprintf( "%s = %s = %.12g\n", p->pd->var_id[k], evaluator_get_string( p->pd->param_expressions[i] ), p->cd->var[k] );
		}
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nFixParam == 0 ) tprintf( "NO fixed parameters.\n" );
		else
		{
			tprintf( "Fixed model parameters (%d):\n", p->pd->nFixParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->calib_type == PPSD ) )
					tprintf( "%s %.12g\n", p->pd->var_id[i], p->cd->var[i] );
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
		if( p->cd->test_func >= 40 )
		{
			phi = 0;
			if( p->cd->check_success && p->cd->obserror > DBL_EPSILON ) success_all = 1;
			else success_all = 0;
			for( k = 0; k < p->od->nTObs; k++ )
			{
				c = f[k];
				f[k] = err = c - p->od->obs_target[k];
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
		if( p->cd->fdebug >= 2 ) tprintf( "\nModel predictions:\n" );
		l = NUM_ANAL_PARAMS_SOURCE * p->cd->num_solutions + ( NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE ); // copy model parameters
		p2 = NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE - 1;
		for( p1 = NUM_ANAL_PARAMS_SOURCE * p->cd->num_solutions; p1 < l; p1++, p2++ )
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
		for( k = 0; k < p->od->nTObs; k++ ) // for computational efficiency performed only for observations with weight > DBL_EPSILON
		{
			if( p->cd->oderiv != -1 ) { k = p->cd->oderiv; }
			if( k >= p->od->nObs )
			{
				// regularization term
				c = evaluator_evaluate( p->rd->regul_expressions[k - p->od->nObs], p->pd->nParam, p->pd->var_id_short, p->cd->var );
			}
			else
			{
				// model predicted calibration target
				i = p->od->obs_well_index[k];
				j = p->od->obs_time_index[k];
				c1 = c2 = p->cd->c_background;
				for( s = 0; s < p->cd->num_solutions; s++ )
				{
					l = NUM_ANAL_PARAMS_SOURCE * ( s + 1 );
					p2 = 0;
					for( p1 = NUM_ANAL_PARAMS_SOURCE * s; p1 < l; p1++, p2++ )
						p->ad->var[p2] = p->cd->var[p1];
					if( p->cd->disp_scaled ) // Scaled dispersivities
					{
						dx = p->ad->var[AX]; dy = p->ad->var[AY]; dz = p->ad->var[AZ];
						x1 = p->wd->x[i] - p->ad->var[SOURCE_X];
						y1 = p->wd->y[i] - p->ad->var[SOURCE_Y];
						z1 = ( p->wd->z1[i] + p->wd->z2[i] ) - p->ad->var[SOURCE_Z];
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
					switch( p->cd->solution_type[s] )
					{
						case POINT:
							c1 += point_source( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							c2 += point_source( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							break;
						case PLANE:
							c1 += rectangle_source( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							c2 += rectangle_source( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							break;
						case PLANE3D:
							c1 += rectangle_source_vz( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							c2 += rectangle_source_vz( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							break;
						default:
						case BOX:
							c1 += box_source( p->wd->x[i], p->wd->y[i], p->wd->z1[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							c2 += box_source( p->wd->x[i], p->wd->y[i], p->wd->z2[i], p->wd->obs_time[i][j], ( void * ) p->ad );
							break;
					}
					if( p->cd->disp_scaled ) // Scaled dispersivities
					{
						p->ad->var[AX] = dx; p->ad->var[AY] = dy; p->ad->var[AZ] = dz;
					}
				}
				c = ( c1 + c2 ) / 2;
			}
			p->od->obs_current[k] = c;
			t = p->od->obs_target[k];
			w = p->od->obs_weight[k];
			min = p->od->obs_min[k];
			max = p->od->obs_max[k];
			if( p->od->obs_log[k] == 0 )
			{
				err = c - t;
				if( p->cd->objfunc_type != SSR )
				{
					if( p->cd->objfunc_type == SSDA )
					{
						err = sqrt( fabs( err ) );
						if( c < t ) err *= -1;
					}
					else err = 0; // SSD0 & SSDX
					if( p->cd->objfunc_type == SSDX ) { dx = max - min; if( p->cd->obsdomain > DBL_EPSILON && p->cd->obsdomain < dx ) dx = p->cd->obsdomain; if( dx > DBL_EPSILON ) { dx /= 10; min += dx; max -= dx; } }
					if( c < min ) err += min - c;
					else if( c > max ) err += c - max;
					// tprintf( "%g %g %g %g\n", err, c, min - c, c - max );
					if( p->cd->objfunc_type == SSDX ) { min = p->od->obs_min[k]; max = p->od->obs_max[k]; }
				}
			}
			else
			{
				if( c < DBL_EPSILON ) c = DBL_EPSILON;
				if( t < DBL_EPSILON ) t = DBL_EPSILON;
				err = log10( c ) - log10( t );
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
				if( !p->cd->compute_phi ) phi += f[i] * f[i];
			}
			if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
		}
	}
	p->success = success_all; // Just in case
	if( fabs( p->cd->obsstep ) > DBL_EPSILON && p->success )
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

void func_levmar( double *x, double *f, int m, int n, void *data ) /* forward run for LevMar */
{
	func_global( x, data, f );
}

void func_dx_levmar( double *x, double *f, double *jac, int m, int n, void *data ) /* forward run for LevMar */
{
	struct opt_data *p = ( struct opt_data * )data;
	double *jacobian;
	int i, j, k;
	if( ( jacobian = ( double * ) malloc( sizeof( double ) * p->pd->nOptParam * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); exit( 1 ); }
	func_dx( x, f, data, jacobian );
	for( k = j = 0; j < p->pd->nOptParam; j++ ) // LEVMAR is using different jacobian order
		for( i = 0; i < p->od->nTObs; i++, k++ )
			jac[i * p->pd->nOptParam + j] = jacobian[k]; // order: obs / param
	free( jacobian );
}

int func_dx( double *x, double *f_x, void *data, double *jacobian ) /* Compute Jacobian using forward numerical derivatives */
{
	struct opt_data *p = ( struct opt_data * )data;
	double *f_xpdx;
	double x_old, dx;
	int i, j, k, old, compute_center = 0, bad_data = 0, ieval;
	if( ( f_xpdx = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
	old = p->cd->compute_phi;
	p->cd->compute_phi = 0;
	if( p->cd->num_proc > 1 && p->cd->solution_type[0] == EXTERNAL ) // Parallel execution of external runs
	{
		ieval = p->cd->neval;
		if( f_x == NULL ) // Model predictions for x are not provided; need to compute
		{
			compute_center = 1;
			if( ( f_x = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 1 ); }
			func_extrn_write( ++ieval, x, data );
		}
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			x_old = x[j];
			if( p->cd->sintrans < DBL_EPSILON ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			x[j] += dx;
			func_extrn_write( ++ieval, x, data );
			x[j] = x_old;
		}
		if( mprun( p->pd->nOptParam + compute_center, data ) < 0 ) // Perform all the runs in parallel
		{
			tprintf( "ERROR: there is a problem with the parallel execution!\n" );
			exit( 1 );
		}
		system( "sleep 2" ); // TODO investigate how much sleep is needed
		ieval -= ( p->pd->nOptParam + compute_center );
		if( compute_center ) func_extrn_read( ++ieval, data, f_x );
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			bad_data = func_extrn_read( ++ieval, data, f_xpdx );
			if( bad_data ) exit( -1 );
			if( p->cd->sintrans < DBL_EPSILON ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			for( i = 0; i < p->od->nTObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f_x[i] ) / dx;
		}
	}
	else
	{
		if( f_x == NULL ) // Model predictions for x are not provided; need to compute
		{
			compute_center = 1;
			if( ( f_x = ( double * ) malloc( sizeof( double ) * p->od->nTObs ) ) == NULL )
			{ tprintf( "Not enough memory!\n" ); return( 1 ); }
			func_global( x, data, f_x );
		}
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			x_old = x[j];
			if( p->cd->sintrans < DBL_EPSILON ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			x[j] += dx;
			func_global( x, data, f_xpdx );
			x[j] = x_old;
			for( i = 0; i < p->od->nTObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f_x[i] ) / dx;
		}
	}
	if( compute_center ) free( f_x );
	free( f_xpdx );
	p->cd->compute_phi = old;
	return GSL_SUCCESS;
}

double func_solver1( double x, double y, double z, double t, void *data ) // Compute for given (x, y, z, t)
{
	int i, j, k, s;
	double c, dx, dy, dz, x1, y1, z1, dist;
	struct calc_data *p = ( struct calc_data * )data;
	struct anal_data ad;
	c = p->c_background;
	k = NUM_ANAL_PARAMS_SOURCE * p->num_solutions + ( NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE );
	j = NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE - 1;
	for( i = NUM_ANAL_PARAMS_SOURCE * p->num_solutions; i < k; i++, j++ )
		ad.var[j] = p->var[i];
	if( p->disp_tied && p->disp_scaled == 0 ) // Tied dispersivities
	{
		if( p->fdebug >= 5 )
		{
			tprintf( "Tied AY %.12g = %.12g / %.12g\n", ad.var[AX] / ad.var[AY], ad.var[AX], ad.var[AY] );
			tprintf( "Tied AZ %.12g = %.12g / %.12g\n", ( ad.var[AX] / ad.var[AY] ) / ad.var[AZ],  ad.var[AX] / ad.var[AY], ad.var[AZ] );
		}
		ad.var[AY] = ad.var[AX] / ad.var[AY];
		ad.var[AZ] = ad.var[AY] / ad.var[AZ];
	}
	for( s = 0; s < p->num_solutions; s++ )
	{
		k = NUM_ANAL_PARAMS_SOURCE * ( s + 1 );
		j = 0;
		for( i = NUM_ANAL_PARAMS_SOURCE * s; i < k; i++, j++ )
			ad.var[j] = p->var[i];
		if( p->disp_scaled ) // Scaled dispersivities
		{
			dx = ad.var[AX]; dy = ad.var[AY]; dz = ad.var[AZ];
			x1 = x - ad.var[SOURCE_X];
			y1 = y - ad.var[SOURCE_Y];
			z1 = z - ad.var[SOURCE_Z];
			dist = sqrt( x1 * x1 + y1 * y1 + z1 * z1 );
			// tprintf( "func_solver1\n" );
			if( p->fdebug >= 5 ) tprintf( "Scaled AX %.12g = %.12g * %.12g\n", ad.var[AX] * dist, ad.var[AX], dist );
			ad.var[AX] *= dist;
			if( p->disp_scaled > 1 && !p->disp_tied ) { ad.var[AY] *= dist; ad.var[AZ] *= dist; }
			else if( p->disp_tied ) { ad.var[AY] = ad.var[AX] / ad.var[AY]; ad.var[AZ] = ad.var[AY] / ad.var[AZ]; };
			if( p->fdebug >= 5 )
			{
				if( p->disp_scaled > 1 && !p->disp_tied ) tprintf( "Transverse dispersivities are scaled!\n" );
				else if( p->disp_tied ) tprintf( "Transverse dispersivities are tied!\n" );
				else tprintf( "Transverse dispersivities are neither tied nor scaled!\n" );
				tprintf( "AY %.12g\n", ad.var[AY] );
				tprintf( "AZ %.12g\n", ad.var[AZ] );
			}
		}
		if( p->fdebug > 6 )
			for( i = 0; i < NUM_ANAL_PARAMS; i++ )
				tprintf( "func_solver1 source #%d parameter #%d %g\n", s + 1, i + 1, ad.var[i] );
		switch( p->solution_type[s] )
		{
			case POINT:
				c += point_source( x, y, z, t, ( void * ) &ad );
				break;
			case PLANE:
				c += rectangle_source( x, y, z, t, ( void * ) &ad );
				break;
			case PLANE3D:
				c += rectangle_source_vz( x, y, z, t, ( void * ) &ad );
				break;
			default:
			case BOX:
				c += box_source( x, y, z, t, ( void * ) &ad );
				break;
		}
		if( p->disp_scaled ) // Scaled dispersivities
		{
			ad.var[AX] = dx; ad.var[AY] = dy; ad.var[AZ] = dz;
		}
	}
	return( c );
}

double func_solver( double x, double y, double z1, double z2, double t, void *data ) // Compute for (x, y, z1, t) and (x, y, z2, t) and average
{
	int i, j, k, s;
	double c1, c2, dx, dy, dz, x1, y1, z3, dist;
	struct calc_data *p = ( struct calc_data * )data;
	struct anal_data ad;
	c1 = c2 = p->c_background;
	k = NUM_ANAL_PARAMS_SOURCE * p->num_solutions + ( NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE );
	j = NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE - 1;
	for( i = NUM_ANAL_PARAMS_SOURCE * p->num_solutions; i < k; i++, j++ )
		ad.var[j] = p->var[i];
	if( p->disp_tied && p->disp_scaled == 0 ) // Tied dispersivities
	{
		if( p->fdebug >= 5 )
		{
			tprintf( "Tied AY %.12g = %.12g / %.12g\n", ad.var[AX] / ad.var[AY], ad.var[AX], ad.var[AY] );
			tprintf( "Tied AZ %.12g = %.12g / %.12g\n", ( ad.var[AX] / ad.var[AY] ) / ad.var[AZ], ad.var[AX] / ad.var[AY], ad.var[AZ] );
		}
		ad.var[AY] = ad.var[AX] / ad.var[AY];
		ad.var[AZ] = ad.var[AY] / ad.var[AZ];
	}
	for( s = 0; s < p->num_solutions; s++ )
	{
		k = NUM_ANAL_PARAMS_SOURCE * ( s + 1 );
		j = 0;
		for( i = NUM_ANAL_PARAMS_SOURCE * s; i < k; i++, j++ )
			ad.var[j] = p->var[i];
		if( p->fdebug >= 6 )
			for( i = 0; i < NUM_ANAL_PARAMS; i++ )
				tprintf( "func_solver source #%d parameter #%d %g\n", s + 1, i + 1, ad.var[i] );
		if( p->disp_scaled ) // Scaled dispersivities
		{
			dx = ad.var[AX]; dy = ad.var[AY]; dz = ad.var[AZ];
			x1 = x - ad.var[SOURCE_X];
			y1 = y - ad.var[SOURCE_Y];
			z3 = ( z1 + z2 ) - ad.var[SOURCE_Z];
			dist = sqrt( x1 * x1 + y1 * y1 + z3 * z3 );
			// tprintf( "func_solver\n" );
			if( p->fdebug >= 5 ) tprintf( "Scaled AX %.12g = %.12g * %.12g\n", ad.var[AX] * dist, ad.var[AX], dist );
			ad.var[AX] *= dist;
			if( p->disp_scaled > 1 && !p->disp_tied ) { ad.var[AY] *= dist; ad.var[AZ] *= dist; }
			else if( p->disp_tied ) { ad.var[AY] = ad.var[AX] / ad.var[AY]; ad.var[AZ] = ad.var[AY] / ad.var[AZ]; };
			if( p->fdebug >= 5 )
			{
				if( p->disp_scaled > 1 && !p->disp_tied ) tprintf( "Transverse dispersivities scaled!\n" );
				else if( p->disp_tied ) tprintf( "Transverse dispersivities tied!\n" );
				else tprintf( "Transverse dispersivities not tied and not scaled!\n" );
				tprintf( "AY %.12g\n", ad.var[AY] );
				tprintf( "AZ %.12g\n", ad.var[AZ] );
			}
		}
		switch( p->solution_type[s] )
		{
			case POINT:
				c1 += point_source( x, y, z1, t, ( void * ) &ad );
				c2 += point_source( x, y, z2, t, ( void * ) &ad );
				break;
			case PLANE:
				c1 += rectangle_source( x, y, z1, t, ( void * ) &ad );
				c2 += rectangle_source( x, y, z2, t, ( void * ) &ad );
				break;
			case PLANE3D:
				c1 += rectangle_source_vz( x, y, z1, t, ( void * ) &ad );
				c2 += rectangle_source_vz( x, y, z2, t, ( void * ) &ad );
				break;
			default:
			case BOX:
				c1 += box_source( x, y, z1, t, ( void * ) &ad );
				c2 += box_source( x, y, z2, t, ( void * ) &ad );
				break;
		}
		if( p->disp_scaled ) // Scaled dispersivities
		{
			ad.var[AX] = dx; ad.var[AY] = dy; ad.var[AZ] = dz;
		}
	}
	return( ( c1 + c2 ) / 2 );
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
		{
			vt[i] = v[i];
			// tprintf( "detrans %s %g -> %g\n", p->pd->var_id[p->pd->var_index[i]], v[i], vt[i] );
		}
	else
		for( i = 0; i < p->pd->nOptParam; i++ )
		{
			k = p->pd->var_index[i];
			vt[i] = ( ( double ) sin( v[i] ) + 1.0 ) / 2.0;
			vt[i] = p->pd->var_min[k] + vt[i] * p->pd->var_range[k];
			// tprintf( "detrans %s %g -> %g\n", p->pd->var_id[p->pd->var_index[i]], v[i], vt[i] );
		}
}
