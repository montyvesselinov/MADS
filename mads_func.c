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
#include "mads.h"
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
int func_extrn( double *x, void *data, double *f );
int func_extrn_write( int ieval, double *x, void *data );
int func_extrn_read( int ieval, void *data, double *f );
int func_extrn_check_read( int ieval, void *data );
int func_intrn( double *x, void *data, double *f );
void func_levmar( double *x, double *f, int m, int n, void *data );
void func_dx_levmar( double *x, double *f, double *jacobian, int m, int n, void *data );
int func_dx( double *x, double *f_x, void *data, double *jacobian );
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
time_t Fdatetime_t( char *filename, int debug );

int func_extrn( double *x, void *data, double *f )
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000];
	double c, t, err, phi = 0.0;
	int i, k, success, success_all = 1, bad_data = 0;
	if( p->cd->num_proc > 1 ) // Parallel execution of a serial job to archive all the intermediate results
	{
		func_extrn_write( p->cd->neval + 1, x, data );
		if( mprun( 1, data ) < 0 ) // Perform one (1) run in parallel
		{
			printf( "ERROR: there is a problem with the parallel execution!\n" );
			exit( 1 );
		}
		bad_data = func_extrn_read( p->cd->neval, data, f ); // p->cd->eval was already incremented in mprun
		if( bad_data ) exit( -1 );
		return GSL_SUCCESS; // DONE
	}
	p->cd->neval++;
	DeTransform( x, p, p->pd->var_current );
	if( p->cd->fdebug >= 3 ) printf( "Model parameters:\n" );
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) p->cd->var[k] = pow( 10, p->pd->var_current[i] );
		else p->cd->var[k] = p->pd->var_current[i];
		if( p->cd->fdebug >= 3 )
			printf( "m %s %.12g\n", p->pd->var_id[k], p->cd->var[k] );
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nParam == p->pd->nOptParam ) printf( "NO fixed parameters.\n" );
		else
		{
			printf( "Fixed parameters:\n" );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->calib_type == PPSD ) )
					printf( "%s %.12g\n", p->pd->var_id[i], p->cd->var[i] );
		}
	}
	if( p->cd->fdebug >= 4 )
	{
		printf( "Objective function: " );
		switch( p->cd->objfunc_type )
		{
			case SSR: printf( "sum of squared residuals" ); break;
			case SSDR: printf( "sum of squared discrepancies and squared residuals" ); break;
			case SSDA: printf( "sum of squared discrepancies and absolute residuals" ); break;
			case SSD0: printf( "sum of squared discrepancies" ); break;
			default: printf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
		}
	}
	for( i = 0; i < p->ed->ntpl; i++ )
		if( par_tpl( p->pd->nParam, p->pd->var_id, p->cd->var, p->ed->fn_tpl[i], p->ed->fn_out[i], p->cd->tpldebug ) == -1 )
			exit( -1 );
	strcpy( buf, "rm -f " );
	for( i = 0; i < p->ed->nins; i++ )
		strcat( buf, p->ed->fn_obs[i] );
	strcat( buf, " >& /dev/null" );
	if( p->cd->tpldebug || p->cd->insdebug ) printf( "\nDelete the expected output files before execution (\'%s\')\n", buf );
	system( buf );
	if( p->cd->tpldebug || p->cd->insdebug ) printf( "Execute external model \'%s\' ... ", p->ed->cmdline );
	system( p->ed->cmdline );
	if( p->cd->tpldebug || p->cd->insdebug ) printf( "done!\n" );
	for( i = 0; i < p->od->nObs; i++ ) p->od->res[i] = -1;
	for( i = 0; i < p->ed->nins; i++ )
		if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, p->od->res, p->ed->fn_ins[i], p->ed->fn_obs[i], p->cd->insdebug ) == -1 )
			exit( -1 );
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( p->od->res[i] < 0 )
		{
			printf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( p->od->res[i] > 1.5 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug )
				printf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], ( int ) p->od->res[i] );
			p->od->obs_current[i] /= p->od->res[i];
		}
	}
	if( bad_data ) exit( -1 );
	if( p->cd->fdebug >= 2 ) printf( "\nModel predictions:\n" );
	for( i = 0; i < p->od->nObs; i++ )
	{
		c = p->od->obs_current[i];
		t = p->od->obs_target[i];
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
				if( c < p->od->obs_min[i] ) err += c - p->od->obs_min[i];
				else if( c > p->od->obs_max[i] ) err += c - p->od->obs_max[i];
			}
		}
		else
		{
			if( c < DBL_EPSILON ) c = DBL_EPSILON;
			if( t < DBL_EPSILON ) t = DBL_EPSILON;
			err = log10( c ) - log10( t );
		}
		f[i] = err * p->od->obs_weight[i];
		if( p->cd->compute_phi ) phi += f[i] * f[i];
		if( p->cd->obserror > 0 )
		{
			if( fabs( c - t ) > p->cd->obserror ) { success = 0; if( p->od->obs_weight[i] > 0.0 ) success_all = 0; }
			else success = 1;
		}
		else // if( p->cd->obsrange )
		{
			if( ( c < p->od->obs_min[i] || c > p->od->obs_max[i] ) ) { success_all = success = 0; if( p->od->obs_weight[i] > 0.0 ) success_all = 0; }
			else success = 1;
		}
		if( p->cd->fdebug >= 2 )
		{
			if( p->od->nObs < 50 || ( i < 20 || i > p->od->nObs - 20 ) )
				printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", p->od->obs_id[i], t, c, err, err * p->od->obs_weight[i], success, p->od->obs_min[i], p->od->obs_max[i] );
			if( p->od->nObs > 50 && i == 21 ) printf( "...\n" );
			if( !p->cd->compute_phi ) phi += f[i] * f[i];
		}
		if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
	}
	if( p->cd->fdebug >= 2 ) printf( "Objective function %g\n", phi );
	if( p->cd->compute_phi ) { p->phi = phi; p->success = success_all; }
	if( p->cd->phi_cutoff > DBL_EPSILON && phi < p->cd->phi_cutoff )
	{
		p->success = 1;
		if( p->cd->fdebug ) printf( "SUCCESS: OF is below predefined cutoff value (%g < %g; func_extrn)!\n", phi, p->cd->phi_cutoff );
	}
	else p->success = 0;
	if( p->cd->check_success ) { p->success = success_all; p->phi = phi; if( p->cd->fdebug && success_all ) printf( "SUCCESS: Model results are within the predefined ranges (func_extrn)!\n" ); }
	return GSL_SUCCESS;
}

int func_extrn_write( int ieval, double *x, void *data ) // Create a series of input files for parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	int i, k;
	DeTransform( x, p, p->pd->var_current );
	if( p->cd->fdebug >= 3 ) printf( "Model parameters (model run = %d):\n", ieval );
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) p->cd->var[k] = pow( 10, p->pd->var_current[i] );
		else p->cd->var[k] = p->pd->var_current[i];
		if( p->cd->fdebug >= 3 )
			printf( "%s %.12g\n", p->pd->var_id[k], p->cd->var[k] );
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nParam == p->pd->nOptParam ) printf( "NO fixed parameters.\n" );
		else
		{
			printf( "Fixed parameters:\n" );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->calib_type == PPSD ) )
					printf( "%s %.12g\n", p->pd->var_id[i], p->cd->var[i] );
		}
	}
	if( p->cd->fdebug >= 4 )
	{
		printf( "Objective function: " );
		switch( p->cd->objfunc_type )
		{
			case SSR: printf( "sum of squared residuals" ); break;
			case SSDR: printf( "sum of squared discrepancies and squared residuals" ); break;
			case SSDA: printf( "sum of squared discrepancies and absolute residuals" ); break;
			case SSD0: printf( "sum of squared discrepancies" ); break;
			default: printf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
		}
	}
	sprintf( dir, "%s_%08d", p->cd->mydir_hosts, ieval ); // Name of directory for parallel runs
	// Delete expected output files in the root directory to prevent the creation of links to these files in the "child" directories
	strcpy( buf, "rm -f " );
	for( i = 0; i < p->ed->nins; i++ )
		strcat( buf, p->ed->fn_obs[i] );
	if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
	if( p->cd->pardebug > 2 ) printf( "Delete the expected output files before execution (\'%s\')\n", buf );
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
	if( p->cd->pardebug > 3 ) printf( "Input files for parallel run #%d are archived!\n", ieval );
	if( p->cd->restart == 0 ) // Do not delete if restart is attempted
	{
		sprintf( buf, "cd ../%s; rm -f ", dir ); // Delete expected output files in the hosts directories
		for( i = 0; i < p->ed->nins; i++ )
			strcat( buf, p->ed->fn_obs[i] );
		if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
		if( p->cd->pardebug > 2 ) printf( "Delete the expected output files before execution (\'%s\')\n", buf );
		system( buf );
	}
	else // Just in case; the restart file should have been already extracted
	{
		sprintf( buf, "unzip -u -: %s ", p->cd->restart_zip_file ); // Archive input files
		for( i = 0; i < p->ed->nins; i++ )
			sprintf( &buf[( int ) strlen( buf )], "../%s/%s ", dir, p->ed->fn_obs[i] );
		if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
		if( p->cd->pardebug > 2 ) printf( "Extract the expected output files before execution (\'%s\')\n", buf );
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
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) printf( "\nWorking directory: ../%s\n", dir );
	if( p->cd->pardebug > 2 )
	{
		sprintf( buf, "cd ../%s; ls -altr ", dir ); // Check directory content
		system( buf );
	}
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) printf( "Execute external model \'%s\' ... ", p->ed->cmdline );
	sprintf( buf, "cd ../%s; %s", dir, p->ed->cmdline );
	system( buf );
	if( p->cd->pardebug || p->cd->tpldebug || p->cd->insdebug ) printf( "done!\n" );
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
		if( Ftest( buf ) == 1 ) { if( p->cd->pardebug ) printf( "File %s does not exist.", buf ); return( 0 ); }
		else if( ins_obs( p->od->nObs, p->od->obs_id, p->od->obs_current, p->od->res, p->ed->fn_ins[i], buf, p->cd->insdebug ) == -1 )
			return( 0 );
	}
	bad_data = 0;
	for( i = 0; i < p->od->nObs; i++ )
	{
		if( p->od->res[i] < 0 )
		{
			if( p->cd->pardebug ) printf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( p->od->res[i] > 1.5 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug || p->cd->pardebug )
				printf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], ( int ) p->od->res[i] );
			p->od->obs_current[i] /= p->od->res[i];
		}
	}
	if( bad_data ) return( 0 );
	if( ( p->cd->time_infile - Fdatetime_t( buf, 0 ) ) > 0 )
	{
		if( p->cd->pardebug ) printf( "File %s is older than the MADS input file.\n", buf );
		if( p->cd->restart == -1 ) return( 1 ); else return( 0 );
	}
	return( 1 );
}

int func_extrn_read( int ieval, void *data, double *f ) // Read a series of output files after parallel execution
{
	struct opt_data *p = ( struct opt_data * )data;
	char buf[1000], dir[500];
	double c, t, err, phi = 0.0;
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
			printf( "ERROR: Observation '\%s\' is not assigned reading the model output files!\n", p->od->obs_id[i] );
			bad_data = 1;
		}
		else if( p->od->res[i] > 1.5 )
		{
			if( p->cd->debug || p->cd->tpldebug || p->cd->insdebug )
				printf( "WARNING: Observation '\%s\' is defined more than once (%d) in the instruction files! Arithmetic average is computed!\n", p->od->obs_id[i], ( int ) p->od->res[i] );
			p->od->obs_current[i] /= p->od->res[i];
		}
	}
	if( bad_data ) return( bad_data );
	sprintf( buf, "zip -u %s ", p->cd->restart_zip_file ); // Archive output files
	for( i = 0; i < p->ed->nins; i++ )
		sprintf( &buf[strlen( buf )], "../%s/%s", dir, p->ed->fn_obs[i] );
	if( p->cd->pardebug <= 3 ) strcat( buf, " >& /dev/null" );
	system( buf );
	if( p->cd->pardebug > 3 ) printf( "Results from parallel run #%d are archived!\n", ieval );
	delete_mprun_dir( dir ); // Delete directory for parallel runs
	if( p->cd->fdebug >= 2 ) printf( "\nModel predictions (model run = %d):\n", ieval );
	for( i = 0; i < p->od->nObs; i++ )
	{
		c = p->od->obs_current[i];
		t = p->od->obs_target[i];
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
				if( c < p->od->obs_min[i] ) err += c - p->od->obs_min[i];
				else if( c > p->od->obs_max[i] ) err += c - p->od->obs_max[i];
			}
		}
		else
		{
			if( c < DBL_EPSILON ) c = DBL_EPSILON;
			if( t < DBL_EPSILON ) t = DBL_EPSILON;
			err = log10( c ) - log10( t );
		}
		f[i] = err * p->od->obs_weight[i];
		if( p->cd->compute_phi ) phi += f[i] * f[i];
		if( p->cd->obserror > 0 )
		{
			if( fabs( c - t ) > p->cd->obserror ) { success = 0; if( p->od->obs_weight[i] > 0.0 ) success_all = 0; }
			else success = 1;
		}
		else // if( p->cd->obsrange )
		{
			if( ( c < p->od->obs_min[i] || c > p->od->obs_max[i] ) ) { success_all = success = 0; if( p->od->obs_weight[i] > 0.0 ) success_all = 0; }
			else success = 1;
		}
		if( p->cd->fdebug >= 2 )
		{
			if( p->od->nObs < 50 || ( i < 20 || i > p->od->nObs - 20 ) )
				printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", p->od->obs_id[i], t, c, err, err * p->od->obs_weight[i], success, p->od->obs_min[i], p->od->obs_max[i] );
			if( p->od->nObs > 50 && i == 21 ) printf( "...\n" );
			if( !p->cd->compute_phi ) phi += f[i] * f[i];
		}
		if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
	}
	if( p->cd->fdebug >= 2 ) printf( "Objective function %g\n", phi );
	if( p->cd->compute_phi ) { p->phi = phi; p->success = success_all; }
	if( p->cd->phi_cutoff > DBL_EPSILON && phi < p->cd->phi_cutoff )
	{
		p->success = 1;
		if( p->cd->fdebug ) printf( "SUCCESS: OF is below predefined cutoff value (%g < %g; func_intrn)!\n", phi, p->cd->phi_cutoff );
	}
	else p->success = 0;
	if( p->cd->check_success ) { p->success = success_all; p->phi = phi; if( p->cd->fdebug && success_all ) printf( "SUCCESS: Model results are within the predefined ranges (func_extrn)!\n" ); }
	return GSL_SUCCESS;
}

int func_intrn( double *x, void *data, double *f ) /* forward run for LM */
{
	int i, j, k, p1, p2, l, s, success, success_all = 1;
	double c, t, c1, c2, err, phi = 0.0;
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
	if( p->cd->fdebug >= 3 ) printf( "Model parameters (%d):\n", p->pd->nOptParam );
	for( i = 0; i < p->pd->nOptParam; i++ )
	{
		k = p->pd->var_index[i];
		if( p->pd->var_log[k] ) p->cd->var[k] = pow( 10, p->pd->var_current[i] );
		else p->cd->var[k] = p->pd->var_current[i];
		if( p->cd->fdebug >= 3 ) printf( "%s %.12g\n", p->pd->var_id[k], p->cd->var[k] );
	}
	if( p->cd->fdebug >= 3 )
	{
		if( p->pd->nParam == p->pd->nOptParam ) printf( "NO fixed parameters.\n" );
		else
		{
			printf( "Fixed parameters (%d):\n", p->pd->nParam - p->pd->nOptParam );
			for( i = 0; i < p->pd->nParam; i++ )
				if( p->pd->var_opt[i] == 0 || ( p->pd->var_opt[i] == 2 && p->cd->calib_type == PPSD ) )
					printf( "%s: %.12g\n", p->pd->var_id[i], p->cd->var[i] );
		}
	}
	if( p->cd->fdebug >= 4 )
	{
		printf( "Objective function: " );
		if( p->cd->solution_type[0] != TEST )
		{
			switch( p->cd->objfunc_type )
			{
				case SSR: printf( "sum of squared residuals" ); break;
				case SSDR: printf( "sum of squared discrepancies and squared residuals" ); break;
				case SSDA: printf( "sum of squared discrepancies and absolute residuals" ); break;
				case SSD0: printf( "sum of squared discrepancies" ); break;
				default: printf( "unknown value; sum of squared residuals assumed" ); p->cd->objfunc_type = SSR; break;
			}
			if( p->cd->compute_phi ) printf( " --- computed\n" );
		}
		else
			printf( "test function\n" );
	}
	// p->cd->compute_phi = 1;
	// if( p->cd->compute_phi ) printf( " --- computed!!!!\n" );
	if( p->cd->solution_type[0] == TEST )
	{
		p->phi = phi = test_problems( p->pd->nOptParam, p->cd->test_func, p->cd->var, p->od->nObs, f );
		if( p->cd->check_success && p->cd->obsrange ) success_all = 0;
		if( p->cd->check_success && p->cd->parerror > 0 )
		{
			success_all = 1;
			for( k = 0; k < p->pd->nOptParam; k++ )
				if( fabs( p->cd->var[k] - p->pd->var_truth[k] ) > p->cd->parerror ) success_all = 0;
			if( p->cd->fdebug >= 4 ) printf( "Test OF %g Success %d\n", phi, success_all );
		}
		if( p->cd->test_func >= 40 )
		{
			phi = 0;
			if( p->cd->check_success && p->cd->obserror > DBL_EPSILON ) success_all = 1;
			else success_all = 0;
			for( k = 0; k < p->od->nObs; k++ )
			{
				c = f[k];
				f[k] = err = c - p->od->obs_target[k];
				if( p->cd->check_success && p->cd->obserror > DBL_EPSILON )
				{
					if( c < p->od->obs_min[k] || c > p->od->obs_max[k] ) { success_all = success = 0; }
					else success = 1;
				}
				// printf( "e %g %g\n", err, p->od->obs_target[k] );
				phi += err * err;
			}
			p->phi = phi;
		}
		else
		{
			if( p->cd->check_success && p->cd->obserror > DBL_EPSILON ) success_all = 0;
			if( p->cd->fdebug >= 5 )
			{
				printf( "Test OF %d %g\n", success_all, phi );
				c = 0;
				for( k = 0; k < p->od->nObs; k++ )
				{
					printf( "%s %g\n", p->od->obs_id[k], f[k] );
					c += f[k] * f[k];
				}
				printf( "Test OF ? %g\n", c );
			}
		}
	}
	else
	{
		if( p->cd->fdebug >= 2 ) printf( "\nModel predictions:\n" );
		for( k = 0; k < p->od->nObs; k++ )
		{
			if( p->cd->oderiv != -1 ) { k = p->cd->oderiv; }
			i = p->od->well_index[k];
			j = p->od->time_index[k];
			c1 = c2 = p->cd->c_background;
			// printf( "HERE intern\n" );
			l = NUM_ANAL_PARAMS_SOURCE * p->cd->num_solutions + ( NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE );
			p2 = NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE - 1;
			for( p1 = NUM_ANAL_PARAMS_SOURCE * p->cd->num_solutions; p1 < l; p1++, p2++ )
				p->ad->var[p2] = p->cd->var[p1];
			if( p->cd->tied ) // Hard wiring; Tied dispersivities
			{
				// printf( "AY: divide %.12g by %.12g\n", p->ad->var[AX], p->ad->var[AY] );
				p->ad->var[AY] = p->ad->var[AX] / p->ad->var[AY];
				// printf( "AZ: divide %.12g by %.12g\n", p->ad->var[AY], p->ad->var[AZ] );
				p->ad->var[AZ] = p->ad->var[AY] / p->ad->var[AZ];
			}
			for( s = 0; s < p->cd->num_solutions; s++ )
			{
				l = NUM_ANAL_PARAMS_SOURCE * ( s + 1 );
				p2 = 0;
				for( p1 = NUM_ANAL_PARAMS_SOURCE * s; p1 < l; p1++, p2++ )
					p->ad->var[p2] = p->cd->var[p1];
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
			}
			c = ( c1 + c2 ) / 2;
			p->od->obs_current[k] = c;
			t = p->wd->obs_target[i][j];
			if( p->wd->obs_log[i][j] == 0 )
			{
				err = c - t;
				if( p->cd->objfunc_type != SSR )
				{
					if( p->cd->objfunc_type == SSD0 ) err = 0;
					if( p->cd->objfunc_type == SSDA )
					{
						err = sqrt( fabs( err ) );
						if( c < t ) err *= -1;
					}
					if( c < p->wd->obs_min[i][j] ) err += c - p->wd->obs_min[i][j];
					else if( c > p->wd->obs_max[i][j] ) err += c - p->wd->obs_max[i][j];
				}
			}
			else
			{
				if( c < DBL_EPSILON ) c = DBL_EPSILON;
				if( t < DBL_EPSILON ) t = DBL_EPSILON;
				err = log10( c ) - log10( t );
			}
			f[k] = err * p->wd->obs_weight[i][j];
			if( p->cd->compute_phi ) phi += f[k] * f[k];
			if( p->cd->obserror > 0 )
			{
				if( fabs( c - t ) > p->cd->obserror ) { success = 0; if( p->od->obs_weight[i] > 0.0 ) success_all = 0; }
				else success = 1;
			}
			else // if( p->cd->obsrange )
			{
				if( ( c < p->od->obs_min[i] || c > p->od->obs_max[i] ) ) { success_all = success = 0; if( p->od->obs_weight[i] > 0.0 ) success_all = 0; }
				else success = 1;
			}
			if( p->cd->fdebug >= 2 )
			{
				if( p->od->nObs < 50 || ( i < 20 || i > p->od->nObs - 20 ) )
					printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", p->od->obs_id[k], t, c, err, err * p->od->obs_weight[k], success, p->od->obs_min[k], p->od->obs_max[k] );
				if( p->od->nObs > 50 && i == 21 ) printf( "...\n" );
				if( !p->cd->compute_phi ) phi += f[i] * f[i];
			}
			if( p->cd->oderiv != -1 ) { return GSL_SUCCESS; }
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
	if( p->cd->fdebug >= 2 ) printf( "Objective function %g\n", phi );
	if( p->cd->compute_phi ) { p->phi = phi; p->success = success_all; }
	if( p->cd->phi_cutoff > DBL_EPSILON && phi < p->cd->phi_cutoff )
	{
		p->success = 1;
		if( p->cd->fdebug ) printf( "SUCCESS: OF is below predefined cutoff value (%g < %g; func_intrn)!\n", phi, p->cd->phi_cutoff );
	}
	else p->success = 0;
	if( p->cd->check_success ) { p->success = success_all; p->phi = phi; if( p->cd->fdebug && success_all ) printf( "SUCCESS: Model results are within the predefined ranges (func_intrn)!\n" ); }
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
	if( ( jacobian = ( double * ) malloc( sizeof( double ) * p->pd->nOptParam * p->od->nObs ) ) == NULL )
	{ printf( "Not enough memory!\n" ); exit( 1 ); }
	func_dx( x, f, data, jacobian );
	for( k = j = 0; j < p->pd->nOptParam; j++ ) // LEVMAR is using different jacobian order
		for( i = 0; i < p->od->nObs; i++, k++ )
			jac[i * p->pd->nOptParam + j] = jacobian[k];
	free( jacobian );
}

int func_dx( double *x, double *f_x, void *data, double *jacobian ) /* Compute Jacobian using forward numerical derivatives */
{
	struct opt_data *p = ( struct opt_data * )data;
	double *f_xpdx;
	double x_old, dx;
	int i, j, k, compute_center = 0, bad_data = 0, ieval;
	ieval = p->cd->neval;
	if( ( f_xpdx = ( double * ) malloc( sizeof( double ) * p->od->nObs ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	p->cd->compute_phi = 0;
	if( p->cd->num_proc > 1 && p->cd->solution_type[0] == EXTERNAL ) // Parallel execution of external runs
	{
		if( f_x == NULL ) // Model predictions for x are not provided; need to compute
		{
			compute_center = 1;
			if( ( f_x = ( double * ) malloc( sizeof( double ) * p->od->nObs ) ) == NULL )
			{ printf( "Not enough memory!\n" ); return( 1 ); }
			func_extrn_write( ++ieval, x, data );
		}
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			x_old = x[j];
			if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			x[j] += dx;
			func_extrn_write( ++ieval, x, data );
			x[j] = x_old;
		}
		if( mprun( p->pd->nOptParam + compute_center, data ) < 0 ) // Perform all the runs in parallel
		{
			printf( "ERROR: there is a problem with the parallel execution!\n" );
			exit( 1 );
		}
		ieval -= ( p->pd->nOptParam + compute_center );
		if( compute_center ) func_extrn_read( ++ieval, data, f_x );
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			bad_data = func_extrn_read( ++ieval, data, f_xpdx );
			if( bad_data ) exit( -1 );
			for( i = 0; i < p->od->nObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f_x[i] ) / dx;
		}
	}
	else
	{
		if( f_x == NULL ) // Model predictions for x are not provided; need to compute
		{
			compute_center = 1;
			if( ( f_x = ( double * ) malloc( sizeof( double ) * p->od->nObs ) ) == NULL )
			{ printf( "Not enough memory!\n" ); return( 1 ); }
			func_global( x, data, f_x );
		}
		for( k = j = 0; j < p->pd->nOptParam; j++ )
		{
			x_old = x[j];
			if( p->cd->sintrans == 0 ) { if( p->pd->var_dx[j] > DBL_EPSILON ) dx = p->pd->var_dx[j]; else dx = p->cd->lindx; }
			else dx = p->cd->sindx;
			x[j] += dx;
			func_global( x, data, f_xpdx );
			x[j] = x_old;
			for( i = 0; i < p->od->nObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f_x[i] ) / dx;
		}
	}
	if( compute_center ) free( f_x );
	free( f_xpdx );
	p->cd->compute_phi = 1;
	return GSL_SUCCESS;
}

double func_solver1( double x, double y, double z, double t, void *data ) // Compute for given (x, y, z, t)
{
	int i, j, k, s;
	double c;
	struct calc_data *p = ( struct calc_data * )data;
	struct anal_data ad;
	c = p->c_background;
	k = NUM_ANAL_PARAMS_SOURCE * p->num_solutions + ( NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE );
	j = NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE - 1;
	for( i = NUM_ANAL_PARAMS_SOURCE * p->num_solutions; i < k; i++, j++ )
		ad.var[j] = p->var[i];
	if( p->tied ) // Hard wiring; Tied dispersivities
	{
		// printf( "AY: divide %.12g by %.12g\n", ad.var[AX], ad.var[AY] );
		ad.var[AY] = ad.var[AX] / ad.var[AY];
		// printf( "AZ: divide %.12g by %.12g\n", ad.var[AY], ad.var[AZ] );
		ad.var[AZ] = ad.var[AY] / ad.var[AZ];
	}
	for( s = 0; s < p->num_solutions; s++ )
	{
		k = NUM_ANAL_PARAMS_SOURCE * ( s + 1 );
		j = 0;
		for( i = NUM_ANAL_PARAMS_SOURCE * s; i < k; i++, j++ )
			ad.var[j] = p->var[i];
		if( p->fdebug )
			for( i = 0; i < NUM_ANAL_PARAMS; i++ )
				printf( "source #%d parameter #%d %g\n", s + 1, i + 1, ad.var[i] );
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
	}
	return( c );
}

double func_solver( double x, double y, double z1, double z2, double t, void *data ) // Compute for (x, y, z1, t) and (x, y, z2, t) and average
{
	int i, j, k, s;
	double c1, c2;
	struct calc_data *p = ( struct calc_data * )data;
	struct anal_data ad;
	c1 = c2 = p->c_background;
	k = NUM_ANAL_PARAMS_SOURCE * p->num_solutions + ( NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE );
	j = NUM_ANAL_PARAMS - NUM_ANAL_PARAMS_SOURCE - 1;
	for( i = NUM_ANAL_PARAMS_SOURCE * p->num_solutions; i < k; i++, j++ )
		ad.var[j] = p->var[i];
	if( p->tied ) // Hard wiring; Tied dispersivities
	{
		// printf( "AY: divide %.12g by %.12g\n", ad.var[AX], ad.var[AY] );
		ad.var[AY] = ad.var[AX] / ad.var[AY];
		// printf( "AZ: divide %.12g by %.12g\n", ad.var[AY], ad.var[AZ] );
		ad.var[AZ] = ad.var[AY] / ad.var[AZ];
	}
	for( s = 0; s < p->num_solutions; s++ )
	{
		k = NUM_ANAL_PARAMS_SOURCE * ( s + 1 );
		j = 0;
		for( i = NUM_ANAL_PARAMS_SOURCE * s; i < k; i++, j++ )
			ad.var[j] = p->var[i];
		if( p->fdebug )
			for( i = 0; i < NUM_ANAL_PARAMS; i++ )
				printf( "source #%d parameter #%d %g\n", s + 1, i + 1, ad.var[i] );
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
			//			printf( "trans %s %g %g %g -> ", p->pd->var_id[p->pd->var_index[i]], v[i], p->pd->var_range[k], p->pd->var_min[k] );
			vt[i] = ( v[i] - p->pd->var_min[k] ) / p->pd->var_range[k];
			vt[i] = asin( ( double ) vt[i] * 2.0 - 1.0 );
			//			printf( "%g\n", vt[i] );
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
			//			printf( "detrans %s %g -> %g\n", p->pd->var_id[p->pd->var_index[i]], v[i], vt[i] );
		}
	else
		for( i = 0; i < p->pd->nOptParam; i++ )
		{
			k = p->pd->var_index[i];
			vt[i] = ( ( double ) sin( v[i] ) + 1.0 ) / 2.0;
			vt[i] = p->pd->var_min[k] + vt[i] * p->pd->var_range[k];
			//			printf( "detrans %s %g -> %g\n", p->pd->var_id[p->pd->var_index[i]], v[i], vt[i] );
		}
}
