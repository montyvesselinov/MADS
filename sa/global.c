// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
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
#include <gsl/gsl_monte.h>
#include <gsl/gsl_integration.h>
#include "do_miser.h"
#include "../mads.h"
#include "global.h"

/* Functions here */
int sa_sobol( struct opt_data *op );
int sa_saltelli( struct opt_data *op );
int sa_moat( struct opt_data *op );
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var );
void ave_sorted( double data[], int n, double *ave, double *ep );

/* Functions elsewhere */
int get_seed( );
int Ftest( char *filename );
FILE *Fwrite( char *filename );
char *Fdatetime( char *filename, int debug );

int sa_sobol( struct opt_data *op )
{
	int i, j, k, count, n_sub, n_obs;
	double *opt_params, **var_a_lhs, **var_b_lhs, *var_a_lhs_local, *var_b_lhs_local;
	char filename[255], buf[255];
	FILE *out, *out2;
	strcpy( op->label, "sobol" );
	double *fhat, *fhat2, gfhat, gfhat2, *phis_full, *t1, *t2;
	double *var_y, **f_a, **f_b, **sens_index, **sens_total, D_hat_t, f_hat_0, f_hat_a, f_hat_b, ep;
	//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_sobol, op->pd->nOptParam );
	n_sub = op->cd->nreal / 2; //  number of samples for subsets a and b; set to half of user specified reals
	n_obs = op->od->nTObs + 1; // nmuber of observation + objective function
	out = out2 = NULL;
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( phis_full = ( double * ) malloc( 2 * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( var_a_lhs_local = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( var_b_lhs_local = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( var_a_lhs = double_matrix( n_sub, op->pd->nOptParam ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( var_b_lhs = double_matrix( n_sub, op->pd->nOptParam ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( var_y = ( double * ) malloc( n_obs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( fhat = ( double * ) malloc( n_obs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( fhat2 = ( double * ) malloc( n_obs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( t1 = ( double * ) malloc( n_obs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( t2 = ( double * ) malloc( n_obs * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( f_a = double_matrix( n_sub, n_obs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( f_b = double_matrix( n_sub, n_obs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( sens_index = double_matrix( op->pd->nOptParam, n_obs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( sens_total = double_matrix( op->pd->nOptParam, n_obs ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nGlobal sensitivity analysis (Sobol) using random sampling:\n" );
	tprintf( "Number of analyzed parameters %d\n", op->pd->nOptParam );
	tprintf( "Number of analyzed observations %d\n", n_obs );
	// Create samples
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	// Create samples
	// Sample A
	tprintf( "Random sampling set 1 (variables %d; realizations %d) using ", op->pd->nOptParam, n_sub );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_a_lhs_local, op, 1 );
	tprintf( "done.\n" );
	// Sample B
	tprintf( "Random sampling set 2 (variables %d; realizations %d) using ", op->pd->nOptParam, n_sub );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_b_lhs_local, op, 1 );
	tprintf( "done.\n" );
	// Copy temp lhs vectors to matrices
	for( count = 0; count < n_sub; count++ )
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			var_a_lhs[count][i] = var_a_lhs_local[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
			var_b_lhs[count][i] = var_b_lhs_local[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
		}
	free( var_a_lhs_local ); free( var_b_lhs_local );
	// Output samples to files
	if( op->cd->mdebug )
	{
		sprintf( filename, "%s.sobol.zip", op->root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "%s \"mv %s.sobol.zip %s.sobol_%s.zip >& /dev/null\"", SHELL, op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "%s \"zip -m %s.sobol.zip %s.sobol_set_* >& /dev/null\"", SHELL, op->root, op->root ); system( buf );
		sprintf( buf, "%s \"mv %s.sobol.zip %s.sobol_%s.zip >& /dev/null\"", SHELL, op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.sobol_set_a", op->root ); out = Fwrite( filename );
		sprintf( filename, "%s.sobol_set_b", op->root ); out2 = Fwrite( filename );
		for( count = 0; count < n_sub; count++ )
		{
			for( k = 0; k < op->pd->nOptParam; k++ )
			{
				fprintf( out, "%.15g ", var_a_lhs[count][k] );
				fprintf( out2, "%.15g ", var_b_lhs[count][k] );
			}
			fprintf( out, "\n" );
			fprintf( out2, "\n" );
		}
		fclose( out );
		fclose( out2 );
		tprintf( "Random sampling sets a and b saved in %s.mcrnd_set_a and %s.mcrnd_set_b\n", op->root, op->root );
	}
	if( op->cd->mdebug > 1 )
	{
		sprintf( filename, "%s.sobol.results", op->root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "%s \"mv %s %s.sobol_%s.results >& /dev/null\"", SHELL, filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
	}
	// Accumulate phis into fhat and fhat2 for total output mean and variance
	tprintf( "Computing model outputs to calculate total output mean and variance ... Sample A ... \n" );
	gfhat = gfhat2 = 0;
	for( j = 0; j < n_obs; j++ )
		fhat[j] = fhat2[j] = 0;
	// Compute sample a phis
	for( count = 0; count < n_sub; count++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = var_a_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		gfhat = fhat[0] += op->phi;
		gfhat2 = fhat2[0] += op->phi * op->phi;
		f_a[count][0] = phis_full[count] = op->phi;
		for( j = 1; j < n_obs; j++ )
		{
			fhat[j] += f_a[count][j] = op->od->res[j - 1];
			fhat2[j] += op->od->res[j - 1] * op->od->res[j - 1];
		}
		if( op->cd->mdebug > 1 )
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
	for( j = 0; j < n_obs; j++ )
		var_y[j] = fhat2[j] / n_sub - ( fhat[j] * fhat[j] / ( ( double ) n_sub * n_sub ) );
	f_hat_a = fhat[0] / n_sub;
	D_hat_t = fhat2[0] / n_sub - f_hat_a * f_hat_a;
	tprintf( "Sample A output mean     (simple) = %g\n", f_hat_a );
	tprintf( "Sample A output variance (simple) = %g\n", D_hat_t );
	// Compute sample b phis
	tprintf( "Computing model outputs to calculate total output mean and variance ... Sample B ... \n" );
	for( j = 0; j < n_obs; j++ )
		fhat[j] = fhat2[j] = 0;
	for( count = 0; count < n_sub; count++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = var_b_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		fhat[0] += op->phi;
		gfhat += op->phi;
		fhat2[0] += op->phi * op->phi;
		gfhat2 += op->phi * op->phi;
		f_b[count][0] = phis_full[n_sub + count] = op->phi;
		for( j = 1; j < n_obs; j++ )
		{
			fhat[j] += f_b[count][j] = op->od->res[j - 1];
			fhat2[j] += op->od->res[j - 1] * op->od->res[j - 1];
		}
		if( op->cd->mdebug > 1 )
		{
			// save to results file
			fprintf( out, "%d : ", n_sub + count + 1 ); // counter
			fprintf( out, "%g :", op->phi );
			for( i = 0; i < op->pd->nParam; i++ )
				if( op->pd->var_opt[i] >= 1 )
					fprintf( out, " %.15g", op->pd->var[i] );
			fprintf( out, "\n" );
			fflush( out );
		}
	}
	if( op->cd->mdebug > 1 )
	{
		tprintf( "Global Sensitivity MC results are saved in %s.sobol.results\n", op->root );
		fclose( out );
	}
	// Calculate total output mean and variance based on sample a
	f_hat_b = fhat[0] / n_sub;
	D_hat_t = fhat2[0] / n_sub - f_hat_b * f_hat_b;
	tprintf( "Sample B output mean     (simple) = %g\n", f_hat_b );
	tprintf( "Sample B output variance (simple) = %g\n", D_hat_t );
	f_hat_0 = gfhat / ( 2 * n_sub );
	D_hat_t = gfhat2 / ( 2 * n_sub ) - f_hat_0 * f_hat_0;
	tprintf( "Total output mean     (simple) = %g\n", f_hat_0 );
	tprintf( "Total output variance (simple) = %g\n", D_hat_t );
	f_hat_0 = gsl_stats_mean( phis_full, 1, n_sub * 2 );
	D_hat_t = gsl_stats_variance( phis_full, 1, n_sub * 2 );
	tprintf( "Total output mean        (gsl) = %g\n", f_hat_0 );
	tprintf( "Total output variance    (gsl) = %g\n", D_hat_t );
	ave_sorted( phis_full, n_sub * 2, &f_hat_0, &ep );
	var_sorted( phis_full, phis_full, n_sub * 2, f_hat_0, ep, &D_hat_t );
	tprintf( "Total output mean         (nr) = %g abs 1st moment = %g\n", f_hat_0, ep );
	tprintf( "Total output variance     (nr) = %g\n", D_hat_t );
	free( phis_full );
	/*
	// Collect matrix of phis for fmat_a
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		tprintf( "Processing parameter %d out of %d ... %s ...\n", i + 1, op->pd->nOptParam, op->pd->var_name[k] );
		for( count = 0; count < n_sub; count++ )
		{
			for( j = 0; j < op->pd->nOptParam; j++ )
			{
				k = op->pd->var_index[j];
				if( i == j ) // then select from sample a
					opt_params[j] = op->pd->var[k] = var_a_lhs[count][j];
				else // else select from sample b
					opt_params[j] = op->pd->var[k] = var_b_lhs[count][j];
			}
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			// Save phi to fmat_a
			fmat_a[i][count] = op->phi;
		}
	}
	*/
	tprintf( "Computing sensitivity indices:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		for( j = 0; j < n_obs; j++ )
			t1[j] = t2[j] = 0;
		k = op->pd->var_index[i];
		tprintf( "Processing parameter %d out of %d ... %s ... ", i + 1, op->pd->nOptParam, op->pd->var_name[k] );
		for( count = 0; count < n_sub; count++ )
		{
			for( j = 0; j < op->pd->nOptParam; j++ )
			{
				k = op->pd->var_index[j];
				if( i == j ) // then select from sample b
					opt_params[j] = op->pd->var[k] = var_b_lhs[count][j];
				else // else select from sample a
					opt_params[j] = op->pd->var[k] = var_a_lhs[count][j];
			}
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			t1[0] += pow( f_b[count][0] - op->phi, 2 ); //t1 = sum in eq. 18
			t2[0] += pow( f_a[count][0] - op->phi, 2 ); //t2 = sum in eq. 19
			for( j = 1; j < n_obs; j++ )
			{
				t1[j] += pow( f_b[count][j] - op->od->res[j - 1], 2 ); //t1 = sum in eq. 18
				t2[j] += pow( f_a[count][j] - op->od->res[j - 1], 2 ); //t2 = sum in eq. 19
			}
		}
		for( j = 0; j < n_obs; j++ )
		{
			sens_index[i][j] = ( double ) 1 - t1[j] / ( 2 * n_sub ) / var_y[j]; //var_y * eq. 18
			sens_total[i][j] = t2[j] / ( 2 * n_sub ) / var_y[j]; // var_y * eq. 19
		}
		k = op->pd->var_index[i];
		tprintf( " %g (total) %g\n", sens_index[i][0], sens_total[i][0] );
	}
	tprintf( "done.\n" );
	// Print sensitivity indices
	tprintf( "\nParameter sensitivity indices for objective function:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		tprintf( "%-39s: %g (total) %g\n", op->pd->var_name[k], sens_index[i][0], sens_total[i][0] );
	}
	if( n_obs > 1 )
	{
		sprintf( filename, "%s.sobol_sens_index", op->root );
		out = Fwrite( filename );
		fprintf( out, "OF: " );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fprintf( out, " %g", ( double ) sens_index[i][0] );
		fprintf( out, "\n" );
		for( j = 1; j < n_obs; j++ )
		{
			fprintf( out, "%s: ", op->od->obs_id[j - 1] );
			for( i = 0; i < op->pd->nOptParam; i++ )
				fprintf( out, " %g", ( double ) sens_index[i][j] );
			fprintf( out, "\n" );
		}
		fclose( out );
		sprintf( filename, "%s.sobol_sens_total", op->root );
		out = Fwrite( filename );
		fprintf( out, "OF: " );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fprintf( out, " %g", ( double ) sens_total[i][0] );
		fprintf( out, "\n" );
		for( j = 1; j < n_obs; j++ )
		{
			fprintf( out, "%s: ", op->od->obs_id[j - 1] );
			for( i = 0; i < op->pd->nOptParam; i++ )
				fprintf( out, " %g", ( double ) sens_total[i][j] );
			fprintf( out, "\n" );
		}
		fclose( out );
	}
	free( opt_params ); free( t1 ); free( t2 ); free( fhat ); free( fhat2 );
	free_matrix( ( void ** ) f_a, n_sub );
	free_matrix( ( void ** ) f_b, n_sub );
	free_matrix( ( void ** ) sens_index, op->pd->nOptParam );
	free_matrix( ( void ** ) sens_total, op->pd->nOptParam );
	free_matrix( ( void ** ) var_a_lhs, n_sub );
	free_matrix( ( void ** ) var_b_lhs, n_sub );
	return( 1 );
}

//Computes the marginal pdf after integrating out the variable associated with special_index
double param_pdf_marginal( struct opt_data *op, double *x, int special_index )
{
	int i, k;
	double c = 1.;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( i != special_index ) c /= ( op->pd->var_max[k] - op->pd->var_min[k] );
	}
	return c;
}

double joint_param_pdf( struct opt_data *op )
{
	int i, k;
	double c = 1.;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		c /= ( op->pd->var_max[k] - op->pd->var_min[k] );
	}
	return c;
}

//Computes the pdf of the parameter associated with special_index at special_value conditioned on all the other parameters
double param_pdf_cond( struct opt_data *op, int special_index, double special_value )
{
	int k;
	double c = 1.;
	k = op->pd->var_index[special_index];
	c /= ( op->pd->var_max[k] - op->pd->var_min[k] );
	return c;
}

//Computes the joint pdf of all the parameters condititioned on the special_index-th parameter having value special_value
double joint_param_pdf_cond( struct opt_data *op, int special_index, double special_value )
{
	int i, k;
	double c = 1.;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( i != special_index ) c /= ( op->pd->var_max[k] - op->pd->var_min[k] );
	}
	return c;
}

double param_pdf( struct saltelli_data *salt )
{
	int k;
	k = salt->op->pd->var_index[salt->special_index];
	return 1. / ( salt->op->pd->var_max[k] - salt->op->pd->var_min[k] );
}

double saltelli_mean( double *x, size_t dim, void *params )
{
	int i, k;
	do_gsl_monte_miser_params *p = ( do_gsl_monte_miser_params * ) params;
	struct saltelli_data *salt = p->func_params;
	struct opt_data *op = salt->op;
	double c;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, x[i] ) : x[i] );
	}
	c = func_solver1( op->wd->x[salt->well_index], op->wd->y[salt->well_index], op->wd->z1[salt->well_index], op->wd->obs_time[salt->well_index][salt->obs_index], op->cd );
	c *= joint_param_pdf( op );
	return c;
}

double saltelli_variance( double *x, size_t dim, void *params )
{
	int i, k;
	do_gsl_monte_miser_params *p = ( do_gsl_monte_miser_params * ) params;
	struct saltelli_data *salt = p->func_params;
	struct opt_data *op = salt->op;
	double c;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, x[i] ) : x[i] );
	}
	c = func_solver1( op->wd->x[salt->well_index], op->wd->y[salt->well_index], op->wd->z1[salt->well_index], op->wd->obs_time[salt->well_index][salt->obs_index], op->cd );
	c -= salt->mean;
	c *= c * joint_param_pdf( op );
	return c;
}

double first_order_sensitivity_integrand_integrand( double *x, size_t dim, void *params )
{
	int i, j, k;
	do_gsl_monte_miser_params *p = ( do_gsl_monte_miser_params * ) params;
	struct saltelli_data *salt = ( struct saltelli_data * ) p->func_params;
	struct opt_data *op = salt->op;
	double c;
	for( i = 0, j = 0; i < op->pd->nOptParam; i++ )
	{
		//need to do an index j here that keeps track of how far along the subindices we are...
		if( i != salt->special_index )
		{
			k = op->pd->var_index[i];
			op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, x[j] ) : x[j] );
			j++;
		}
	}
	k = op->pd->var_index[salt->special_index];
	op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, salt->special_value ) : salt->special_value );
	c = func_solver1( op->wd->x[salt->well_index], op->wd->y[salt->well_index], op->wd->z1[salt->well_index], op->wd->obs_time[salt->well_index][salt->obs_index], op->cd );
	c *= joint_param_pdf_cond( op, salt->special_index, salt->special_value );
	return c;
}

double first_order_sensitivity_integrand( double x, void *params )
{
	struct saltelli_data *salt = ( struct saltelli_data * ) params;
	struct opt_data *op = salt->op;
	double *lower_bounds;
	double *upper_bounds;
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_monte_function F;
	do_gsl_monte_miser_state *s;
	do_gsl_monte_miser_params p;
	int i, j, k;
	double cond_mean, err;
	lower_bounds = ( double * ) malloc( ( op->pd->nOptParam - 1 ) * sizeof( double ) );
	upper_bounds = ( double * ) malloc( ( op->pd->nOptParam - 1 ) * sizeof( double ) );
	for( i = 0, j = 0; i < op->pd->nOptParam; i++ )
	{
		if( i != salt->special_index )
		{
			k = op->pd->var_index[i];
			lower_bounds[j] = op->pd->var_min[k];
			upper_bounds[j] = op->pd->var_max[k];
			j++;
		}
	}
	F.f = &first_order_sensitivity_integrand_integrand;
	F.dim = op->pd->nOptParam - 1;
	F.params = &p;
	T = gsl_rng_default;
	r = gsl_rng_alloc( T );
	s = do_gsl_monte_miser_alloc( op->pd->nOptParam - 1 );
	p.func_params = ( void * ) salt;
	salt->special_value = x;
	//compute the mean
	do_gsl_monte_miser_integrate( &F, lower_bounds, upper_bounds, op->pd->nOptParam - 1, lround( pow( op->cd->nreal, ( op->pd->nOptParam - 1. ) / op->pd->nOptParam ) ), r, s, &cond_mean, &err );
	do_gsl_monte_miser_free( s );
	gsl_rng_free( r );
	free( lower_bounds );
	free( upper_bounds );
	return ( cond_mean - salt->mean ) * ( cond_mean - salt->mean ) * param_pdf( salt );
}

double first_order_sensitivity( struct saltelli_data *salt )
{
	struct opt_data *op = salt->op;
	gsl_integration_glfixed_table *table;
	gsl_function F;
	double si;
	int i = salt->special_index;
	table = gsl_integration_glfixed_table_alloc( lround( pow( salt->num_evals, 1. / op->pd->nOptParam ) ) );
	F.function = &first_order_sensitivity_integrand;
	F.params = ( void * ) salt;
	si = gsl_integration_glfixed( &F, salt->lower_bounds[i], salt->upper_bounds[i], table );
	//si -= salt->mean * salt->mean;
	si /= salt->variance;
	gsl_integration_glfixed_table_free( table );
	return si;
}

double total_effect_integrand_integrand( double x, void *params )
{
	struct saltelli_data *salt = ( struct saltelli_data * ) params;
	struct opt_data *op = salt->op;
	int k;
	double c;
	k = op->pd->var_index[salt->special_index];
	op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, x ) : x );
	c = func_solver1( op->wd->x[salt->well_index], op->wd->y[salt->well_index], op->wd->z1[salt->well_index], op->wd->obs_time[salt->well_index][salt->obs_index], op->cd );
	c -= salt->cond_mean;
	c *= c * param_pdf_cond( op, salt->special_index, salt->special_value );
	return c;
}

double total_effect_integrand_cond_mean( double x, void *params )
{
	struct saltelli_data *salt = ( struct saltelli_data * ) params;
	struct opt_data *op = salt->op;
	int k;
	double c;
	k = op->pd->var_index[salt->special_index];
	op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, x ) : x );
	c = func_solver1( op->wd->x[salt->well_index], op->wd->y[salt->well_index], op->wd->z1[salt->well_index], op->wd->obs_time[salt->well_index][salt->obs_index], op->cd );
	c *= param_pdf_cond( op, salt->special_index, salt->special_value );
	return c;
}

double total_effect_integrand( double *x, size_t dim, void *params )
{
	do_gsl_monte_miser_params *p = ( do_gsl_monte_miser_params * ) params;
	struct saltelli_data *salt = ( struct saltelli_data * ) p->func_params;
	struct opt_data *op = salt->op;
	gsl_integration_glfixed_table *table;
	gsl_function F;
	double cond_var;
	int i, j, k;
	for( i = 0, j = 0; i < op->pd->nOptParam; i++ )
	{
		if( i != salt->special_index )
		{
			k = op->pd->var_index[i];
			op->cd->var[k] = ( op->pd->var_log[k] ? pow( 10, x[j] ) : x[j] );
			j++;
		}
	}
	table = gsl_integration_glfixed_table_alloc( lround( pow( salt->num_evals, 1. / op->pd->nOptParam ) ) );
	F.params = ( void * ) salt;
	i = salt->special_index;
	F.function = &total_effect_integrand_cond_mean;
	salt->cond_mean = gsl_integration_glfixed( &F, salt->lower_bounds[i], salt->upper_bounds[i], table );
	F.function = &total_effect_integrand_integrand;
	cond_var = gsl_integration_glfixed( &F, salt->lower_bounds[i], salt->upper_bounds[i], table );
	gsl_integration_glfixed_table_free( table );
	return cond_var * param_pdf_marginal( op, x, salt->special_index );
}

double total_effect( struct saltelli_data *salt )
{
	struct opt_data *op = salt->op;
	double *lower_bounds;
	double *upper_bounds;
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_monte_function F;
	do_gsl_monte_miser_state *s;
	do_gsl_monte_miser_params p;
	int i, j, k;
	double expected_variance, err;
	lower_bounds = ( double * ) malloc( ( op->pd->nOptParam - 1 ) * sizeof( double ) );
	upper_bounds = ( double * ) malloc( ( op->pd->nOptParam - 1 ) * sizeof( double ) );
	for( i = 0, j = 0; i < op->pd->nOptParam; i++ )
	{
		if( i != salt->special_index )
		{
			k = op->pd->var_index[i];
			lower_bounds[j] = op->pd->var_min[k];
			upper_bounds[j] = op->pd->var_max[k];
			j++;
		}
	}
	F.f = &total_effect_integrand;
	F.dim = op->pd->nOptParam - 1;
	F.params = &p;
	T = gsl_rng_default;
	r = gsl_rng_alloc( T );
	s = do_gsl_monte_miser_alloc( op->pd->nOptParam - 1 );
	p.func_params = ( void * ) salt;
	//compute the mean
	do_gsl_monte_miser_integrate( &F, lower_bounds, upper_bounds, op->pd->nOptParam - 1, lround( pow( op->cd->nreal, ( op->pd->nOptParam - 1. ) / op->pd->nOptParam ) ), r, s, &expected_variance, &err );
	do_gsl_monte_miser_free( s );
	gsl_rng_free( r );
	free( lower_bounds );
	free( upper_bounds );
	//printf("cond_mean: %g, param_pdf: %g, sv: %g\n", cond_mean, param_pdf( salt ), salt->special_value );
	//return cond_mean * cond_mean * param_pdf( salt );
	return expected_variance / salt->variance;
}

int sa_saltelli( struct opt_data *op )
{
	struct saltelli_data salt;
	double *lower_bounds;
	double *upper_bounds;
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_monte_function F;
	do_gsl_monte_miser_state *s;
	do_gsl_monte_miser_params p;
	int i, j, k, n;
	double err;
	double mean, variance;
	double si, ti;//first order sensitivity and total effect
	salt.op = op;
	if( salt.op->pd->nOptParam < 2 ) { fprintf( stderr, "You must have at least 2 \"opt\" params to do saltelli.\n" ); mads_quits( op->root ); }
	salt.num_evals = op->cd->nreal;
	lower_bounds = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) );
	upper_bounds = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) );
	salt.lower_bounds = lower_bounds;
	salt.upper_bounds = upper_bounds;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		lower_bounds[i] = op->pd->var_min[k];
		upper_bounds[i] = op->pd->var_max[k];
	}
	F.dim = op->pd->nOptParam;
	F.params = &p;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc( T );
	s = do_gsl_monte_miser_alloc( op->pd->nOptParam );
	p.func_params = ( void * ) &salt;
	for( j = 0; j < op->wd->nW; j++ )//loop through the wells
	{
		for( n = 0; n < op->wd->nWellObs[j]; n++ )//loop through the observations
		{
			if( op->wd->obs_weight[j][n] < 0 )//making the weight less than zero is a trick to "flag" the observation
			{
				tprintf( "%s at t=%g:\n", op->wd->id[j], op->wd->obs_time[j][n] );
				salt.well_index = j;
				salt.obs_index = n;
				F.f = &saltelli_mean;
				//compute the mean
				do_gsl_monte_miser_integrate( &F, lower_bounds, upper_bounds, op->pd->nOptParam, op->cd->nreal, r, s, &mean, &err );
				salt.mean = mean;
				F.f = &saltelli_variance;
				//compute the variance
				do_gsl_monte_miser_integrate( &F, lower_bounds, upper_bounds, op->pd->nOptParam, op->cd->nreal, r, s, &variance, &err );
				salt.variance = variance;
				tprintf( "mean: %g\nvariance: %g\n", mean, variance );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					salt.special_index = i;
					si = first_order_sensitivity( &salt );
					ti = total_effect( &salt );
					k = op->pd->var_index[i];
					tprintf( "%-39s: %g (total) %g\n", op->pd->var_name[k], si, ti );
				}
				printf( "\n" );
			}
		}
	}
	gsl_rng_free( r );
	do_gsl_monte_miser_free( s );
	free( lower_bounds );
	free( upper_bounds );
	return 1;
}

int sa_moat( struct opt_data *op )
{
	tprintf( "MOAT\n" );
	return( 1 );
}

// Modified from Numerical Recipes in C: The Art of Scientific Computing (ISBN 0-521-43108-5)
// corrected three-pass algorithm to minimize roundoff error in variance
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var )
{
	int j;
	double *dev2;
	if( ( dev2 = ( double * ) malloc( n * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return; }
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
	free( dev2 );
}

void ave_sorted( double data[], int n, double *ave, double *ep )
{
	int j;
	double s, *dev;
	if( ( dev = ( double * ) malloc( n * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return; }
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
	free( dev );
}
