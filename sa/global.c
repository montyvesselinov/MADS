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
#include "../mads.h"

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
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.sobol.zip %s.sobol_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.sobol.zip %s.sobol_set_* >& /dev/null", op->root, op->root ); system( buf );
		sprintf( buf, "mv %s.sobol.zip %s.sobol_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
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
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.sobol_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
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
			f_a[count][j] = fhat[j] = op->od->res[j - 1];
			fhat2[j] = op->od->res[j - 1] * op->od->res[j - 1];
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
		gfhat = fhat[0] += op->phi;
		gfhat2 = fhat2[0] += op->phi * op->phi;
		f_b[count][0] = phis_full[n_sub + count] = op->phi;
		for( j = 1; j < n_obs; j++ )
		{
			f_b[count][j] = fhat[j] = op->od->res[j - 1];
			fhat2[j] = op->od->res[j - 1] * op->od->res[j - 1];
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

int sa_saltelli( struct opt_data *op )
{
	int num_opt_params;
	int num_samples;
	int num_samples_per_param;
	int i, j, k, n, m;
	double *func_evals;
	double *opt_params;
	double var_value;
	double mean;
	double variance;
	double ep;//used by ave_sorted and var_sorted
	num_opt_params = op->pd->nOptParam;
	//Round up on the number of samples per param
	num_samples_per_param = ceil( pow( op->cd->nreal, 1. / num_opt_params ) );
	num_samples = pow( num_samples_per_param, num_opt_params );
	if( num_samples != op->cd->nreal ) { tprintf( "You requested %d samples. You got %d samples. Life is cruel.\n", op->cd->nreal, num_samples ); }
	if( ( func_evals = ( double * ) malloc( num_samples * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return 0; }
	if( ( opt_params = ( double * ) malloc( num_opt_params * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return 0; }
	//do all the function evaluations
	for( i = 0; i < num_samples; i++ )
	{
		n = i;
		//map i into num_opt_params dimensional indicies
		for( j = 0; j < num_opt_params; j++ )
		{
			k = op->pd->var_index[j];
			m = n % num_samples_per_param;
			n = ( n - m ) / num_opt_params;
			var_value = op->pd->var_min[k] + ( m + .5 ) / num_samples_per_param * op->pd->var_range[k];
			opt_params[j] = op->pd->var[k] = var_value;
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		func_evals[i] = op->phi;
	}
	//compute the mean and variance
	ave_sorted( func_evals, num_samples, &mean, &ep );
	mean /= ( num_samples - 1 );
	var_sorted( func_evals, func_evals, num_samples, mean, ep, &variance );
	//compute the first order sensitivities
	//compute the total sensitivities
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
