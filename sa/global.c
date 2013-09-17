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

/* Functions elsewhere */
int get_seed( );
int Ftest( char *filename );
FILE *Fwrite( char *filename );
char *Fdatetime( char *filename, int debug );

int sa_sobol( struct opt_data *op )
{
	int i, j, k, count;
	double *opt_params, *var_a_lhs, *var_b_lhs;
	char filename[255], buf[255];
	FILE *out, *out2;
	struct gsa_data gs;
	strcpy( op->label, "sobol" );
	double fhat, fhat2, *phis_full, *phis_half;
	int n_sub; //! number of samples for subsets a and b
	//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_sobol, op->pd->nOptParam );
	n_sub = op->cd->nreal / 2;	// set to half of user specified reals
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store op->cd->nreal phis
	if( ( phis_full = ( double * ) malloc( op->cd->nreal * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store m_sub phis
	if( ( phis_half = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store random sample a
	if( ( var_a_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Sample a phis
	if( ( gs.f_a = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Sample b phis
	if( ( gs.f_b = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store random sample b
	if( ( var_b_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// matrices to store lhs samples
	gs.var_a_lhs = double_matrix( n_sub, op->pd->nOptParam );
	gs.var_b_lhs = double_matrix( n_sub, op->pd->nOptParam );
	// Matrices to store phis with different combinations of parameters from samples a and b
	if( ( gs.fmat_a = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL )
	{ tprintf( "Error creating 3D matrix\n" ); return( 0 ); }
	if( ( gs.fmat_b = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL )
	{ tprintf( "Error creating 3D matrix\n" ); return( 0 ); }
	// Vector of variances for individual component contribution
	if( ( gs.D_hat = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Vector of variances for total component contribution
	if( ( gs.D_hat_n = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nGlobal sensitivity analysis (Sobol) using random sampling:\n" );
	// Create samples
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	tprintf( "Random sampling set 1 (variables %d; realizations %d) using ", op->pd->nOptParam, op->cd->nreal );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_a_lhs, op, 1 );
	tprintf( "done.\n" );
	tprintf( "Random sampling set 2 (variables %d; realizations %d) using ", op->pd->nOptParam, op->cd->nreal );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_b_lhs, op, 1 );
	tprintf( "done.\n" );
	// Create samples using Sobol's quasi-random sequence
	/*		for( count = 0; count < n_sub; count++ )
			{
				double v[ op->pd->nOptParam ];
				gsl_qrng_get( q, v);
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					gs.var_a_lhs[count][i] = v[i] * op->pd->var_range[k] + op->pd->var_min[k];
				}
			}

			for( count = 0; count < n_sub; count++ )
			{
				double v[ op->pd->nOptParam ];
				gsl_qrng_get( q, v);
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					gs.var_b_lhs[count][i] = v[i] * op->pd->var_range[k] + op->pd->var_min[k];
				}
			}*/
	// Copy temp lhs vectors to matrices
	for( count = 0; count < n_sub; count++ )
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			gs.var_a_lhs[count][i] = var_a_lhs[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
			gs.var_b_lhs[count][i] = var_b_lhs[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
		}
	free( var_a_lhs );
	free( var_b_lhs );
	// Output samples to files
	if( op->cd->mdebug )
	{
		sprintf( filename, "%s.sobol.zip", op->root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.sobol.zip %s.sobol_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.sobol.zip %s.sobol_set_* >& /dev/null", op->root, op->root ); system( buf );
		sprintf( buf, "mv %s.sobol.zip %s.sobol_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.sobol_set_a", op->root ); out = Fwrite( filename );
		sprintf( filename, "%s.sobol_set_b", op->root ); out2 = Fwrite( filename );
		for( count = 0; count < n_sub; count ++ )
		{
			for( k = 0; k < op->pd->nOptParam; k++ )
			{
				fprintf( out, "%.15g ", gs.var_a_lhs[count][k] );
				fprintf( out2, "%.15g ", gs.var_b_lhs[count][k] );
			}
			fprintf( out, "\n" );
			fprintf( out2, "\n" );
		}
		fclose( out );
		fclose( out2 );
		tprintf( "Random sampling sets a and b saved in %s.mcrnd_set_a and %s.mcrnd_set_b\n", op->root, op->root );
	}
	sprintf( filename, "%s.sobol.results", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.sobol_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	out = Fwrite( filename );
	// Accumulate phis into fhat and fhat2 for total output mean and variance
	fhat = fhat2 = 0;
	tprintf( "Computing phis to calculate total output mean and variance...\n" );
	// Compute sample a phis
	for( count = 0; count < n_sub; count ++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_a_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		fhat += op->phi;
		fhat2 += pow( op->phi, 2 );
		// Save sample a phis
		gs.f_a[count] = op->phi;
		phis_full[count] = op->phi;
		// save to results file
		fprintf( out, "%d : ", count + 1 ); // counter
		fprintf( out, "%g :", op->phi );
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] >= 1 )
				fprintf( out, " %.15g", op->pd->var[i] );
		fprintf( out, "\n" );
		fflush( out );
	}
	// Compute sample b phis
	for( count = 0; count < n_sub; count ++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_b_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		fhat += op->phi;
		fhat2 += pow( op->phi, 2 );
		// Save sample b phis
		gs.f_b[count] = op->phi;
		phis_full[ n_sub + count ] = op->phi;
		// save to results file
		fprintf( out, "%d : ", n_sub + count ); // counter
		fprintf( out, "%g :", op->phi );
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] >= 1 )
				fprintf( out, " %.15g", op->pd->var[i] );
		fprintf( out, "\n" );
		fflush( out );
	}
	fclose( out );
	tprintf( "Global Sensitivity MC results are saved in %s.sobol.results\n", op->root );
	// Calculate total output mean and variance based on sample a
	gs.f_hat_0 = fhat / ( 2 * n_sub );
	gs.D_hat_t = fhat2 / ( 2 * n_sub ) - gs.f_hat_0 * gs.f_hat_0;
	tprintf( "Total output mean = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gsl_stats_mean( phis_full, 1, op->cd->nreal );
	gs.D_hat_t = gsl_stats_variance( phis_full, 1, op->cd->nreal );
	tprintf( "Total output mean = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gs.D_hat_t = 0.0;
	ave_sorted( phis_full, op->cd->nreal, &gs.f_hat_0, &gs.ep );
	tprintf( "Total output mean = %g abs 1st moment = %g\n", gs.f_hat_0, gs.ep );
	var_sorted( phis_full, phis_full, op->cd->nreal, gs.f_hat_0, gs.ep, &gs.D_hat_t );
	tprintf( "Total output variance = %g\n", gs.D_hat_t );
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
		 	tprintf( "Total output variance = %g\n", gs.D_hat_t );
			gs.D_hat_t = gsl_stats_variance( phis_full, 1, op->cd->nreal );
		 	tprintf( "Total output variance = %g\n", gs.D_hat_t );
	 */		free( phis_full );
	// Collect matrix of phis for fmat_a
	tprintf( "Computing phis for calculation of individual output variances:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		tprintf( "Parameter %d...\n", i + 1 );
		for( count = 0; count < n_sub; count ++ )
		{
			for( j = 0; j < op->pd->nOptParam; j++ )
			{
				k = op->pd->var_index[j];
				if( i == j ) // then select from sample a
					opt_params[j] = op->pd->var[k] = gs.var_a_lhs[count][j];
				else // else select from sample b
					opt_params[j] = op->pd->var[k] = gs.var_b_lhs[count][j];
			}
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			// Save phi to fmat_a
			gs.fmat_a[i][count] = op->phi;
		}
	}
	// Collect matrix of phis for fmat_b
	tprintf( "Computing phis for calculation of individual plus interaction output variances:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		tprintf( "Parameter %d...\n", i + 1 );
		for( count = 0; count < n_sub; count ++ )
		{
			for( j = 0; j < op->pd->nOptParam; j++ )
			{
				k = op->pd->var_index[j];
				if( i == j ) // then select from sample b
					opt_params[j] = op->pd->var[k] = gs.var_b_lhs[count][j];
				else // else select from sample a
					opt_params[j] = op->pd->var[k] = gs.var_a_lhs[count][j];
			}
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			// Save phi to fmat_b
			gs.fmat_b[i][count] = op->phi;
		}
	}
	tprintf( "done.\n" );
	// Calculate individual and interaction output variances
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		fhat2 = 0;
		for( j = 0; j < n_sub; j++ )
		{
			fhat2 += ( gs.f_a[j] * gs.fmat_a[i][j] );
			phis_half[ j ] = ( gs.f_a[j] * gs.fmat_a[i][j] );
		}
		gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		gs.D_hat[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		gs.D_hat[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_a[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		var_sorted( gs.f_a, gs.fmat_a[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat[i] );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		//gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		fhat2 = 0;
		for( j = 0; j < n_sub; j++ )
		{
			fhat2 += ( gs.f_a[j] * gs.fmat_b[i][j] );
			phis_half[ j ] = ( gs.f_a[j] * gs.fmat_b[i][j] );
		}
		gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		gs.D_hat_n[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		gs.D_hat_n[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_b[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		var_sorted( gs.f_a, gs.fmat_b[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_n[i] );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		//gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
	}
	// Print sensitivity indices
	tprintf( "\nParameter sensitivity indices:\n" );
	tprintf( "parameter individual interaction\n" );
	for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( "%d %g %g\n", i + 1, gs.D_hat[i] / gs.D_hat_t, 1 - ( gs.D_hat_n[i] / gs.D_hat_t ) );
	tprintf( "\n" );
	free( opt_params ); free( phis_half ); free( gs.f_a ); free( gs.f_b ); free( gs.D_hat ); free( gs.D_hat_n );
	free_matrix( ( void ** ) gs.var_a_lhs, n_sub );
	free_matrix( ( void ** ) gs.var_b_lhs, n_sub );
	free_matrix( ( void ** ) gs.fmat_a, op->pd->nOptParam );
	free_matrix( ( void ** ) gs.fmat_b, op->pd->nOptParam );
	return( 1 );
}

int sa_saltelli( struct opt_data *op )
{
	int i, j, k, count;
	double *opt_params, *var_a_lhs, *var_b_lhs;
	char filename[255], buf[255];
	FILE *out, *out2;
	struct gsa_data gs;
	strcpy( op->label, "saltelli" );
	double fhat, fhat2, *phis_full, *phis_half;
	int n_sub; //! number of samples for subsets a and b
	//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_saltelli, op->pd->nOptParam );
	n_sub = op->cd->nreal / 2;	// set to half of user specified reals
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store op->cd->nreal phis
	if( ( phis_full = ( double * ) malloc( op->cd->nreal * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store m_sub phis
	if( ( phis_half = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store random sample a
	if( ( var_a_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Sample a phis
	if( ( gs.f_a = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Sample b phis
	if( ( gs.f_b = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Temporary variable to store random sample b
	if( ( var_b_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// matrices to store lhs samples
	gs.var_a_lhs = double_matrix( n_sub, op->pd->nOptParam );
	gs.var_b_lhs = double_matrix( n_sub, op->pd->nOptParam );
	// Matrices to store phis with different combinations of parameters from samples a and b
	if( ( gs.fmat_a = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL )
	{ tprintf( "Error creating 3D matrix\n" ); return( 0 ); }
	if( ( gs.fmat_b = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL )
	{ tprintf( "Error creating 3D matrix\n" ); return( 0 ); }
	// Vector of variances for individual component contribution
	if( ( gs.D_hat = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	// Vector of variances for total component contribution
	if( ( gs.D_hat_n = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nGlobal sensitivity analysis (Saltelli) using random sampling:\n" );
	// Create samples
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	tprintf( "Random sampling set 1 (variables %d; realizations %d) using ", op->pd->nOptParam, op->cd->nreal );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_a_lhs, op, 1 );
	tprintf( "done.\n" );
	tprintf( "Random sampling set 2 (variables %d; realizations %d) using ", op->pd->nOptParam, op->cd->nreal );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_b_lhs, op, 1 );
	tprintf( "done.\n" );
	// Create samples using Saltelli's quasi-random sequence
	/*		for( count = 0; count < n_sub; count++ )
			{
				double v[ op->pd->nOptParam ];
				gsl_qrng_get( q, v);
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					gs.var_a_lhs[count][i] = v[i] * op->pd->var_range[k] + op->pd->var_min[k];
				}
			}

			for( count = 0; count < n_sub; count++ )
			{
				double v[ op->pd->nOptParam ];
				gsl_qrng_get( q, v);
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					gs.var_b_lhs[count][i] = v[i] * op->pd->var_range[k] + op->pd->var_min[k];
				}
			}*/
	// Copy temp lhs vectors to matrices
	for( count = 0; count < n_sub; count++ )
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			gs.var_a_lhs[count][i] = var_a_lhs[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
			gs.var_b_lhs[count][i] = var_b_lhs[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
		}
	free( var_a_lhs );
	free( var_b_lhs );
	// Output samples to files
	if( op->cd->mdebug )
	{
		sprintf( filename, "%s.saltelli.zip", op->root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.saltelli.zip %s.saltelli_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.saltelli.zip %s.saltelli_set_* >& /dev/null", op->root, op->root ); system( buf );
		sprintf( buf, "mv %s.saltelli.zip %s.saltelli_%s.zip >& /dev/null", op->root, op->root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.saltelli_set_a", op->root ); out = Fwrite( filename );
		sprintf( filename, "%s.saltelli_set_b", op->root ); out2 = Fwrite( filename );
		for( count = 0; count < n_sub; count ++ )
		{
			for( k = 0; k < op->pd->nOptParam; k++ )
			{
				fprintf( out, "%.15g ", gs.var_a_lhs[count][k] );
				fprintf( out2, "%.15g ", gs.var_b_lhs[count][k] );
			}
			fprintf( out, "\n" );
			fprintf( out2, "\n" );
		}
		fclose( out );
		fclose( out2 );
		tprintf( "Random sampling sets a and b saved in %s.mcrnd_set_a and %s.mcrnd_set_b\n", op->root, op->root );
	}
	sprintf( filename, "%s.saltelli.results", op->root );
	if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.saltelli_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
	out = Fwrite( filename );
	// Accumulate phis into fhat and fhat2 for total output mean and variance
	fhat = fhat2 = 0;
	tprintf( "Computing phis to calculate total output mean and variance...\n" );
	// Compute sample a phis
	for( count = 0; count < n_sub; count ++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_a_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		fhat += op->phi;
		fhat2 += pow( op->phi, 2 );
		// Save sample a phis
		gs.f_a[count] = op->phi;
		phis_full[count] = op->phi;
		// save to results file
		fprintf( out, "%d : ", count + 1 ); // counter
		fprintf( out, "%g :", op->phi );
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] >= 1 )
				fprintf( out, " %.15g", op->pd->var[i] );
		fprintf( out, "\n" );
		fflush( out );
	}
	// Compute sample b phis
	for( count = 0; count < n_sub; count ++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_b_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		fhat += op->phi;
		fhat2 += pow( op->phi, 2 );
		// Save sample b phis
		gs.f_b[count] = op->phi;
		phis_full[ n_sub + count ] = op->phi;
		// save to results file
		fprintf( out, "%d : ", n_sub + count ); // counter
		fprintf( out, "%g :", op->phi );
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] >= 1 )
				fprintf( out, " %.15g", op->pd->var[i] );
		fprintf( out, "\n" );
		fflush( out );
	}
	fclose( out );
	tprintf( "Global Sensitivity MC results are saved in %s.saltelli.results\n", op->root );
	// Calculate total output mean and variance based on sample a
	gs.f_hat_0 = fhat / ( 2 * n_sub );
	gs.D_hat_t = fhat2 / ( 2 * n_sub ) - gs.f_hat_0;
	tprintf( "Total output mean = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gsl_stats_mean( phis_full, 1, op->cd->nreal );
	gs.D_hat_t = gsl_stats_variance( phis_full, 1, op->cd->nreal );
	tprintf( "Total output mean = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gs.D_hat_t = 0.0;
	ave_sorted( phis_full, op->cd->nreal, &gs.f_hat_0, &gs.ep );
	tprintf( "Total output mean = %g abs 1st moment = %g\n", gs.f_hat_0, gs.ep );
	var_sorted( phis_full, phis_full, op->cd->nreal, gs.f_hat_0, gs.ep, &gs.D_hat_t );
	tprintf( "Total output variance = %g\n", gs.D_hat_t );
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
		 	tprintf( "Total output variance = %g\n", gs.D_hat_t );
			gs.D_hat_t = gsl_stats_variance( phis_full, 1, op->cd->nreal );
		 	tprintf( "Total output variance = %g\n", gs.D_hat_t );
	 */		free( phis_full );
	// Collect matrix of phis for fmat_a
	tprintf( "Computing phis for calculation of individual output variances:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		tprintf( "Parameter %d...\n", i + 1 );
		for( count = 0; count < n_sub; count ++ )
		{
			for( j = 0; j < op->pd->nOptParam; j++ )
			{
				k = op->pd->var_index[j];
				if( i == j ) // then select from sample a
					opt_params[j] = op->pd->var[k] = gs.var_a_lhs[count][j];
				else // else select from sample b
					opt_params[j] = op->pd->var[k] = gs.var_b_lhs[count][j];
			}
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			// Save phi to fmat_a
			gs.fmat_a[i][count] = op->phi;
		}
	}
	// Collect matrix of phis for fmat_b
	tprintf( "Computing phis for calculation of individual plus interaction output variances:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		tprintf( "Parameter %d...\n", i + 1 );
		for( count = 0; count < n_sub; count ++ )
		{
			for( j = 0; j < op->pd->nOptParam; j++ )
			{
				k = op->pd->var_index[j];
				if( i == j ) // then select from sample b
					opt_params[j] = op->pd->var[k] = gs.var_b_lhs[count][j];
				else // else select from sample a
					opt_params[j] = op->pd->var[k] = gs.var_a_lhs[count][j];
			}
			Transform( opt_params, op, opt_params );
			func_global( opt_params, op, op->od->res );
			// Save phi to fmat_b
			gs.fmat_b[i][count] = op->phi;
		}
	}
	tprintf( "done.\n" );
	// Calculate individual and interaction output variances
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		fhat2 = 0;
		for( j = 0; j < n_sub; j++ )
		{
			fhat2 += ( gs.f_a[j] * gs.fmat_a[i][j] );
			phis_half[ j ] = ( gs.f_a[j] * gs.fmat_a[i][j] );
		}
		gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		gs.D_hat[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		gs.D_hat[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_a[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		var_sorted( gs.f_a, gs.fmat_a[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat[i] );
		tprintf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
		//gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		fhat2 = 0;
		for( j = 0; j < n_sub; j++ )
		{
			fhat2 += ( gs.f_a[j] * gs.fmat_b[i][j] );
			phis_half[ j ] = ( gs.f_a[j] * gs.fmat_b[i][j] );
		}
		gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		gs.D_hat_n[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		gs.D_hat_n[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_b[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		var_sorted( gs.f_a, gs.fmat_b[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_n[i] );
		tprintf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
		//gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
	}
	// Print sensitivity indices
	tprintf( "\nParameter sensitivity indices:\n" );
	tprintf( "parameter individual interaction\n" );
	for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( "%d %g %g\n", i + 1, gs.D_hat[i] / gs.D_hat_t, 1 - ( gs.D_hat_n[i] / gs.D_hat_t ) );
	tprintf( "\n" );
	free( opt_params ); free( phis_half ); free( gs.f_a ); free( gs.f_b ); free( gs.D_hat ); free( gs.D_hat_n );
	free_matrix( ( void ** ) gs.var_a_lhs, n_sub );
	free_matrix( ( void ** ) gs.var_b_lhs, n_sub );
	free_matrix( ( void ** ) gs.fmat_a, op->pd->nOptParam );
	free_matrix( ( void ** ) gs.fmat_b, op->pd->nOptParam );
	return( 1 );
}

int sa_moat( struct opt_data *op )
{
	tprintf( "MOAT\n" );
	return( 1 );
}
