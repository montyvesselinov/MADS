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

int sa_sobol_dh( struct opt_data *op )
{
	int i, j, k, count;
	double *opt_params, *var_a_lhs, *var_b_lhs;
	char filename[255], buf[255];
	FILE *out, *out2;
	struct gsa_data gs;
	strcpy( op->label, "sobol" );
	double fhat, fhat2, gfhat, gfhat2, *phis_full, *phis_half, tmp;
	int n_sub; //! number of samples for subsets a and b
	//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_sobol, op->pd->nOptParam );
	n_sub = op->cd->nreal / 2;	// set to half of user specified reals
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store op->cd->nreal phis
	if( ( phis_full = ( double * ) malloc( op->cd->nreal * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store m_sub phis
	if( ( phis_half = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store random sample a
	if( ( var_a_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Sample a phis
	if( ( var_b_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // matrices to store lhs samples
	if( ( gs.f_a = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Sample b phis
	if( ( gs.f_b = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store random sample b
	gs.var_a_lhs = double_matrix( n_sub, op->pd->nOptParam );
	gs.var_b_lhs = double_matrix( n_sub, op->pd->nOptParam ); // Matrices to store phis with different combinations of parameters from samples a and b
	if( ( gs.fmat_a = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( gs.fmat_b = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Vector of variances for individual component contribution
	if( ( gs.D_hat_a = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Vector of variances for total component contribution
	if( ( gs.D_hat_b = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nGlobal sensitivity analysis (Sobol) using random sampling:\n" );
	// Create samples
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	// Create samples
	// Sample A
	tprintf( "Random sampling set 1 (variables %d; realizations %d) using ", op->pd->nOptParam, op->cd->nreal );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_a_lhs, op, 1 );
	tprintf( "done.\n" );
	// Sample B
	tprintf( "Random sampling set 2 (variables %d; realizations %d) using ", op->pd->nOptParam, op->cd->nreal );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_b_lhs, op, 1 );
	tprintf( "done.\n" );
	// Create samples using Sobol's quasi-random sequence
	/*		for( count = 0; count < n_sub; count++ )
			{
				double v[ op->pd->nOptParam ];
				gsl_qrng_get( q, v );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					gs.var_a_lhs[count][i] = v[i] * op->pd->var_range[k] + op->pd->var_min[k];
				}
			}

			for( count = 0; count < n_sub; count++ )
			{
				double v[ op->pd->nOptParam ];
				gsl_qrng_get( q, v );
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
	free( var_a_lhs ); free( var_b_lhs );
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
	tprintf( "Computing phis to calculate total output mean and variance...\n" );
	gfhat = gfhat2 = fhat = fhat2 = 0;
	// Compute sample a phis
	for( count = 0; count < n_sub; count++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_a_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		gfhat = fhat += op->phi;
		gfhat2 = fhat2 += op->phi * op->phi;
		// Save sample a phis
		gs.f_a[count] = phis_full[count] = op->phi;
		// save to results file
		fprintf( out, "%d : ", count + 1 ); // counter
		fprintf( out, "%g :", op->phi );
		for( i = 0; i < op->pd->nParam; i++ )
			if( op->pd->var_opt[i] >= 1 )
				fprintf( out, " %.15g", op->pd->var[i] );
		fprintf( out, "\n" );
		fflush( out );
	}
	gs.f_hat_a = fhat / n_sub;
	gs.D_hat_t = fhat2 / n_sub - gs.f_hat_a * gs.f_hat_a;
	tprintf( "Sample A output mean     (simple) = %g\n", gs.f_hat_a );
	tprintf( "Sample A output variance (simple) = %g\n", gs.D_hat_t );
	// Compute sample b phis
	fhat = fhat2 = 0;
	for( count = 0; count < n_sub; count++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_b_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		gfhat = fhat += op->phi;
		gfhat2 = fhat2 += op->phi * op->phi;
		// Save sample b phis
		gs.f_b[count] = phis_full[n_sub + count] = op->phi;
		// save to results file
		fprintf( out, "%d : ", n_sub + count + 1 ); // counter
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
	gs.f_hat_b = fhat / n_sub;
	gs.D_hat_t = fhat2 / n_sub - gs.f_hat_b * gs.f_hat_b;
	tprintf( "Sample B output mean     (simple) = %g\n", gs.f_hat_b );
	tprintf( "Sample B output variance (simple) = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gfhat / ( 2 * n_sub );
	gs.D_hat_t = gfhat2 / ( 2 * n_sub ) - gs.f_hat_0 * gs.f_hat_0;
	tprintf( "Total output mean     (simple) = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance (simple) = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gsl_stats_mean( phis_full, 1, n_sub * 2 );
	gs.D_hat_t = gsl_stats_variance( phis_full, 1, n_sub * 2 );
	tprintf( "Total output mean        (gsl) = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance    (gsl) = %g\n", gs.D_hat_t );
	ave_sorted( phis_full, n_sub * 2, &gs.f_hat_0, &gs.ep );
	var_sorted( phis_full, phis_full, n_sub * 2, gs.f_hat_0, gs.ep, &gs.D_hat_t );
	tprintf( "Total output mean         (nr) = %g abs 1st moment = %g\n", gs.f_hat_0, gs.ep );
	tprintf( "Total output variance     (nr) = %g\n", gs.D_hat_t );
	/*
	// Subtract f_hat_0 from phis and recalculate total output variance
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
	*/
	free( phis_full );
	// Collect matrix of phis for fmat_a
	tprintf( "Computing phis for calculation of individual output variances:\n" );
	// gfhat = gfhat2 = 0;
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		tprintf( "Parameter %d...\n", i + 1 );
		for( count = 0; count < n_sub; count++ )
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
			gfhat += op->phi;
			gfhat2 += op->phi * op->phi;
		}
	}
	// Collect matrix of phis for fmat_b
	tprintf( "Computing phis for calculation of individual plus interaction output variances:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		tprintf( "Parameter %d...\n", i + 1 );
		for( count = 0; count < n_sub; count++ )
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
			gfhat += op->phi;
			gfhat2 += op->phi * op->phi;
		}
	}
	tprintf( "done.\n" );
	// gs.f_hat_0 = gfhat / ( 2 * n_sub * ( op->pd->nOptParam + 1 ) );
	// gs.D_hat_t = gfhat2 / ( 2 * n_sub * ( op->pd->nOptParam + 1 ) ) - gs.f_hat_0 * gs.f_hat_0;
	tprintf( "Total output mean     (simple) = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance (simple) = %g\n", gs.D_hat_t );
	// Calculate individual and interaction output variances
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		fhat2 = 0;
		for( count = 0; count < n_sub; count++ )
		{
			// tmp = gs.f_a[count] * gs.fmat_a[i][count];
			tmp = ( gs.f_a[count] - gs.f_hat_a ) * ( gs.fmat_a[i][count] - gs.f_hat_a );
			fhat2 += tmp;
			phis_half[count] = tmp;
		}
		// gs.D_hat_a[i] = ( fhat2 / n_sub ) - gs.f_hat_a * gs.f_hat_a;
		gs.D_hat_a[i] = ( fhat2 / n_sub );
		tprintf( "hat{D}_a %d (simple) %g", i + 1, gs.D_hat_a[i] );
		// gs.D_hat_a[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_a[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( " (gsl) %g", gs.D_hat_a[i] );
		// var_sorted( gs.f_a, gs.fmat_a[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_a[i] );
		tprintf( " (nr) %g\n", gs.D_hat_a[i] );
		fhat2 = 0;
		for( count = 0; count < n_sub; count++ )
		{
			// tmp = gs.f_b[count] * gs.fmat_b[i][count];
			tmp = ( gs.f_b[count] - gs.f_hat_0 ) * ( gs.fmat_b[i][count] - gs.f_hat_0 );
			fhat2 += tmp;
			phis_half[count] = tmp;
		}
		// gs.D_hat_b[i] = ( fhat2 / n_sub ) - gs.f_hat_b * gs.f_hat_b;
		gs.D_hat_b[i] = ( fhat2 / n_sub );
		tprintf( "hat{D}_b %d (simple) %g", i + 1, gs.D_hat_b[i] );
		// gs.D_hat_b[i] = gsl_stats_covariance_m( gs.f_b, 1, gs.fmat_b[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( " (gsl) %g", i, gs.D_hat_b[i] );
		// var_sorted( gs.f_b, gs.fmat_b[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_b[i] );
		tprintf( " (nr) %g\n", gs.D_hat_b[i] );
	}
	tprintf( "Total output variance (simple) = %g\n", gs.D_hat_t );
	// Print sensitivity indices
	tprintf( "\nParameter sensitivity indices:\n" );
	tprintf( "parameter SI (individual sensitvity index) ST (total sensitvity index)\n" );
	for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( "%d %g %g\n", i + 1, ( double ) gs.D_hat_a[i] / gs.D_hat_t, ( double ) 1 - ( gs.D_hat_b[i] / gs.D_hat_t ) );
	tprintf( "\n" );
	free( opt_params ); free( phis_half ); free( gs.f_a ); free( gs.f_b ); free( gs.D_hat_a ); free( gs.D_hat_b );
	free_matrix( ( void ** ) gs.var_a_lhs, n_sub );
	free_matrix( ( void ** ) gs.var_b_lhs, n_sub );
	free_matrix( ( void ** ) gs.fmat_a, op->pd->nOptParam );
	free_matrix( ( void ** ) gs.fmat_b, op->pd->nOptParam );
	return( 1 );
}

int sa_sobol( struct opt_data *op )
{
	int i, j, k, count;
	double *opt_params, *var_a_lhs, *var_b_lhs;
	char filename[255], buf[255];
	FILE *out, *out2;
	struct gsa_data gs;
	strcpy( op->label, "sobol" );
	double fhat, fhat2, gfhat, gfhat2, *phis_full, *phis_half, t1, t2;
	double var_y;
	int n_sub; //! number of samples for subsets a and b
	//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_sobol, op->pd->nOptParam );
	n_sub = op->cd->nreal / 2;	// set to half of user specified reals
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store op->cd->nreal phis
	if( ( phis_full = ( double * ) malloc( 2 * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store m_sub phis
	if( ( phis_half = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store random sample a
	if( ( var_a_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Sample a phis
	if( ( gs.f_a = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Sample b phis
	if( ( gs.f_b = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Temporary variable to store random sample b
	if( ( var_b_lhs = ( double * ) malloc( op->pd->nOptParam * n_sub * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // matrices to store lhs samples
	gs.var_a_lhs = double_matrix( n_sub, op->pd->nOptParam );
	gs.var_b_lhs = double_matrix( n_sub, op->pd->nOptParam ); // Matrices to store phis with different combinations of parameters from samples a and b
	if( ( gs.fmat_a = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( ( gs.fmat_b = double_matrix( op->pd->nOptParam, n_sub ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Vector of variances for individual component contribution
	if( ( gs.D_hat_a = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); } // Vector of variances for total component contribution
	if( ( gs.D_hat_b = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\nGlobal sensitivity analysis (Sobol) using random sampling:\n" );
	// Create samples
	if( op->cd->seed < 0 ) { op->cd->seed *= -1; tprintf( "Imported seed: %d\n", op->cd->seed ); }
	else if( op->cd->seed == 0 ) { tprintf( "New " ); op->cd->seed_init = op->cd->seed = get_seed(); }
	else tprintf( "Current seed: %d\n", op->cd->seed );
	// Create samples
	// Sample A
	tprintf( "Random sampling set 1 (variables %d; realizations %d) using ", op->pd->nOptParam, n_sub );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_a_lhs, op, 1 );
	tprintf( "done.\n" );
	// Sample B
	tprintf( "Random sampling set 2 (variables %d; realizations %d) using ", op->pd->nOptParam, n_sub );
	sampling( op->pd->nOptParam, n_sub, &op->cd->seed, var_b_lhs, op, 1 );
	tprintf( "done.\n" );
	// Copy temp lhs vectors to matrices
	for( count = 0; count < n_sub; count++ )
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			gs.var_a_lhs[count][i] = var_a_lhs[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
			gs.var_b_lhs[count][i] = var_b_lhs[i + count * op->pd->nOptParam] * op->pd->var_range[k] + op->pd->var_min[k];
		}
	free( var_a_lhs ); free( var_b_lhs );
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
	if( op->cd->mdebug > 1 )
	{
		sprintf( filename, "%s.sobol.results", op->root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.sobol_%s.results >& /dev/null", filename, op->root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
	}
	// Accumulate phis into fhat and fhat2 for total output mean and variance
	tprintf( "Computing model outputs to calculate total output mean and variance ... Sample A ... \n" );
	gfhat = gfhat2 = fhat = fhat2 = 0;
	// Compute sample a phis
	for( count = 0; count < n_sub; count++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_a_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		gfhat = fhat += op->phi;
		gfhat2 = fhat2 += op->phi * op->phi;
		// Save sample a phis
		gs.f_a[count] = phis_full[count] = op->phi;
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
	var_y = fhat2 / n_sub - ( fhat * fhat / ( ( double ) n_sub * n_sub ) );
	gs.f_hat_a = fhat / n_sub;
	gs.D_hat_t = fhat2 / n_sub - gs.f_hat_a * gs.f_hat_a;
	tprintf( "Sample A output mean     (simple) = %g\n", gs.f_hat_a );
	tprintf( "Sample A output variance (simple) = %g\n", gs.D_hat_t );
	// Compute sample b phis
	tprintf( "Computing model outputs to calculate total output mean and variance ... Sample B ... \n" );
	fhat = fhat2 = 0;
	for( count = 0; count < n_sub; count++ )
	{
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			opt_params[i] = op->pd->var[k] = gs.var_b_lhs[count][i];
		}
		Transform( opt_params, op, opt_params );
		func_global( opt_params, op, op->od->res );
		// Sum phi and phi^2
		gfhat = fhat += op->phi;
		gfhat2 = fhat2 += op->phi * op->phi;
		// Save sample b phis
		gs.f_b[count] = phis_full[n_sub + count] = op->phi;
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
	gs.f_hat_b = fhat / n_sub;
	gs.D_hat_t = fhat2 / n_sub - gs.f_hat_b * gs.f_hat_b;
	tprintf( "Sample B output mean     (simple) = %g\n", gs.f_hat_b );
	tprintf( "Sample B output variance (simple) = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gfhat / ( 2 * n_sub );
	gs.D_hat_t = gfhat2 / ( 2 * n_sub ) - gs.f_hat_0 * gs.f_hat_0;
	tprintf( "Total output mean     (simple) = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance (simple) = %g\n", gs.D_hat_t );
	gs.f_hat_0 = gsl_stats_mean( phis_full, 1, n_sub * 2 );
	gs.D_hat_t = gsl_stats_variance( phis_full, 1, n_sub * 2 );
	tprintf( "Total output mean        (gsl) = %g\n", gs.f_hat_0 );
	tprintf( "Total output variance    (gsl) = %g\n", gs.D_hat_t );
	ave_sorted( phis_full, n_sub * 2, &gs.f_hat_0, &gs.ep );
	var_sorted( phis_full, phis_full, n_sub * 2, gs.f_hat_0, gs.ep, &gs.D_hat_t );
	tprintf( "Total output mean         (nr) = %g abs 1st moment = %g\n", gs.f_hat_0, gs.ep );
	tprintf( "Total output variance     (nr) = %g\n", gs.D_hat_t );
	free( phis_full );
	// Collect matrix of phis for fmat_a
	tprintf( "Computing phis for calculation of individual output variances:\n" );
	// gfhat = gfhat2 = 0;
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
		k = op->pd->var_index[i];
		tprintf( "Processing parameter %d out of %d ... %s ...\n", i + 1, op->pd->nOptParam, op->pd->var_name[k] );
		for( count = 0; count < n_sub; count++ )
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
		fhat = fhat2 = 0;
		t1 = 0;
		t2 = 0;
		for( count = 0; count < n_sub; count++ )
		{
			// t1 = gs.f_a[count] * gs.fmat_a[i][count]
			/*t1 = ( gs.fmat_b[i][count] - gs.f_a[count] );
			t2 = gs.f_b[count] * t1;
			fhat += t2;
			fhat2 += t1 * t1;
			*/
			// phis_half[count] = t2;
			// These next two sums are based on Eqs. 18 and 19 in "Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index" by Saltelli, et al.
			t1 += pow( gs.f_b[count] - gs.fmat_b[i][count], 2 ); //t1 = sum in eq. 18
			t2 += pow( gs.f_a[count] - gs.fmat_b[i][count], 2 ); //t2 = sum in eq. 19
			// tprintf("f_b[%d]: %g\tf_a[%d]: %g\tfmat_b[%d][%d]: %g\n", count, gs.f_b[count], count, gs.f_a[count], i, count, gs.fmat_b[i][count]);
			//t1 += gs.f_b[count] * ( gs.fmat_a[i][count] - gs.f_a[count] );
			//t2 += pow( gs.f_a[i] - gs.fmat_a[i][count], 2 );
		}
		// gs.D_hat_a[i] = ( fhat2 / n_sub ) - gs.f_hat_a * gs.f_hat_a;
		//gs.D_hat_a[i] = fhat / n_sub;
		//gs.D_hat_b[i] = fhat2 / ( 2 * n_sub );
		gs.D_hat_a[i] = var_y - t1 / ( 2 * n_sub ); //var_y * eq. 18
		gs.D_hat_b[i] = t2 / ( 2 * n_sub ); // var_y * eq. 19
		//gs.D_hat_a[i] = t1 / n_sub;
		//gs.D_hat_b[i] = t2 / ( 2 * n_sub );
		tprintf( "hat{D}_a %d (simple) %g", i + 1, gs.D_hat_a[i] );
		// gs.D_hat_a[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_a[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
		tprintf( " (gsl) %g", gs.D_hat_a[i] );
		// var_sorted( gs.f_a, gs.fmat_a[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_a[i] );
		tprintf( " (nr) %g\n", gs.D_hat_a[i] );
	}
	// Print sensitivity indices
	tprintf( "\nParameter sensitivity indices:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		tprintf( "%-39s: %g (total) %g\n", op->pd->var_name[k], ( double ) gs.D_hat_a[i] / var_y, gs.D_hat_b[i] / var_y );
	}
	tprintf( "\n" );
	free( opt_params ); free( phis_half ); free( gs.f_a ); free( gs.f_b ); free( gs.D_hat_a ); free( gs.D_hat_b );
	free_matrix( ( void ** ) gs.var_a_lhs, n_sub );
	free_matrix( ( void ** ) gs.var_b_lhs, n_sub );
	free_matrix( ( void ** ) gs.fmat_a, op->pd->nOptParam );
	free_matrix( ( void ** ) gs.fmat_b, op->pd->nOptParam );
	return( 1 );
}

int sa_saltelli( struct opt_data *op )
{
	return 0;
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
