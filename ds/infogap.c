// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
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
int infogap_obs( struct opt_data *op );
int infogap( struct opt_data *op );

/* Functions elsewhere */
int igrnd( struct opt_data *op );
int count_lines( char *filename );
int count_cols( char *filename, int row );

int infogap_obs( struct opt_data *op )
{
	int i, j, ig_index, status, success = 0, count, neval_total, njac_total, no_memory = 0;
	double *opt_params, phi_min = HUGE_VAL;
	int ( *optimize_func )( struct opt_data * op );
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL ) no_memory = 1;
	if( no_memory ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	tprintf( "\n\nInfo-gap analysis: observation step %g observation domain %g\nInfo-gap search: ", op->cd->obsstep, op->cd->obsdomain );
	if( op->cd->obsstep > DBL_EPSILON ) tprintf( "maximum\n" ); else tprintf( "minimum\n" );
	tprintf( "Number of predictions %d\n", op->preds->nTObs );
	for( i = 0; i < op->preds->nTObs; i++ )
	{
		// op->preds->obs_best are updated in mads_func.c
		if( op->cd->obsstep > DBL_EPSILON ) op->preds->obs_best[i] = -HUGE_VAL; // max search
		else op->preds->obs_best[i] = HUGE_VAL; // min search
		j = op->preds->obs_index[i];
		op->od->obs_weight[j] *= -1;
	}
	ig_index = op->preds->obs_index[0]; // first prediction is applied only
	tprintf( "Info-gap observation:\n" );
	tprintf( "%-20s: info-gap target %12g weight %12g range %12g - %12g\n", op->od->obs_id[ig_index], op->od->obs_target[ig_index], op->od->obs_weight[ig_index], op->od->obs_min[ig_index], op->od->obs_max[ig_index] );
	if( op->cd->obsstep > DBL_EPSILON ) { op->od->obs_target[ig_index] = op->od->obs_min[ig_index]; op->od->obs_min[ig_index] -= op->cd->obsstep / 2; } // obsstep is positive
	else { op->od->obs_target[ig_index] = op->od->obs_max[ig_index]; op->od->obs_max[ig_index] -= op->cd->obsstep / 2; } // obsstep is negative
	if( strncasecmp( op->cd->opt_method, "lm", 2 ) == 0 ) optimize_func = optimize_lm; // Define optimization method: LM
	else optimize_func = optimize_pso; // Define optimization method: PSO
	neval_total = njac_total = count = 0;
	while( 1 )
	{
		tprintf( "\nInfo-gap analysis #%d\n", ++count );
		tprintf( "%-20s: info-gap target %12g weight %12g range %12g - %12g\n", op->od->obs_id[ig_index], op->od->obs_target[ig_index], op->od->obs_weight[ig_index], op->od->obs_min[ig_index], op->od->obs_max[ig_index] );
		op->cd->neval = op->cd->njac = 0;
		if( op->cd->calib_type == IGRND ) status = igrnd( op );
		else status = optimize_func( op );
		neval_total += op->cd->neval;
		njac_total += op->cd->njac;
		if( !status ) break;
		if( op->success )
		{
			for( i = 0; i < op->pd->nOptParam; i++ ) opt_params[i] = op->pd->var[op->pd->var_index[i]];
			for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_best[i] = op->od->obs_current[i];
			phi_min = op->phi;
			success = 1;
		}
		tprintf( "Intermediate info-gap results for model predictions:\n" );
		tprintf( "%-20s: info-gap target %12g weight %12g range %12g - %12g\n", op->od->obs_id[ig_index], op->od->obs_target[ig_index], op->od->obs_weight[ig_index], op->od->obs_min[ig_index], op->od->obs_max[ig_index] );
		for( i = 0; i < op->preds->nTObs; i++ )
		{
			j = op->preds->obs_index[i];
			if( op->cd->obsstep > DBL_EPSILON ) tprintf( "%-20s: Current info-gap max %12g Observation step %g Observation domain %g Success %d\n", op->od->obs_id[j], op->preds->obs_best[i], op->cd->obsstep, op->cd->obsdomain, op->success );
			else                           tprintf( "%-20s: Current info-gap min %12g Observation step %g Observation domain %g Success %d\n", op->od->obs_id[j], op->preds->obs_best[i], op->cd->obsstep, op->cd->obsdomain, op->success );
		}
		if( !op->success ) break;
		if( op->cd->debug ) print_results( op, 1 );
		save_results( 1, "infogap", op, op->gd );
		op->od->obs_target[ig_index] += op->cd->obsstep;
		if( op->cd->obsstep > DBL_EPSILON ) // max search
		{
			if( op->od->obs_target[ig_index] > op->od->obs_max[ig_index] ) break;
			if( fabs( op->preds->obs_best[0] - op->od->obs_max[ig_index] ) < DBL_EPSILON ) break;
			op->od->obs_min[ig_index] += op->cd->obsstep;
			j = ( int )( ( double )( op->preds->obs_best[0] - op->od->obs_min[ig_index] + op->cd->obsstep / 2 ) / op->cd->obsstep + 1 );
			op->od->obs_target[ig_index] += op->cd->obsstep * j;
			op->od->obs_min[ig_index] += op->cd->obsstep * j;
			if( op->od->obs_target[ig_index] > op->od->obs_max[ig_index] ) op->od->obs_target[ig_index] = op->od->obs_max[ig_index];
			if( op->od->obs_min[ig_index] > op->od->obs_max[ig_index] ) op->od->obs_min[ig_index] = op->od->obs_max[ig_index];
		}
		else // min search
		{
			if( op->od->obs_target[ig_index] < op->od->obs_min[ig_index] ) break;
			if( fabs( op->preds->obs_best[0] - op->od->obs_min[ig_index] ) < DBL_EPSILON ) break;
			op->od->obs_max[ig_index] += op->cd->obsstep; // obsstep is negative
			j = ( int )( ( double )( op->od->obs_max[ig_index] - op->preds->obs_best[0] - op->cd->obsstep / 2 ) / -op->cd->obsstep + 1 ); // obsstep is negative
			op->od->obs_target[ig_index] += op->cd->obsstep * j;
			op->od->obs_max[ig_index] += op->cd->obsstep * j;
			if( op->od->obs_target[ig_index] < op->od->obs_min[ig_index] ) op->od->obs_target[ig_index] = op->od->obs_min[ig_index];
			if( op->od->obs_max[ig_index] < op->od->obs_min[ig_index] ) op->od->obs_max[ig_index] = op->od->obs_min[ig_index];
		}
	}
	op->cd->neval = neval_total; // provide the correct number of total evaluations
	op->cd->njac = njac_total; // provide the correct number of total evaluations
	if( success )
	{
		for( i = 0; i < op->pd->nOptParam; i++ ) op->cd->var[i] = op->pd->var[op->pd->var_index[i]] = op->pd->var_current[i] = op->pd->var_best[i] = opt_params[i];
		for( i = 0; i < op->od->nTObs; i++ ) op->od->obs_current[i] = op->od->obs_best[i];
		op->phi = phi_min;
		op->success = success;
	}
	else tprintf( "\nWARNING: Info-gap analysis failed to find acceptable solutions!\n" );
	tprintf( "\nInfo-gap results for model predictions:\n" );
	for( i = 0; i < op->preds->nTObs; i++ )
	{
		j = op->preds->obs_index[i];
		if( op->cd->obsstep > DBL_EPSILON ) tprintf( "%-20s: Info-gap max %12g Observation step %g Observation domain %g\n", op->od->obs_id[j], op->preds->obs_best[i], op->cd->obsstep, op->cd->obsdomain ); // max search
		else                           tprintf( "%-20s: Info-gap min %12g Observation step %g Observation domain %g\n", op->od->obs_id[j], op->preds->obs_best[i], op->cd->obsstep, op->cd->obsdomain ); // min search
		op->od->obs_target[j] = op->preds->obs_target[i];
		op->od->obs_min[j] = op->preds->obs_min[i];
		op->od->obs_max[j] = op->preds->obs_max[i];
		op->od->obs_weight[j] *= -1;
	}
	tprintf( "\n" );
	free( opt_params );
	print_results( op, 1 );
	save_results( 1, "", op, op->gd );
	return( 1 );
}

int infogap( struct opt_data *op )
{
	FILE *fl, *outfl;
	double *opt_params, of, maxof;
	char buf[255], filename[255];
	int i, j, k, n, npar, nrow, ncol, *nPreds, col;
	gsl_matrix *ig_mat; //! info gap matrix for sorting
	gsl_permutation *p;
	nPreds = &op->preds->nTObs; // Set pointer to nObs for convenience
	if( op->cd->infile[0] == 0 ) { tprintf( "\nInfile must be specified for infogap run\n" ); return( 0 );}
	nrow = count_lines( op->cd->infile ); nrow--; // Determine number of parameter sets in file
	npar = count_cols( op->cd->infile, 2 ); npar = npar - 2; // Determine number of parameter sets in file
	if( npar != op->pd->nOptParam ) { tprintf( "Number of optimization parameters in %s does not match input file\n", op->cd->infile ); return( 0 ); } // Make sure MADS input file and PSSA file agree
	tprintf( "\n%s contains %d parameters and %d parameter sets\n", op->cd->infile, npar, nrow );
	ncol = npar + *nPreds + 1; // Number of columns for ig_mat = #pars + #preds + #ofs
	ig_mat = gsl_matrix_alloc( nrow, ncol );
	p = gsl_permutation_alloc( nrow );
	fl = fopen( op->cd->infile, "r" );
	if( fl == NULL ) { tprintf( "\nError opening %s\n", op->cd->infile ); return( 0 ); }
	tprintf( "Computing predictions for %s...", op->cd->infile );
	if( ( opt_params = ( double * ) malloc( npar * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
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
		if( outfl == NULL ) { tprintf( "\nError opening %s\n", filename ); return( 0 ); }
		fprintf( outfl, " %-12s", op->preds->obs_id[k] );
		fprintf( outfl, " OFmax OF" );
		for( i = 0; i < npar; i++ )
			fprintf( outfl, " (%-12s)", op->pd->var_name[i] );
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
		tprintf( "Done\n" );
		tprintf( "Results written to %s\n\n", filename );
	}
	gsl_matrix_free( ig_mat );
	return( 1 );
}
