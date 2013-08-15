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
int glue( struct opt_data *op );

/* Functions elsewhere */
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );

int glue( struct opt_data *op )
{
	FILE *in, *out;
	double *phi, **preds, phi_temp, *percentile, *pred_temp, *sum;
	char buf[200], filename[80], pred_id[100][30];
	int num_lines = 0, j;
	gsl_matrix *glue_mat; // matrix for sorting predictions
	gsl_permutation *p1;
	// Open postpua output file
	if( op->cd->infile[0] == 0 ) { tprintf( "\nInfile (results file from postpua run) must be specified for glue run\n" ); return( 0 );}
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
	tprintf( "\nNumber of solutions with phi <= %g: %d\n", op->cd->phi_cutoff, num_lines );
	tprintf( "\nPerforming GLUE analysis for %s...", op->cd->infile );
	// Read in data
	rewind( in );
	fscanf( in, "%*s %*s %[^\n]s", buf ); // Skip first part of header ("Number OF")
	int i = 0;
	// Read in names of predictions (e.g. combination of well names and times)
	while( sscanf( buf, " %s %[^\n]s", pred_id[i], buf ) > 1 ) { i++; }
	int num_preds = i + 1;
	// Allocate memory for phis and predictions
	if( ( phi = ( double * ) malloc( num_lines * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	//phi = ( double * ) malloc( num_lines * sizeof( double ) );
	preds = double_matrix( num_lines, num_preds );
	// Collect acceptable solutions
	glue_mat = gsl_matrix_alloc( num_lines, num_preds + 1 );
	int phi_index = num_preds; // phi_index indicates column of phis in glue_mat
	i = 0;
	// tprintf( "\n\nAcceptable lines from %s:\n", op->cd->infile );
	while( fgets( buf, sizeof buf, in ) != NULL )
	{
		sscanf( buf, "%*d %lf %[^\n]s", &phi_temp, buf );
		if( phi_temp <= op->cd->phi_cutoff )
		{
			// tprintf( "%lf %s\n", phi_temp, buf );
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
	// tprintf( "\nglue_mat:\n" );
	// gsl_matrix_fprintf( stdout, glue_mat, "%g" );
	// Calculate weighted percentile of each phi; note: low phis imply high percentile
	p1 = gsl_permutation_alloc( num_lines );
	percentile = ( double * ) malloc( num_lines * sizeof( double ) );
	pred_temp = ( double * ) malloc( num_lines * sizeof( double ) );
	sum = ( double * ) malloc( num_lines * sizeof( double ) );
	double p05, p95; // 5th and 95th percentiles
	tprintf( "\n\nprediction p05 p95\n" );
	int count;
	for( i = 0; i < num_preds; i++ )
	{
		gsl_vector_view column = gsl_matrix_column( glue_mat, i );
		gsl_sort_vector_index( p1, &column.vector );
		sum[0] = gsl_matrix_get( glue_mat, gsl_permutation_get( p1, 0 ), num_preds );
		pred_temp[0] = gsl_matrix_get( glue_mat, gsl_permutation_get( p1, 0 ), i );
		// tprintf( "\nSample sum prediction:\n" );
		// tprintf( "0 %g %g\n", sum[0], pred_temp[0] );
		// Collect summation of weights and ordered predictions
		for( j = 1; j < num_lines; j++ )
		{
			sum[j] = sum[j - 1] + gsl_matrix_get( glue_mat, gsl_permutation_get( p1, j ), num_preds );
			pred_temp[j] = gsl_matrix_get( glue_mat, gsl_permutation_get( p1, j ), i );
			// tprintf( "%d %g %g\n", j, sum[j], pred_temp[j] );
		}
		// tprintf( "\nno prediction percentile:\n" );
		for( j = 0; j < num_lines; j++ ) { percentile[j] = ( 1.0 / sum[num_lines - 1] ) * ( sum[j] - pred_temp[j] / 2.0 ); /*printf( "%d %g %g\n", j+1, pred_temp[j], percentile[j]);*/ }
		if( percentile[0] > 0.05 ) p05 = pred_temp[0];
		else if( percentile[num_lines - 1] < 0.05 ) p05 = pred_temp[num_lines - 1];
		else
		{
			count = 0;
			for( j = 1; j < num_lines; j++ ) { if( percentile[j] < 0.05 ) count++; else {break;} }
			p05 = pred_temp[j - 1] + ( ( 0.05 - percentile[j - 1] ) / ( percentile[j] - percentile[j - 1] ) ) * ( pred_temp[j] - pred_temp[j - 1] );
		}
		// tprintf( "\n%d\n", j );
		// tprintf( "\n%g %g %g %g %g\n", pred_temp[j], pred_temp[j-1], percentile[j], percentile[j-1], p05 );
		if( percentile[0] > 0.95 ) p95 = pred_temp[0];
		else if( percentile[num_lines - 1] < 0.95 ) p95 = pred_temp[num_lines - 1];
		else
		{
			count = 0;
			for( j = 1; j < num_lines; j++ ) { if( percentile[j] < 0.95 ) count++; else {break;}  }
			p95 = pred_temp[j - 1] + ( ( 0.95 - percentile[j - 1] ) / ( percentile[j] - percentile[j - 1] ) ) * ( pred_temp[j] - pred_temp[j - 1] );
		}
		// tprintf( "\n%d\n", j );
		// tprintf( "\n%g %g %g %g %g\n", pred_temp[j], pred_temp[j-1], percentile[j], percentile[j-1], p95 );
		tprintf( "%d %g %g\n", i + 1, p05, p95 );
		// tprintf( "%g ", gsl_interp_eval( pred_interp, pred_temp, percentile, 0.95, accelerator ) );
		// tprintf( " %g\n", gsl_interp_eval( pred_interp, pred_temp, percentile, 0.05, accelerator ) );
	}
	tprintf( "\n" );
	fclose( out );
	tprintf( "Done.\n" );
	tprintf( "Results written to %s\n\n", filename );
	gsl_matrix_free( glue_mat );
	return( 1 );
}
