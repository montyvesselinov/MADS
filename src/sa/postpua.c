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
int postpua( struct opt_data *op );

/* Functions elsewhere */
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );

int postpua( struct opt_data *op )
{
	FILE *in, *out;
	double *opt_params, of;
	char buf[255], filename[255];
	int i, n;
	op->od = op->preds;
	if( op->cd->infile[0] == 0 ) { tprintf( "\nInfile (results file from abagus run) must be specified for postpua run\n" ); return( 0 );}
	in = Fread( op->cd->infile );
	// Create postpua output file
	sprintf( filename, "%s.pua", op->root );
	out = Fwrite( filename );
	tprintf( "\nComputing predictions for %s...", op->cd->infile );
	if( ( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return( 0 ); }
	fgets( buf, sizeof buf, in ); // Skip header
	fprintf( out, "Number       OF           " );
	for( i = 0; i < op->od->nTObs; i++ )
		fprintf( out, " %-12s", op->od->obs_id[i] );
	fprintf( out, "\n" );
	while( fscanf( in, "%d %lf", &n, &of ) > 0 )
	{
		fprintf( out, "%-12d %-12lf ", n, of );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fscanf( in, "%lf", &opt_params[i] );
		fscanf( in, " \n" );
		func_global( opt_params, op, op->od->res );
		for( i = 0; i < op->od->nTObs; i++ )
			fprintf( out, " %-12g", op->od->obs_current[i] );
		fprintf( out, "\n" );
	}
	fclose( in );
	fclose( out );
	tprintf( "Done.\n" );
	tprintf( "Results written to %s\n\n", filename );
	return( 1 );
}
