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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define iswhite(c) ((c)== ' ' || (c)=='\t' || (c)=='\n')

char **char_matrix( int maxCols, int maxRows );
float **float_matrix( int maxCols, int maxRows );
double **double_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
void zero_double_matrix( double **matrix, int maxCols, int maxRows );
void *malloc_check( const char *what, size_t n );
char *white_trim( char *x );
void white_skip( char **s );

char **char_matrix( int maxCols, int maxRows )
{
	char **matrix;
	int i;
	if( ( matrix = ( char ** ) malloc( maxCols * sizeof( char * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if( ( matrix[i] = ( char * ) malloc( maxRows * sizeof( char ) ) ) == NULL )
		{
			for( i--; i >= 0; i-- )
				free( matrix[i] );
			free( matrix );
			return( NULL );
		}
	return( matrix );
}

float **float_matrix( int maxCols, int maxRows )
{
	float **matrix;
	int i;
	if( ( matrix = ( float ** ) malloc( maxCols * sizeof( float * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if( ( matrix[i] = ( float * ) malloc( maxRows * sizeof( float ) ) ) == NULL )
		{
			for( i--; i >= 0; i-- )
				free( matrix[i] );
			free( matrix );
			return( NULL );
		}
	return( matrix );
}

double **double_matrix( int maxCols, int maxRows )
{
	double **matrix;
	int i;
	if( ( matrix = ( double ** ) malloc( maxCols * sizeof( double * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if( ( matrix[i] = ( double * ) malloc( maxRows * sizeof( double ) ) ) == NULL )
		{
			for( i--; i >= 0; i-- )
				free( matrix[i] );
			free( matrix );
			return( NULL );
		}
	return( matrix );
}

void free_matrix( void **matrix, int maxCols )
{
	int i;
	for( i = 0; i < maxCols; i++ )
		free( matrix[i] );
	free( matrix );
}

void zero_double_matrix( double **matrix, int maxCols, int maxRows )
{
	int i;
	for( i = 0; i < maxCols; i++ )
		memset( matrix[i], 0, maxRows * sizeof( double ) );
}

void *malloc_check( const char *what, size_t n )
{
	void *p = malloc( n );
	if( p == NULL )
	{
		fprintf( stderr, "Cannot allocate %zu bytes to %s\n", n, what );
		exit( 2 );
	}
	return p;
}

void white_skip( char **s ) // remove the white space at the beginning of the string
{
	while( iswhite( **s ) )( *s )++;
}

char *white_trim( char *x ) // remove the white space at the end of the string
{
	char *y;
	if( !x ) return( x );
	y = x + strlen( x ) - 1;
	while( y >= x && iswhite( *y ) ) *y-- = 0;
	return x;
}
