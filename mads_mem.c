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

void white_skip( char **s )
{
	while( iswhite( **s ) )( *s )++;
}

char *white_trim( char *x )
{
	char *y;
	if( !x ) return( x );
	y = x + strlen( x ) - 1;
	while( y >= x && iswhite( *y ) ) *y-- = 0;
	return x;
}
