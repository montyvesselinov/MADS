#include <stdlib.h>
#include <string.h>
#include <stdio.h>

char **char_matrix( int maxCols, int maxRows );
float **float_matrix( int maxCols, int maxRows );
double **double_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
void zero_double_matrix( double **matrix, int maxCols, int maxRows );
void *malloc_check( const char *what, size_t n );
int count_lines( char *filename );
int count_cols( char *filename, int row );

char **char_matrix( int maxCols, int maxRows )
{
	char **matrix;
	int i;
	if(( matrix = ( char ** ) malloc( maxCols * sizeof( char * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if(( matrix[i] = ( char * ) malloc( maxRows * sizeof( char ) ) ) == NULL )
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
	if(( matrix = ( float ** ) malloc( maxCols * sizeof( float * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if(( matrix[i] = ( float * ) malloc( maxRows * sizeof( float ) ) ) == NULL )
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
	if(( matrix = ( double ** ) malloc( maxCols * sizeof( double * ) ) ) == NULL )
		return( NULL );
	for( i = 0; i < maxCols; i++ )
		if(( matrix[i] = ( double * ) malloc( maxRows * sizeof( double ) ) ) == NULL )
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

int count_lines( char *filename )
{
	int nol=0;
	FILE *fl;
	char buf[1000];

	fl = fopen( filename, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", filename ); exit( 0 ); }
	while( (fgets( buf, sizeof buf, fl )) != NULL) nol++;
	fclose( fl );
	return nol;
}

//! Count number of columns at row
int count_cols( char *filename, int row )
{
	int ncol=0, i, n=0;
	FILE *fl;
	char buf[1000], entry[16], *ln;

	fl = fopen( filename, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", filename ); exit( 0 ); }
	for( i = 1; i < row; i++ ) ln = fgets( buf, sizeof buf, fl );
	while( sscanf( ln, "%10s%n", entry, &n ) == 1 )
	{
		ncol++;
		ln += n;
	}
	fclose( fl );
	return ncol;
}
