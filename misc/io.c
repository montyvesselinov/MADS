#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

/* Functions here */
int Ftest( char *filename );
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );
FILE *Fappend( char *filename );
char *Fdatetime( char *filename, int debug );
time_t Fdatetime_t( char *filename, int debug );

int Ftest( char *filename )
{
	return( access( filename, R_OK ) );
}

FILE *Fread( char *filename )
{
	FILE *in;
	if( ( in = fopen( filename, "rb" ) ) == NULL )
	{
		printf( "ERROR: The file %s could not opened to read!\n", filename );
		exit( 0 );
	}
	return( in );
}

FILE *Fwrite( char *filename )
{
	FILE *out;
	if( ( out = fopen( filename, "w" ) ) == NULL )
	{
		printf( "ERROR: The file %s could not opened to write!\n", filename );
		exit( 0 );
	}
	return( out );
}

FILE *Fappend( char *filename )
{
	FILE *out;
	if( ( out = fopen( filename, "a" ) ) == NULL )
	{
		printf( "ERROR: The file %s could not opened to write!\n", filename );
		exit( 0 );
	}
	return( out );
}

char *Fdatetime( char *filename, int debug )
{
	struct tm *ptr_ts;
	struct stat b;
	char *datetime;
	datetime = malloc( 16 * sizeof( char ) );
	if( Ftest( filename ) != 0 ) { datetime[0] = 0; if( debug ) printf( "File %s: does not exist\n", filename ); }
	else if( !stat( filename, &b ) )
	{
		ptr_ts = localtime( &b.st_mtime );
		sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
		if( debug ) printf( "File %s: last modified at %s\n", filename, datetime );
	}
	else { datetime[0] = 0; if( debug ) printf( "File %s: cannot display the time\n", filename ); }
	return( datetime );
}

time_t Fdatetime_t( char *filename, int debug )
{
	struct tm *ptr_ts;
	struct stat b;
	char *datetime;
	datetime = malloc( 16 * sizeof( char ) );
	if( Ftest( filename ) != 0 ) { if( debug ) printf( "File %s: does not exist\n", filename ); return ( 0 ); }
	else if( !stat( filename, &b ) )
	{
		ptr_ts = localtime( &b.st_mtime );
		sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
		if( debug ) printf( "File %s: last modified at %s\n", filename, datetime );
		return( b.st_mtime );
	}
	else { if( debug ) printf( "File %s: cannot display the time\n", filename ); return( 0 ); }
}
