#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

/* Functions here */
int Ftest( char *filename );
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );

int Ftest( char *filename )
{
	return( access( filename, R_OK ) );
}

FILE *Fread( char *filename )
{
	FILE *in;
	if(( in = fopen( filename, "rb" ) ) == NULL )
	{
		printf( "ERROR: The file %s could not opened to read!\n", filename );
		exit( 0 );
	}
	return( in );
}

FILE *Fwrite( char *filename )
{
	FILE *out;
	if(( out = fopen( filename, "w" ) ) == NULL )
	{
		printf( "ERROR: The file %s could not opened to write!\n", filename );
		exit( 0 );
	}
	return( out );
}
