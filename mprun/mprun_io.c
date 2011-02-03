#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include "../mads.h"

/* Functions elsewhere */
char **char_matrix( int maxCols, int maxRows );
int Ftest( char *filename );
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );
void *malloc_check( const char *what, size_t n );

int create_mprun_dirs( int nDir, char **dirs )
{
	char *cwd, *p, *mdir, buf[1000];
	time_t t;
	int i, r;
	cwd = getenv( "PWD" );
	p = strrchr( cwd, '/' );
//	printf( "%s %s\n", cwd, p );
	mdir = &p[1];
	p[0] = 0;
//	printf( "%s %s\n", cwd, mdir );
	t = time( NULL );
	for( i = 0; i < nDir; i++ )
	{
		sprintf( dirs[i], "%s_%ld_p%d", mdir, t, i );
		sprintf( buf, "mkdir ../%s", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); exit( -1 ); }
		sprintf( buf, "cd ../%s; ln -sf ../%s/* .; cd ../%s", dirs[i], mdir, mdir );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); exit( -1 ); }
	}
	return( 0 );
}

int delete_mprun_dirs( int nDir, char **dirs )
{
	char buf[1000];
	int i, r;
	for( i = 0; i < nDir; i++ )
	{
		sprintf( buf, "rm -fR ../%s", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	}
	return( 0 );
}
