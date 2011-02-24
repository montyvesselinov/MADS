#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../mads.h"

int Ftest( char *filename );

char *dir_hosts( void *data )
{
	struct opt_data *p = ( struct opt_data * )data;
	char *dir;
	dir = ( char * ) malloc(( strlen( p->cd->mydir ) + strlen( p->root ) + 3 ) * sizeof( char ) );
	sprintf( dir, "%s_%s", p->cd->mydir, p->root );
	return( dir );
}

int create_mprun_dir( char *dir )
{
	char buf[1000];
	int r;
	sprintf( buf, "mkdir ../%s >& /dev/null", dir );
	r = system( buf );
	if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	sprintf( buf, "ln -s $PWD/* ../%s >& /dev/null", dir );
	r = system( buf );
	if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	return( 0 );
}

int delete_mprun_dir( char *dirs )
{
	char buf[1000];
	int r;
	sprintf( buf, "rm -fR ../%s >& /dev/null", dirs );
	r = system( buf );
	if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	return( 0 );
}

int create_mprun_dirs( int nDir, char **dirs )
{
	char buf[1000];
	int i, r;
	for( i = 0; i < nDir; i++ )
	{
		sprintf( buf, "mkdir ../%s  >& /dev/null", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
		sprintf( buf, "ln -s $cwd/* ../%s", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	}
	return( 0 );
}

int delete_mprun_dirs( int nDir, char **dirs )
{
	char buf[1000];
	int i, r;
	for( i = 0; i < nDir; i++ )
	{
		sprintf( buf, "rm -fR ../%s >& /dev/null", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	}
	return( 0 );
}
