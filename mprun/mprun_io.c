#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int create_mprun_dir( char *dir )
{
	char buf[1000];
	int r;
	sprintf( buf, "mkdir ../%s", dir );
	r = system( buf );
	if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	sprintf( buf, "ln -sf $PWD/* ../%s", dir );
	r = system( buf );
	if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	return( 0 );
}

int delete_mprun_dir( char *dirs )
{
	char buf[1000];
	int r;
	sprintf( buf, "rm -fR ../%s", dirs );
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
		sprintf( buf, "mkdir ../%s", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
		sprintf( buf, "ln -sf $cwd/* ../%s", dirs[i] );
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
		sprintf( buf, "rm -fR ../%s", dirs[i] );
		r = system( buf );
		if( r == -1 || r == 127 ) { printf( "ERROR: System call failed!\n" ); return( -1 ); }
	}
	return( 0 );
}
