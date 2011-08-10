// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dylan Harp, dharp@lanl.gov
//
// http://www.ees.lanl.gov/staff/monty/codes/mads
//
// LA-CC-10-055; LA-CC-11-035
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../mads.h"

int Ftest( char *filename );

char *dir_hosts( void *data, char *timedate_stamp )
{
	struct opt_data *p = ( struct opt_data * )data;
	char *dir;
	dir = ( char * ) malloc( ( strlen( p->cd->mydir ) + strlen( p->root ) + strlen( timedate_stamp ) + 255 ) * sizeof( char ) );
	sprintf( dir, "%s_%s_%s", p->cd->mydir, p->root, timedate_stamp );
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
