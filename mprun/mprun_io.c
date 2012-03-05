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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../mads.h"

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
