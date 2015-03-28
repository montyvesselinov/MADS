// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
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
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include "../mads.h"

char *dir_hosts( void *data, char *timedate_stamp )
{
	struct opt_data *p = ( struct opt_data * ) data;
	char *dir;
	dir = ( char * ) malloc( ( strlen( p->cd->mydir ) + strlen( p->root ) + strlen( timedate_stamp ) + 255 ) * sizeof( char ) );
	sprintf( dir, "%s_%s_%06d_%s", p->cd->mydir, p->root, getpid(), timedate_stamp );
	return ( dir );
}

int create_mprun_dir( char *dir )
{
	char buf[1000];
	char buf2[1000];
	char cwd[1000];
	sprintf( buf, "../%s", dir );
	mkdir( buf, S_IRWXU ); // mprun directory
	DIR *dp;
	struct dirent *ep;
	getcwd( cwd, 1000 );
	// tprintf( "Working directory %s\n", cwd );
	char *root_dot = strrchr( cwd, '/' );
	char *mydir = &root_dot[1];
	// tprintf( "Working directory %s\n", mydir );
	dp = opendir( "./" ); // working directory
	if( dp != NULL )
	{
		while( ( ep = readdir( dp ) ) )
		{
			sprintf( buf, "../%s/%s", dir, ep->d_name );
			sprintf( buf2, "../%s/%s", mydir, ep->d_name );
			// tprintf( "Sym link %s\n", buf2 );
			symlink( buf2, buf );
		}
		( void ) closedir( dp );
	}
	else
		tprintf( "Couldn't open the working directory\n" );
	return ( 0 );
}

int delete_mprun_dir( char *dir )
{
	char buf[1000];
	sprintf( buf, "../%s", dir );
	// tprintf( "Delete directory %s\n", dir );
	DIR *dp;
	struct dirent *ep;
	dp = opendir( buf );
	if( dp != NULL )
	{
		while( ( ep = readdir( dp ) ) )
		{
			sprintf( buf, "../%s/%s", dir, ep->d_name );
			unlink( buf );
		}
		( void ) closedir( dp );
		sprintf( buf, "../%s", dir );
		rmdir( buf );
	}
	else
		tprintf( "Couldn't delete the directory %s\n", dir );
	return ( 0 );
}

int create_mprun_dirs( int nDir, char **dirs )
{
	int i;
	for( i = 0; i < nDir; i++ )
		create_mprun_dir( dirs[i] );
	return ( 0 );
}

int delete_mprun_dirs( int nDir, char **dirs )
{
	int i;
	for( i = 0; i < nDir; i++ )
		delete_mprun_dir( dirs[i] );
	return ( 0 );
}
