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
#include <stddef.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "../mads.h"

/* Functions here */
int Ftest( char *filename );
int Ftestread( char *filename );
FILE *Fread( char *filename );
void Fwrite( struct io_output_object *output_obj, struct opt_data *op );
FILE *Fappend( char *filename );
char *Fdatetime( char *filename, int debug );
time_t Fdatetime_t( char *filename, int debug );
void removeChars( char *str, char *garbage );
/* Functions elsewhere */
void tprintf( char const *fmt, ... );

int Ftest( char *filename )
{
	return( access( filename, R_OK ) );
}

int Ftestread( char *filename )
{
	FILE *in;
	if( ( in = fopen( filename, "rb" ) ) == NULL ) return( 1 );
	fclose( in );
	return( 0 );
}

FILE *Fread( char *filename )
{
	FILE *in;
	if( ( in = fopen( filename, "rb" ) ) == NULL )
	{
		tprintf( "ERROR: The file %s could not opened to read!\n", filename );
		exit( 0 );
	}
	return( in );
}

void Fwrite( struct io_output_object *output_obj, struct opt_data *op )
{
        FILE *hold_file_ptr;

	if( ( output_obj->outfile = fopen( output_obj->filename, "w" ) ) == NULL )
	{
		tprintf( "ERROR: The file %s could not opened to write!\n", output_obj->filename );
		exit( 0 );
	}

	// Test to see if pointers have been initialized
 	if(op->wd)
	  {
	    if(op->wd->x)
	      output_obj->use_netcdf = 1 ;
	    else
	      output_obj->use_netcdf = 0 ;
	  }
	else
	  output_obj->use_netcdf = 0 ;

	if( output_obj->use_netcdf )
	  {
	    create_cdf_file(output_obj, op->wd) ;
	    // I want to call write_problem but not output to the text file
	    
	    hold_file_ptr = output_obj->outfile ;

	    // This directions text output to nowhere
	    output_obj->outfile = fopen("/dev/null","w") ;

	    write_problem_text(output_obj, op) ;

	    output_obj->outfile = hold_file_ptr ;
	  }
}

void Fclose( struct io_output_object *output_obj )
{
  fclose( output_obj->outfile ) ;

  if( output_obj->use_netcdf )
    close_cdf_file(output_obj) ;
		      
}

FILE *Fappend( char *filename )
{
	FILE *out;
	if( ( out = fopen( filename, "a" ) ) == NULL )
	{
		tprintf( "ERROR: The file %s could not opened to write!\n", filename );
		exit( 0 );
	}
	return( out );
}

char *Fdatetime( char *filename, int debug )
{
	struct tm *ptr_ts;
	struct stat b;
	char *datetime;
	datetime = ( char * ) malloc( 16 * sizeof( char ) );
	if( Ftest( filename ) != 0 ) { datetime[0] = 0; if( debug ) tprintf( "File %s: does not exist\n", filename ); }
	else if( !stat( filename, &b ) )
	{
		ptr_ts = localtime( &b.st_mtime );
		sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
		if( debug )	printf( "File %s: last modified at %s\n", filename, datetime );
	}
	else { datetime[0] = 0; if( debug )	printf( "File %s: cannot display the time\n", filename ); }
	return( datetime );
}

time_t Fdatetime_t( char *filename, int debug )
{
	struct tm *ptr_ts;
	struct stat b;
	char *datetime;
	if( Ftest( filename ) != 0 ) { if( debug ) tprintf( "File %s: does not exist\n", filename ); return ( 0 ); }
	else if( !stat( filename, &b ) )
	{
		ptr_ts = localtime( &b.st_mtime );
		datetime = ( char * ) malloc( 16 * sizeof( char ) );
		sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
		if( debug )	tprintf( "File %s: last modified at %s\n", filename, datetime );
		free( datetime );
		return( b.st_mtime );
	}
	else { if( debug ) printf( "File %s: cannot display the time\n", filename ); return( 0 ); }
}

void removeChars( char *str, char *garbage )
{
	char *src, *dst, *ch;
	for( ch = garbage; *ch != '\0'; ch++ )
	{
		for( src = dst = str; *src != '\0'; src++ )
		{
			*dst = *src;
			if( *dst != *ch ) dst++;
		}
		*dst = '\0';
	}
}
