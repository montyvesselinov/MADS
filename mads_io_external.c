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

#define _GNU_SOURCE
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<ctype.h>
#include<math.h>

#include "mads.h"

/* Functions here */
int load_pst( char *filename, struct opt_data *op );
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
/* Functions elsewhere */
char *white_trim( char *x );
void white_skip( char **s );
char **char_matrix( int maxCols, int maxRows );

int load_pst( char *filename, struct opt_data *op )
{
	FILE *in;
	double d;
	char code[20], buf[1000];
	int i, k, npar_groups, nobs_groups;
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od;
	struct extrn_data *ed;
	cd = op->cd;
	pd = op->pd;
	od = op->od;
	ed = op->ed;
	op->gd->min_t = op-> gd->time = 0;
	if( ( in = fopen( filename, "r" ) ) == NULL )
	{
		printf( "PEST control file %s cannot be opened to read problem data!\n", filename );
		return( -1 );
	}
	cd->opt_method = ( char * ) malloc( 50 * sizeof( char ) );
	cd->solution_id = ( char * ) malloc( 50 * sizeof( char ) );
	strcpy( cd->solution_id, "extertnal" );
	cd->solution_type = -1;
	for( i = 0; i < 4; i++ ) // skip 4 lines
		fgets( buf, 1000, in );
	sscanf( buf, "%d %d %d %*d %d", &( *pd ).nParam, &( *od ).nObs, &npar_groups, &nobs_groups );
	printf( "Parameters = %d (groups %d)\n", pd->nParam, npar_groups );
	printf( "Observations = %d (groups %d)\n", od->nObs, nobs_groups );
	fgets( buf, 1000, in );
	sscanf( buf, "%d %d", &( *ed ).ntpl, &( *ed ).nins );
	printf( "Number of template files = %d\nNumber of instruction files = %d\n", ( *ed ).ntpl, ( *ed ).nins );
	pd->var_id = char_matrix( ( *pd ).nParam, 50 );
	pd->var = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_current = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_best = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	cd->var = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc( ( *pd ).nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc( ( *pd ).nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	printf( "Parameters = %d:\n", pd->nParam );
	for( i = 0; i < 6; i++ ) // skip 6 lines
		fgets( buf, 1000, in );
	for( i = 0; i < npar_groups; i++ )
		fgets( buf, 1000, in );
	fgets( buf, 1000, in );
	pd->nFlgParam = 0;
	pd->nOptParam = 0;
	for( i = 0; i < pd->nParam; i++ )
	{
		fscanf( in, "%s %s %*s %lf %lf %lf %*s %*f %*f %*f\n", pd->var_id[i], code, &pd->var[i], &pd->var_min[i], &pd->var_max[i] );
		printf( "%-26s: init %15.12g min %12g max %12g\n", pd->var_id[i], pd->var[i], ( *pd ).var_min[i], ( *pd ).var_max[i] );
		if( strcmp( code, "fixed" ) == 0 ) pd->var_opt[i] = 0; else { pd->nOptParam++; pd->var_opt[i] = 1; }
		if( strcmp( code, "log" ) == 0 ) pd->var_log[i] = 1; else pd->var_log[i] = 0;
		if( ( *pd ).var_log[i] == 1 )
		{
			pd->var[i] = log10( pd->var[i] );
			pd->var_min[i] = log10( pd->var_min[i] );
			pd->var_max[i] = log10( pd->var_max[i] );
		}
		pd->var_range[i] = pd->var_max[i] - pd->var_min[i];
		pd->var_dx[i] = pd->var_range[i] / 10;
	}
	pd->var_index = ( int * ) malloc( ( *pd ).nOptParam * sizeof( int ) );
	printf( "Optimized parameters = %d\n", pd->nOptParam );
	for( k = i = 0; i < ( *pd ).nParam; i++ )
		if( ( *pd ).var_opt[i] == 1 )
		{
			if( ( *pd ).var_log[i] == 1 ) d = log10( pd->var[i] ); else d = pd->var[i];
			printf( "%-26s: init %15.12g min %12g max %12g\n", pd->var_id[i], d, ( *pd ).var_min[i], ( *pd ).var_max[i] );
			( *pd ).var_index[k++] = i;
		}
	fgets( buf, 1000, in ); // skip line
	for( i = 0; i < nobs_groups; i++ )
		fgets( buf, 1000, in );
	fgets( buf, 1000, in ); // skip line
	od->obs_id = char_matrix( ( *od ).nObs, 50 );
	od->obs_target = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_weight = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_min = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_max = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_current = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_best = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->res = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
	od->obs_log = ( int * ) malloc( ( *od ).nObs * sizeof( int ) );
	for( i = 0; i < od->nObs; i++ )
		fscanf( in, "%s %lf %lf %*s\n", od->obs_id[i], &od->obs_target[i], &od->obs_weight[i] );
	printf( "Calibration targets = %d\n", ( *od ).nObs );
	for( i = 0; i < od->nObs; i++ )
	{
		if( od->nObs < 50 || ( i < 20 || i > od->nObs - 20 ) )
			printf( "%-13s: value %15.12g weight %g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i] );
		if( od->nObs > 50 && i == 21 ) printf( "...\n" );
		od->obs_min[i] = 0; od->obs_max[i] = od->obs_target[i] * 2;
		od->obs_log[i] = 0;
	}
	fgets( buf, 1000, in ); // skip line
	ed->cmdline = ( char * ) malloc( 80 * sizeof( char ) );
	fgets( ed->cmdline, 80, in );
	ed->cmdline[strlen( ed->cmdline ) - 1] = 0;
	printf( "Execution command: %s\n", ed->cmdline );
	printf( "External files:\n" );
	ed->fn_ins = char_matrix( ed->nins, 80 );
	ed->fn_obs = char_matrix( ed->nins, 80 );
	ed->fn_tpl = char_matrix( ed->ntpl, 80 );
	ed->fn_out = char_matrix( ed->ntpl, 80 );
	fgets( buf, 1000, in ); // skip line
	for( i = 0; i < ed->ntpl; i++ )
		fscanf( in, "%s %s\n", ed->fn_tpl[i], ed->fn_out[i] );
	printf( "- to provide current model parameters:\n" );
	for( i = 0; i < ed->ntpl; i++ )
		printf( "%s -> %s\n", ed->fn_tpl[i], ed->fn_out[i] );
	for( i = 0; i < ed->nins; i++ )
		fscanf( in, "%s %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	printf( "- to read current model predictions:\n" );
	for( i = 0; i < ed->nins; i++ )
		printf( "%s <- %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	fclose( in );
	printf( "\n" );
	return( 1 );
}

int check_ins_obs( int nobs, char **obs_id, double *check, char *fn_in_i, int debug )
{
	FILE *infile_inst;
	char *separator = " \t\n";
	char *word_inst, *word_search, token_obs[2], token_search[2], dummy_var[6], buf_inst[1000], *pnt_inst;
	int i, c, bad_data = 0;
	if( ( infile_inst = fopen( fn_in_i, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_i );
		return( -1 );
	}
	if( debug ) printf( "\nChecking instruction file \'%s\'.\n", fn_in_i );
	fgets( buf_inst, 1000, infile_inst );
	for( c = 0, word_inst = strtok( buf_inst, separator ); word_inst; c++, word_inst = strtok( NULL, separator ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word_inst );
			if( strcasestr( word_inst, "pif" ) )
			{
				if( debug ) printf( "PEST Instruction file\n" );
				token_obs[0] = '!';
			}
			else if( strcasestr( word_inst, "instruction" ) )
			{
				if( debug ) printf( "MADS Instruction file; user-specified search/variable tokens are expected\n" );
			}
			else
			{
				if( debug ) printf( "MADS Instruction file\n" );
				rewind( infile_inst );
				token_search[0] = '@';
				token_obs[0] = '!';
				break;
			}
		}
		else if( c == 1 ) // second entry; "search" token
		{
			white_trim( word_inst );
			if( debug ) printf( "Search token %s\n", word_inst );
			token_search[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				printf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_search[0] == 0 ) token_search[0] = '@';
		}
		else if( c == 2 ) // third entry; "variable" token
		{
			white_trim( word_inst );
			if( debug ) printf( "Variable token %s\n", word_inst );
			token_obs[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				printf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_obs[0] == 0 ) token_obs[0] = '!';
			break;
		}
	}
	token_search[1] = token_obs[1] = 0;
	dummy_var[0] = token_obs[0];
	dummy_var[1] = 0;
	strcat( dummy_var, "dum" );
	dummy_var[4] = token_obs[0];
	dummy_var[5] = 0;
	if( debug )
	{
		printf( "Search separator: %s\n", token_search );
		printf( "Observation separator: %s\n", token_obs );
		printf( "Dummy observation: %s\n", dummy_var );
	}
	while( 1 )
	{
		fgets( buf_inst, 1000, infile_inst );
		if( feof( infile_inst ) ) break;
		white_trim( buf_inst );
		if( debug ) printf( "\nCurrent instruction line: %s\n", buf_inst );
		pnt_inst = &buf_inst[0];
		word_inst = strtok_r( buf_inst, separator, &pnt_inst );
		if( buf_inst[0] == 'l' ) // skip lines in the "data" file
		{
			sscanf( &word_inst[1], "%d", &c );
			if( debug ) printf( "Skip %d lines\n", c );
			word_inst = strtok_r( NULL, separator, &pnt_inst );
		}
		for( ; word_inst; word_inst = strtok_r( NULL, separator, &pnt_inst ) )
		{
			white_skip( &word_inst ); white_trim( word_inst );
			if( debug ) printf( "TEMPLETE word \'%s\' : ", word_inst ); fflush( stdout );
			if( word_inst[0] == token_search[0] ) // search for keyword
			{
				word_search = strtok( word_inst, token_search );
				white_skip( &word_search ); white_trim( word_search );
				if( debug ) printf( "Search for keyword \'%s\' in the data file ...\n", word_search );
			}
			else if( strncmp( word_inst, dummy_var, 5 ) == 0 ) // dummy variable
			{
				if( debug ) printf( "Skip dummy data!\n" );
			}
			else if( word_inst[0] == 'w' ) // white space
			{
				if( debug ) printf( "Skip white space!\n" );
			}
			else if( word_inst[0] == token_obs[0] ) // observation variable
			{
				c = 0;
				if( strlen( word_inst ) == 1 ) word_inst = strtok_r( NULL, separator, &pnt_inst );
				else word_inst = &word_inst[1];
				if( word_inst[strlen( word_inst ) - 1] == token_obs[0] ) word_inst[strlen( word_inst ) - 1] = 0;
				else strtok_r( NULL, separator, &pnt_inst );
				if( debug ) printf( "Observation keyword \'%s\' ... ", word_inst );
				white_skip( &word_inst );
				white_trim( word_inst );
				for( i = 0; i < nobs; i++ )
				{
					if( strcmp( word_inst, obs_id[i] ) == 0 )
					{
						if( check[i] < 0 ) { check[i] = 1; }
						else { check[i] += 1; }
						if( debug ) printf( "\'%s\' detected %g times\n", obs_id[i], check[i] );
						break;
					}
				}
				if( nobs == i )
				{
					printf( "\nERROR: Observation keyword \'%s\' does not match any of observation variables!\n", word_inst );
					bad_data = 1;
				}
			}
			else
			{
				printf( "ERROR: Instruction file %s does not follow the expected format!\n", fn_in_i );
				printf( "White space (w), search (%s) or observation (%s) tokens are expected!\n", token_search, token_obs );
				bad_data = 1;
				break;
			}
		}
	}
	fclose( infile_inst );
	if( bad_data ) return( -1 );
	else return( 0 );
}


int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_i, char *fn_in_d, int debug )
{
	FILE *infile_inst, *infile_data;
	char *separator = " \t\n";
	char *word_inst, *word_data, *word_search, token_search[2], token_obs[2], dummy_var[6], buf_data[1000], buf_inst[1000], *pnt_inst, *pnt_data;
	int i, c, bad_data = 0;
	double v;
	if( ( infile_inst = fopen( fn_in_i, "r" ) ) == NULL )
	{
		printf( "\nERROR: File %s cannot be opened to read template data!\n", fn_in_i );
		return( -1 );
	}
	if( ( infile_data = fopen( fn_in_d, "r" ) ) == NULL )
	{
		printf( "\nERROR: File %s cannot be opened to read the model-predicted observations!\n", fn_in_d );
		return( -1 );
	}
	if( debug ) printf( "\nReading output file \'%s\' obtained from external model execution using instruction file \'%s\'.\n", fn_in_d, fn_in_i );
	fgets( buf_inst, 1000, infile_inst );
	for( c = 0, word_inst = strtok( buf_inst, separator ); word_inst; c++, word_inst = strtok( NULL, separator ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word_inst );
			if( strcasestr( word_inst, "pif" ) )
			{
				if( debug ) printf( "PEST Instruction file\n" );
				token_obs[0] = '!';
			}
			else if( strcasestr( word_inst, "instruction" ) )
			{
				if( debug ) printf( "MADS Instruction file; user-specified search/variable tokens are expected\n" );
			}
			else
			{
				if( debug ) printf( "MADS Instruction file\n" );
				rewind( infile_inst );
				token_search[0] = '@';
				token_obs[0] = '!';
				break;
			}
		}
		else if( c == 1 ) // second entry; "search" token
		{
			white_trim( word_inst );
			token_search[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				printf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_search[0] == 0 ) token_search[0] = '@';
		}
		else if( c == 2 ) // third entry; "variable" token
		{
			white_trim( word_inst );
			token_obs[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				printf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_obs[0] == 0 ) token_obs[0] = '!';
			break;
		}
	}
	token_search[1] = token_obs[1] = 0;
	dummy_var[0] = token_obs[0];
	dummy_var[1] = 0;
	strcat( dummy_var, "dum" );
	dummy_var[4] = token_obs[0];
	dummy_var[5] = 0;
	if( debug )
	{
		printf( "Search separator: %s\n", token_search );
		printf( "Observation separator: %s\n", token_obs );
		printf( "Dummy observation: %s\n", dummy_var );
	}
	buf_data[0] = 0; word_data = pnt_data = NULL;
	while( 1 )
	{
		fgets( buf_inst, 1000, infile_inst );
		if( feof( infile_inst ) ) break;
		white_trim( buf_inst );
		if( debug ) printf( "\nCurrent instruction line: %s\n", buf_inst );
		pnt_inst = &buf_inst[0];
		word_inst = strtok_r( buf_inst, separator, &pnt_inst );
		if( buf_inst[0] == 'l' ) // skip lines in the "data" file
		{
			sscanf( &word_inst[1], "%d", &c );
			word_inst = strtok_r( NULL, separator, &pnt_inst );
			if( debug ) printf( "Skip %d lines\n", c );
			for( i = 0; i < c; i++ )
				fgets( buf_data, 1000, infile_data );
			if( feof( infile_data ) ) { printf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
			white_trim( buf_data );
			pnt_data = &buf_data[0];
			word_data = NULL;
		}
		if( pnt_data == NULL ) // if there was no "l" (skip line) command, read the next "data" line
		{
			fgets( buf_data, 1000, infile_data );
			white_trim( buf_data );
			pnt_data = &buf_data[0];
			word_data = NULL;
		}
		if( debug ) printf( "Current location in model output file: \'%s\' Remaining line: \'%s\'", word_data, pnt_data );
		if( debug ) { if( pnt_data == NULL ) printf( "\n" ); else { if( pnt_data[strlen( pnt_data ) - 2] != '\n' ) printf( "\n" ); } }
		c = 0;
		for( ; word_inst; word_inst = strtok_r( NULL, separator, &pnt_inst ) )
		{
			white_skip( &word_inst ); white_trim( word_inst );
			if( debug ) printf( "TEMPLETE word \'%s\' : ", word_inst ); fflush( stdout );
			if( word_inst[0] == token_search[0] ) // search for keyword
			{
				word_search = strtok( word_inst, token_search );
				white_skip( &word_search ); white_trim( word_search );
				if( debug ) printf( "Search for keyword \'%s\' in the data file ...\n", word_search );
				bad_data = 1;
				while( !feof( infile_data ) )
				{
					if( ( pnt_data = strstr( pnt_data, word_search ) ) != NULL )
					{
						if( debug ) printf( "Matching line found in the data file: \'%s\' Location \'%s\'\n", buf_data, pnt_data );
						bad_data = 0;
						break;
					}
					fgets( buf_data, 1000, infile_data );
					if( feof( infile_data ) ) { printf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
					white_trim( buf_data );
					pnt_data = &buf_data[0];
					word_data = NULL;
				}
				if( bad_data == 1 )
				{
					printf( "\nERROR: Search keyword \'%s\' cannot be found in the data file \'%s\'!\n", word_inst, fn_in_d );
					return( -1 );
				}
			}
			else if( strncmp( word_inst, dummy_var, 5 ) == 0 ) // dummy variable
			{
				if( debug ) printf( "Skip dummy data!\n" );
				word_data = strtok_r( NULL, separator, &pnt_data );
				if( debug ) printf( "Current data location: \'%s\' Remaining line: \'%s\'\n", word_data, pnt_data );
			}
			else if( word_inst[0] == 'w' ) // white space
			{
				if( debug ) printf( "Skip white space!\n" );
				if( pnt_data[0] != ' ' )
				{
					if( word_data == NULL ) word_data = strtok_r( NULL, separator, &pnt_data );
					word_data = strtok_r( NULL, separator, &pnt_data );
				}
				else
				{
					word_data = strtok_r( NULL, separator, &pnt_data );
				}
				if( debug ) printf( "Current data location: \'%s\' Remaining line: \'%s\'\n", word_data, pnt_data );
			}
			else if( word_inst[0] == token_obs[0] ) // observation variable
			{
				c++;
				if( word_data == NULL || c > 1 ) word_data = strtok_r( NULL, separator, &pnt_data );
				if( strlen( word_inst ) == 1 ) word_inst = strtok_r( NULL, separator, &pnt_inst );
				else word_inst = &word_inst[1];
				if( word_inst[strlen( word_inst ) - 1] == token_obs[0] ) word_inst[strlen( word_inst ) - 1] = 0;
				else strtok_r( NULL, separator, &pnt_inst );
				if( debug ) printf( "Observation keyword \'%s\' & data field \'%s\' ... ", word_inst, word_data );
				white_skip( &word_inst );
				white_trim( word_inst );
				for( i = 0; i < nobs; i++ )
				{
					if( strcmp( word_inst, obs_id[i] ) == 0 )
					{
						sscanf( word_data, "%lf", &v );
						if( check[i] < 0 ) { obs[i] = v; check[i] = 1; }
						else { obs[i] += v; check[i] += 1; }
						if( debug ) printf( "\'%s\'=%g\n", obs_id[i], obs[i] );
						break;
					}
				}
				if( nobs == i )
				{
					printf( "\nERROR: Observation keyword \'%s\' does not match any of observation variables!\n", word_inst );
					bad_data = 1;
				}
			}
			else
			{
				printf( "ERROR: Instruction file %s does not follow the expected format!\n", fn_in_i );
				printf( "White space (w), search (%s) or observation (%s) tokens are expected!\n", token_search, token_obs );
				bad_data = 1;
				break;
			}
		}
	}
	fclose( infile_data ); fclose( infile_inst );
	if( bad_data ) return( -1 );
	else return( 0 );
}

int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug )
{
	FILE *in;
	char *sep = " \t\n"; // White spaces
	char *word, token[2], buf[1000];
	int i, l, c, start = 0, bad_data = 0;
	if( ( in = fopen( fn_in_t, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_t );
		return( -1 );
	}
	if( debug ) printf( "\nChecking the template file \'%s\'.\n", fn_in_t );
	fgets( buf, 1000, in );
	for( c = 0, word = strtok( buf, sep ); word; c++, word = strtok( NULL, sep ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word );
			if( strstr( word, "ptf" ) )
			{
				if( debug ) printf( "PEST Template file\n" );
			}
			else if( strcasestr( word, "template" ) )
			{
				if( debug ) printf( "MADS Template file; user-specified parameter token is expected\n" );
			}
			else
			{
				if( debug ) printf( "MADS Template file\n" );
				rewind( in );
				token[0] = '#'; // default tokes
				break; // quit the loop; done
			}
		}
		if( c == 1 ) // second entry in the case of PEST Template file
		{
			white_trim( word );
			if( strlen( word ) > 1 )
				printf( "WARNING: expecting a single character as parameter keyword separator on the first line of template file (\'%s\'; assumed \'%s\')\n", word, token );
			token[0] = word[0];
			if( token[0] == 0 ) token[0] = '#';
			break;
		}
	}
	token[1] = 0;
	if( debug ) printf( "Parameter separator: %s\n", token );
	while( !feof( in ) )
	{
		fgets( buf, 1000, in );
		if( feof( in ) ) break; // fgets does not produce EOF after reading the last line ...
		l = strlen( buf );
		buf[l - 1] = 0;
		if( buf[0] == token[0] ) start = 0; else start = 1;
		for( c = 0, word = strtok( buf, token ); word; c++, word = strtok( NULL, token ) ) // separation between the tokens is expected; e.g. "# a   # space # b  #"
		{
//			printf( "%d %s\n", c, word );
			if( c % 2 == start )
			{
				if( debug ) printf( "Parameter keyword \'%s\' ", word );
				l = strlen( word );
				white_skip( &word );
				white_trim( word );
				for( i = 0; i < npar; i++ )
				{
					if( strcmp( word, par_id[i] ) == 0 )
					{
						if( debug ) printf( "will be replaced with the value of model parameter \'%s\'\n", par_id[i] );
						if( par[i] < 0 ) par[i] = 1;
						else par[i] += 1;
						break;
					}
				}
				if( i == npar )
				{
					if( debug ) printf( "ERROR: does not match defined model parameters!!!\n" );
					else printf( "\nERROR: Parameter keyword \'%s\' in template file \'%s\' does not match defined model parameters!\n", word, fn_in_t );
					bad_data = 1;
				}
			}
		}
	}
	fclose( in );
	if( bad_data == 1 ) return( -1 );
	else return( 0 );
}

int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug )
{
	FILE *in, *out;
	char *sep = " \t\n";
	char *word, token[2], number[80], buf[1000];
	int i, l, l2, c, start, space = 0, bad_data = 0;
	if( ( in = fopen( fn_in_t, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_t );
		return( -1 );
	}
	sprintf( buf, "rm %s >& /dev/null", fn_out ); system( buf );
	if( ( out = fopen( fn_out, "w" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to write data!\n", fn_out );
		return( -1 );
	}
	if( debug ) printf( "\nCreating model input file \'%s\' for external model execution using template file \'%s\'.\n", fn_out, fn_in_t );
	fgets( buf, 1000, in );
	for( c = 0, word = strtok( buf, sep ); word; c++, word = strtok( NULL, sep ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word );
			if( strcasestr( word, "ptf" ) )
			{
				if( debug ) printf( "PEST Template file\n" );
			}
			else if( strcasestr( word, "template" ) )
			{
				if( debug ) printf( "MADS Template file; user-specified parameter token is expected\n" );
			}
			else
			{
				if( debug ) printf( "MADS Template file\n" );
				rewind( in );
				token[0] = '#';
				break; // quit the loop; done
			}
		}
		if( c == 1 ) // second entry in the case of PEST Template file
		{
			white_trim( word );
			if( strlen( word ) > 1 )
				printf( "WARNING: expecting a single character as parameter keyword separator on the first line of template file (\'%s\'; assumed \'%s\')\n", word, token );
			token[0] = word[0];
			if( token[0] == 0 ) token[0] = '#';
			break;
		}
	}
	token[1] = 0;
	if( debug ) printf( "Parameter separator: %s\n", token );
	while( !feof( in ) )
	{
		fgets( buf, 1000, in );
		if( feof( in ) ) break; // fgets does not produce EOF after reading the last line ...
		l = strlen( buf );
		buf[l - 1] = 0; // remove 'new line' character
		if( buf[0] == token[0] ) start = 0; else start = 1; // if first character is a token it will be not considered a separator
		space = 0;
		for( c = 0, word = strtok( buf, token ); word; c++, word = strtok( NULL, token ) ) // separation between the tokens is expected; e.g. "# a   # space # b  #"
		{
			if( c % 2 == start )
			{
				if( debug ) printf( "Parameter keyword \'%s\' ", word );
				l = strlen( word );
				white_skip( &word );
				white_trim( word );
				for( i = 0; i < npar; i++ )
				{
					if( strcmp( word, par_id[i] ) == 0 )
					{
						sprintf( number, "%-*.*g", l, l - 5, par[i] );
						l2 = strlen( number );
						if( l2 > l ) number[l] = 0;
						if( space ) fprintf( out, " %s", number );
						else { space = 1; fprintf( out, "%s", number ); }
						if( debug ) printf( "replaced with \'%s\' (parameter \'%s\')\n", number, par_id[i] );
						break;
					}
				}
				if( i == npar )
				{
					if( debug ) printf( "ERROR: does not match defined model parameters!!!\n" );
					else printf( "\nERROR: Parameter keyword \'%s\' in template file \'%s\' does not match defined model parameters!\n", word, fn_in_t );
					bad_data = 1;
				}
			}
			else
			{
				if( space ) fprintf( out, " %s", word );
				else { space = 1; fprintf( out, "%s", word ); }
			}
		}
		fprintf( out, "\n" );
	}
	fclose( in ); fclose( out );
	if( bad_data == 1 ) return( -1 );
	else return( 0 );
}
