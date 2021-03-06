// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
// http://www.ees.lanl.gov/staff/monty/codes/mads
// http://gitlab.com/monty/mads
//
// Licensing: GPLv3: http://www.gnu.org/licenses/gpl-3.0.html
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<ctype.h>
#include<math.h>

#include "mads.h"

/* Functions here */
int load_pst( char *filename, struct opt_data *op );
int check_ins_obs( int nobs, char **obs_id, int *obs_count, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, int *obs_count, char *fn_in_t, char *fn_in_d, int debug );
int check_par_tpl( int npar, char **par_id, int *par_count, char *fn_in_t, int debug );
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
	int i, j, k, npar_groups, nobs_groups, bad_data = 0;
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od;
	struct extrn_data *ed;
	cd = op->cd;
	pd = op->pd;
	od = op->od;
	ed = op->ed;
	pd->nParam = pd->nFlgParam = pd->nOptParam = 0;
	od->nTObs = od->nCObs = od->nObs = 0;
	ed->ntpl = ed->nins = 0;
	bad_data = 0;
	op->gd->min_t = op-> gd->time = 0;
	if( ( in = fopen( filename, "r" ) ) == NULL )
	{
		tprintf( "PEST control file %s cannot be opened to read problem data!\n", filename );
		return( -1 );
	}
	cd->opt_method = ( char * ) malloc( 50 * sizeof( char ) );
	cd->solution_id = ( char * ) malloc( 50 * sizeof( char ) );
	cd->solution_type = ( int * ) malloc( 1 * sizeof( int ) );
	strcpy( cd->solution_id, "external" );
	cd->num_sources = 1;
	cd->solution_type = ( int * ) malloc( sizeof( int ) );
	cd->solution_type[0] = EXTERNAL;
	for( i = 0; i < 4; i++ ) // skip 4 lines
		fgets( buf, 1000, in );
	sscanf( buf, "%d %d %d %*d %d", &pd->nParam, &od->nTObs, &npar_groups, &nobs_groups );
	tprintf( "Parameters = %d (groups %d)\n", pd->nParam, npar_groups );
	tprintf( "Observations = %d (groups %d)\n", od->nTObs, nobs_groups );
	od->nObs = od->nCObs = od->nTObs;
	fgets( buf, 1000, in );
	sscanf( buf, "%d %d", &ed->ntpl, &ed->nins );
	tprintf( "Number of template files = %d\nNumber of instruction files = %d\n", ed->ntpl, ed->nins );
	pd->var_name = char_matrix( pd->nParam, 50 );
	pd->var_id = char_matrix( pd->nParam, 50 );
	pd->var = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_current = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_best = ( double * ) malloc( pd->nParam * sizeof( double ) );
	cd->var = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc( pd->nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc( pd->nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc( pd->nParam * sizeof( double ) );
	tprintf( "Parameters = %d:\n", pd->nParam );
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
		strcpy( pd->var_name[i], pd->var_id[i] );
		tprintf( "%-27s: init %15.12g min %12g max %12g\n", pd->var_name[i], pd->var[i], pd->var_min[i], pd->var_max[i] );
		if( strcmp( code, "fixed" ) == 0 ) pd->var_opt[i] = 0; else { pd->nOptParam++; pd->var_opt[i] = 1; }
		if( strcmp( code, "log" ) == 0 ) pd->var_log[i] = 1; else pd->var_log[i] = 0;
		if( pd->var_log[i] == 1 )
		{
			pd->var[i] = log10( pd->var[i] );
			pd->var_min[i] = log10( pd->var_min[i] );
			pd->var_max[i] = log10( pd->var_max[i] );
		}
		pd->var_range[i] = pd->var_max[i] - pd->var_min[i];
		pd->var_dx[i] = pd->var_range[i] / 10;
	}
	pd->var_index = ( int * ) malloc( pd->nOptParam * sizeof( int ) );
	tprintf( "Optimized parameters = %d\n", pd->nOptParam );
	for( k = i = 0; i < pd->nParam; i++ )
		if( pd->var_opt[i] == 1 )
		{
			if( pd->var_log[i] == 1 ) d = log10( pd->var[i] ); else d = pd->var[i];
			tprintf( "%-27s: init %15.12g min %12g max %12g\n", pd->var_name[i], d, pd->var_min[i], pd->var_max[i] );
			pd->var_index[k++] = i;
		}
	for( i = 0; i < pd->nParam; i++ )
		for( j = i + 1; j < pd->nParam; j++ )
			if( strcmp( pd->var_name[i], pd->var_name[j] ) == 0 )
			{
				tprintf( "ERROR: Parameter names #%i (%s) and #%i (%s) are identical!\n", i + 1, pd->var_name[i], j + 1, pd->var_name[j] );
				bad_data = 1;
			}
	if( bad_data ) return( 0 );
	fgets( buf, 1000, in ); // skip line
	for( i = 0; i < nobs_groups; i++ )
		fgets( buf, 1000, in );
	fgets( buf, 1000, in ); // skip line
	od->obs_id = char_matrix( od->nTObs, 50 );
	od->obs_target = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_weight = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_min = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_max = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_current = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_best = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->res = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_log = ( int * ) malloc( od->nTObs * sizeof( int ) );
	for( i = 0; i < od->nTObs; i++ )
		fscanf( in, "%s %lf %lf %*s\n", od->obs_id[i], &od->obs_target[i], &od->obs_weight[i] );
	tprintf( "Calibration targets = %d\n", od->nTObs );
	for( i = 0; i < od->nTObs; i++ )
	{
		if( od->nTObs < 50 || ( i < 20 || i > od->nTObs - 20 ) )
			tprintf( "%-13s: value %15.12g weight %g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i] );
		if( od->nTObs > 50 && i == 21 ) tprintf( "...\n" );
		od->obs_min[i] = 0; od->obs_max[i] = od->obs_target[i] * 2;
		od->obs_log[i] = 0;
	}
	if( od->nObs < 10000 || cd->problem_type == CHECK || cd->debug > 10 )
	{
		tprintf( "Checking for duplicate observations ... \n" );
		if( od->nObs >= 10000 ) tprintf( "WARNING: The number of observations is large (%d); this may take a long time ... \n", od->nObs );
		for( i = 0; i < od->nTObs; i++ )
			for( j = i + 1; j < od->nTObs; j++ )
				if( strcmp( od->obs_id[i], od->obs_id[j] ) == 0 )
				{
					tprintf( "ERROR: Observation names #%i (%s) and #%i (%s) are identical!\n", i + 1, od->obs_id[i], j + 1, od->obs_id[j] );
					bad_data = 1;
				}
	}
	if( bad_data ) return( 0 );
	fgets( buf, 1000, in ); // skip line
	ed->cmdline = ( char * ) malloc( 255 * sizeof( char ) );
	fgets( ed->cmdline, 255, in );
	ed->cmdline[strlen( ed->cmdline ) - 1] = 0;
	tprintf( "Execution command: %s\n", ed->cmdline );
	tprintf( "External files:\n" );
	ed->fn_ins = char_matrix( ed->nins, 255 );
	ed->fn_obs = char_matrix( ed->nins, 255 );
	ed->fn_tpl = char_matrix( ed->ntpl, 255 );
	ed->fn_out = char_matrix( ed->ntpl, 255 );
	fgets( buf, 1000, in ); // skip line
	for( i = 0; i < ed->ntpl; i++ )
		fscanf( in, "%s %s\n", ed->fn_tpl[i], ed->fn_out[i] );
	tprintf( "- to provide current model parameters:\n" );
	for( i = 0; i < ed->ntpl; i++ )
		tprintf( "%s -> %s\n", ed->fn_tpl[i], ed->fn_out[i] );
	for( i = 0; i < ed->nins; i++ )
		fscanf( in, "%s %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	tprintf( "- to read current model predictions:\n" );
	for( i = 0; i < ed->nins; i++ )
		tprintf( "%s <- %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	fclose( in );
	tprintf( "\n" );
	return( 1 );
}

int check_ins_obs( int nobs, char **obs_id, int *obs_count, char *fn_in_i, int debug )
{
	FILE *infile_inst;
	char *separator = " \t\n";
	char *word_inst, *word_search, token_obs[2], token_search[2], comment[2], dummy_var[6], buf_inst[1000], *pnt_inst;
	int i, c, bad_data = 0;
	if( debug ) tprintf( "\nChecking instruction file \'%s\'.\n", fn_in_i );
	if( ( infile_inst = fopen( fn_in_i, "r" ) ) == NULL )
	{
		tprintf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_i );
		return( -1 );
	}
	fgets( buf_inst, 1000, infile_inst );
	if( debug ) tprintf( "\nFirst instruction line: %s\n", buf_inst );
	pnt_inst = &buf_inst[0];
	for( c = 0, word_inst = strtok_r( buf_inst, separator, &pnt_inst ); word_inst; c++, word_inst = strtok_r( NULL, separator, &pnt_inst ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word_inst );
			if( strcasestr( word_inst, "pif" ) )
			{
				if( debug ) tprintf( "PEST Instruction file\n" );
				token_search[0] = '@'; // just in case
				token_obs[0] = '!';
				comment[0] = 0;
			}
			else if( strcasestr( word_inst, "instruction" ) )
			{
				if( debug ) tprintf( "MADS Instruction file; user-specified search/variable tokens are expected\n" );
				token_search[0] = '@'; // just in case
				token_obs[0] = '!';
				comment[0] = '#';
			}
			else
			{
				if( debug ) tprintf( "MADS Instruction file\n" );
				rewind( infile_inst );
				token_search[0] = '@';
				token_obs[0] = '!';
				comment[0] = '#';
				break;
			}
		}
		else if( c == 1 ) // second entry; "search" token
		{
			white_trim( word_inst );
			if( debug > 1 ) tprintf( "Search token %s\n", word_inst );
			token_search[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				tprintf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_search[0] == 0 ) token_search[0] = '@';
		}
		else if( c == 2 ) // third entry; "variable" token
		{
			white_trim( word_inst );
			if( debug > 1 ) tprintf( "Variable token %s\n", word_inst );
			token_obs[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				tprintf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_obs[0] == 0 ) token_obs[0] = '!';
		}
		else if( c == 3 ) // third entry; "comment" token
		{
			white_trim( word_inst );
			if( debug > 1 ) tprintf( "Comment token %s\n", word_inst );
			comment[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				tprintf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( comment[0] == 0 ) comment[0] = '#';
		}
	}
	token_search[1] = token_obs[1] = 0;
	dummy_var[0] = token_obs[0];
	dummy_var[1] = 0;
	strcat( dummy_var, "dum" );
	dummy_var[4] = token_obs[0];
	dummy_var[5] = 0;
	token_obs[1] = token_search[1] = comment[1] = 0;
	if( debug )
	{
		tprintf( "Search separator: %s\n", token_search );
		tprintf( "Observation separator: %s\n", token_obs );
		tprintf( "Dummy observation: %s\n", dummy_var );
		if( comment[0] ) tprintf( "Comment: %s\n", comment );
	}
	while( !feof( infile_inst ) ) // IMPORTANT: strtok below modifies buf_inst by adding '\0's; if needed strcpy buf_inst
	{
		if( fgets( buf_inst, 1000, infile_inst ) == NULL ) { if( debug > 1 ) tprintf( "END of instruction file.\n" ); break; }
		pnt_inst = &buf_inst[0];
		word_inst = 0;
		white_trim( pnt_inst ); white_skip( &pnt_inst );
		if( debug ) tprintf( "\nCurrent instruction line: %s\n", pnt_inst );
		if( comment[0] && pnt_inst[0] == comment[0] ) { if( debug > 1 ) tprintf( "Comment; skip this line.\n" ); continue; } // Instruction line is a comment
		if( strlen( pnt_inst ) == 0 ) { if( debug ) tprintf( "Empty line; will be skipped.\n" ); continue; } // Empty line
		if( pnt_inst[0] == 'l' ) // skip lines in the "data" file
		{
			sscanf( &pnt_inst[1], "%d", &c );
			if( debug > 1 ) tprintf( "Skip %d lines\n", c );
			word_inst = strtok_r( NULL, separator, &pnt_inst ); // skip l command
		}
		while( 1 )
		{
			if( pnt_inst[0] == token_search[0] ) // search for keyword
			{
				if( debug ) tprintf( "KEYWORD search " );
				word_search = strtok_r( NULL, token_search, &pnt_inst ); // read search keyword
				if( debug ) tprintf( "\'%s\' in the data file ...\n", word_search );
			}
			else
			{
				word_inst = strtok_r( NULL, separator, &pnt_inst ); // read TEMPLETE word
				if( debug > 1 ) tprintf( "Current location in instruction input file: => \'%s\' <= \'%s\'\n", word_inst, pnt_inst );
				white_trim( word_inst );
				if( debug ) tprintf( "INSTRUCTION word \'%s\' : ", word_inst );
				if( strncmp( word_inst, dummy_var, 5 ) == 0 ) // dummy variable
				{
					if( debug ) tprintf( "Skip dummy data!\n" );
				}
				else if( word_inst[0] == 'w' ) // white space
				{
					if( debug ) tprintf( "Skip white space!\n" );
				}
				else if( word_inst[0] == token_obs[0] ) // observation variable
				{
					c = 0;
					if( strlen( word_inst ) == 1 ) word_inst = strtok_r( NULL, separator, &pnt_inst );
					else word_inst = &word_inst[1];
					if( word_inst[strlen( word_inst ) - 1] == token_obs[0] ) word_inst[strlen( word_inst ) - 1] = 0;
					else strtok_r( NULL, separator, &pnt_inst );
					if( debug ) tprintf( "Observation keyword \'%s\' ... ", word_inst );
					white_skip( &word_inst );
					white_trim( word_inst );
					for( i = 0; i < nobs; i++ )
					{
						if( strcmp( word_inst, obs_id[i] ) == 0 )
						{
							obs_count[i]++;
							if( debug ) tprintf( "\'%s\' detected %d times\n", obs_id[i], obs_count[i] );
							break;
						}
					}
					if( nobs == i )
					{
						tprintf( "\nERROR: Observation keyword \'%s\' does not match any of observation variables!\n", word_inst );
						bad_data = 1;
					}
				}
				else if( comment[0] && word_inst[0] == comment[0] ) // comment
				{
					if( debug ) tprintf( "Comment. Skip rest of the instruction line!\n" );
					break;
				}
				else
				{
					tprintf( "\nERROR: Instruction file %s does not follow the expected format!\n", fn_in_i );
					tprintf( "White space (w), search (%s) or observation (%s) tokens are expected!\n", token_search, token_obs );
					bad_data = 1;
					break;
				}
			}
			if( pnt_inst == NULL || strlen( pnt_inst ) == 0 ) break;
		}
	}
	fclose( infile_inst );
	if( bad_data ) return( -1 );
	else return( 0 );
}

int ins_obs( int nobs, char **obs_id, double *obs, int *obs_count, char *fn_in_i, char *fn_in_d, int debug )
{
	FILE *infile_inst, *infile_data;
	char *separator = " \t\n";
	char *word_inst, *word_data, *word_search, token_search[2], token_obs[2], comment[2], dummy_var[6], buf_data[1000], buf_inst[1000], *pnt_inst, *pnt_data;
	int i, c, bad_data = 0, sl;
	double v;
	if( ( infile_inst = fopen( fn_in_i, "r" ) ) == NULL )
	{
		tprintf( "\nERROR: File %s cannot be opened to read template data!\n", fn_in_i );
		return( -1 );
	}
	if( ( infile_data = fopen( fn_in_d, "r" ) ) == NULL )
	{
		tprintf( "\nERROR: File %s cannot be opened to read the model-predicted observations!\n", fn_in_d );
		return( -1 );
	}
	if( debug ) tprintf( "\nReading output file \'%s\' obtained from external model execution using instruction file \'%s\'.\n", fn_in_d, fn_in_i );
	fgets( buf_inst, 1000, infile_inst );
	if( debug > 1 ) tprintf( "First instruction line: %s\n", buf_inst );
	pnt_inst = &buf_inst[0];
	for( c = 0, word_inst = strtok_r( buf_inst, separator, &pnt_inst ); word_inst; c++, word_inst = strtok_r( NULL, separator, &pnt_inst ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word_inst );
			if( strcasestr( word_inst, "pif" ) )
			{
				if( debug > 1 ) tprintf( "PEST Instruction file\n" );
				token_search[0] = '@'; // just in case
				token_obs[0] = '!';
				comment[0] = 0;
			}
			else if( strcasestr( word_inst, "instruction" ) )
			{
				if( debug > 1 ) tprintf( "MADS Instruction file; user-specified search/variable tokens are expected\n" );
				token_search[0] = '@'; // just in case
				token_obs[0] = '!';
				comment[0] = '#';
			}
			else
			{
				if( debug > 1 ) tprintf( "MADS Instruction file\n" );
				rewind( infile_inst );
				token_search[0] = '@';
				token_obs[0] = '!';
				comment[0] = '#';
				break;
			}
		}
		else if( c == 1 ) // second entry; "search" token
		{
			white_trim( word_inst );
			token_search[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				tprintf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_search[0] == 0 ) token_search[0] = '@';
		}
		else if( c == 2 ) // third entry; "variable" token
		{
			white_trim( word_inst );
			token_obs[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				tprintf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( token_obs[0] == 0 ) token_obs[0] = '!';
		}
		else if( c == 3 ) // third entry; "comment" token
		{
			white_trim( word_inst );
			comment[0] = word_inst[0];
			if( strlen( word_inst ) > 1 )
				tprintf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_inst, token_search );
			if( comment[0] == 0 ) comment[0] = '#';
			break;
		}
	}
	token_search[1] = token_obs[1] = 0;
	dummy_var[0] = token_obs[0];
	dummy_var[1] = 0;
	strcat( dummy_var, "dum" );
	dummy_var[4] = token_obs[0];
	dummy_var[5] = 0;
	token_obs[1] = token_search[1] = comment[1] = 0;
	if( debug > 1 )
	{
		tprintf( "Search separator: %s\n", token_search );
		tprintf( "Observation separator: %s\n", token_obs );
		tprintf( "Dummy observation: %s\n", dummy_var );
		if( comment[0] ) tprintf( "Comment: %s\n", comment );
	}
	buf_data[0] = 0; word_data = NULL;
	while( !feof( infile_inst ) )
	{
		if( fgets( buf_inst, 1000, infile_inst ) == NULL ) { if( debug > 1 ) tprintf( "END of instruction file.\n" ); break; }
		pnt_inst = &buf_inst[0];
		word_inst = 0;
		white_trim( pnt_inst ); white_skip( &pnt_inst );
		if( debug > 1 ) tprintf( "\n\nCurrent instruction line: %s\n", pnt_inst );
		if( comment[0] && pnt_inst[0] == comment[0] ) { if( debug > 1 ) tprintf( "Comment; skip this line.\n" ); continue; } // Instruction line is a comment
		if( strlen( pnt_inst ) == 0 ) { if( debug ) tprintf( "Empty line; will be skipped.\n" ); continue; }
		pnt_data = NULL;
		if( pnt_inst[0] == 'l' ) // skip lines in the "data" file
		{
			sscanf( &pnt_inst[1], "%d", &c );
			if( debug > 1 ) tprintf( "Skip %d lines\n", c );
			for( i = 0; i < c; i++ )
				if( fgets( buf_data, 1000, infile_data ) == NULL ) { tprintf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
			word_inst = strtok_r( NULL, separator, &pnt_inst ); // skip l command
			if( feof( infile_data ) ) { tprintf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
			white_trim( buf_data );
			pnt_data = &buf_data[0];
			word_data = NULL;
		}
		if( pnt_data == NULL ) // if there was no "l" (skip line) command, read the next "data" line
		{
			if( debug > 1 ) tprintf( "Read the next \'data\' line (there was no \'l\' (skip line) command)\n" );
			fgets( buf_data, 1000, infile_data );
			white_trim( buf_data );
			pnt_data = &buf_data[0];
			word_data = NULL;
		}
		if( debug > 1 ) tprintf( "Current location in model output file: => \'%s\' <= \'%s\'\n", word_data, pnt_data );
		if( debug ) { if( pnt_data != NULL ) { if( pnt_data[strlen( pnt_data ) - 2] != '\n' ) {} } }
		c = 0;
		while( 1 )
		{
			if( pnt_inst[0] == token_search[0] ) // search for keyword
			{
				if( debug > 1 ) tprintf( "KEYWORD search " );
				word_search = strtok_r( NULL, token_search, &pnt_inst ); // read search keyword
				if( debug > 1 ) tprintf( "\'%s\' in the data file ...\n", word_search );
				bad_data = 1;
				while( !feof( infile_data ) )
				{
					if( ( pnt_data = strstr( pnt_data, word_search ) ) != NULL )
					{
						pnt_data += strlen( word_search );
						if( debug > 1 ) tprintf( "Matching data file location \'=>%s<=%s\'\n", word_search, pnt_data );
						bad_data = 0;
						break;
					}
					if( fgets( buf_data, 1000, infile_data ) == NULL ) { tprintf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
					white_trim( buf_data );
					pnt_data = &buf_data[0];
					word_data = NULL; // Force reading
				}
				if( bad_data == 1 )
				{
					tprintf( "\nERROR: Search keyword \'%s\' cannot be found in the data file \'%s\'!\n", word_search, fn_in_d );
					return( -1 );
				}
			}
			else // no keyword search
			{
				word_inst = strtok_r( NULL, separator, &pnt_inst ); // read TEMPLETE word
				if( debug > 1 ) tprintf( "Current location in instruction input file: => \'%s\' <= \'%s\'\n", word_inst, pnt_inst );
				white_trim( word_inst );
				if( debug > 1 ) tprintf( "INSTRUCTION word \'%s\' : ", word_inst );
				if( strncmp( word_inst, dummy_var, 5 ) == 0 ) // dummy variable
				{
					if( debug > 1 ) tprintf( "Skip dummy data!\n" );
					if( word_data == NULL ) word_data = strtok_r( NULL, separator, &pnt_data );
					word_data = strtok_r( NULL, separator, &pnt_data );
					if( debug > 1 ) tprintf( "Current location in model output file: => \'%s\' <= \'%s\'\n", word_data, pnt_data );
				}
				else if( word_inst[0] == 'w' ) // white space
				{
					if( debug > 1 ) tprintf( "Skip white space!\n" );
					if( !iswhite( pnt_data[0] ) )
					{
						if( word_data == NULL ) word_data = strtok_r( NULL, separator, &pnt_data );
						word_data = strtok_r( NULL, separator, &pnt_data );
					}
					else
					{
						word_data = strtok_r( NULL, separator, &pnt_data );
					}
					if( debug > 1 ) tprintf( "Current location in model output file: => \'%s\' <= \'%s\'\n", word_data, pnt_data );
				}
				else if( word_inst[0] == token_obs[0] ) // observation variable
				{
					if( debug ) tprintf( "Observation variable\n" );
					c++;
					if( word_data == NULL || c > 1 ) word_data = strtok_r( NULL, separator, &pnt_data );
					if( strlen( word_inst ) == 1 ) word_inst = strtok_r( NULL, separator, &pnt_inst );
					else word_inst = &word_inst[1];
					sl = strlen( word_inst );
					if( word_inst[sl - 1] == token_obs[0] ) word_inst[sl - 1] = 0;
					else strtok_r( NULL, separator, &pnt_inst );
					white_skip( &word_inst );
					white_trim( word_inst );
					if( debug ) tprintf( "Observation keyword \'%s\' & data field \'%s\' ... ", word_inst, word_data );
					if( word_data == NULL || strlen( word_data ) == 0 )
					{
						tprintf( "ERROR: Mismatch between the instruction file \'%s\' and the data file \'%s\'!\n", fn_in_i, fn_in_d );
						tprintf( "INSTRUCTION word \'%s\'\n", word_inst );
						tprintf( "Current location in instruction input file: => \'%s\' <= \'%s\'\n", word_inst, pnt_inst );
						tprintf( "Current location in model output file: => \'%s\' <= \'%s\'\n", word_data, pnt_data );
						bad_data = 1;
						break;
					}
					for( i = 0; i < nobs; i++ )
					{
						if( strcmp( word_inst, obs_id[i] ) == 0 )
						{
							sscanf( word_data, "%lf", &v );
							if( obs_count[i] == 0 ) { obs[i] = v; obs_count[i] = 1; }
							else { obs[i] += v; obs_count[i]++; }
							if( debug ) tprintf( "\'%s\'=%d\n", obs_id[i], obs[i] );
							break;
						}
					}
					if( nobs == i )
					{
						tprintf( "\nERROR: Observation keyword \'%s\' does not match any of observation variables!\n", word_inst );
						bad_data = 1;
					}
				}
				else if( comment[0] && word_inst[0] == comment[0] ) // comment
				{
					if( debug > 1 ) tprintf( "Comment. Skip rest of the instruction line!\n" );
					break;
				}
				else
				{
					tprintf( "\nERROR: Instruction file %s does not follow the expected format!\n", fn_in_i );
					tprintf( "White space (w), search (%s) or observation (%s) tokens are expected!\n", token_search, token_obs );
					bad_data = 1;
					break;
				}
			}
			if( pnt_inst == NULL || strlen( pnt_inst ) == 0 ) break;
		}
	}
	fclose( infile_data ); fclose( infile_inst );
	if( bad_data ) return( -1 );
	else return( 0 );
}

int check_par_tpl( int npar, char **par_id, int *par_count, char *fn_in_t, int debug )
{
	FILE *in;
	char *sep = " \t\n"; // White spaces
	char *word, token[2], buf[1000], *pnt_inst;
	int i, l, c, start = 0, bad_data = 0;
	if( ( in = fopen( fn_in_t, "r" ) ) == NULL )
	{
		tprintf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_t );
		return( -1 );
	}
	if( debug ) tprintf( "\nChecking the template file \'%s\'.\n", fn_in_t );
	fgets( buf, 1000, in );
	pnt_inst = &buf[0];
	for( c = 0, word = strtok_r( buf, sep, &pnt_inst ); word; c++, word = strtok_r( NULL, sep, &pnt_inst ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word );
			if( strstr( word, "ptf" ) )
			{
				if( debug ) tprintf( "PEST Template file\n" );
			}
			else if( strcasestr( word, "template" ) )
			{
				if( debug ) tprintf( "MADS Template file; user-specified parameter token is expected\n" );
			}
			else
			{
				if( debug ) tprintf( "MADS Template file\n" );
				rewind( in );
				token[0] = '#'; // default tokes
				break; // quit the loop; done
			}
		}
		if( c == 1 ) // second entry in the case of PEST Template file
		{
			white_trim( word );
			if( strlen( word ) > 1 )
				tprintf( "WARNING: expecting a single character as parameter keyword separator on the first line of template file (\'%s\'; assumed \'%s\')\n", word, token );
			token[0] = word[0];
			if( token[0] == 0 ) token[0] = '#';
			break;
		}
	}
	token[1] = 0;
	if( debug ) tprintf( "Parameter separator: %s\n", token );
	while( !feof( in ) )
	{
		if( fgets( buf, 1000, in ) == NULL ) { if( debug > 1 ) tprintf( "END of template file.\n" ); break; }
		l = strlen( buf );
		buf[l - 1] = 0;
		if( buf[0] == token[0] ) start = 0; else start = 1;
		pnt_inst = &buf[0];
		for( c = 0, word = strtok_r( buf, token, &pnt_inst ); word; c++, word = strtok_r( NULL, token, &pnt_inst ) ) // separation between the tokens is expected; e.g. space # b  #"
		{
//			tprintf( "%d %s\n", c, word );
			if( c % 2 == start )
			{
				if( debug ) tprintf( "Parameter keyword \'%s\' ", word );
				l = strlen( word );
				white_skip( &word );
				white_trim( word );
				for( i = 0; i < npar; i++ )
				{
					if( strcmp( word, par_id[i] ) == 0 )
					{
						if( debug ) tprintf( "will be replaced with the value of model parameter \'%s\'\n", par_id[i] );
						if( par_count[i] < 0 ) par_count[i] = 1;
						else par_count[i] += 1;
						break;
					}
				}
				if( i == npar )
				{
					if( debug ) tprintf( "ERROR: does not match defined model parameters!!!\n" );
					else tprintf( "\nERROR: Parameter keyword \'%s\' in template file \'%s\' does not match defined model parameters!\n", word, fn_in_t );
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
	char *word, token[2], number[80], buf[1000], *pnt_inst;
	word = ( char * ) malloc( 1000 * sizeof( char ) );
	int i, l, l2, c, start, space = 0, bad_data = 0, preserve;
	if( ( in = fopen( fn_in_t, "r" ) ) == NULL )
	{
		tprintf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_t );
		return( -1 );
	}
	if( debug ) tprintf( "Remove files for model inputs: %s\n", fn_out );
	remove( fn_out );
	if( ( out = fopen( fn_out, "w" ) ) == NULL )
	{
		tprintf( "\n\nERROR: File %s cannot be opened to write data!\n", fn_out );
		return( -1 );
	}
	if( debug ) tprintf( "\nCreating model input file \'%s\' for external model execution using template file \'%s\'.\n", fn_out, fn_in_t );
	fgets( buf, 1000, in );
	pnt_inst = &buf[0];
	for( c = 0, word = strtok_r( buf, sep, &pnt_inst ); word; c++, word = strtok_r( NULL, sep, &pnt_inst ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word );
			if( strcasestr( word, "ptf" ) )
			{
				if( debug ) tprintf( "PEST Template file\n" );
			}
			else if( strcasestr( word, "template" ) )
			{
				if( debug ) tprintf( "MADS Template file; user-specified parameter token is expected\n" );
			}
			else
			{
				if( debug ) tprintf( "MADS Template file\n" );
				rewind( in );
				token[0] = '#';
				break; // quit the loop; done
			}
		}
		if( c == 1 ) // second entry in the case of PEST Template file
		{
			white_trim( word );
			if( strlen( word ) > 1 )
				tprintf( "WARNING: expecting a single character as parameter keyword separator on the first line of template file (\'%s\'; assumed \'%s\')\n", word, token );
			token[0] = word[0];
			if( token[0] == 0 ) token[0] = '#';
			break;
		}
	}
	token[1] = 0;
	if( debug > 1 ) tprintf( "Parameter separator: %s\n", token );
	while( !feof( in ) )
	{
		if( fgets( buf, 1000, in ) == NULL ) { if( debug > 1 ) tprintf( "END of template file.\n" ); break; }
		l = strlen( buf );
		buf[l - 1] = 0; // remove 'new line' character
		if( buf[0] == token[0] ) start = 0; else start = 1; // if first character is a token it will be not considered a separator
		space = 0;
		pnt_inst = &buf[0];
		for( c = 0, word = strtok_r( buf, token, &pnt_inst ); word; c++, word = strtok_r( NULL, token, &pnt_inst ) ) // separation between the tokens is expected; e.g. "# a   # space # b  #"
		{
			if( c % 2 == start )
			{
				if( debug ) tprintf( "Parameter keyword \'%s\' ", word );
				l = strlen( word );
				white_skip( &word );
				white_trim( word );
				l2 = strlen( word );
				if( l > ( l2 + 2 ) ) preserve = 1;
				else preserve = 0;
				for( i = 0; i < npar; i++ )
				{
					if( strcmp( word, par_id[i] ) == 0 )
					{
						if( preserve == 1 )
						{
							if( par[i] > 0 ) sprintf( number, "%.*g", l - 1, par[i] );
							else sprintf( number, "%.*g", l - 2, par[i] );
							l2 = strlen( number );
							if( l2 > l ) tprintf( "WARNING: The parameter does not fit the requested field (%s length %d > %d)!\n", number, l2, l );
						}
						else
							sprintf( number, "%.15g", par[i] );
						if( space ) fprintf( out, " %s", number );
						else { space = 0; fprintf( out, "%s", number ); } // TODO originally was space = 1
						if( debug ) tprintf( "is replaced with \'%s\'\n", number );
						break;
					}
				}
				if( i == npar )
				{
					if( debug )	tprintf( "ERROR: does not match defined model parameters!!!\n" );
					else tprintf( "\nERROR: Parameter keyword \'%s\' in template file \'%s\' does not match defined model parameters!\n", word, fn_in_t );
					bad_data = 1;
				}
			}
			else
			{
				if( space ) fprintf( out, " %s", word );
				else { space = 0; fprintf( out, "%s", word ); }  // TODO originally was space = 1
			}
		}
		fprintf( out, "\n" );
	}
	fclose( in ); fclose( out );
	if( bad_data == 1 ) return( -1 );
	else return( 0 );
}
