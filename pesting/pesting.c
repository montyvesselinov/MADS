#define _GNU_SOURCE
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<ctype.h>
#include<math.h>

#include "../mads.h"

#define iswhite(c) ((c)== ' ' || (c)=='\t' || (c)=='\n')

/* Functions here */
int load_pst( char *filename, struct calc_data *cd, struct param_data *pd, struct obs_data *od, struct extrn_data *extrn );
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
static char *white_trim( char *x );
static void white_skip( char **s );
char **char_matrix( int maxCols, int maxRows );

int load_pst( char *filename, struct calc_data *cd, struct param_data *pd, struct obs_data *od, struct extrn_data *extrn )
{
	FILE *in;
	double d;
	char code[20], buf[1000];
	int i, k, npar_groups, nobs_groups;
	if(( in = fopen( filename, "r" ) ) == NULL )
	{
		printf( "PEST control file %s cannot be opened to read problem data!\n", filename );
		return( -1 );
	}
	cd->opt_method = ( char * ) malloc( 50 * sizeof( char ) );
	cd->solution_id = ( char * ) malloc( 50 * sizeof( char ) );
	strcpy( cd->solution_id, "extertnal" );
	cd->solution_type = -1;
	fgets( buf, 1000, in );
	fgets( buf, 1000, in );
	fgets( buf, 1000, in );
	fgets( buf, 1000, in );
	sscanf( buf, "%d %d %d %*d %d", &( *pd ).nParam, &( *od ).nObs, &npar_groups, &nobs_groups );
	printf( "Parameters = %d (groups %d)\n", pd->nParam, npar_groups );
	printf( "Observations = %d (groups %d)\n", od->nObs, nobs_groups );
	fgets( buf, 1000, in );
	sscanf( buf, "%d %d", &( *extrn ).ntpl, &( *extrn ).nins );
	printf( "Number of template files = %d\nNumber of instruction files = %d\n", ( *extrn ).ntpl, ( *extrn ).nins );
	pd->var_id = char_matrix(( *pd ).nParam, 50 );
	pd->var = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_current = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	cd->var = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc(( *pd ).nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc(( *pd ).nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	printf( "Parameters = %d:\n", pd->nParam );
	for( i = 0; i < 6; i++ )
		fgets( buf, 1000, in );
	for( i = 0; i < npar_groups; i++ )
		fgets( buf, 1000, in );
	fgets( buf, 1000, in );
	pd->nOptParam = 0;
	for( i = 0; i < pd->nParam; i++ )
	{
		fscanf( in, "%s %s %*s %lf %lf %lf %*s %*f %*f %*f\n", pd->var_id[i], code, &pd->var[i], &pd->var_min[i], &pd->var_max[i] );
		printf( "%-26s: init %15.12g min %12g max %12g\n", pd->var_id[i], pd->var[i], ( *pd ).var_min[i], ( *pd ).var_max[i] );
		if( strcmp( code, "fixed" ) == 0 ) pd->var_opt[i] = 0; else { pd->nOptParam++; pd->var_opt[i] = 1; }
		if( strcmp( code, "log" ) == 0 ) pd->var_log[i] = 1; else pd->var_log[i] = 0;
		if(( *pd ).var_log[i] == 1 )
		{
			pd->var[i] = log10( pd->var[i] );
			pd->var_min[i] = log10( pd->var_min[i] );
			pd->var_max[i] = log10( pd->var_max[i] );
		}
		pd->var_range[i] = pd->var_max[i] - pd->var_min[i];
		pd->var_dx[i] = pd->var_range[i] / 10;
	}
	pd->var_index = ( int * ) malloc(( *pd ).nOptParam * sizeof( int ) );
	printf( "Optimized parameters = %d\n", pd->nOptParam );
	for( k = i = 0; i < ( *pd ).nParam; i++ )
		if(( *pd ).var_opt[i] == 1 )
		{
			if(( *pd ).var_log[i] == 1 ) d = log10( pd->var[i] ); else d = pd->var[i];
			printf( "%-26s: init %15.12g min %12g max %12g\n", pd->var_id[i], d, ( *pd ).var_min[i], ( *pd ).var_max[i] );
			( *pd ).var_index[k++] = i;
		}
	fgets( buf, 1000, in );
	for( i = 0; i < nobs_groups; i++ )
		fgets( buf, 1000, in );
	fgets( buf, 1000, in );
	od->obs_id = char_matrix(( *od ).nObs, 50 );
	od->obs_target = ( double * ) malloc(( *od ).nObs * sizeof( double ) );
	od->obs_weight = ( double * ) malloc(( *od ).nObs * sizeof( double ) );
	od->obs_min = ( double * ) malloc(( *od ).nObs * sizeof( double ) );
	od->obs_max = ( double * ) malloc(( *od ).nObs * sizeof( double ) );
	od->obs_current = ( double * ) malloc(( *od ).nObs * sizeof( double ) );
	od->res = ( double * ) malloc(( *od ).nObs * sizeof( double ) );
	od->obs_log = ( int * ) malloc(( *od ).nObs * sizeof( int ) );
	for( i = 0; i < od->nObs; i++ )
		fscanf( in, "%s %lf %lf %*s\n", od->obs_id[i], &od->obs_target[i], &od->obs_weight[i] );
	printf( "Calibration targets = %d\n", ( *od ).nObs );
	for( i = 0; i < od->nObs; i++ )
	{
		printf( "%-13s: value %15.12g weight %g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i] );
		od->obs_min[i] = 0; od->obs_max[i] = od->obs_target[i] * 2;
		od->obs_log[i] = 0;
	}
	fgets( buf, 1000, in );
	extrn->cmdline = ( char * ) malloc( 80 * sizeof( char ) );
	fgets( extrn->cmdline, 80, in );
	extrn->cmdline[strlen( extrn->cmdline ) - 1] = 0;
	printf( "Execution command: %s\n", extrn->cmdline );
	printf( "External files:\n" );
	extrn->fn_ins = char_matrix( extrn->nins, 80 );
	extrn->fn_obs = char_matrix( extrn->nins, 80 );
	extrn->fn_tpl = char_matrix( extrn->ntpl, 80 );
	extrn->fn_out = char_matrix( extrn->ntpl, 80 );
	fgets( buf, 1000, in );
	for( i = 0; i < extrn->ntpl; i++ )
		fscanf( in, "%s %s\n", extrn->fn_tpl[i], extrn->fn_out[i] );
	printf( "- to provide current model parameters:\n" );
	for( i = 0; i < extrn->ntpl; i++ )
		printf( "%s -> %s\n", extrn->fn_tpl[i], extrn->fn_out[i] );
	for( i = 0; i < extrn->nins; i++ )
		fscanf( in, "%s %s\n", extrn->fn_ins[i], extrn->fn_obs[i] );
	printf( "- to read current model predictions:\n" );
	for( i = 0; i < extrn->nins; i++ )
		printf( "%s <- %s\n", extrn->fn_ins[i], extrn->fn_obs[i] );
	fclose( in );
	printf( "\n" );
	return( 0 );
}

int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_i, int debug )
{
	FILE *in_i;
	char *sep = " \t\n";
	char *sep_s = "!";
	char *word_i, token[2], buf_i[1000], *p_i;
	int i, c, bad_data = 0;
	if(( in_i = fopen( fn_in_i, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_i );
		return( -1 );
	}
	if( debug ) printf( "\nChecking instruction file \'%s\' ... \n", fn_in_i );
	fgets( buf_i, 1000, in_i );
	for( c = 0, word_i = strtok( buf_i, sep ); word_i; c++, word_i = strtok( NULL, sep ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word_i );
			if( strcasestr( word_i, "ins" ) )
			{
				if( debug ) printf( "PEST Instruction file\n" );
			}
			else
			{
				if( debug ) printf( "MADS Instruction file\n" );
				rewind( in_i );
				token[0] = '@';
				break;
			}
		}
		if( c == 1 ) // second entry
		{
			white_trim( word_i );
			token[0] = word_i[0];
			if( strlen( word_i ) > 1 )
				printf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_i, token );
			if( token[0] == 0 ) token[0] = '@';
			break;
		}
	}
	token[1] = 0;
	if( debug )
	{
		printf( "Search separator: %s\n", token );
		printf( "Parameter separator: %s\n", sep_s );
	}
	while( !feof( in_i ) )
	{
		fgets( buf_i, 1000, in_i );
		white_trim( buf_i );
		if( debug ) printf( "\nCurrent instruction line: %s\n", buf_i );
		if( buf_i[0] == 'l' )
		{
			word_i = strtok_r( buf_i, sep, &p_i );
			sscanf( &word_i[1], "%d", &c );
			if( debug ) printf( "Skip %d lines\n", c );
		}
		else if( buf_i[0] == token[0] )
		{
			word_i = strtok_r( buf_i, token, &p_i );
			white_trim( word_i );
			white_skip( &word_i );
			if( debug ) printf( "Search for  keyword \'%s\' in the data file ...", word_i );
		}
		for( word_i = strtok_r( NULL, sep, &p_i ); word_i; word_i = strtok_r( NULL, sep, &p_i ) )
		{
			white_skip( &word_i );
			white_trim( word_i );
			if( debug ) printf( "Current template keyword %s : ", word_i );
			if( strncmp( word_i, "!dum!", 5 ) == 0 )
			{
				if( debug ) printf( "Skip dummy data! \n" );
			}
			else if( word_i[0] == 'w' )
			{
				if( debug ) printf( "Skip white space! " );
			}
			else if( word_i[0] == sep_s[0] )
			{
				c = 0;
				if( strlen( word_i ) == 1 ) word_i = strtok_r( NULL, sep, &p_i );
				else word_i = &word_i[1];
				if( word_i[strlen( word_i )-1] == sep_s[0] ) word_i[strlen( word_i )-1] = 0;
				else strtok_r( NULL, sep, &p_i );
				if( debug ) printf( "Observation keyword \'%s\' ", word_i );
				white_skip( &word_i );
				white_trim( word_i );
				for( i = 0; i < nobs; i++ )
				{
					if( strcmp( word_i, obs_id[i] ) == 0 )
					{
						if( debug ) printf( "will be applied to read observation \'%s\'", obs_id[i] );
						if( obs[i] < 0 ) obs[i] = 1;
						else obs[i] += 1;
						break;
					}
				}
				if( nobs == i )
				{
					printf( "\nERROR: Observation keyword \'%s\' in instruction file \'%s\' does not match any of observation variables!\n", word_i, fn_in_i );
					bad_data = 1;
				}
			}
			if( debug ) printf( "\n" );
		}
	}
	fclose( in_i );
	if( bad_data ) return( -1 );
	else return( 0 );
}

int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_i, char *fn_in_d, int debug )
{
	FILE *in_i, *in_d;
	char *sep = " \t\n";
	char *sep_s = "!";
	char *word_i, *word_d, token[2], buf_d[1000], buf_i[1000], *p_i, *p_d;
	int i, c, bad_data = 0;
	double v;
	if(( in_i = fopen( fn_in_i, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_i );
		return( -1 );
	}
	if(( in_d = fopen( fn_in_d, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read the model-predicted observations!\n", fn_in_d );
		return( -1 );
	}
	if( debug ) printf( "\nReading output file \'%s\' obtained from external model execution using instruction file \'%s\'.\n", fn_in_d, fn_in_i );
	fgets( buf_i, 1000, in_i );
	for( c = 0, word_i = strtok( buf_i, sep ); word_i; c++, word_i = strtok( NULL, sep ) )
	{
		if( c == 0 ) // first entry
		{
			white_trim( word_i );
			if( strcasestr( word_i, "ins" ) )
			{
				if( debug ) printf( "PEST Instruction file\n" );
			}
			else
			{
				if( debug ) printf( "MADS Instruction file\n" );
				rewind( in_i );
				token[0] = '@';
				break;
			}
		}
		if( c == 1 ) // second entry
		{
			white_trim( word_i );
			token[0] = word_i[0];
			if( strlen( word_i ) > 1 )
				printf( "WARNING: expecting a single character as search separator on the first line of instruction file (\'%s\'; assumed \'%s\')\n", word_i, token );
			if( token[0] == 0 ) token[0] = '@';
			break;
		}
	}
	token[1] = 0;
	if( debug )
	{
		printf( "Search separator: %s\n", token );
		printf( "Parameter separator: %s\n", sep_s );
	}
	while( 1 )
	{
		fgets( buf_i, 1000, in_i );
		if( feof( in_i ) ) break;
		white_trim( buf_i );
		if( debug ) printf( "\nCurrent instruction line: %s\n", buf_i );
		if( buf_i[0] == 'l' )
		{
			word_i = strtok_r( buf_i, sep, &p_i );
			sscanf( &word_i[1], "%d", &c );
			if( debug ) printf( "Skip %d lines\n", c );
			for( i = 0; i < c; i++ )
				fgets( buf_d, 1000, in_d );
			if( feof( in_d ) ) { printf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
			white_trim( buf_d );
			p_d = &buf_d[0];
			word_d = strtok_r( p_d, sep, &p_d );
		}
		else if( buf_i[0] == token[0] )
		{
			word_i = strtok_r( buf_i, token, &p_i );
			white_trim( word_i );
			white_skip( &word_i );
			if( debug ) printf( "Search for  keyword \'%s\' in the data file ...", word_i );
			bad_data = 1;
			while( !feof( in_d ) )
			{
				fgets( buf_d, 1000, in_d );
				if( feof( in_d ) ) { printf( "\nERROR: Model output file \'%s\' is incomplete or instruction file \'%s\' is inaccurate!\n       Model output file \'%s\' ended before instruction file \'%s\' is completely processed!\n", fn_in_d, fn_in_i, fn_in_d, fn_in_i ); break; }
				if(( p_d = strstr( buf_d, word_i ) ) != NULL )
				{
					white_trim( buf_d );
					if( debug ) printf( "\nMatching line found in the data file: %s\n", buf_d );
					word_d = strtok_r( p_d, sep, &p_d );
					bad_data = 0;
					break;
				}
			}
			if( bad_data == 1 )
			{
				printf( "\nERROR: Search keyword \'%s\' cannot be found in the data file \'%s\'!\n", word_i, fn_in_d );
				return( -1 );
			}
		}
		if( debug ) printf( "Current location in model output file: \'%s\' %s\n", word_d, p_d );
		for( word_i = strtok_r( NULL, sep, &p_i ); word_i; word_i = strtok_r( NULL, sep, &p_i ) )
		{
			white_skip( &word_i );
			white_trim( word_i );
			if( debug ) printf( "Template word %s : ", word_i );
			if( strncmp( word_i, "!dum!", 5 ) == 0 )
			{
				if( debug ) printf( "Skip dummy data! \n" );
				word_d = strtok_r( NULL, sep, &p_d );
				if( debug ) printf( "Current data location: \'%s\' %s\n", word_d, p_d );
			}
			else if( word_i[0] == 'w' )
			{
				if( debug ) printf( "Skip white space! " );
				word_d = strtok_r( NULL, sep, &p_d );
				if( debug ) printf( "Current data location: \'%s\' %s\n", word_d, p_d );
			}
			else if( word_i[0] == sep_s[0] )
			{
				c = 0;
				if( strlen( word_i ) == 1 ) word_i = strtok_r( NULL, sep, &p_i );
				else word_i = &word_i[1];
				if( word_i[strlen( word_i )-1] == sep_s[0] ) word_i[strlen( word_i )-1] = 0;
				else strtok_r( NULL, sep, &p_i );
				if( debug ) printf( "Observation keyword \'%s\' & data field \'%s\' ... ", word_i, word_d );
				white_skip( &word_i );
				white_trim( word_i );
				for( i = 0; i < nobs; i++ )
				{
					if( strcmp( word_i, obs_id[i] ) == 0 )
					{
						sscanf( word_d, "%lf", &v );
						if( check[i] < 0 ) { obs[i] = v; check[i] = 1; }
						else { obs[i] += v; check[i] += 1; }
						if( debug ) printf( "\'%s\'=%g\n", obs_id[i], obs[i] );
						break;
					}
				}
				if( nobs == i )
				{
					printf( "\nERROR: Observation keyword \'%s\' does not match any of observation variables!\n", word_i );
					bad_data = 1;
				}
			}
		}
	}
	fclose( in_d ); fclose( in_i );
	if( bad_data ) return( -1 );
	else return( 0 );
}

int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug )
{
	FILE *in;
	char *sep = " \t\n"; // White spaces
	char *word, token[2], buf[1000], *p;
	int i, l, c, bad_data = 0;
	if(( in = fopen( fn_in_t, "r" ) ) == NULL )
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
			else
			{
				if( debug ) printf( "MADS Template file\n" );
				rewind( in );
				token[0] = '#';
				break;
			}
		}
		if( c == 1 ) // second entry
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
		buf[l-1] = 0;
		p = &buf[0];
		for( c = 0, word = strtok( buf, token ); word; c++, word = strtok( NULL, token ) )
		{
			if( buf[0] == token[0] || c % 2 == 1 )
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
	char *word, token[2], number[80], buf[1000], *p;
	int i, l, c, bad_data = 0;
	if(( in = fopen( fn_in_t, "r" ) ) == NULL )
	{
		printf( "\n\nERROR: File %s cannot be opened to read template data!\n", fn_in_t );
		return( -1 );
	}
	if(( out = fopen( fn_out, "w" ) ) == NULL )
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
			else
			{
				if( debug ) printf( "MADS Template file\n" );
				rewind( in );
				token[0] = '#';
				break;
			}
		}
		if( c == 1 ) // second entry
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
		buf[l-1] = 0;
		p = &buf[0];
		for( c = 0, word = strtok( buf, token ); word; c++, word = strtok( NULL, token ) )
		{
			if( buf[0] == token[0] || c % 2 == 1 )
			{
				if( debug ) printf( "Parameter keyword \'%s\' ", word );
				l = strlen( word );
				white_skip( &word );
				white_trim( word );
				for( i = 0; i < npar; i++ )
				{
					if( strcmp( word, par_id[i] ) == 0 )
					{
						sprintf( number, "%-30.25g", par[i] );
						number[l] = 0;
						fprintf( out, "%s", number );
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
			else fprintf( out, " %s", word );
		}
		fprintf( out, "\n" );
	}
	fclose( in ); fclose( out );
	if( bad_data == 1 ) return( -1 );
	else return( 0 );
}

static void white_skip( char **s )
{
	while( iswhite( **s ) )( *s )++;
}

static char *white_trim( char *x )
{
	char *y;
	if( !x ) return( x );
	y = x + strlen( x ) - 1;
	while( y >= x && iswhite( *y ) ) *y-- = 0;
	return x;
}
