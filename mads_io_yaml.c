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


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

struct opt_data *gop;

#include "mads.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <yaml.h>
#include <glib.h>

#define C2P(c)          ((gpointer) ((long) (c)))
#define P2C(p)          ((gchar) ((long) (p)))
enum storage_flags { VAR, VAL, SEQ }; // "Store as" switch

/* Functions here */
int load_yaml_problem( char *filename, int argn, char *argv[], struct opt_data *op );
void yaml_parse_layer( yaml_parser_t *parser, GNode *data );
gboolean gnode_tree_dump( GNode *n, gpointer data );
void gnode_tree_dump_classes( GNode *n, gpointer data );
void gnode_tree_parse_classes( GNode *node, gpointer data );
int load_ymal_wells( GNode *node, gpointer data );
int load_ymal_grid( GNode *node, gpointer data );
int load_ymal_params( GNode *node, gpointer data );

/* Functions in mads_io */
int parse_cmd( char *buf, struct calc_data *cd );
int load_problem( char *filename, int argn, char *argv[], struct opt_data *op );
int save_problem( char *filename, struct opt_data *op );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc2( char *filename, char *filename2, struct opt_data *op );
void compute_btc( char *filename, struct opt_data *op );
static char *strsave( const char *s, const char *lim );
char **shellpath( void );
void freeshellpath( char *shellpath[] );
unsigned maxpathlen( char *path[], const char *base );
void execvepath( char *path[], const char *base, char *const argv[], char *const envp[] );
int count_lines( char *filename );
int count_cols( char *filename, int row );
char *timestamp(); // create time stamp
char *datestamp(); // create date stamp
char *str_replace( char *orig, char *rep, char *with ); // replace all string occurrences

/* Functions elsewhere */
char **char_matrix( int maxCols, int maxRows );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
int set_test_problems( struct opt_data *op );
void *malloc_check( const char *what, size_t n );
int Ftest( char *filename );
FILE *Fread( char *filename );


int load_yaml_problem( char *filename, int argn, char *argv[], struct opt_data *op )
{
	FILE *infile;
	gop = op;
	GNode *gnode_data = g_node_new( filename );
	yaml_parser_t parser;
	if( ( infile = fopen( filename, "rb" ) ) == NULL )
	{
		tprintf( "File \'%s\' cannot be opened to read problem information!\n", filename );
		tprintf( "ERROR: Input file is needed!\n\n" );
		return( -1 );
	}
	else
	{
		yaml_parser_initialize( &parser );
		yaml_parser_set_input_file( &parser, infile );
		yaml_parse_layer( &parser, gnode_data ); // Recursive parsing into GNODE data
		yaml_parser_delete( &parser ); // Destroy YAML parser
		fclose( infile );
		tprintf( "YAML/GNODE Tree:\n" );
		g_node_traverse( gnode_data, G_PRE_ORDER, G_TRAVERSE_ALL, -1, gnode_tree_dump, NULL );
		g_node_traverse( gnode_data, G_PRE_ORDER, G_TRAVERSE_NON_LEAVES, -1, gnode_tree_dump, NULL );
		g_node_children_foreach( gnode_data, G_TRAVERSE_ALL, ( GNodeForeachFunc )gnode_tree_dump_classes, NULL );
		tprintf( "Tree Depth %d\n", g_node_depth( gnode_data ) );
		tprintf( "Tree Max Height %d\n", g_node_max_height( gnode_data ) );
		tprintf( "Tree Number of Classes %d\n", g_node_n_children( gnode_data ) );
		tprintf( "Tree Number of Data Sets %d\n", g_node_n_nodes( gnode_data, G_TRAVERSE_NON_LEAVES ) );
		g_node_children_foreach( gnode_data, G_TRAVERSE_ALL, ( GNodeForeachFunc )gnode_tree_parse_classes, NULL );
		g_node_destroy( gnode_data ); // Destroy GNODE data
	}
	return( 1 );
}

void yaml_parse_layer( yaml_parser_t *parser, GNode *data )
{
	GNode *last_leaf = data;
	yaml_event_t event;
	int storage = VAR; // mapping cannot start with VAL definition without VAR key
	while( 1 )
	{
		yaml_parser_parse( parser, &event );
		// Parse value either as a new leaf in the mapping or as a leaf value (one of them, in case it's a sequence)
		if( event.type == YAML_SCALAR_EVENT )
		{
			if( storage ) g_node_append_data( last_leaf, g_strdup( ( gchar * ) event.data.scalar.value ) ); // if sequence or val
			else last_leaf = g_node_append( data, g_node_new( g_strdup( ( gchar * ) event.data.scalar.value ) ) ); // if var
			storage ^= VAL; // Flip VAR/VAL switch for the next event
		}
		// Sequence - all the following scalars will be appended to the last_leaf
		else if( event.type == YAML_SEQUENCE_START_EVENT ) { storage = SEQ; }
		else if( event.type == YAML_SEQUENCE_END_EVENT ) {  storage = VAR; }
		// depth += 1
		else if( event.type == YAML_MAPPING_START_EVENT )
		{
			yaml_parse_layer( parser, last_leaf );
			storage ^= VAL; // Flip VAR/VAL, without touching SEQ
			// storage = VAR; // Var should be expected ...
		}
		// depth -= 1
		else if( event.type == YAML_MAPPING_END_EVENT || event.type == YAML_STREAM_END_EVENT ) { yaml_event_delete( &event ); break; } // Quit; yaml_event_delete(&event) is needed
		yaml_event_delete( &event );
	}
}

/*
void yaml_parse_layer(yaml_parser_t *parser, GNode *data)
{
	GNode *last_leaf = data;
	GNode *last_leaf2 = data;
	yaml_event_t event;
	int storage = VAR; // mapping cannot start with VAL definition without VAR key
	int event_type_last = YAML_STREAM_START_EVENT, mapping = 0;
	char buf[50];
	while (1)
	{
		yaml_parser_parse(parser, &event);
		// Parse value either as a new leaf in the mapping or as a leaf value (one of them, in case it's a sequence)
		if (event.type == YAML_SCALAR_EVENT) {
			if (storage) g_node_append_data(last_leaf, g_strdup((gchar*) event.data.scalar.value)); // if sequence or val
			else last_leaf = g_node_append(data, g_node_new(g_strdup((gchar*) event.data.scalar.value))); // if var
			storage ^= VAL; // Flip VAR/VAL switch for the next event
			event_type_last = event.type;
		}
		// Sequence - all the following scalars will be appended to the last_leaf
		else if (event.type == YAML_SEQUENCE_START_EVENT) { mapping = 0; event_type_last = event.type; storage = SEQ; }
		else if (event.type == YAML_SEQUENCE_END_EVENT) { event_type_last = event.type; storage = VAR; }
		// depth += 1
		else if (event.type == YAML_MAPPING_START_EVENT) {
			if( event_type_last == YAML_SEQUENCE_START_EVENT || event_type_last == YAML_MAPPING_END_EVENT )
			{
				sprintf( buf, "Map %i", mapping );
				last_leaf2 = g_node_append(last_leaf, g_node_new(g_strdup( (gchar*) buf ) ));
				mapping++;
			}
			yaml_parse_layer(parser, last_leaf2);
			event_type_last = YAML_MAPPING_END_EVENT;
			storage ^= VAL; // Flip VAR/VAL, without touching SEQ
			storage = VAR; // Var should be expected ...
		}
		// depth -= 1
		else if ( event.type == YAML_MAPPING_END_EVENT || event.type == YAML_STREAM_END_EVENT ) { yaml_event_delete(&event); break; } // Quit; yaml_event_delete(&event) is needed
		yaml_event_delete(&event);
	}
}
 */

gboolean gnode_tree_dump( GNode *node, gpointer data )
{
	int i = g_node_depth( node );
	while( --i ) printf( "    " );
	printf( "%s\n", ( char * ) node->data );
	return( 0 );
}

void gnode_tree_parse_classes( GNode *node, gpointer data )
{
	int found = 0;
	int i = g_node_depth( node );
	printf( "Process Class:" );
	while( --i ) printf( " " );
	printf( "%s ", ( char * ) node->data );
	if( strcasestr( ( char * ) node->data, "Wells" ) ) { tprintf( " ... process wells ..." ); load_ymal_wells( node, data ); found = 1; }
	if( strcasestr( ( char * ) node->data, "Sources" ) ) { tprintf( " ... process sources ..." ); found = 1; }
	if( strcasestr( ( char * ) node->data, "Parameters" ) ) { tprintf( " ... process parameters ..." ); load_ymal_params( node, data );  found = 1; }
	if( strcasestr( ( char * ) node->data, "Grid" ) ) { tprintf( " ... process parameters ..." ); load_ymal_grid( node, data ); found = 1; }
	if( found == 0 ) tprintf( " WARNING: ... this data set will be not processed ..." );
	tprintf( "\n" );
}

void gnode_tree_dump_classes( GNode *node, gpointer data )
{
	int i = g_node_depth( node );
	printf( "- Class:" );
	while( --i ) printf( " " );
	printf( "%s -> %i components\n", ( char * ) node->data, g_node_n_children( node ) );
}

int load_ymal_params( GNode *node, gpointer data )
{
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od, *preds;
	GNode *node_key, *node_value, *node_par;
	cd = gop->cd;
	od = gop->od;
	pd = gop->pd;
	preds = gop->preds;
	int i, j, k, m, bad_data = 0, include_predictions;
	cd->debug = 5;
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	pd->nParam = g_node_n_children( node );
	tprintf( "Number of parameters: %i\n", pd->nParam );
	pd->var_id = char_matrix( pd->nParam, 50 );
	pd->var = ( double * ) malloc( pd->nParam * sizeof( double ) );
	cd->var = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc( pd->nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc( pd->nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc( pd->nParam * sizeof( double ) );
	pd->param_expressions_index = ( int * ) malloc( pd->nParam * sizeof( int ) );
	pd->param_expressions = ( void ** ) malloc( pd->nParam * sizeof( void * ) );
	pd->nOptParam = pd->nFlgParam = 0;
	for( i = 0; i < pd->nParam; i++ )
	{
		node_par = g_node_nth_child( node, i );
		tprintf( "\nParameter %s\n", ( char * ) node_par->data );
		strcpy( pd->var_id[i], node_par->data );
		pd->var_opt[i] = 1; pd->var_log[i] = 0;
		for( k = 0; k < g_node_n_children( node_par ); k++ )  // Number of parameter arguments
		{
			node_key = g_node_nth_child( node_par, k );
			node_value = g_node_nth_child( node_key, 0 );
			if( cd->debug > 1 ) tprintf( "Key %s = %s\n", ( char * ) node_key->data, ( char * ) node_value->data );
			if( strcasestr( ( char * ) node_key->data, "longname" ) ) strcpy( pd->var_id[i], ( char * ) node_value->data );
			if( strcasestr( ( char * ) node_key->data, "log" ) ) if( strcasestr( ( char * ) node_value->data, "yes" ) || strcasestr( ( char * ) node_value->data, "1" ) ) pd->var_log[i] = 1;
			if( strcasestr( ( char * ) node_key->data, "type" ) )
			{
				if( strcasestr( ( char * ) node_value->data, "opt" ) || strcasestr( ( char * ) node_value->data, "1" ) ) pd->var_opt[i] = 1;
				else if( strcasestr( ( char * ) node_value->data, "flag" ) || strcasestr( ( char * ) node_value->data, "2" ) ) pd->var_opt[i] = 2;
			}
			if( strcasestr( ( char * ) node_key->data, "init" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var[i] );
			if( strcasestr( ( char * ) node_key->data, "step" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var_dx[i] );
			if( strcasestr( ( char * ) node_key->data, "min" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var_min[i] );
			if( strcasestr( ( char * ) node_key->data, "max" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var_max[i] );
		}
		if( cd->debug ) tprintf( "%-26s: ", pd->var_id[i] );
		cd->var[i] = pd->var[i];
		if( cd->debug ) tprintf( "init %9g opt %1d log %1d step %7g min %9g max %9g\n", pd->var[i], pd->var_opt[i], pd->var_log[i], pd->var_dx[i], pd->var_min[i], pd->var_max[i] );
		if( pd->var_opt[i] == 1 ) pd->nOptParam++;
		else if( pd->var_opt[i] == 2 )
		{
			pd->nFlgParam++;
			if( cd->calib_type != PPSD ) pd->nOptParam++;
		}
		if( pd->var_opt[i] >= 1 )
		{
			if( pd->var_max[i] < pd->var[i] || pd->var_min[i] > pd->var[i] )
			{
				tprintf( "ERROR: Parameter initial value is outside the specified min/max range! " );
				tprintf( "Parameter %s: %g min %g max %g\n", pd->var_id[i], pd->var[i], pd->var_min[i], pd->var_max[i] );
				bad_data = 1;
			}
			if( pd->var_max[i] < pd->var_min[i] )
			{
				tprintf( "ERROR: Parameter min/max range is not correctly specified! " );
				tprintf( "Parameter %s: min %g max %g\n", pd->var_id[i], pd->var_min[i], pd->var_max[i] );
				bad_data = 1;
			}
			if( cd->plogtrans == 1 ) pd->var_log[i] = 1;
			else if( cd->plogtrans == 0 ) pd->var_log[i] = 0;
			if( pd->var_log[i] == 1 )
			{
				if( pd->var_min[i] < 0 || pd->var[i] < 0 )
				{
					tprintf( "ERROR: Parameter cannot be log transformed (negative values)!\n" );
					tprintf( "Parameter %s: min %g max %g\n", pd->var_id[i], pd->var_min[i], pd->var_max[i] );
					if( cd->plogtrans ) { pd->var_log[i] = 0; pd->var_range[i] = pd->var_max[i] - pd->var_min[i]; continue; }
					else bad_data = 1;
				}
				double d;
				if( pd->var_dx[i] < 2 ) d = ( pd->var_max[i] - pd->var_min[i] ) / pd->var_dx[i];
				if( pd->var[i] < DBL_EPSILON ) pd->var[i] = DBL_EPSILON;
				if( pd->var_min[i] < DBL_EPSILON ) pd->var_min[i] = DBL_EPSILON;
				pd->var[i] = log10( pd->var[i] );
				pd->var_min[i] = log10( pd->var_min[i] );
				pd->var_max[i] = log10( pd->var_max[i] );
				if( pd->var_dx[i] < 2 ) pd->var_dx[i] = ( pd->var_max[i] - pd->var_min[i] ) / d;
				else pd->var_dx[i] = log10( pd->var_dx[i] );
			}
			pd->var_range[i] = pd->var_max[i] - pd->var_min[i];
			if( pd->var_dx[i] > DBL_EPSILON ) cd->pardx = 1; // discretization is ON
		}
	}
	return( 1 );
}

int load_ymal_wells( GNode *node, gpointer data )
{
	struct well_data *wd;
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od, *preds;
	GNode *node_key, *node_value, *node_well, *node_key2, *node_value2, *node_obs;
	wd = gop->wd;
	cd = gop->cd;
	od = gop->od;
	pd = gop->pd;
	preds = gop->preds;
	int i, j, k, m, bad_data = 0, include_predictions;
	cd->debug = 5;
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	wd->nW = g_node_n_children( node );
	tprintf( "Number of wells: %i\n", wd->nW );
	wd->id = char_matrix( wd->nW, 40 );
	wd->x = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->y = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->z1 = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->z2 = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->xa = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->ya = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->za1 = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->za2 = ( double * ) malloc( wd->nW * sizeof( double ) );
	wd->nWellObs = ( int * ) malloc( wd->nW * sizeof( int ) );
	wd->obs_target = ( double ** ) malloc( wd->nW * sizeof( double * ) );
	wd->obs_log = ( int ** ) malloc( wd->nW * sizeof( int * ) );
	wd->obs_time = ( double ** ) malloc( wd->nW * sizeof( double * ) );
	wd->obs_weight = ( double ** ) malloc( wd->nW * sizeof( double * ) );
	wd->obs_min = ( double ** ) malloc( wd->nW * sizeof( double * ) );
	wd->obs_max = ( double ** ) malloc( wd->nW * sizeof( double * ) );
	// od->nObs = preds->nTObs = 0;
	if( cd->debug ) tprintf( "\nObservation data:\n" );
	for( i = 0; i < wd->nW; i++ ) // Number of wells loop
	{
		node_well = g_node_nth_child( node, i );
		tprintf( "\nWell %s\n", ( char * ) node_well->data );
		strcpy( wd->id[i], node_well->data );
		for( k = 0; k < g_node_n_children( node_well ); k++ )  // Number of well parameters
		{
			node_key = g_node_nth_child( node_well, k );
			node_value = g_node_nth_child( node_key, 0 );
			if( cd->debug > 1 ) tprintf( "Key %s = %s\n", ( char * ) node_key->data, ( char * ) node_value->data );
			if( strcasestr( ( char * ) node_key->data, "x" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->x[i] );
			if( strcasestr( ( char * ) node_key->data, "y" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->y[i] );
			if( strcasestr( ( char * ) node_key->data, "z0" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->z1[i] );
			if( strcasestr( ( char * ) node_key->data, "z1" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->z2[i] );
			if( strcasestr( ( char * ) node_key->data, "obs" ) )
			{
				wd->nWellObs[i] = g_node_n_children( node_key );
				if( cd->debug ) tprintf( "Well %-6s x %8g y %8g z0 %6g z1 %6g nObs %2i ", wd->id[i], wd->x[i], wd->y[i], wd->z1[i], wd->z2[i], wd->nWellObs[i] );
				if( wd->nWellObs[i] > 0 )
				{
					wd->obs_target[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) );
					wd->obs_time[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) );
					wd->obs_log[i] = ( int * ) malloc( wd->nWellObs[i] * sizeof( int ) );
					wd->obs_weight[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) );
					wd->obs_min[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) );
					wd->obs_max[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) );
					for( j = 0; j < wd->nWellObs[i]; j++ )
					{
						node_obs = g_node_nth_child( node_key, j );
						wd->obs_min[i][j] = -1e6; wd->obs_max[i][j] = 1e6; wd->obs_weight[i][j] = 1; wd->obs_log[i][j] = 0;
						for( m = 0; m < g_node_n_children( node_obs ); m++ )  // Number of well parameters
						{
							node_key2 = g_node_nth_child( node_obs, m );
							node_value2 = g_node_nth_child( node_key2, 0 );
							if( cd->debug > 1 ) tprintf( "Key %s = %s\n", ( char * ) node_key2->data, ( char * ) node_value2->data );
							if( strcasestr( ( char * ) node_key2->data, "t" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_time[i][j] );
							if( strcasestr( ( char * ) node_key2->data, "c" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_target[i][j] );
							if( strcasestr( ( char * ) node_key2->data, "weight" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_weight[i][j] );
							if( strcasestr( ( char * ) node_key2->data, "log" ) ) sscanf( ( char * ) node_value2->data, "%i", &wd->obs_log[i][j] );
							if( strcasestr( ( char * ) node_key2->data, "min" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_min[i][j] );
							if( strcasestr( ( char * ) node_key2->data, "max" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_max[i][j] );
						}
						if( cd->obsdomain > DBL_EPSILON && wd->obs_weight[i][j] > DBL_EPSILON ) { wd->obs_min[i][j] = wd->obs_target[i][j] - cd->obsdomain; wd->obs_max[i][j] = wd->obs_target[i][j] + cd->obsdomain; }
						if( cd->ologtrans == 1 ) wd->obs_log[i][j] = 1;
						else if( cd->ologtrans == 0 ) wd->obs_log[i][j] = 0;
						if( cd->oweight == 1 ) wd->obs_weight[i][j] = 1;
						else if( cd->oweight == 0 ) wd->obs_weight[i][j] = 0;
						else if( cd->oweight == 2 ) { if( fabs( wd->obs_target[i][j] ) > DBL_EPSILON ) wd->obs_weight[i][j] = ( double ) 1.0 / wd->obs_target[i][j]; else wd->obs_weight[i][j] = HUGE_VAL; }
						if( cd->debug )
							tprintf( "t %5g c %5g weight %7g log %1d acceptable range: min %5g max %5g\n", wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_weight[i][j], wd->obs_log[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
						if( wd->obs_max[i][j] < wd->obs_target[i][j] || wd->obs_min[i][j] > wd->obs_target[i][j] )
						{
							tprintf( "ERROR: Observation target is outside the specified min/max range! " );
							tprintf( "Observation %s(%g): %g min %g max %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
							bad_data = 1;
						}
						if( wd->obs_max[i][j] <= wd->obs_min[i][j] )
						{
							tprintf( "ERROR: Calibration range is not correctly specified! " );
							tprintf( "Observation %s(%g): min %g max %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
							bad_data = 1;
						}
						if( wd->obs_weight[i][j] > DBL_EPSILON ) od->nObs++;
						if( wd->obs_weight[i][j] < -DBL_EPSILON ) { preds->nTObs++; if( include_predictions ) od->nObs++; } // Predictions have negative weights
						if( j + 1 < wd->nWellObs[i] ) if( cd->debug ) tprintf( "\t\t\t\t\t\t\t      " );
					}
				}
			}
		}
	}
	od->nCObs = od->nObs;
	for( i = 0; i < wd->nW; i++ )
	{
		if( wd->nWellObs[i] <= 0 )
			tprintf( "WARNING: Well %s has no observations!\n", wd->id[i] );
		for( j = 0; j < wd->nWellObs[i]; j++ )
			if( wd->obs_time[i][j] < DBL_EPSILON )
				tprintf( "WARNING: Observation #%d time for well %s is too small (%g); potential error in the input file %s!\n", j + 1, wd->id[i], wd->obs_time[i][j], gop->filename );
		for( j = i + 1; j < wd->nW; j++ )
			if( strcmp( wd->id[i], wd->id[j] ) == 0 )
				tprintf( "WARNING: Well names #%i (%s) and #%i (%s) are identical!\n", i + 1, wd->id[i], j + 1, wd->id[j] );
	}
	if( od->nObs == 0 )
	{
		if( cd->problem_type != FORWARD && cd->problem_type != MONTECARLO )
		{ tprintf( "\nERROR: Number of calibration targets is equal to zero!\n\n" ); bad_data = 1; }
		else tprintf( "\nWARNING: Number of calibration targets is equal to zero!\n\n" );
	}
	if( bad_data ) return( 0 );
	if( cd->debug ) tprintf( "\n" );
	return( 1 );
}

int load_ymal_grid( GNode *node, gpointer data )
{
	struct grid_data *gd;
	struct calc_data *cd;
	GNode *node_key, *node_value;
	gd = gop->gd;
	cd = gop->cd;
	int i;
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	for( i = 0; i < g_node_n_children( node ); i++ )
	{
		node_key = g_node_nth_child( node, i );
		node_value = g_node_nth_child( node_key, 0 );
		if( cd->debug > 1 ) tprintf( "Key %s = %s\n", ( char * ) node_key->data, ( char * ) node_value->data );
		if( strcasestr( ( char * ) node_key->data, "time" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->time );
		if( strcasestr( ( char * ) node_key->data, "xcount" ) ) sscanf( ( char * ) node_value->data, "%i", &gd->nx );
		if( strcasestr( ( char * ) node_key->data, "ycount" ) ) sscanf( ( char * ) node_value->data, "%i", &gd->ny );
		if( strcasestr( ( char * ) node_key->data, "zcount" ) ) sscanf( ( char * ) node_value->data, "%i", &gd->nz );
		if( strcasestr( ( char * ) node_key->data, "xmin" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_x );
		if( strcasestr( ( char * ) node_key->data, "ymin" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_y );
		if( strcasestr( ( char * ) node_key->data, "zmin" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_z );
		if( strcasestr( ( char * ) node_key->data, "xmax" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_x );
		if( strcasestr( ( char * ) node_key->data, "ymax" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_y );
		if( strcasestr( ( char * ) node_key->data, "zmax" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_z );
	}
	if( cd->debug )
	{
		tprintf( "\nGrid Time: %g\n", gd->time );
		tprintf( "Grid lines: %i %i %i\n", gd->nx, gd->ny, gd->nz );
		tprintf( "Grid Minimums: %g %g %g\n", gd->min_x, gd->min_y, gd->min_z );
		tprintf( "Grid Maximums: %g %g %g\n", gd->max_x, gd->max_y, gd->max_z );
	}
	if( gd->nx == 1 ) gd->dx = 0;
	else gd->dx = ( gd->max_x - gd->min_x ) / ( gd->nx - 1 );
	if( gd->ny == 1 ) gd->dy = 0;
	gd->dy = ( gd->max_y - gd->min_y ) / ( gd->ny - 1 );
	// if(gd->nz == 1 ) gd->dz = gd->max_z - gd->min_z ); // In this way compute_grid computed for min_z
	if( gd->nz == 1 ) gd->dz = 0;
	else gd->dz = ( gd->max_z - gd->min_z ) / ( gd->nz - 1 );
	return( 1 );
}
