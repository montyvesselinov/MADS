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

// this file is loaded only if YAML is defined
// #ifdef YAML is not needed

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "mads.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
// YAML includes
#include <yaml.h>
#include <glib.h>
#ifdef MATHEVAL
#include <matheval.h>
#endif

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
int load_ymal_sources( GNode *node, gpointer data );
int load_ymal_params( GNode *node, gpointer data, int num_keys, char **keywords, int *keyindex );
int load_ymal_observations( GNode *node, gpointer data );
int load_ymal_regularizations( GNode *node, gpointer data );
int load_ymal_time( GNode *node, gpointer data );
int load_ymal_problem( GNode *node, gpointer data );
int load_ymal_solution( GNode *node, gpointer data );
int load_ymal_command( GNode *node, gpointer data );
int load_ymal_templates( GNode *node, gpointer data );
int load_ymal_instructions( GNode *node, gpointer data );
static gboolean g_node_find_func( GNode *node, gpointer data );
gpointer g_node_find_key( GNode *gnode_data, char **key );
void set_param_arrays( int num_param, struct opt_data *op );

/* Functions in mads_io */
int set_param_id( struct opt_data *op );
int set_param_names( struct opt_data *op );
int parse_cmd( char *buf, struct calc_data *cd );
int set_optimized_params( struct opt_data *op );
int map_obs( struct opt_data *op );
int map_well_obs( struct opt_data *op );
char **shellpath( void );
void freeshellpath( char *shellpath[] );

/* Functions elsewhere */
int set_test_problems( struct opt_data *op );
char **char_matrix( int maxCols, int maxRows );

int load_yaml_problem( char *filename, int argn, char *argv[], struct opt_data *op )
{
	FILE *infile;
	struct calc_data *cd;
	struct obs_data *od;
	struct aquifer_data *qd;
	GNode *gnode_data = g_node_new( filename );
	yaml_parser_t parser;
	char buf[1000];
	char **keywords;
	int *keyindex;
	gpointer key_pointer;
	int i, k, ier, num_keys = 0;
	keyindex = NULL;
	keywords = NULL;
	ier = 1;
	qd = op->qd;
	cd = op->cd;
	od = op->od;
	op->gd->min_t = op->gd->time = 0;
	op->od->include_predictions = 1;
	if( op->cd->problem_type == INFOGAP ) op->od->include_predictions = 0;
	if( fabs( op->cd->obsstep ) > DBL_EPSILON ) op->od->include_predictions = 1;
	if( ( infile = fopen( filename, "rb" ) ) == NULL )
	{
		tprintf( "File \'%s\' cannot be opened to read problem information!\n", filename );
		tprintf( "ERROR: Input file is needed!\n\n" );
		return( -1 );
	}
	yaml_parser_initialize( &parser );
	yaml_parser_set_input_file( &parser, infile );
	yaml_parse_layer( &parser, gnode_data ); // Recursive parsing into GNODE data
	yaml_parser_delete( &parser ); // Destroy YAML parser
	fclose( infile );
	if( cd->debug > 5 )
	{
		tprintf( "YAML/GNODE Complete Tree:\n" );
		g_node_traverse( gnode_data, G_PRE_ORDER, G_TRAVERSE_ALL, -1, gnode_tree_dump, NULL );
		g_node_traverse( gnode_data, G_PRE_ORDER, G_TRAVERSE_NON_LEAVES, -1, gnode_tree_dump, NULL );
		tprintf( "Tree Depth %d\n", g_node_depth( gnode_data ) );
		tprintf( "Tree Max Height %d\n", g_node_max_height( gnode_data ) );
		tprintf( "Tree Number of Data Sets %d\n", g_node_n_nodes( gnode_data, G_TRAVERSE_NON_LEAVES ) );
	}
	tprintf( "Number of YAML classes %d\n", g_node_n_children( gnode_data ) );
	tprintf( "YAML Classes:\n" );
	g_node_children_foreach( gnode_data, G_TRAVERSE_ALL, ( GNodeForeachFunc )gnode_tree_dump_classes, NULL );
	tprintf( "Process YAML Classes ...\n" );
	// Problem
	key_pointer = g_node_find_key( gnode_data, ( char ** ) "Problem" );
	buf[0] = 0;
	if( key_pointer != NULL )
	{
		tprintf( "Process Problem ... " );
		GNode *node_key, *node_value;
		if( cd->debug > 2 )tprintf( "Number of components: %i ... Components:", g_node_n_children( key_pointer ) );
		for( i = 0; i < g_node_n_children( key_pointer ); i++ )
		{
			node_key = g_node_nth_child( key_pointer, i );
			if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
			strcat( buf, " " ); strcat( buf, ( char * ) node_key->data );
			if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
			{
				if( cd->debug > 1 ) tprintf( " = %s", ( char * ) node_value->data );
				strcat( buf, "=" ); strcat( buf, ( char * ) node_value->data );
			}
		}
		tprintf( "\n" );
	}
	else tprintf( "WARNING: Problem class not found!\n" );
	if( ier != 1 ) return( -1 );
	// Parse commands
	for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
	cd->solution_type = ( int * ) malloc( sizeof( int ) );
	if( parse_cmd( buf, cd ) == -1 ) return( -1 );
	od->include_predictions = 1;
	if( cd->problem_type == INFOGAP ) od->include_predictions = 0;
	if( fabs( cd->obsstep ) > DBL_EPSILON ) od->include_predictions = 1;
	// Solution
	key_pointer = g_node_find_key( gnode_data, ( char ** ) "Solution" );
	if( key_pointer != NULL ) { tprintf( "Process Solution ... " ); load_ymal_solution( ( GNode * ) key_pointer, ( void * ) op ); }
	else tprintf( "WARNING: Solution class not found!\n" );
	if( cd->solution_type[0] != EXTERNAL )
	{
		// Sources
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Sources" );
		if( key_pointer != NULL )
		{
			tprintf( "Process Sources ... " );
			ier = load_ymal_sources( key_pointer, ( void * ) op );
			if( cd->num_sources > 1 ) tprintf( "\nModels:" );
			else tprintf( "Model: " );
			for( i = 0; i < cd->num_sources; i++ )
			{
				if( cd->num_sources > 1 ) tprintf( " (%d) ", i + 1 );
				switch( cd->solution_type[i] )
				{
					case EXTERNAL: { tprintf( "external" ); break; }
					case POINT: { tprintf( "internal point contaminant source" ); break; }
					case PLANE: { tprintf( "internal rectangular contaminant source" ); break; }
					case GAUSSIAN2D: { tprintf( "internal planar (2d) gaussian contaminant source" ); break; }
					case GAUSSIAN3D: { tprintf( "internal spatial (3d) gaussian contaminant source" ); break; }
					case PLANE3D: { tprintf( "internal rectangular contaminant source with vertical flow component" ); break; }
					case BOX: { tprintf( "internal box contaminant source" ); break; }
					case TEST: { tprintf( "internal test optimization problem #%d: ", cd->test_func ); set_test_problems( op ); break; }
					default: tprintf( "WARNING! UNDEFINED model type!" ); break;
				}
			}
			if( cd->c_background ) tprintf( " | background concentration = %g", cd->c_background );
			tprintf( "\n" );
		}
		else tprintf( "WARNING: Sources class not found!\n" );
	}
	if( ier != 1 ) return( -1 );
	// Parameters
	key_pointer = g_node_find_key( gnode_data, ( char ** ) "Parameters" );
	if( key_pointer != NULL )
	{
		tprintf( "Process Parameters ... " );
		if( cd->solution_type[0] != EXTERNAL )
		{
			if( cd->num_sources > 0 )
			{
				num_keys = cd->num_aquifer_params;
				keyindex = ( int * ) malloc( cd->num_aquifer_params * sizeof( int ) );
				k = cd->num_sources * cd->num_source_params;
				for( i = 0; i < cd->num_aquifer_params; i++ ) keyindex[i] = k + i;
				keywords = qd->param_id;
			}
			ier = load_ymal_params( key_pointer, op, num_keys, keywords, keyindex );
			if( num_keys && keyindex != NULL ) free( keyindex );
		}
		else ier = load_ymal_params( key_pointer, op, 0, NULL, NULL );
	}
	else tprintf( "WARNING: Wells class not found!\n" );
	if( ier != 1 ) return( -1 );
	if( cd->solution_type[0] != EXTERNAL )
	{
		// Wells
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Wells" );
		if( key_pointer != NULL ) { tprintf( "Process Wells ... " ); ier = load_ymal_wells( ( GNode * ) key_pointer, ( void * ) op ); map_well_obs( op ); }
		else tprintf( "WARNING: Wells class not found!\n" );
	}
	else
	{
		// Observations
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Observations" );
		if( key_pointer != NULL ) { tprintf( "Process Observations ... " ); ier = load_ymal_observations( ( GNode * ) key_pointer, ( void * ) op ); }
		else tprintf( "WARNING: Observations class not found!\n" );
	}
	if( ier != 1 ) return( -1 );
	// Regularization
	key_pointer = g_node_find_key( gnode_data, ( char ** ) "Regularizations" );
	if( key_pointer != NULL ) { tprintf( "Process Regularizations ... " ); ier = load_ymal_regularizations( ( GNode * ) key_pointer, ( void * ) op ); }
	else tprintf( "WARNING: Regularization class not found!\n" );
	if( ier != 1 ) return( -1 );
	if( cd->solution_type[0] != EXTERNAL )
	{
		// Grid
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Grid" );
		if( key_pointer != NULL ) { tprintf( "Process Grid ... " ); ier = load_ymal_grid( ( GNode * ) key_pointer, ( void * ) op ); }
		else tprintf( "WARNING: Grid class not found!\n" );
		if( ier != 1 ) return( -1 );
		// Time
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Time" );
		if( key_pointer != NULL ) { tprintf( "Process Time ... " ); ier = load_ymal_time( ( GNode * ) key_pointer, ( void * ) op ); }
		else tprintf( "WARNING: Time class not found!\n" );
	}
	else
	{
		// Execution Command
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Command" );
		if( key_pointer != NULL ) { tprintf( "Process Execution Command ... " ); ier = load_ymal_command( ( GNode * ) key_pointer, ( void * ) op ); }
		else tprintf( "WARNING: Command class not found!\n" );
		if( ier != 1 ) return( -1 );
		// Templates
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Templates" );
		if( key_pointer != NULL ) { tprintf( "Process Templates and Output files ... " ); ier = load_ymal_templates( ( GNode * ) key_pointer, ( void * ) op ); }
		else tprintf( "WARNING: Templates class not found!\n" );
		if( ier != 1 ) return( -1 );
		// Instructions
		key_pointer = g_node_find_key( gnode_data, ( char ** ) "Instructions" );
		if( key_pointer != NULL ) { tprintf( "Process Instructions and Input files ... " ); ier = load_ymal_instructions( ( GNode * ) key_pointer, ( void * ) op ); }
		else tprintf( "WARNING: Instructions class not found!\n" );
		if( ier != 1 ) return( -1 );
	}
	// tprintf( "Process More YAML Classes ...\n" );
	// g_node_children_foreach( gnode_data, G_TRAVERSE_ALL, ( GNodeForeachFunc )gnode_tree_parse_classes, (void *) op );
	g_node_destroy( gnode_data ); // Destroy GNODE data
	if( !set_optimized_params( op ) ) return( -1 );
	if( op->rd->nRegul > 0 ) map_obs( op ); // add regularizations to the observations
	tprintf( "Number of regularization terms = %d\n", op->rd->nRegul );
	tprintf( "Number of predictions = %d\n", op->preds->nTObs );
	return( ier );
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

gboolean gnode_tree_dump( GNode *node, gpointer data )
{
	int i = g_node_depth( node );
	while( --i ) tprintf( "    " );
	tprintf( "%s\n", ( char * ) node->data );
	return( 0 );
}

void gnode_tree_dump_classes( GNode *node, gpointer data )
{
	int i = g_node_depth( node );
	while( --i ) tprintf( " " );
	tprintf( "%s -> %i components\n", ( char * ) node->data, g_node_n_children( node ) );
}

void gnode_tree_parse_classes( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * )data;
	struct calc_data *cd;
	struct aquifer_data *qd;
	int found = 0;
	int k, num_keys = 0, c, i = g_node_depth( node );
	char **keywords;
	int *keyindex;
	cd = op->cd;
	qd = op->qd;
	keyindex = NULL;
	keywords = NULL;
	tprintf( "Processing Class:" );
	while( --i ) tprintf( " " );
	tprintf( "%s ", ( char * ) node->data );
	if( !strcasecmp( ( char * ) node->data, "Problem" ) ) { tprintf( " ... process solution ... " ); load_ymal_problem( node, op ); found = 1; }
	if( !strcasecmp( ( char * ) node->data, "Solution" ) ) { tprintf( " ... process solution ... " ); load_ymal_solution( node, op ); found = 1; }
	if( !strcasecmp( ( char * ) node->data, "Sources" ) )
	{
		tprintf( " ... process sources ... " );
		load_ymal_sources( node, op );
		if( cd->num_sources > 1 ) tprintf( "\nModels:" );
		else tprintf( "Model: " );
		for( c = 0; c < cd->num_sources; c++ )
		{
			if( cd->num_sources > 1 ) tprintf( " (%d) ", c + 1 );
			switch( cd->solution_type[c] )
			{
				case EXTERNAL: { tprintf( "external" ); break; }
				case POINT: { tprintf( "internal point contaminant source" ); break; }
				case PLANE: { tprintf( "internal rectangular contaminant source" ); break; }
				case GAUSSIAN2D: { tprintf( "internal planar (2d) gaussian contaminant source" ); break; }
				case GAUSSIAN3D: { tprintf( "internal spatial (3d) gaussian contaminant source" ); break; }
				case PLANE3D: { tprintf( "internal rectangular contaminant source with vertical flow component" ); break; }
				case BOX: { tprintf( "internal box contaminant source" ); break; }
				case TEST: { tprintf( "internal test optimization problem #%d: ", cd->test_func ); set_test_problems( op ); break; }
				default: tprintf( "WARNING! UNDEFINED model type!" ); break;
			}
		}
		if( cd->c_background ) tprintf( " | background concentration = %g", cd->c_background );
		tprintf( "\n" );
		found = 1;
	}
	if( !strcasecmp( ( char * ) node->data, "Parameters" ) )
	{
		tprintf( " ... process parameters ... " );
		if( cd->num_sources > 0 || cd->solution_type[0] != EXTERNAL || cd->solution_type[0] != TEST )
		{
			num_keys = cd->num_aquifer_params;
			keyindex = ( int * ) malloc( cd->num_aquifer_params * sizeof( int ) );
			k = cd->num_sources * cd->num_source_params;
			for( i = 0; i < cd->num_aquifer_params; i++ ) keyindex[i] = k + i;
			keywords = qd->param_id;
		}
		load_ymal_params( node, op, num_keys, keywords, keyindex );
		if( num_keys && keyindex != NULL ) free( keyindex );
		found = 1;
	}
	if( !strcasecmp( ( char * ) node->data, "Wells" ) ) { tprintf( " ... process wells ... " ); load_ymal_wells( node, op ); found = 1; }
	if( !strcasecmp( ( char * ) node->data, "Grid" ) ) { tprintf( " ... process grid ... " ); load_ymal_grid( node, op ); found = 1; }
	if( !strcasecmp( ( char * ) node->data, "Time" ) ) { tprintf( " ... process breakthrough time data ... " ); load_ymal_time( node, op ); found = 1; }
	if( found == 0 ) tprintf( " WARNING: this data set will be not processed!\n" );
}

void set_param_arrays( int num_param, struct opt_data *op )
{
	struct calc_data *cd;
	struct param_data *pd;
	int i;
	pd = op->pd;
	cd = op->cd;
	pd->var_name = char_matrix( num_param, 50 );
	pd->var_id = char_matrix( num_param, 10 );
	for( i = 0; i < num_param; i++ )
		pd->var_name[i][0] = pd->var_id[i][0] = 0;
	pd->var = ( double * ) malloc( num_param * sizeof( double ) );
	cd->var = ( double * ) malloc( num_param * sizeof( double ) );
	pd->var_opt = ( int * ) malloc( num_param * sizeof( int ) );
	pd->var_log = ( int * ) malloc( num_param * sizeof( int ) );
	pd->var_dx = ( double * ) malloc( num_param * sizeof( double ) );
	pd->var_min = ( double * ) malloc( num_param * sizeof( double ) );
	pd->var_max = ( double * ) malloc( num_param * sizeof( double ) );
	pd->var_range = ( double * ) malloc( num_param * sizeof( double ) );
	pd->param_expressions_index = ( int * ) malloc( num_param * sizeof( int ) );
	pd->param_expression = ( void ** ) malloc( num_param * sizeof( void * ) );
	pd->nOptParam = pd->nFlgParam = pd->nExpParam = 0;
}

int load_ymal_sources( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct param_data *pd;
	struct source_data *sd;
	GNode *node_par;
	cd = op->cd;
	pd = op->pd;
	sd = op->sd;
	int i, k, bad_data = 0, *keyindex;
	if( cd->debug > 1 ) tprintf( "\nClass %s\n", ( char * ) node->data );
	cd->num_sources = g_node_n_children( node );
	tprintf( "Number of sources: %i\n", cd->num_sources );
	set_param_id( op ); // set analytical parameter id's
	pd->nParam = cd->num_sources * cd->num_source_params + cd->num_aquifer_params;
	set_param_arrays( pd->nParam, op );
	set_param_names( op );
	cd->solution_type = ( int * ) malloc( cd->num_sources * sizeof( int ) );
	keyindex = ( int * ) malloc( cd->num_source_params * sizeof( int ) );
	for( i = 0; i < cd->num_sources; i++ )
	{
		node_par = g_node_nth_child( node, i );
		tprintf( "Source #%d type: %s ... ", i + 1, ( char * ) node_par->data );
		if( !strncasecmp( ( char * ) node_par->data, "poi", 3 ) ) cd->solution_type[i] = POINT;
		if( !strncasecmp( ( char * ) node_par->data, "gau", 3 ) ) { if( strcasestr( ( char * ) node_par->data, "2" ) ) cd->solution_type[i] = GAUSSIAN2D; else cd->solution_type[i] = GAUSSIAN3D; }
		if( !strncasecmp( ( char * ) node_par->data, "rec", 3 ) ) { if( strcasestr( ( char * ) node_par->data, "ver" ) ) cd->solution_type[i] = PLANE3D; else cd->solution_type[i] = PLANE; }
		if( !strncasecmp( ( char * ) node_par->data, "box", 3 ) ) cd->solution_type[i] = BOX;
		for( k = 0; k < cd->num_source_params; k++ ) keyindex[k] = i * cd->num_source_params + k;
		node_par = g_node_nth_child( node, 0 );
		load_ymal_params( node_par, op, cd->num_source_params, sd->param_id, keyindex );
	}
	if( keyindex != NULL ) free( keyindex );
	if( bad_data ) return( 0 );
	else return( 1 );
}

int load_ymal_params( GNode *node, gpointer data, int num_keys, char **keywords, int *keyindex )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct param_data *pd;
	struct regul_data *rd;
	char **expvar_names;
	int expvar_count;
	GNode *node_key, *node_value, *node_par;
	cd = op->cd;
	pd = op->pd;
	rd = op->rd;
	rd->nRegul = 0;
	int i, index, k, num_param, internal, found, bad_data = 0;
	if( num_keys <= 0 ) internal = 0;
	else internal = 1;
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	num_param = g_node_n_children( node );
	tprintf( "Number of parameters: %i\n", num_param );
	if( !internal )
	{
		pd->nParam = num_param;
		set_param_arrays( pd->nParam, op );
	}
	else
	{
		if( num_param != num_keys )
			tprintf( "WARNING: The number of provided parameters (%d) is different than the number of expected parameters (%d)\n", num_param, num_keys );
		k = cd->num_source_params * ( cd->num_sources - 1 );
		pd->var[k + TSCALE_DISP] = 2; pd->var[k + TSCALE_ADV] = 0; pd->var[k + TSCALE_REACT] = 0;
		pd->var_opt[k + TSCALE_DISP] = 0; pd->var_opt[k + TSCALE_ADV] = 0; pd->var_opt[k + TSCALE_REACT] = 0;
		pd->var_log[k + TSCALE_DISP] = 0; pd->var_log[k + TSCALE_ADV] = 0; pd->var_log[k + TSCALE_REACT] = 0;
		pd->var_dx[k + TSCALE_DISP] = 0.1; pd->var_dx[k + TSCALE_ADV] = 0.1; pd->var_dx[k + TSCALE_REACT] = 0.1;
		pd->var_min[k + TSCALE_DISP] = 0; pd->var_min[k + TSCALE_ADV] = 0; pd->var_min[k + TSCALE_REACT] = 0;
		pd->var_max[k + TSCALE_DISP] = 10; pd->var_max[k + TSCALE_ADV] = 10; pd->var_max[k + TSCALE_REACT] = 10;
	}
	for( i = 0; i < num_param; i++ )
	{
		node_par = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Parameter %s", ( char * ) node_par->data );
		if( internal )
		{
			found = 0;
			if( cd->debug > 1 ) tprintf( " ..." );
			for( k = 0; k < num_keys; k++ )
				if( !strcasecmp( ( char * ) node_par->data, keywords[k] ) )
				{
					index = keyindex[k];
					strcpy( pd->var_id[index], ( char * ) node_par->data );
					if( cd->debug > 1 ) tprintf( " found -- indices %d -> %d\n", k, keyindex[k] );
					found = 1;
					break;
				}
			if( !found ) { tprintf( "WARNING: parameter name \'%s\' is not recognized; parameter values are ignored!\n", node_par->data ); continue; }
		}
		else
		{
			if( cd->debug > 1 ) tprintf( "\n" );
			index = i;
			strcpy( pd->var_id[index], node_par->data );
		}
		pd->var_min[index] = -HUGE_VAL; pd->var_max[index] = HUGE_VAL; pd->var[index] = 0; pd->var_opt[index] = 1; pd->var_log[index] = 0;
		for( k = 0; k < g_node_n_children( node_par ); k++ )  // Number of parameter arguments
		{
			node_key = g_node_nth_child( node_par, k );
			node_value = g_node_nth_child( node_key, 0 );
			if( cd->debug > 2 ) tprintf( "Key %s = %s\n", ( char * ) node_key->data, ( char * ) node_value->data );
			if( !strcasecmp( ( char * ) node_key->data, "longname" ) ) strcpy( pd->var_name[index], ( char * ) node_value->data );
			if( !strcasecmp( ( char * ) node_key->data, "log" ) ) if( !strcasecmp( ( char * ) node_value->data, "yes" ) || !strcasecmp( ( char * ) node_value->data, "1" ) ) pd->var_log[index] = 1;
			if( !strcasecmp( ( char * ) node_key->data, "type" ) )
			{
				if( !strcasecmp( ( char * ) node_value->data, "null" ) || !strcasecmp( ( char * ) node_value->data, "0" ) ) pd->var_opt[index] = 0;
				else if( !strcasecmp( ( char * ) node_value->data, "flag" ) || !strcasecmp( ( char * ) node_value->data, "2" ) ) pd->var_opt[index] = 2;
			}
			if( !strcasecmp( ( char * ) node_key->data, "init" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var[index] );
			if( !strcasecmp( ( char * ) node_key->data, "step" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var_dx[index] );
			if( !strcasecmp( ( char * ) node_key->data, "min" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var_min[index] );
			if( !strcasecmp( ( char * ) node_key->data, "max" ) ) sscanf( ( char * ) node_value->data, "%lf", &pd->var_max[index] );
			if( !strcasecmp( ( char * ) node_key->data, "exp" ) )
			{
				pd->var_opt[i] = pd->var_log[i] = 0;
#ifdef MATHEVAL
				pd->param_expressions_index[pd->nExpParam] = index;
				pd->param_expression[pd->nExpParam] = evaluator_create( ( char * ) node_value->data );
				assert( pd->param_expression[pd->nExpParam] );
				evaluator_get_variables( pd->param_expression[pd->nExpParam], &expvar_names, &expvar_count );
#else
				expvar_count = 0;
#endif
			}
		}
		cd->var[index] = pd->var[index];
		if( pd->var_name[index][0] == 0 ) strcpy( pd->var_name[index], pd->var_id[index] );
		if( cd->debug )
		{
			tprintf( "%-26s :%-6s: ", pd->var_name[index], pd->var_id[index] );
			if( pd->var_opt[index] > -1 )
				tprintf( "init %9g opt %1d log %1d step %7g min %9g max %9g\n", pd->var[index], pd->var_opt[index], pd->var_log[index], pd->var_dx[index], pd->var_min[index], pd->var_max[index] );
			else
			{
				if( expvar_count > 0 )
				{
					if( cd->debug )
					{
						tprintf( " -> variables:" );
						int j;
						for( j = 0; j < expvar_count; j++ )
							tprintf( " %s", expvar_names[j] );
						tprintf( "\n" );
					}
					pd->var_opt[index] = -1;
					pd->nExpParam++;
				}
				else
				{
#ifdef MATHEVAL
					pd->var[index] = cd->var[index] = evaluator_evaluate_x( pd->param_expression[pd->nExpParam], 0 );
					if( cd->debug ) tprintf( " = %g (NO variables; fixed parameter)\n", pd->var[index] );
#else
					tprintf( " MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
				}
			}
		}
		if( pd->var_opt[index] == 1 ) pd->nOptParam++;
		else if( pd->var_opt[index] == 2 )
		{
			pd->nFlgParam++;
			if( cd->calib_type != PPSD ) pd->nOptParam++;
		}
		if( pd->var_opt[index] >= 1 )
		{
			if( pd->var_max[index] < pd->var[index] || pd->var_min[index] > pd->var[index] )
			{
				tprintf( "ERROR: Parameter initial value is outside the specified min/max range! " );
				tprintf( "Parameter %s: %g min %g max %g\n", pd->var_name[index], pd->var[index], pd->var_min[index], pd->var_max[index] );
				bad_data = 1;
			}
			if( pd->var_max[index] < pd->var_min[index] )
			{
				tprintf( "ERROR: Parameter min/max range is not correctly specified! " );
				tprintf( "Parameter %s: min %g max %g\n", pd->var_name[index], pd->var_min[index], pd->var_max[index] );
				bad_data = 1;
			}
			if( cd->plogtrans == 1 ) pd->var_log[index] = 1;
			else if( cd->plogtrans == 0 ) pd->var_log[index] = 0;
			if( pd->var_log[index] == 1 )
			{
				if( pd->var_min[index] < 0 || pd->var[index] < 0 )
				{
					tprintf( "ERROR: Parameter cannot be log transformed (negative values)!\n" );
					tprintf( "Parameter %s: min %g max %g\n", pd->var_name[index], pd->var_min[index], pd->var_max[index] );
					if( cd->plogtrans ) { pd->var_log[index] = 0; pd->var_range[index] = pd->var_max[index] - pd->var_min[index]; continue; }
					else bad_data = 1;
				}
				double d = 0;
				if( pd->var_dx[index] < 2 ) d = ( pd->var_max[index] - pd->var_min[index] ) / pd->var_dx[index];
				if( pd->var[index] < DBL_EPSILON ) pd->var[index] = DBL_EPSILON;
				if( pd->var_min[index] < DBL_EPSILON ) pd->var_min[index] = DBL_EPSILON;
				pd->var[index] = log10( pd->var[index] );
				pd->var_min[index] = log10( pd->var_min[index] );
				pd->var_max[index] = log10( pd->var_max[index] );
				if( pd->var_dx[index] < 2 ) pd->var_dx[index] = ( pd->var_max[index] - pd->var_min[index] ) / d;
				else pd->var_dx[index] = log10( pd->var_dx[index] );
			}
			pd->var_range[index] = pd->var_max[index] - pd->var_min[index];
			if( pd->var_dx[index] > DBL_EPSILON ) cd->pardx = 1; // discretization is ON
		}
	}
	if( bad_data ) return( 0 );
	else return( 1 );
}

int load_ymal_regularizations( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct param_data *pd;
	struct regul_data *rd;
	struct obs_data *od;
	GNode *node_key, *node_value, *node_regul;
	char **expvar_names;
	int expvar_count;
	int i, k, bad_data = 0;
	cd = op->cd;
	pd = op->pd;
	rd = op->rd;
	od = op->od;
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	rd->nRegul = g_node_n_children( node );
	tprintf( "Number of regularization terms = %d\n", rd->nRegul );
	rd->regul_expression = ( void ** ) malloc( rd->nRegul * sizeof( void * ) );
	rd->regul_id = char_matrix( rd->nRegul, 10 );
	rd->regul_target = ( double * ) malloc( rd->nRegul * sizeof( double ) );
	rd->regul_weight = ( double * ) malloc( rd->nRegul * sizeof( double ) );
	rd->regul_min = ( double * ) malloc( rd->nRegul * sizeof( double ) );
	rd->regul_max = ( double * ) malloc( rd->nRegul * sizeof( double ) );
	rd->regul_log = ( int * ) malloc( rd->nRegul * sizeof( int ) );
#ifdef MATHEVAL
	rd->regul_nMap = pd->nParam + od->nObs;
	rd->regul_map_id = ( char ** ) malloc( ( rd->regul_nMap ) * sizeof( char * ) );
	rd->regul_map_val = ( double * ) malloc( ( rd->regul_nMap + rd->nRegul ) * sizeof( double ) );
#endif
#ifndef MATHEVAL
	tprintf( "WARNING: MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
	for( i = 0; i < rd->nRegul; i++ )
	{
		node_regul = g_node_nth_child( node, i );
		strcpy( rd->regul_id[i], node_regul->data );
		rd->regul_min[i] = -HUGE_VAL; rd->regul_max[i] = HUGE_VAL; rd->regul_weight[i] = 1; rd->regul_log[i] = 0;
		for( k = 0; k < g_node_n_children( node_regul ); k++ )  // Number of regulization components
		{
			node_key = g_node_nth_child( node_regul, k );
			if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
			if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
			{
				if( cd->debug > 1 ) tprintf( "=%s\n", ( char * ) node_value->data );
			}
			else { if( cd->debug > 1 ) tprintf( " WARNING: No data\n" ); continue; }
			if( !strcasecmp( ( char * ) node_key->data, "target" ) ) sscanf( ( char * ) node_value->data, "%lf", &rd->regul_target[i] );
			if( !strcasecmp( ( char * ) node_key->data, "weight" ) ) sscanf( ( char * ) node_value->data, "%lf", &rd->regul_weight[i] );
			if( !strcasecmp( ( char * ) node_key->data, "log" ) ) if( !strcasecmp( ( char * ) node_value->data, "yes" ) || !strcasecmp( ( char * ) node_value->data, "1" ) ) rd->regul_log[i] = 1;
			if( !strcasecmp( ( char * ) node_key->data, "max" ) ) sscanf( ( char * ) node_value->data, "%lf", &rd->regul_max[i] );
			if( !strcasecmp( ( char * ) node_key->data, "min" ) ) sscanf( ( char * ) node_value->data, "%lf", &rd->regul_min[i] );
			if( !strcasecmp( ( char * ) node_key->data, "equation" ) )
			{
#ifdef MATHEVAL
				rd->regul_expression[i] = evaluator_create( ( char * ) node_value->data );
				assert( rd->regul_expression[i] );
				evaluator_get_variables( rd->regul_expression[i], &expvar_names, &expvar_count );
#else
				expvar_count = 0;
#endif
			}
		}
		if( cd->debug ) tprintf( "%-12s: target %g weight %g log %i min %g max %g : equation %s", rd->regul_id[i], rd->regul_target[i], rd->regul_weight[i], rd->regul_log[i], rd->regul_min[i], rd->regul_max[i], evaluator_get_string( rd->regul_expression[i] ) );
		if( !( rd->regul_weight[i] > DBL_EPSILON ) ) tprintf( " WARNING Weight <= 0 " );
		if( expvar_count > 0 )
		{
			int j, l1, status;
			if( cd->debug )
			{
				tprintf( " -> variables:" );
				for( j = 0; j < expvar_count; j++ )
					tprintf( " %s", expvar_names[j] );
				tprintf( "\n" );
			}
			for( j = 0; j < expvar_count; j++ )
			{
				l1 = strlen( expvar_names[j] );
				status = 0;
				for( k = 0; k < pd->nParam; k++ )
					if( !strncmp( expvar_names[j], pd->var_id[k], l1 ) ) { status = 1; break; }
				for( k = 0; k < od->nObs; k++ )
					if( !strncmp( expvar_names[j], od->obs_id[k], l1 ) ) { status = 1; break; }
#ifdef MATHEVAL
				if( status == 0 ) { tprintf( "ERROR: parameter name \'%s\' in regularization term \'%s\' is not defined!\n", expvar_names[j], evaluator_get_string( rd->regul_expression[i] ) ); bad_data = 1; }
#endif
			}
		}
#ifdef MATHEVAL
		else { tprintf( "ERROR: no variables\n" ); bad_data = 1; }
#endif
	}
#ifdef MATHEVAL
	for( k = 0; k < pd->nParam; k++ ) { rd->regul_map_id[k] = pd->var_id[k]; rd->regul_map_val[k] = cd->var[k]; }
	for( i = pd->nParam, k = 0; k < od->nObs; k++, i++ ) { rd->regul_map_id[i] = od->obs_id[k]; rd->regul_map_val[i] = od->obs_current[k] = od->obs_target[k]; }
	free( cd->var );
	cd->var = &rd->regul_map_val[0];
	free( od->obs_current );
	od->obs_current = &rd->regul_map_val[pd->nParam];
	// for( k = 0; k < rd->regul_nMap; k++ ) { tprintf( "%s %g\n", rd->regul_map_id[k], rd->regul_map_val[k] ); }
#endif
	if( cd->debug )
	{
#ifdef MATHEVAL
		tprintf( "Regularization expressions evaluated (initial values applied for parameters; target values applied for observations).\n" );
#endif
		float result;
		for( i = 0; i < rd->nRegul; i++ )
		{
			tprintf( "%-12s= ", rd->regul_id[i] );
#ifdef MATHEVAL
			tprintf( "%s", evaluator_get_string( rd->regul_expression[i] ) );
			result = evaluator_evaluate( rd->regul_expression[i], rd->regul_nMap, rd->regul_map_id, rd->regul_map_val );
			tprintf( " = %g\n", result );
#else
			tprintf( "MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
		}
	}
	if( bad_data ) return ( 0 );
	else return( 1 );
}

int load_ymal_observations( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct obs_data *preds;
	struct obs_data *od;
	GNode *node_key, *node_value, *node_obs;
	int i, k, j, bad_data = 0;
	cd = op->cd;
	od = op->od;
	preds = op->preds;
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	od->nObs = g_node_n_children( node );
	tprintf( "Number of observations = %d\n", od->nObs );
	od->obs_id = char_matrix( od->nObs, 50 );
	od->obs_target = ( double * ) malloc( od->nObs * sizeof( double ) );
	od->obs_weight = ( double * ) malloc( od->nObs * sizeof( double ) );
	od->obs_min = ( double * ) malloc( od->nObs * sizeof( double ) );
	od->obs_max = ( double * ) malloc( od->nObs * sizeof( double ) );
	od->obs_log = ( int * ) malloc( od->nObs * sizeof( int ) );
	od->obs_current = ( double * ) malloc( od->nObs * sizeof( double ) );
	od->obs_best = ( double * ) malloc( od->nObs * sizeof( double ) );
	od->res = ( double * ) malloc( od->nObs * sizeof( double ) );
	preds->nTObs = 0; // TODO INFOGAP and GLUE analysis for external problems
	for( i = 0; i < od->nObs; i++ )
	{
		node_obs = g_node_nth_child( node, i );
		strcpy( od->obs_id[i], ( char * ) node_obs->data );
		if( cd->debug > 1 ) tprintf( "\nObservation %s", ( char * ) node_obs->data );
		if( ( node_value = g_node_nth_child( node_obs, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( " = %s\n", ( char * ) node_value->data );
			sscanf( ( char * ) node_value->data, "%lf", &od->obs_target[i] );
		}
		else tprintf( "\n" );
		od->obs_min[i] = -HUGE_VAL; od->obs_max[i] = HUGE_VAL; od->obs_weight[i] = 1; od->obs_log[i] = 0;
		for( k = 0; k < g_node_n_children( node_obs ); k++ )  // Number of regulization components
		{
			node_key = g_node_nth_child( node_obs, k );
			if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
			if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
			{
				if( cd->debug > 1 ) tprintf( "=%s\n", ( char * ) node_value->data );
			}
			else { if( cd->debug > 1 ) tprintf( " WARNING: No data\n" ); continue; }
			if( !strcasecmp( ( char * ) node_key->data, "target" ) ) sscanf( ( char * ) node_value->data, "%lf", &od->obs_target[i] );
			if( !strcasecmp( ( char * ) node_key->data, "weight" ) ) sscanf( ( char * ) node_value->data, "%lf", &od->obs_weight[i] );
			if( !strcasecmp( ( char * ) node_key->data, "log" ) ) if( !strcasecmp( ( char * ) node_value->data, "yes" ) || !strcasecmp( ( char * ) node_value->data, "1" ) ) od->obs_log[i] = 1;
			if( !strcasecmp( ( char * ) node_key->data, "max" ) ) sscanf( ( char * ) node_value->data, "%lf", &od->obs_max[i] );
			if( !strcasecmp( ( char * ) node_key->data, "min" ) ) sscanf( ( char * ) node_value->data, "%lf", &od->obs_min[i] );
		}
		if( cd->debug ) tprintf( "%-12s: target %g weight %g log %i min %g max %g", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
		if( cd->obsdomain > DBL_EPSILON && &od->obs_weight[i] > 0 ) { od->obs_min[i] = od->obs_target[i] - cd->obsdomain; od->obs_max[i] = od->obs_target[i] + cd->obsdomain; }
		if( od->obs_max[i] < od->obs_target[i] || od->obs_min[i] > od->obs_target[i] )
		{
			tprintf( "ERROR: Observation target is outside the specified min/max range! " );
			tprintf( "Observation %s: %g min %g max %g\n", od->obs_id[i], od->obs_target[i], od->obs_min[i], od->obs_max[i] );
			bad_data = 1;
		}
		if( od->obs_max[i] <= od->obs_min[i] )
		{
			tprintf( "ERROR: Calibration range is not correctly specified! " );
			tprintf( "Observation %s: min %g max %g\n", od->obs_id[i], od->obs_min[i], od->obs_max[i] );
			bad_data = 1;
		}
		if( cd->ologtrans == 1 ) od->obs_log[i] = 1;
		else if( cd->ologtrans == 0 ) od->obs_log[i] = 0;
		if( cd->oweight == 1 ) od->obs_weight[i] = 1;
		else if( cd->oweight == 0 ) od->obs_weight[i] = 0;
		else if( cd->oweight == 2 ) { if( fabs( od->obs_target[i] ) > DBL_EPSILON ) od->obs_weight[i] = ( double ) 1.0 / od->obs_target[i]; else od->obs_weight[i] = HUGE_VAL; }
		if( od->obs_weight[i] > DBL_EPSILON ) od->nCObs++;
		if( od->obs_weight[i] < -DBL_EPSILON ) { preds->nTObs++; if( od->include_predictions ) od->nCObs++; } // Predictions have negative weights
	}
	tprintf( "Number of calibration targets = %d\n", od->nCObs );
	tprintf( "Number of predictions = %d\n", preds->nTObs );
	if( bad_data ) return( 0 );
	if( cd->debug )
	{
		tprintf( "\n" );
		for( i = 0; i < od->nObs; i++ )
		{
			if( cd->debug > 10 || od->nObs <= 50 || ( i < 20 || i > od->nObs - 20 ) )
				tprintf( "%-20s: %15g weight %7g log %1d acceptable range: min %15g max %15g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
			if( ( !( cd->debug > 10 ) || od->nObs > 50 ) && i == 21 ) tprintf( "...\n" );
		}
		for( i = od->nObs; i < od->nTObs; i++ )
			tprintf( "%-20s: %15g weight %7g log %1d acceptable range: min %15g max %15g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
	}
	for( i = 0; i < od->nObs; i++ )
		for( j = i + 1; j < od->nObs; j++ )
			if( strcmp( od->obs_id[i], od->obs_id[j] ) == 0 )
			{
				tprintf( "ERROR: Observation names #%i (%s) and #%i (%s) are identical!\n", i + 1, od->obs_id[i], j + 1, od->obs_id[j] );
				bad_data = 1;
			}
	if( bad_data ) return ( 0 );
	else return( 1 );
}

int load_ymal_wells( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct well_data *wd;
	struct calc_data *cd;
	struct obs_data *od, *preds;
	GNode *node_key, *node_value, *node_well, *node_key2, *node_value2, *node_obs;
	wd = op->wd;
	cd = op->cd;
	od = op->od;
	preds = op->preds;
	int i, j, k, m, bad_data = 0;
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
	od->nObs = preds->nTObs = 0;
	for( i = 0; i < wd->nW; i++ ) // Number of wells loop
	{
		node_well = g_node_nth_child( node, i );
		strcpy( wd->id[i], node_well->data );
		for( k = 0; k < g_node_n_children( node_well ); k++ )  // Number of well parameters
		{
			node_key = g_node_nth_child( node_well, k );
			if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
			if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
			{
				if( cd->debug > 1 ) tprintf( "=%s\n", ( char * ) node_value->data );
			}
			else { if( cd->debug > 1 ) tprintf( " WARNING: No data\n" ); continue; }
			if( !strcasecmp( ( char * ) node_key->data, "x" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->x[i] );
			if( !strcasecmp( ( char * ) node_key->data, "y" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->y[i] );
			if( !strcasecmp( ( char * ) node_key->data, "z0" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->z1[i] );
			if( !strcasecmp( ( char * ) node_key->data, "z1" ) ) sscanf( ( char * ) node_value->data, "%lf", &wd->z2[i] );
			if( !strcasecmp( ( char * ) node_key->data, "obs" ) )
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
							if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key2->data );
							if( ( node_value2 = g_node_nth_child( node_key2, 0 ) ) != NULL )
							{
								if( cd->debug > 1 ) tprintf( "=%s\n", ( char * ) node_value2->data );
							}
							else { if( cd->debug > 1 ) tprintf( " WARNING: No data\n" ); continue; }
							if( !strcasecmp( ( char * ) node_key2->data, "t" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_time[i][j] );
							if( !strcasecmp( ( char * ) node_key2->data, "c" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_target[i][j] );
							if( !strcasecmp( ( char * ) node_key2->data, "weight" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_weight[i][j] );
							if( !strcasecmp( ( char * ) node_key2->data, "log" ) ) if( !strcasecmp( ( char * ) node_value2->data, "yes" ) || !strcasecmp( ( char * ) node_value2->data, "1" ) ) wd->obs_log[i][j] = 1;
							if( !strcasecmp( ( char * ) node_key2->data, "min" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_min[i][j] );
							if( !strcasecmp( ( char * ) node_key2->data, "max" ) ) sscanf( ( char * ) node_value2->data, "%lf", &wd->obs_max[i][j] );
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
						if( wd->obs_weight[i][j] < -DBL_EPSILON ) { preds->nTObs++; if( od->include_predictions ) od->nObs++; } // Predictions have negative weights
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
				tprintf( "WARNING: Observation #%d time for well %s is too small (%g); potential error in the input file %s!\n", j + 1, wd->id[i], wd->obs_time[i][j], op->filename );
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
	struct opt_data *op = ( struct opt_data * ) data;
	struct grid_data *gd;
	struct calc_data *cd;
	GNode *node_key, *node_value;
	gd = op->gd;
	cd = op->cd;
	int i;
	tprintf( "Number of grid components: %i\n", g_node_n_children( node ) );
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	for( i = 0; i < g_node_n_children( node ); i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( "=%s\n", ( char * ) node_value->data );
		}
		else { if( cd->debug > 1 ) tprintf( " WARNING: No data\n" ); continue; }
		if( !strcasecmp( ( char * ) node_key->data, "time" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->time );
		if( !strcasecmp( ( char * ) node_key->data, "xcount" ) ) sscanf( ( char * ) node_value->data, "%i", &gd->nx );
		if( !strcasecmp( ( char * ) node_key->data, "ycount" ) ) sscanf( ( char * ) node_value->data, "%i", &gd->ny );
		if( !strcasecmp( ( char * ) node_key->data, "zcount" ) ) sscanf( ( char * ) node_value->data, "%i", &gd->nz );
		if( !strcasecmp( ( char * ) node_key->data, "xmin" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_x );
		if( !strcasecmp( ( char * ) node_key->data, "ymin" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_y );
		if( !strcasecmp( ( char * ) node_key->data, "zmin" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_z );
		if( !strcasecmp( ( char * ) node_key->data, "xmax" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_x );
		if( !strcasecmp( ( char * ) node_key->data, "ymax" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_y );
		if( !strcasecmp( ( char * ) node_key->data, "zmax" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_z );
	}
	if( cd->debug )
	{
		tprintf( "Grid Time: %g\n", gd->time );
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

int load_ymal_time( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct grid_data *gd;
	struct calc_data *cd;
	GNode *node_key, *node_value;
	gd = op->gd;
	cd = op->cd;
	int i;
	tprintf( "Number of time components: %i\n", g_node_n_children( node ) );
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	for( i = 0; i < g_node_n_children( node ); i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( "=%s\n", ( char * ) node_value->data );
		}
		else { if( cd->debug > 1 ) tprintf( " WARNING: No data\n" ); continue; }
		if( !strcasecmp( ( char * ) node_key->data, "start" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->min_t );
		if( !strcasecmp( ( char * ) node_key->data, "end" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->max_t );
		if( !strcasecmp( ( char * ) node_key->data, "step" ) ) sscanf( ( char * ) node_value->data, "%lf", &gd->dt );
	}
	if( cd->debug ) tprintf( "Breakthrough-curve time window: start %g end %g step %g\n", gd->min_t, gd->max_t, gd->dt );
	return( 1 );
}

int load_ymal_problem( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	GNode *node_key, *node_value;
	cd = op->cd;
	int i;
	tprintf( "Number of problem components: %i\nComponents: ", g_node_n_children( node ) );
	if( cd->debug > 1 ) tprintf( "\n%s\n", ( char * ) node->data );
	for( i = 0; i < g_node_n_children( node ); i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( "=%s", ( char * ) node_value->data );
		}
		tprintf( "\n" );
	}
	return( 1 );
}

int load_ymal_solution( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	GNode *node_key, *node_value;
	cd = op->cd;
	int i;
	tprintf( "Number of solution components: %i\n", g_node_n_children( node ) );
	if( cd->debug > 1 ) tprintf( "Components:\n%s\n", ( char * ) node->data );
	for( i = 0; i < g_node_n_children( node ); i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		cd->solution_type[0] = EXTERNAL;
		if( !strncasecmp( ( char * ) node_key->data, "internal", 8 ) ) cd->solution_type[0] = POINT;
		if( !strncasecmp( ( char * ) node_key->data, "external", 8 ) ) cd->solution_type[0] = EXTERNAL;
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( "=%s", ( char * ) node_value->data );
		}
		if( cd->debug > 1 ) tprintf( " [%d]", cd->solution_type[0] );
		tprintf( "\n" );
	}
	return( 1 );
}

int load_ymal_command( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct extrn_data *ed;
	GNode *node_key, *node_value;
	cd = op->cd;
	ed = op->ed;
	int i, k, bad_data = 0;
	char buf[1000], *file, **path, exec[1000];
	tprintf( "Number of command components: %i\n", g_node_n_children( node ) );
	if( cd->debug > 1 ) tprintf( "Components:\n%s\n", ( char * ) node->data );
	for( i = 0; i < g_node_n_children( node ); i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( "=%s", ( char * ) node_value->data );
		}
		tprintf( "\n" );
		ed->cmdline = ( char * ) malloc( 80 * sizeof( char ) );
		strcpy( ed->cmdline, node_key->data );
		if( cd->debug > 1 ) tprintf( "Execution command: %s\n", ed->cmdline );
		buf[0] = 0;
		file = &buf[0];
		strcpy( file, ed->cmdline );
		file = strsep( &file, " \t" );
		k = 1;
		if( access( file, X_OK ) == -1 )
		{
			k = 0;
			if( file[0] == '~' )
			{
				sprintf( exec, "%s/%s", getenv( "HOME" ), &file[1] );
				if( access( exec, X_OK ) == 0 )
					k = 1;
			}
			if( k == 0 )
			{
				path = shellpath();
				for( i = 0; k == 0 && path[i]; i++ )
				{
					if( cd->debug > 2 ) tprintf( "%s\n", path[i] );
					sprintf( exec, "%s/%s", path[i], file );
					if( access( exec, X_OK ) == 0 )
						k = 1;
				}
			}
		}
		if( k == 0 )
		{
			tprintf( "ERROR: Program \'%s\' does not exist or cannot be executed!\n", file );
			bad_data = 1;
		}
	}
	if( bad_data ) return( 0 );
	return( 1 );
}

int load_ymal_templates( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct extrn_data *ed;
	GNode *node_key, *node_key2, *node_value;
	cd = op->cd;
	ed = op->ed;
	int i, j, bad_data = 0;
	ed->ntpl = g_node_n_children( node );
	tprintf( "Number of template components: %i\n", ed->ntpl );
	ed->fn_tpl = char_matrix( ed->ntpl, 80 );
	ed->fn_out = char_matrix( ed->ntpl, 80 );
	if( cd->debug > 1 ) tprintf( "Components:\n%s\n", ( char * ) node->data );
	for( i = 0; i < ed->ntpl; i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( " = %s", ( char * ) node_value->data );
		}
		else { if( cd->debug > 1 ) { tprintf( " WARNING: No data\n" ); bad_data = 1; } else tprintf( "\n" ); continue; }
		for( j = 0; j < g_node_n_children( node_key ); j++ )
		{
			node_key2 = g_node_nth_child( node_key, j );
			if( cd->debug > 1 ) tprintf( "\nKey %s", ( char * ) node_key2->data );
			if( ( node_value = g_node_nth_child( node_key2, 0 ) ) != NULL )
			{
				if( cd->debug > 1 ) tprintf( " = %s", ( char * ) node_value->data );
			}
			else { if( cd->debug > 1 ) { tprintf( " WARNING: No data\n" ); bad_data = 1; } else tprintf( "\n" ); continue; }
			if( !strncasecmp( ( char * ) node_key2->data, "tpl", 3 ) ) strcpy( ed->fn_tpl[i], ( char * )  node_value->data );
			if( !strncasecmp( ( char * ) node_key2->data, "write", 5 ) ) strcpy( ed->fn_out[i], ( char * ) node_value->data );
		}
	}
	if( cd->debug > 1 )
	{
		tprintf( "\nExternal files to provide current model parameters:\n" );
		for( i = 0; i < ed->ntpl; i++ )
			tprintf( "%s -> %s\n", ed->fn_tpl[i], ed->fn_out[i] );
	}
	if( bad_data ) return( 0 );
	return( 1 );
}

int load_ymal_instructions( GNode *node, gpointer data )
{
	struct opt_data *op = ( struct opt_data * ) data;
	struct calc_data *cd;
	struct extrn_data *ed;
	GNode *node_key, *node_key2, *node_value;
	cd = op->cd;
	ed = op->ed;
	int i, j, bad_data = 0;
	ed->nins = g_node_n_children( node );
	tprintf( "Number of instruction components: %i\n", ed->nins );
	ed->fn_ins = char_matrix( ed->nins, 80 );
	ed->fn_obs = char_matrix( ed->nins, 80 );
	if( cd->debug > 1 ) tprintf( "Components:\n%s\n", ( char * ) node->data );
	for( i = 0; i < ed->nins; i++ )
	{
		node_key = g_node_nth_child( node, i );
		if( cd->debug > 1 ) tprintf( "Key %s", ( char * ) node_key->data );
		if( ( node_value = g_node_nth_child( node_key, 0 ) ) != NULL )
		{
			if( cd->debug > 1 ) tprintf( " = %s", ( char * ) node_value->data );
		}
		else { if( cd->debug > 1 ) { tprintf( " WARNING: No data\n" ); bad_data = 1; } else tprintf( "\n" ); continue; }
		for( j = 0; j < g_node_n_children( node_key ); j++ )
		{
			node_key2 = g_node_nth_child( node_key, j );
			if( cd->debug > 1 ) tprintf( "\nKey %s", ( char * ) node_key2->data );
			if( ( node_value = g_node_nth_child( node_key2, 0 ) ) != NULL )
			{
				if( cd->debug > 1 ) tprintf( " = %s", ( char * ) node_value->data );
			}
			else { if( cd->debug > 1 ) { tprintf( " WARNING: No data\n" ); bad_data = 1; } else tprintf( "\n" ); continue; }
			if( !strncasecmp( ( char * ) node_key2->data, "ins", 3 ) ) strcpy( ed->fn_ins[i], ( char * )  node_value->data );
			if( !strncasecmp( ( char * ) node_key2->data, "read", 5 ) ) strcpy( ed->fn_obs[i], ( char * ) node_value->data );
		}
	}
	if( cd->debug > 1 )
	{
		tprintf( "\nExternal files to read current model predictions:\n" );
		for( i = 0; i < ed->nins; i++ )
			tprintf( "%s <- %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	}
	if( bad_data ) return( 0 );
	return( 1 );
}

gpointer g_node_find_key( GNode *gnode_data, char **key )
{
	gpointer search_pointer[2];
	search_pointer[0] = ( char ** ) key;
	search_pointer[1] = NULL;
	g_node_traverse( gnode_data, G_LEVEL_ORDER, G_TRAVERSE_MASK, -1, g_node_find_func, search_pointer );
	return( search_pointer[1] );
}

static gboolean g_node_find_func( GNode *node, gpointer data )
{
	register gpointer *d = data;
	// tprintf( "Search %s %s\n", ( char * ) node->data, (gchar *)d[0] );
	if( !node->data || strcmp( ( gchar * )d[0], ( gchar * )node->data ) ) return FALSE;
	// tprintf( "Found %s ", ( char * ) node->data );
	d[1] = node;
	return TRUE;
}
