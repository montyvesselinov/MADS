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
// XML
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <glib.h>

/* Functions here */
int load_xml_problem( char *filename, int argn, char *argv[], struct opt_data *op );
static void print_xml_elements( xmlNode *xml_node);
static void parse_xml( xmlNodePtr xml_node, GNode *gdata, int debug );

/* Functions in mads_io_yaml */
gboolean gnode_tree_dump( GNode *node, gpointer data );
void gnode_tree_dump_classes( GNode *node, gpointer data );
int parse_gnode_classes( GNode *node,  int argn, char *argv[], struct opt_data *op );

int load_xml_problem( char *filename, int argn, char *argv[], struct opt_data *op )
{
	FILE *infile;
	xmlDocPtr xmldoc; // pointer to XML document
	xmlNodePtr xmlcur; // pointer to current node within document
	struct calc_data *cd;
	GNode *gnode_data;
	int ier = 1;
	cd = op->cd;
	op->preds->nTObs = 0;
	op->gd->min_t = op->gd->time = 0;
	op->od->include_predictions = 1;
	if( op->cd->problem_type == INFOGAP ) op->od->include_predictions = 0;
	if( fabs( op->cd->obsstep ) > DBL_EPSILON ) op->od->include_predictions = 1;
	if( ( infile = fopen( filename, "rb" ) ) == NULL )
	{
		tprintf( "ERROR: XML file \'%s\' cannot be opened to read problem information!\n", filename );
		mads_quits( 0 );
	}
	fclose( infile );
	xmldoc = xmlParseFile( filename ); // attempt to parse xml file and get pointer to parsed document
	if( xmldoc == NULL )
	{
		tprintf( "ERROR: XML file \'%s\' is empty.\n", filename );
		mads_quits( 0 );
	}
	xmlcur = xmlDocGetRootElement( xmldoc ); // Start at the root of the XML document
	if( xmlcur == NULL )
	{
		tprintf( "ERROR: XML file \'%s\' cannot be parsed.\n", filename );
		xmlFreeDoc( xmldoc );
		mads_quits( 0 );
	}
	print_xml_elements( xmlcur );
	if( xmlStrcmp( xmlcur->name, ( xmlChar * ) "mads" ) )
	{
		tprintf( "ERROR: XML file \'%s\' is of a wrong type (root node should be named \'mads\')\n", filename );
		xmlFreeDoc( xmldoc );
		return 0;
	}
	gnode_data = g_node_new( filename );
	if( !gnode_data ) { tprintf( "ERROR: GNode array cannot be created.\n" ); mads_quits( 0 ); }
	parse_xml( xmlcur, gnode_data, ( int )( op->cd->debug > 5 ) );  // Recursive parsing into GNODE data
	xmlFreeDoc( xmldoc ); // Destroy XML
	if( cd->debug > 5 )
	{
		tprintf( "XML/GNODE Complete Tree:\n" );
		g_node_traverse( gnode_data, G_PRE_ORDER, G_TRAVERSE_ALL, -1, gnode_tree_dump, NULL );
		g_node_traverse( gnode_data, G_PRE_ORDER, G_TRAVERSE_NON_LEAVES, -1, gnode_tree_dump, NULL );
		tprintf( "Tree Depth %d\n", g_node_depth( gnode_data ) );
		tprintf( "Tree Max Height %d\n", g_node_max_height( gnode_data ) );
		tprintf( "Tree Number of Data Sets %d\n", g_node_n_nodes( gnode_data, G_TRAVERSE_NON_LEAVES ) );
	}
	tprintf( "Number of XML classes %d\n", g_node_n_children( gnode_data ) );
	tprintf( "XML Classes:\n" );
	g_node_children_foreach( gnode_data, G_TRAVERSE_ALL, ( GNodeForeachFunc )gnode_tree_dump_classes, NULL );
	tprintf( "Process XML Classes ...\n" );
	ier = parse_gnode_classes( gnode_data, argn, argv, op );
	return( ier );
}

/*
GNode *parse_templates_xml( const char *filename )
{
	xmlDocPtr doc; // pointer to XML document
	xmlNodePtr cur; // pointer to current node within document
	GNode *gnode_data; // Return value, all templates stored as children here.
	GNode *tnode_prev = 0;
	doc = xmlParseFile( filename ); // attempt to parse xml file and get pointer to parsed document
	if( doc == NULL )
	{
		tprintf( "ERROR: XML file did not parsed successfully.\n" );
		mads_quits( 0 );
	}

	cur = xmlDocGetRootElement( doc ); // Start at the root of the XML document
	if( cur == NULL )
	{
		tprintf( "ERROR: Empty XML File.\n" );
		xmlFreeDoc( doc );
		mads_quits( 0 );
	}
	print_element_names( cur );

	while( cur != NULL )
	{
		if( !ignore_xml_node( cur ) )
		{
			if( ( !xmlStrcasecmp( cur->name, ( xmlChar * ) "mads" ) ) )
			{
				GNode *tnode;
				tnode = parseTemplate( cur );
				 Parse template, and check if parse failed.
				if( tnode )
				{
					g_node_insert_after( gnode_data, tnode_prev, tnode );
					tnode_prev = tnode;
				}
				else
				{
					tprintf( "Parsing a template failed." );
					xmlFreeDoc( doc );
					return 0;
				}
			}
			else
			{
				fprintf( stderr, "Warning: Unknown templates child: %s\n", cur->name );
				fprintf( stderr, "Continuing to parse templates file ...\n" );
			}
		}
		cur = cur->next;
	}

	xmlFreeDoc( doc );
	return gnode_data;
}
*/

static void print_xml_elements( xmlNode *xml_node )
{
	xmlNode *cur_node = NULL;
	xmlChar *content;
	int count;

	for( cur_node = xml_node; cur_node; cur_node = cur_node->next )
	{
		if( cur_node->type == XML_ELEMENT_NODE )
		{
			tprintf( "Node Name: %s", cur_node->name );
			content = xmlNodeGetContent( cur_node );
			count = xmlChildElementCount( cur_node );
			if( content != NULL && count == 0 ) tprintf( " Content: %s\n", content );
			else { tprintf( "\n" ); print_xml_elements(cur_node->children); }
		}
	}
}

static void parse_xml( xmlNodePtr xml_node, GNode *gdata, int debug )
{
	GNode *last_leaf = gdata;
	xmlNode *cur_node = NULL;
	xmlChar *content;
	int count;

	for( cur_node = (xmlNode *) xml_node; cur_node; cur_node = cur_node->next )
	{
		if( cur_node->type == XML_ELEMENT_NODE )
		{
			last_leaf = g_node_append( gdata, g_node_new( g_strdup( ( gchar * ) cur_node->name ) ) );
			if( debug ) tprintf( "Node Name: %s", cur_node->name );
			content = xmlNodeGetContent( cur_node );
			count = xmlChildElementCount( cur_node );
			if( content != NULL && count == 0 )
			{
				g_node_append_data( last_leaf, g_strdup( ( gchar * ) content ) );
				if( debug ) tprintf( " Content: %s\n", content );
			}
			else { tprintf( "\n" ); parse_xml( cur_node->children, last_leaf, debug ); }
		}
	}
}
