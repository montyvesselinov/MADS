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
#include "xml-template.h"
#include <glib.h>
#ifdef MATHEVAL
#include <matheval.h>
#endif

#define C2P(c)          ((gpointer) ((long) (c)))
#define P2C(p)          ((gchar) ((long) (p)))
enum storage_flags { VAR, VAL, SEQ }; // "Store as" switch

/* Functions here */
int load_xml_problem( char *filename, int argn, char *argv[], struct opt_data *op );
// void xml_parse_layer( yaml_parser_t *parser, GNode *data, int debug );
GNode *parse_templates_xml( const char *filename );

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
		tprintf( "File \'%s\' cannot be opened to read problem information!\n", filename );
		tprintf( "ERROR: Input file is needed!\n\n" );
		return( -1 );
	}
	fclose( infile );
	gnode_data = parse_templates_xml( filename );
	/*
		xmldoc = xmlParseFile( filename ); // attempt to parse xml file and get pointer to parsed document
		if( xmldoc == NULL )
		{
			tprintf( "File \'%s\' cannot be opened to read problem information!\n", filename );
			tprintf( "ERROR: XML document not parsed successfully. \n" );
			return( -1 );
		}
		xmlcur = xmlDocGetRootElement( xmldoc ); // Start at the root of the XML document
		if( xmlcur == NULL )
		{
			tprintf( "File \'%s\' cannot be parsed!\n", filename );
			tprintf( "ERROR: XML document not parsed successfully. \n" );
			xmlFreeDoc( xmldoc );
			return( -1 );
		}
		xml_parse_layer( &xmldoc, gnode_data, ( int )( op->cd->debug > 5 ) );  // Recursive parsing into GNODE data
		xmlFreeDoc( xmldoc ); // Destroy XML
	*/
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
void xml_parse_layer( xmlDocPtr *xmldoc, GNode *data, int debug )
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
			if( storage ) { g_node_append_data( last_leaf, g_strdup( ( gchar * ) event.data.scalar.value ) ); if( debug ) tprintf( "Data: %s\n", event.data.scalar.value ); } // if sequence or val
			else { last_leaf = g_node_append( data, g_node_new( g_strdup( ( gchar * ) event.data.scalar.value ) ) ); if( debug ) tprintf( "Key : %s\n", event.data.scalar.value ); } // if var
			storage ^= VAL; // Flip VAR/VAL switch for the next event
		}
		// Sequence - all the following scalars will be appended to the last_leaf
		else if( event.type == YAML_SEQUENCE_START_EVENT ) { storage = SEQ; if( debug ) tprintf( "Sequence start\n" ); }
		else if( event.type == YAML_SEQUENCE_END_EVENT ) {  storage = VAR; if( debug ) tprintf( "Sequence end\n" ); }
		// depth += 1
		else if( event.type == YAML_MAPPING_START_EVENT )
		{
			yaml_parse_layer( parser, last_leaf, debug );
			storage ^= VAL; // Flip VAR/VAL, without touching SEQ
			// storage = VAR; // Var should be expected ...
		}
		// depth -= 1
		else if( event.type == YAML_MAPPING_END_EVENT || event.type == YAML_STREAM_END_EVENT ) { yaml_event_delete( &event ); break; } // Quit; yaml_event_delete(&event) is needed
		yaml_event_delete( &event );
	}
}

*/
/* Private (static) headers. */
static GNode *parseTemplate( xmlNodePtr cur );
static void set_field_attributes( xmlNodePtr node, FieldType *tfield );
static gboolean ignore_xml_node( xmlNodePtr cur );
static void parser_walk_children( xmlNodePtr cur,
								  GNode *parent, GNode *tnode_prev );
static GNode *new_parsed_field( xmlNodePtr xmlnode );
static gboolean field_type_match( xmlNodePtr node,
								  FieldTypeIdentifier type );
static gboolean parse_operator( xmlNodePtr xmlnode, FieldType *tfield );
static gboolean parse_decimal( xmlNodePtr xmlnode, FieldType *tfield, GNode *tnode );
static gboolean parse_sequence( xmlNodePtr xmlnode, FieldType *tfield, GNode *tnode );
static gboolean operator_type_match( xmlNodePtr node, FieldOperatorIdentifier type );

/*! \brief  Convert an XML file into an internal representation of
 *          the templates.
 * \param filename  Name of the XML file to parse.
 * \return  An internal tree of FieldTypes.
 */
GNode *parse_templates_xml( const char *filename )
{
	xmlDocPtr doc; /* pointer to XML document */
	xmlNodePtr cur; /* pointer to current node within document */
	GNode *templates; /* Return value, all templates stored as children here. */
	GNode *tnode_prev = 0;
	doc = xmlParseFile( filename ); /* attempt to parse xml file and get pointer to parsed document */
	if( doc == NULL )
	{
		fprintf( stderr, "Document not parsed successfully. \n" );
		return 0;
	}
	/* Start at the root of the XML document. */
	cur = xmlDocGetRootElement( doc );
	if( cur == NULL )
	{
		fprintf( stderr, "empty document\n" );
		xmlFreeDoc( doc );
		return 0;
	}
	templates = g_node_new( 0 );
	if( !templates ) mads_quits( 0 );
	/* Check if root is of type "templates". */
	if( xmlStrcmp( cur->name, ( xmlChar * ) "templates" ) )
	{
		fprintf( stderr, "document of the wrong type, root node != templates\n" );
		xmlFreeDoc( doc );
		return 0;
	}
	cur = cur->xmlChildrenNode;
	while( cur != NULL )
	{
		if( !ignore_xml_node( cur ) )
		{
			if( ( !xmlStrcasecmp( cur->name, ( xmlChar * ) "template" ) ) )
			{
				GNode *tnode;
				tnode = parseTemplate( cur );
				/* Parse template, and check if parse failed. */
				if( tnode )
				{
					g_node_insert_after( templates, tnode_prev, tnode );
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
				fprintf( stderr, "Warning: Unkown templates child: %s\n", cur->name );
				fprintf( stderr, "Continuing to parse templates file ...\n" );
			}
		}
		cur = cur->next;
	}
	xmlFreeDoc( doc );
	return templates;
}

/*! \brief Parse a template section of the XML document
 * \param cur  Current position in the XML structure.
 * \return  Tree containing the internal template definition.
 */
GNode *parseTemplate( xmlNodePtr cur )
{
	GNode *tnode;
	FieldType *tfield;
	tnode = create_field( FieldTypeUInt32, FieldOperatorCopy );
	if( !tnode ) { tprintf( "Error allocating memory." ); mads_quits( 0 ); }
	tfield = ( FieldType * ) tnode->data;
	/* Get name, id, etc. */
	set_field_attributes( cur, tfield );
	cur = cur->xmlChildrenNode;
	parser_walk_children( cur, tnode, 0 );
	return tnode;
}

/*! \brief  Build up the parse tree from XML data,
 *          starting at an arbitrary node.
 *
 * \param cur  Current node in an xml tree.
 * \param parent  This parse tree level's parent.
 * \param tnode_prev  Previous node in parse tree at this level.
 */
void parser_walk_children( xmlNodePtr cur,
						   GNode *parent, GNode *tnode_prev )
{
	while( cur != NULL )
	{
		if( !ignore_xml_node( cur ) )
		{
			GNode *tnode;
			tnode = new_parsed_field( cur );
			if( tnode )
			{
				/* DBG0("adding tnode"); */
				g_node_insert_after( parent, tnode_prev, tnode );
				tnode_prev = tnode;
			}
		}
		cur = cur->next;
	}
}


/*! \brief  Create a new field in the parse tree based on an XML node.
 * \param  The XML node which /should/ be a field.
 * \return  The new GNode* containing a FieldType.
 *          NULL if something went wrong.
 */
GNode *new_parsed_field( xmlNodePtr xmlnode )
{
	gboolean found = FALSE;  /* Field type found. */
	gboolean valid = FALSE;  /* Field type valid. */
	GNode *tnode;
	FieldType *tfield;
	tnode = create_field( FieldTypeInvalid, FieldOperatorNone );
	if( !tnode )
	{
		return 0;
	}
	tfield = ( FieldType * ) tnode->data;
	/* Try the integers. */
	if( !found )
	{
		found = TRUE;
		if( field_type_match( xmlnode, FieldTypeUInt32 ) )
		{
			tfield->type = FieldTypeUInt32;
		}
		else if( field_type_match( xmlnode, FieldTypeUInt64 ) )
		{
			tfield->type = FieldTypeUInt64;
		}
		else if( field_type_match( xmlnode, FieldTypeInt32 ) )
		{
			tfield->type = FieldTypeInt32;
		}
		else if( field_type_match( xmlnode, FieldTypeInt64 ) )
		{
			tfield->type = FieldTypeInt64;
		}
		else
		{
			found = FALSE;
		}
		if( found )
		{
			/* check if field has valid operators */
			valid = parse_operator( xmlnode, tfield );
		}
	}
	/* Try decimal. */
	if( !found && field_type_match( xmlnode, FieldTypeDecimal ) )
	{
		found = TRUE;
		valid = parse_decimal( xmlnode, tfield, tnode );
	}
	/* Try string */
	if( !found && field_type_match( xmlnode, FieldTypeAsciiString ) )
	{
		tfield->type = FieldTypeAsciiString;
		found = TRUE;
		valid = parse_operator( xmlnode, tfield );
	}
	else if( !found && field_type_match( xmlnode, FieldTypeUnicodeString ) )
	{
		tfield->type = FieldTypeUnicodeString;
		found = TRUE;
		valid = parse_operator( xmlnode, tfield );
	}
	/* Try bytevector */
	if( !found && field_type_match( xmlnode, FieldTypeByteVector ) )
	{
		tfield->type = FieldTypeByteVector;
		found = TRUE;
		valid = parse_operator( xmlnode, tfield );
	}
	/* Try group */
	if( !found && field_type_match( xmlnode, FieldTypeGroup ) )
	{
		tfield->type = FieldTypeGroup;
		found = TRUE;
		valid = TRUE;
		parser_walk_children( xmlnode->xmlChildrenNode, tnode, 0 );
	}
	/* Try sequence */
	if( !found && field_type_match( xmlnode, FieldTypeSequence ) )
	{
		found = TRUE;
		valid = parse_sequence( xmlnode, tfield, tnode );
	}
	/* If we retrieved built up the field,
	 * fill its attributes.
	 */
	if( found && valid )
	{
		set_field_attributes( xmlnode, tfield );
	}
	else
	{
		/* TODO: Free parse tree. */
		tnode = 0;
		if( !valid )
		{
			tprintf( "Field %s could not be parsed.", xmlnode->name );
		}
		else
		{
			tprintf( "Unknown field type %s.", xmlnode->name );
		}
	}
	return tnode;
}


/*! \brief  Fill in a decimal field in the parse tree.
 * \param  The XML node which /should/ be a decimal field.
 * \param  A pointer to the template within the parse tree.
 * \return  True if sucessfully parsed
 */
gboolean parse_decimal( xmlNodePtr xmlnode, FieldType *tfield, GNode *tnode )
{
	GNode *exptNode;
	GNode *mantNode;
	FieldType *expt;
	FieldType *mant;
	tfield->type = FieldTypeDecimal;
	/* Add exponent and mantissa. */
	exptNode = create_field( FieldTypeInt32, FieldOperatorNone );
	if( !exptNode )  { tprintf( "Error creating exponent field." ); mads_quits( 0 ); }
	g_node_insert_after( tnode, 0,        exptNode );
	mantNode = create_field( FieldTypeInt64, FieldOperatorNone );
	if( !mantNode )  { tprintf( "Error creating mantissa field." ); mads_quits( 0 ); }
	g_node_insert_after( tnode, exptNode, mantNode );
	/* Get exponent and mantissa operators and values (if given) */
	xmlnode = xmlnode->xmlChildrenNode;
	while( xmlnode != NULL )
	{
		if( !ignore_xml_node( xmlnode ) )
		{
			if( 0 == xmlStrcasecmp( xmlnode->name, ( xmlChar * )"exponent" ) )
			{
				expt = ( FieldType * ) exptNode->data;
				if( !parse_operator( xmlnode, expt ) ) return FALSE;
			}
			else if( 0 == xmlStrcasecmp( xmlnode->name, ( xmlChar * )"mantissa" ) )
			{
				mant = ( FieldType * ) mantNode->data;
				if( !parse_operator( xmlnode, mant ) ) return FALSE;
			}
			else
			{
				tprintf( "Unknown decimal subfield %s.", xmlnode->name );
				return FALSE;
			}
		}
		xmlnode = xmlnode->next;
	}
	return TRUE;
}


/*! \brief  Fill in a sequence field in the parse tree.
 * \param  The XML node which /should/ be a sequence field.
 * \param  A pointer to the template within the parse tree.
 * \return  True if sucessfully parsed
 */
gboolean parse_sequence( xmlNodePtr xmlnode, FieldType *tfield, GNode *tnode )
{
	GNode *group_tnode;
	tfield->type = FieldTypeSequence;
	group_tnode = create_field( FieldTypeGroup, FieldOperatorNone );
	if( !group_tnode )
	{
		return FALSE;
	}
	g_node_insert_after( tnode, 0, group_tnode );
	/* Descend the tree on the group. */
	parser_walk_children( xmlnode->xmlChildrenNode, group_tnode, 0 );
	return TRUE;
}


/*! \brief  Fill in a field in the parse tree with operator info.
 * \param  The XML node which /should/ be a field.
 * \param  A pointer to the field within the parse tree.
 * \return  True if sucessfully parsed
 */
gboolean parse_operator( xmlNodePtr xmlnode, FieldType *tfield )
{
	xmlChar *prop;
	const xmlChar *name;
	name = xmlnode->name;
	/* loop through field to find operators */
	xmlnode = xmlnode->xmlChildrenNode;
	while( xmlnode != NULL )
	{
		if( !ignore_xml_node( xmlnode ) )
		{
			if( operator_type_match( xmlnode, FieldOperatorConstant ) )
			{
				tfield->op = FieldOperatorConstant;
			}
			else if( operator_type_match( xmlnode, FieldOperatorDefault ) )
			{
				tfield->op = FieldOperatorDefault;
			}
			else if( operator_type_match( xmlnode, FieldOperatorCopy ) )
			{
				tfield->op = FieldOperatorCopy;
			}
			else if( operator_type_match( xmlnode, FieldOperatorIncrement ) )
			{
				tfield->op = FieldOperatorIncrement;
			}
			else if( operator_type_match( xmlnode, FieldOperatorDelta ) )
			{
				tfield->op = FieldOperatorDelta;
			}
			else if( operator_type_match( xmlnode, FieldOperatorTail ) )
			{
				tfield->op = FieldOperatorTail;
			}
			else
			{
				tprintf( "Invalid operator (%s) for field %s", xmlnode->name, name );
				return FALSE;
			}
			/* get value of operator if given */
			prop = xmlGetProp( xmlnode, ( xmlChar * )"value" );
			if( prop != NULL )
			{
				tfield->value = prop;
			}
			else
			{
				tfield->value = NULL;
			}
		}
		xmlnode = xmlnode->next;
	}
	return TRUE;
}


/*! \brief  Set standard attributes of a field.
 *
 * Specifically, set 'id', 'name', and 'presence'.
 * \param node  XML node to query for attributes.
 * \param tfield  Return value. Must already be allocated.
 */
void set_field_attributes( xmlNodePtr xmlnode, FieldType *tfield )
{
	const xmlChar *str;
	/* Name. */
	str = xmlGetProp( xmlnode, ( xmlChar * ) "name" );
	if( str )
	{
		tfield->name = g_strdup( ( char * )str );
	}
	xmlFree( ( void * )str );
	/* Identifier number. */
	str = xmlGetProp( xmlnode, ( xmlChar * ) "id" );
	if( str )
	{
		tfield->id = atoi( ( char * )str );
	}
	xmlFree( ( void * )str );
	/* Presence. */
	str = xmlGetProp( xmlnode, ( xmlChar * ) "presence" );
	if( str )
	{
		if( xmlStrcasecmp( str, ( xmlChar * ) "optional" ) )
		{
			tfield->mandatory = FALSE;
		}
		else if( xmlStrcasecmp( str, ( xmlChar * ) "mandatory" ) )
		{
			tfield->mandatory = TRUE;
		}
		else
		{
			tprintf( "Error, bad presence option '%s'.", ( char * ) str );
		}
	}
	xmlFree( ( void * )str );
}

/*! \brief  Check if this XML node should be ignored.
 *
 * Ignored nodes include comments and text between XML tags.
 *
 * \param xmlnode  Current XML node.
 * \return  TRUE iff the node is something we ignore.
 */
gboolean ignore_xml_node( xmlNodePtr xmlnode )
{
	return ( 0 == xmlStrcasecmp( xmlnode->name, ( xmlChar * ) "text" ) ||
			 0 == xmlStrcasecmp( xmlnode->name, ( xmlChar * ) "comment" ) );
}

/*! \brief  Check if a string from XML.
 * \param node  Node in the XML tree. Its name is to be used in comparison.
 * \param type  Type whose name you want to check against.
 * \return  TRUE iff the node name matches.
 */
gboolean field_type_match( xmlNodePtr node, FieldTypeIdentifier type )
{
	const xmlChar *str1 = node->name;
	const char *str2 = field_typename( type );
	if( type == FieldTypeAsciiString )
	{
		if( 0 == xmlStrcasecmp( str1, ( xmlChar * ) "string" ) )
		{
			/* check if charset is ascii or unicode */
			xmlChar *prop = xmlGetProp( node, ( xmlChar * ) "charset" );
			if( prop == NULL )
			{
				/* Assume ascii if not given */
				return TRUE;
			}
			if( 0 == xmlStrcasecmp( prop, ( xmlChar * ) "ascii" ) )
			{
				return TRUE;
			}
		}
		return FALSE;
	}
	else if( type == FieldTypeUnicodeString )
	{
		if( 0 == xmlStrcasecmp( str1, ( xmlChar * ) "string" ) )
		{
			/* check if charset is ascii or unicode */
			xmlChar *prop = xmlGetProp( node, ( xmlChar * )"charset" );
			if( prop == NULL )
			{
				/* Assume ascii if not given */
				return FALSE;
			}
			if( 0 == xmlStrcasecmp( prop, ( xmlChar * ) "unicode" ) )
			{
				return TRUE;
			}
		}
		return FALSE;
	}
	else
	{
		return ( 0 == xmlStrcasecmp( str1, ( xmlChar * ) str2 ) );
	}
}

/*! \brief  Check if a string from XML matches the given Operator's type string
 * \param node  Node in the XML tree. Its name is to be used in comparison.
 * \param type  Type whose name you want to check against.
 * \return  TRUE iff the node name matches.
 */
gboolean operator_type_match( xmlNodePtr node, FieldOperatorIdentifier type )
{
	const xmlChar *str1 = node->name;
	const char *str2 = operator_typename( type );
	return ( 0 == xmlStrcasecmp( str1, ( xmlChar * ) str2 ) );
}

