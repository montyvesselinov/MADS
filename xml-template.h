/*
 * xml-template.h
 *
 *  Created on: Jul 29, 2014
 *      Author: monty
 */

#ifndef XML_TEMPLATE_H_
#define XML_TEMPLATE_H_
#include <glib.h>

enum field_type_identifier_enum
{
	FieldTypeUInt32,
	FieldTypeUInt64,
	FieldTypeInt32,
	FieldTypeInt64,
	FieldTypeDecimal,
	FieldTypeAsciiString,
	FieldTypeUnicodeString,
	FieldTypeByteVector,
	FieldTypeGroup,
	FieldTypeSequence,
	FieldTypeEnumLimit,
	FieldTypeInvalid
};
typedef enum field_type_identifier_enum FieldTypeIdentifier;

enum field_operator_identifier_enum
{
	FieldOperatorNone,
	FieldOperatorConstant,
	FieldOperatorDefault,
	FieldOperatorCopy,
	FieldOperatorIncrement,
	FieldOperatorDelta,
	FieldOperatorTail,
	FieldOperatorEnumLimit
};
typedef enum field_operator_identifier_enum FieldOperatorIdentifier;

typedef void FieldValue;

/*! \brief  Hold data relevant to a template definition.
 */
struct field_type_struct
{
	char *name;
	/*! \brief  Template id. Typed for the hash lookup. */
	gint id;
	gboolean mandatory;
	FieldTypeIdentifier type;
	FieldOperatorIdentifier op;
	FieldValue *value;
};
typedef struct field_type_struct FieldType;

void add_templates( GNode *tmpl );
const gchar *field_typename( FieldTypeIdentifier type );
const gchar *operator_typename( FieldOperatorIdentifier type );
GNode *create_field( FieldTypeIdentifier type, FieldOperatorIdentifier op );
GNode *find_template( guint32 id );

#endif /* XML_TEMPLATE_H_ */
