// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
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
/*
 *      Author: Daniel O'Malley <omalled@lanl.gov>
 */

#include <stdlib.h>
#include <stdio.h>

#include "pqueue.h"

/* Util macros */
#define LEFT(x) (2 * (x) + 1)
#define RIGHT(x) (2 * (x) + 2)
#define PARENT(x) ((x-1) / 2)

void pqueue_heapify( PQueue *q, size_t idx );

/**
* Allocates memory for a new Priority Queue structure .

* 'cmp' function:
*   returns 0 if d1 and d2 have the same priorities
*   returns [negative value] if d1 have a smaller priority than d2
*   returns [positive value] if d1 have a greater priority than d2
*/
PQueue *pqueue_new( int ( *cmp )( const void *d1, const void *d2 ),
					size_t capacity )
{
	PQueue *res = NULL;
	NP_CHECK( cmp );
	res = malloc( sizeof( *res ) );
	NP_CHECK( res );
	res->cmp = cmp;
	/* The inner representation of data inside the queue is an array of void* */
	res->data = malloc( capacity * sizeof( *( res->data ) ) );
	NP_CHECK( res->data );
	res->size = 0;
	res->capacity = capacity;
	return ( res );
}

/**
* De-allocates memory for a given Priority Queue structure .
*/
void pqueue_delete( PQueue *q )
{
	if( NULL == q )
	{
		DEBUG( "Priority Queue is already NULL. Nothing to free." );
		return;
	}
	free( q->data );
	free( q );
}

/**
* Adds a new element to the Priority Queue .
*/
void pqueue_enqueue( PQueue *q, const void *data )
{
	size_t i;
	void *tmp = NULL;
	NP_CHECK( q );
	if( q->size >= q->capacity )
	{
		DEBUG( "Priority Queue is full. Cannot add another element ." );
		return;
	}
	/* Adds element last */
	q->data[q->size] = ( void * ) data;
	i = q->size;
	q->size++;
	/* The new element is swapped with its parent as long as its
	precedence is higher */
	while( i > 0 && q->cmp( q->data[i], q->data[PARENT( i )] ) > 0 )
	{
		tmp = q->data[i];
		q->data[i] = q->data[PARENT( i )];
		q->data[PARENT( i )] = tmp;
		i = PARENT( i );
	}
}

/**
* Returns the element with the biggest priority from the queue .
*/
void *pqueue_dequeue( PQueue *q )
{
	void *data = NULL;
	NP_CHECK( q );
	if( q->size < 1 ) {         /* Priority Queue is empty */         DEBUG( "Priority Queue underflow . Cannot remove another element ." );         return NULL;     }     data = q->data[0];
	q->data[0] = q->data[q->size - 1];
	q->size--;
	/* Restore heap property */
	pqueue_heapify( q, 0 );
	return ( data );
}

/**
* Turn an "almost-heap" into a heap .
*/
void pqueue_heapify( PQueue *q, size_t idx )
{
	/* left index, right index, largest */
	void *tmp = NULL;
	size_t l_idx, r_idx, lrg_idx;
	NP_CHECK( q );
	l_idx = LEFT( idx );
	r_idx = RIGHT( idx );
	/* Left child exists, compare left child with its parent */
	if( l_idx < q->size && q->cmp( q->data[l_idx], q->data[idx] ) > 0 )
	{
		lrg_idx = l_idx;
	}
	else
	{
		lrg_idx = idx;
	}
	/* Right child exists, compare right child with the largest element */
	if( r_idx < q->size && q->cmp( q->data[r_idx], q->data[lrg_idx] ) > 0 )
	{
		lrg_idx = r_idx;
	}
	/* At this point largest element was determined */
	if( lrg_idx != idx )
	{
		/* Swap between the index at the largest element */
		tmp = q->data[lrg_idx];
		q->data[lrg_idx] = q->data[idx];
		q->data[idx] = tmp;
		/* Heapify again */
		pqueue_heapify( q, lrg_idx );
	}
}
