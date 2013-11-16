// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov/
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
/*
 *      Author: Daniel O'Malley <omalled@lanl.gov>
 */

#ifndef __PQUEUE__H__
#define __PQUEUE__H__

/**
* Debugging macro .
*
* Checks for a NULL pointer, and prints the error message, source file and
* line via 'stderr' .
* If the check fails the program exits with error code (-1) .
*/
#define NP_CHECK(ptr) \
	{ \
		if (NULL == (ptr)) { \
			fprintf(stderr, "%s:%d NULL POINTER: %s n", \
					__FILE__, __LINE__, #ptr); \
			exit(-1); \
		} \
	} \
	 
#define DEBUG(msg) fprintf(stderr, "%s:%d %s", __FILE__, __LINE__, (msg))

/**
* Priority Queue Structure
*/
typedef struct PQueue_s
{
	/* The actual size of heap at a certain time */
	size_t size;
	/* The amount of allocated memory for the heap */
	size_t capacity;
	/* An array of (void*), the actual max-heap */
	void **data;
	/* A pointer to a comparator function, used to prioritize elements */
	int ( *cmp )( const void *d1, const void *d2 );
} PQueue;

/** Allocates memory for a new Priority Queue .
Needs a pointer to a comparator function, thus establishing priorities .
*/
PQueue *pqueue_new( int ( *cmp )( const void *d1, const void *d2 ),
					size_t capacity );

/** De-allocates memory for a given Priority Queue */
void pqueue_delete( PQueue *q );

/** Add an element inside the Priority Queue */
void pqueue_enqueue( PQueue *q, const void *data );

/** Removes the element with the greatest priority from within the Queue */
void *pqueue_dequeue( PQueue *q );

#endif
