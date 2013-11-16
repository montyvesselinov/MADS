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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "astable.h"
#include "pqueue.h"

#define DEFAULT_PERCENTILE 0.0001
#define DEFAULT_ABSERR 1e-6

struct interval
{
	double left;
	double right;
	double func_left;
	double func_right;
	double cdf_left;
	double cdf_right;
};

void copy_interpolant_params( struct interpolant_params *from, struct interpolant_params *to )
{
	to->alpha = from->alpha;
	to->beta = from->beta;
	to->acc = from->acc;
	to->spline = from->spline;
	to->ymax = from->ymax;
	to->ymin = from->ymin;
	to->size = from->size;
	to->limiting_coefficient = from->limiting_coefficient;
	to->cdf_or_pdf = from->cdf_or_pdf;
}

void symmetric_astable_pdf_interp( double x, double alpha, double gamma, double lambda, double *val )
{
	const int num_interps = 19;
	static struct interpolant_params ips[19];
	static int setup = 0;
	struct interpolant_params *ip;
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	int i;
	double alpha1, alpha2;
	double interp1, interp2;
	if( alpha < 0.1 )
	{
		fprintf( stderr, "value of alpha is too small for symmetric_astable_pdf_interp\n" );
		exit( 1 );
	}
	else
	{
		if( !setup )
		{
			pthread_mutex_lock( &mutex );
			if( !setup )
			{
				setup = 1;
				for( i = 0; i < num_interps; i++ )
				{
					ip = automate_interpolant( 2. * ( i + 1. ) / ( num_interps + 1. ), 0., DEFAULT_PERCENTILE, DEFAULT_ABSERR, INTERP_PDF );
					copy_interpolant_params( ip, ips + i );
				}
			}
			pthread_mutex_unlock( &mutex );
		}
		if( alpha == 2. ) *val = astable_pdf( x, 2., 0., gamma, lambda );
		else
		{
			i = floor( alpha * ( num_interps + 1. ) / 2. ) - 1;
			alpha1 = 2. * ( i + 1. ) / ( num_interps + 1. );
			alpha2 = 2. * ( i + 2. ) / ( num_interps + 1. );
			interp1 = interpolate( x, gamma, lambda, ips + i );
			if( alpha2 < 2. - .5 / ( num_interps + 1. ) ) interp2 = interpolate( x, gamma, lambda, ips + i + 1 );
			else interp2 = astable_pdf( x, 2., 0., gamma, lambda );
			*val = interp2 * ( alpha1 - alpha ) + interp1 * ( alpha - alpha2 );
			*val /= ( alpha1 - alpha2 );
		}
	}
}

void symmetric_astable_cdf_interp( double x, double alpha, double gamma, double lambda, double *val )
{
	const int num_interps = 19;
	static struct interpolant_params ips[19];
	static int setup = 0;
	struct interpolant_params *ip;
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	int i;
	double alpha1, alpha2;
	double interp1, interp2;
	if( alpha < 0.1 )
	{
		fprintf( stderr, "value of alpha is too small for symmetric_astable_cdf_interp\n" );
		exit( 1 );
	}
	else
	{
		if( !setup )
		{
			pthread_mutex_lock( &mutex );
			if( !setup )
			{
				setup = 1;
				for( i = 0; i < num_interps; i++ )
				{
					ip = automate_interpolant( 2. * ( i + 1. ) / ( num_interps + 1. ), 0., DEFAULT_PERCENTILE, DEFAULT_ABSERR, INTERP_CDF );
					copy_interpolant_params( ip, ips + i );
				}
			}
			pthread_mutex_unlock( &mutex );
		}
		if( alpha == 2. ) *val = astable_cdf( x, 2., 0., gamma, lambda );
		else
		{
			i = floor( alpha * ( num_interps + 1. ) / 2. ) - 1;
			alpha1 = 2. * ( i + 1. ) / ( num_interps + 1. );
			alpha2 = 2. * ( i + 2. ) / ( num_interps + 1. );
			interp1 = interpolate( x, gamma, lambda, ips + i );
			if( alpha2 < 2. - .5 / ( num_interps + 1. ) ) interp2 = interpolate( x, gamma, lambda, ips + i + 1 );
			else interp2 = astable_cdf( x, 2., 0., gamma, lambda );
			*val = interp2 * ( alpha1 - alpha ) + interp1 * ( alpha - alpha2 );
			*val /= ( alpha1 - alpha2 );
		}
	}
}

struct interpolant_params *automate_interpolant( double alpha, double beta, double percentile, double abserr, int CDF_OR_PDF )
{
	double left;
	double right;
	struct interpolant_params *ip;

	right = 1;
	while( standard_astable_cdf( right, alpha, beta ) < 1. - percentile )
	{
		right *= 2;
	}
	left = -1;
	while( standard_astable_cdf( left, alpha, beta ) > percentile )
	{
		left *= 2;
	}

	//the 100 + 10. / sqrt(abserr) gives a very rough approximation of the number of node points needed to acheive the desired accuracy
	//the sqrt arises because the truncation error associated with linear interpolation is O(h^2)
	//the approximation will likely be problematic for small alpha and beta near +-1
	ip = setup_interpolant( alpha, beta, left, right, 100 + 10. / sqrt( abserr ), CDF_OR_PDF );

	return ip;
}

double interpolate( double x, double gamma, double lambda, struct interpolant_params *ip )
{
	double retval;
	double normalized_x;
	//TODO: find the expressions for when \beta=-1,x\to\infty and \beta=1,x\to-\infty
	normalized_x = ( x - gamma ) / lambda;
	if( normalized_x <= ip->spline->interp->xmax && normalized_x >= ip->spline->interp->xmin )
	{
		retval = gsl_spline_eval( ip->spline, normalized_x, ip->acc ); // / lambda;
		if( ip->cdf_or_pdf == INTERP_PDF ) retval /= lambda;
	}
	else if( normalized_x > ip->spline->interp->xmax )
	{
		if( ip->cdf_or_pdf == INTERP_PDF )
		{
			if( ip->beta > -1. )
			{
				retval = ip->limiting_coefficient * ( 1 + ip->beta ) * ( ip->alpha * pow( lambda, ip->alpha ) * pow( x, -1 - ip->alpha ) );
			}
			else if( ip->alpha < 1 && x > gamma + lambda * tan( M_PI * ip->alpha / 2. ) ) //beta == -1 in this case
			{
				retval = 0.;
			}
			else
			{
				fprintf( stderr, "Warning: You are asking for a value outside the interpolation range for which no asymptotic estimate is used. Returning value at right end point.\n" );
				retval = ip->ymax;
			}
		}
		else
		{
			//TODO: fix this code (account for beta, etc)
			retval =  1. - ( 1. - ip->ymax ) * pow( ip->spline->interp->xmax / normalized_x, ip->alpha );
		}
	}
	else
	{
		if( ip->cdf_or_pdf == INTERP_PDF )
		{
			if( ip->beta < 1. )
			{
				retval = ip->limiting_coefficient * ( 1 + ip->beta ) * ( ip->alpha * pow( lambda, ip->alpha ) * pow( -x, -1 - ip->alpha ) );
			}
			else if( ip->alpha < 1 && x < gamma - lambda * tan( M_PI * ip->alpha / 2. ) ) //beta == 1 in this case
			{
				retval = 0.;
			}
			else
			{
				fprintf( stderr, "Warning: You are asking for a value outside the interpolation range for which no asymptotic estimate is used. Returning value at left end point.\n" );
				retval = ip->ymin;
			}
		}
		else
		{
			//TODO: fix this code (account for beta, etc)
			//return astable_cdf(x, ip->alpha, ip->beta, gamma, lambda);
			retval = ip->ymin * pow( ip->spline->interp->xmin / normalized_x, ip->alpha );
		}
	}
	return retval;
}

int interval_comp( const void *d1, const void *d2 )
{
	const struct interval *i1;
	const struct interval *i2;
	double difference;
	i1 = ( const struct interval * )d1;
	i2 = ( const struct interval * )d2;
	difference = i1->cdf_right + i2->cdf_left - ( i2->cdf_right + i1->cdf_left );
	if( difference > 0 ) return 1;
	else if( difference > 0 ) return -1;
	else return 0;
}

int interval_comp_x( const void *d1, const void *d2 )
{
	const struct interval *i1;
	const struct interval *i2;
	double difference;
	i1 = ( const struct interval * )d1;
	i2 = ( const struct interval * )d2;
	difference = i2->left - i1->left;
	if( difference < 0 ) return -1;
	else if( difference > 0 ) return 1;
	else return 0;
}

void get_sorted_vals( double alpha, double beta, double left, double right, int N, int CDF_OR_PDF, double *x_sorted, double *y_sorted )
{
	struct interval *worst_interval;
	struct interval *new_interval[2];
	struct interval *current_interval;
	double func_left;
	double func_right;
	double cdf_left;
	double cdf_right;
	double mid;
	double func_mid;
	double cdf_mid;
	PQueue *q;
	PQueue *q_x;
	int i;
	if( left > -1 ) left = -1.; //We need the lower end point to be < 0 in order for the limiting behavior code to work.
	if( right < 1 ) right = 1.; //We need the upper end point to be > 0 in order for the limiting behavior code to work.
	q = pqueue_new( &interval_comp, N );
	worst_interval = ( struct interval * )malloc( sizeof( struct interval ) );
	worst_interval->left = left;
	worst_interval->right = right;
	if( CDF_OR_PDF == INTERP_CDF )
	{
		func_left = standard_astable_cdf( left, alpha, beta );
		func_right = standard_astable_cdf( right, alpha, beta );
		cdf_left = func_left;
		cdf_right = func_right;
	}
	else
	{
		func_left = standard_astable_pdf( left, alpha, beta );
		func_right = standard_astable_pdf( right, alpha, beta );
		cdf_left = standard_astable_cdf( left, alpha, beta );
		cdf_right = standard_astable_cdf( right, alpha, beta );
	}
	worst_interval->func_right = func_right;
	worst_interval->func_left = func_left;
	worst_interval->cdf_left = cdf_left;
	worst_interval->cdf_right = cdf_right;
	pqueue_enqueue( q, worst_interval );
	for( i = 0; i < N - 1; i++ )
	{
		worst_interval = ( struct interval * )pqueue_dequeue( q );
		left = worst_interval->left;
		right = worst_interval->right;
		mid = .5 * ( right + left );
		func_left = worst_interval->func_left;
		func_right = worst_interval->func_right;
		cdf_left = worst_interval->cdf_left;
		cdf_right = worst_interval->cdf_right;
		func_mid = ( CDF_OR_PDF == INTERP_CDF ? standard_astable_cdf( mid, alpha, beta ) : standard_astable_pdf( mid, alpha, beta ) );
		cdf_mid = ( CDF_OR_PDF == INTERP_CDF ? func_mid : standard_astable_cdf( mid, alpha, beta ) );
		new_interval[0] = ( struct interval * )malloc( sizeof( struct interval ) );
		new_interval[0]->left = left;
		new_interval[0]->right = mid;
		new_interval[0]->func_left = func_left;
		new_interval[0]->func_right = func_mid;
		new_interval[0]->cdf_left = cdf_left;
		new_interval[0]->cdf_right = cdf_mid;
		pqueue_enqueue( q, new_interval[0] );
		new_interval[1] = ( struct interval * )malloc( sizeof( struct interval ) );
		new_interval[1]->left = mid;
		new_interval[1]->right = right;
		new_interval[1]->func_left = func_mid;
		new_interval[1]->func_right = func_right;
		new_interval[1]->cdf_left = cdf_mid;
		new_interval[1]->cdf_right = cdf_right;
		pqueue_enqueue( q, new_interval[1] );
		free( worst_interval );
	}
	/* Now sort them in increasing values of the left interval point. */
	q_x = pqueue_new( &interval_comp_x, N );
	i = 0;
	/* Take them off the old queue and add them to the new queue. */
	current_interval = NULL;
	while( q->size > 0 )
	{
		i++;
		current_interval = ( struct interval * )pqueue_dequeue( q );
		pqueue_enqueue( q_x, current_interval );
	}
	i = 0;
	/* Put them on the new queue. */
	while( q_x->size > 0 )
	{
		if( i > 0 ) free( current_interval );
		current_interval = ( struct interval * )pqueue_dequeue( q_x );
		x_sorted[i] = current_interval->left;
		y_sorted[i] = current_interval->func_left;
		i++;
	}
	x_sorted[i] = current_interval->right;
	y_sorted[i] = current_interval->func_right;
	free( current_interval );
	pqueue_delete( q );
	pqueue_delete( q_x );
}

struct interpolant_params *setup_interpolant( double alpha, double beta, double left, double right, int N, int CDF_OR_PDF )
{
	struct interpolant_params *ip;
	double *x_sorted;
	double *y_sorted;

	x_sorted = ( double * )malloc( ( N + 1 ) * sizeof( double ) );
	y_sorted = ( double * )malloc( ( N + 1 ) * sizeof( double ) );

	get_sorted_vals( alpha, beta, left, right, N, CDF_OR_PDF, x_sorted, y_sorted );

	ip = ( struct interpolant_params * )malloc( sizeof( struct interpolant_params ) );
	ip->alpha = alpha;
	ip->beta = beta;
	ip->ymax = y_sorted[N];
	ip->ymin = y_sorted[0];
	ip->size = N + 1;
	ip->acc = gsl_interp_accel_alloc();
	ip->spline = gsl_spline_alloc( gsl_interp_linear, N + 1 );
	gsl_spline_init( ip->spline, x_sorted, y_sorted, N + 1 );
	ip->limiting_coefficient = sin( M_PI *alpha / 2. ) *gsl_sf_gamma( alpha ) / M_PI;
	ip->cdf_or_pdf = CDF_OR_PDF;

	free( x_sorted );
	free( y_sorted );

	return ip;
}

void free_interpolant_params( struct interpolant_params *ip )
{
	gsl_spline_free( ip->spline );
	gsl_interp_accel_free( ip->acc );
	free( ip );
}
