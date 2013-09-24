/*
 * interpolation.c
 *
 *  Created on: Jul 24, 2013
 *      Author: Daniel O'Malley <omalled@lanl.gov>
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "astable.h"
#include "pqueue.h"

struct interval
{
	double left;
	double right;
	double func_left;
	double func_right;
	double cdf_left;
	double cdf_right;
};

struct interpolant_params global_ip[1024];//Global memory! Seriously?!?!?!? Yep!!!!!

void astable_cdf_interp( double x, double alpha, double beta, double gamma, double lambda, double *val)
{
	static int num_in_use = 0;
	static int last_use = 0;

	struct interpolant_params *ip;
	int i;

	if(alpha == 2)
	{
		double a;
		a = astable_cdf( x, alpha, beta, gamma, lambda );
		val[0] = a;
		return;
	}

	i = 0;
	while( i < num_in_use && ( global_ip[( i + last_use ) % num_in_use].alpha != alpha || global_ip[( i + last_use ) % num_in_use].beta != beta ) )
	{
		i++;
	}
	if( i == num_in_use )
	{
		if( num_in_use == 1024 )
		{
			fprintf( stderr, "Too many alpha, beta combinations in use\n" );
			exit(1);
		}
		else
		{
			ip = automate_interpolant( alpha, beta, 0.00001, 1e-6, INTERP_CDF );
			global_ip[num_in_use].alpha = alpha;
			global_ip[num_in_use].beta = beta;
			global_ip[num_in_use].acc = ip->acc;
			global_ip[num_in_use].spline = ip->spline;
			global_ip[num_in_use].ymax = ip->ymax;
			global_ip[num_in_use].ymin = ip->ymin;
			global_ip[num_in_use].size = ip->size;
			global_ip[num_in_use].limiting_coefficient = ip->limiting_coefficient;
			global_ip[num_in_use].cdf_or_pdf = ip->cdf_or_pdf;
			num_in_use++;
		}
	}

	last_use = i;

	val[0] = interpolate( x, gamma, lambda, global_ip + i );
	return;

	//return interpolate( x, gamma, lambda, &global_ip[i] );
}

struct interpolant_params *automate_interpolant(double alpha, double beta, double percentile, double abserr, int CDF_OR_PDF)
{
	double left;
	double right;
	struct interpolant_params *ip;

	right = 1;
	while(standard_astable_cdf(right, alpha, beta) < 1. - percentile)
	{
		right *= 2;
	}
	left = -1;
	while(standard_astable_cdf(left, alpha, beta) > percentile)
	{
		left *= 2;
	}

	//the 100 + 10. / sqrt(abserr) gives a very rough approximation of the number of node points needed to acheive the desired accuracy
	//the sqrt arises because the truncation error associated with linear interpolation is O(h^2)
	//the approximation will likely be problematic for small alpha and beta near +-1
	ip = setup_interpolant(alpha, beta, left, right, 100 + 10. / sqrt(abserr), CDF_OR_PDF);

	return ip;
}

double interpolate(double x, double gamma, double lambda, struct interpolant_params *ip)
{
	double normalized_x;

	//TODO: find the expressions for when \beta=-1,x\to\infty and \beta=1,x\to-\infty
	normalized_x = (x - gamma) / lambda;
	if(normalized_x <= ip->spline->interp->xmax && normalized_x >= ip->spline->interp->xmin)
	{
		double a;
		a = gsl_spline_eval(ip->spline, normalized_x, ip->acc);// / lambda;
		if(ip->cdf_or_pdf == INTERP_PDF) a = a / lambda;
		return a;
		//return gsl_spline_eval(ip->spline, (x - gamma) / lambda, ip->acc) / lambda;
	}
	else if(normalized_x > ip->spline->interp->xmax)
	{
		if(ip->cdf_or_pdf == INTERP_PDF)
		{
			if(ip->beta > -1.)
			{
				return ip->limiting_coefficient * (1 + ip->beta) * (ip->cdf_or_pdf == INTERP_CDF ? pow(gamma / x, ip->alpha) : ip->alpha * pow(gamma, ip->alpha) * pow(x, -1 - ip->alpha));
			}
			else if(ip->alpha < 1 && x > gamma + lambda * tan(M_PI * ip->alpha / 2.))//beta == -1 in this case
			{
				return 0.;
			}
			else
			{
				fprintf(stderr, "Warning: You are asking for a value outside the interpolation range for which no asymptotic estimate is used. Returning value at right end point.\n");
				return ip->ymax;
			}
		}
		else
		{
			//TODO: fix this code (account for beta, etc)
			return 1. - (1. - ip->ymax) * pow(ip->spline->interp->xmax / normalized_x, ip->alpha);
		}
	}
	else
	{
		if(ip->cdf_or_pdf == INTERP_PDF)
		{
			if(ip->beta < 1.)
			{
				return ip->limiting_coefficient * (1 + ip->beta) * (ip->cdf_or_pdf == INTERP_CDF ? pow(gamma / -x, ip->alpha) : ip->alpha * pow(gamma, ip->alpha) * pow(-x, -1 - ip->alpha));
			}
			else if(ip->alpha < 1 && x < gamma - lambda * tan(M_PI * ip->alpha / 2.))//beta == 1 in this case
			{
				return 0.;
			}
			else
			{
				fprintf(stderr, "Warning: You are asking for a value outside the interpolation range for which no asymptotic estimate is used. Returning value at left end point.\n");
				return ip->ymin;
			}
		}
		else
		{
			//TODO: fix this code (account for beta, etc)
			//return astable_cdf(x, ip->alpha, ip->beta, gamma, lambda);
			return ip->ymin * pow(ip->spline->interp->xmin / normalized_x, ip->alpha);
			//return 0.;
		}
	}
}

int interval_comp(const void *d1, const void *d2)
{
	const struct interval *i1;
	const struct interval *i2;
	double difference;

	i1 = (const struct interval *)d1;
	i2 = (const struct interval *)d2;

	difference = i1->cdf_right + i2->cdf_left - (i2->cdf_right + i1->cdf_left);

	if(difference > 0) return 1;
	else if(difference > 0) return -1;
	else return 0;
}

int interval_comp_x(const void *d1, const void *d2)
{
	const struct interval *i1;
	const struct interval *i2;
	double difference;

	i1 = (const struct interval *)d1;
	i2 = (const struct interval *)d2;

	difference = i2->left - i1->left;

	if(difference < 0) return -1;
	else if(difference > 0) return 1;
	else return 0;
}

struct interpolant_params *setup_interpolant(double alpha, double beta, double left, double right, int N, int CDF_OR_PDF)
{
	struct interpolant_params *ip;
	struct interval *worst_interval;
	struct interval *new_interval[2];
	struct interval *current_interval;
	double *x_sorted;
	double *y_sorted;
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

	if(left > -1) left = -1.;//We need the lower end point to be < 0 in order for the limiting behavior code to work.
	if(right < 1) right = 1.;//We need the upper end point to be > 0 in order for the limiting behavior code to work.

	q = pqueue_new(&interval_comp, N);

	x_sorted = (double *)malloc((N + 1) * sizeof(double));
	y_sorted = (double *)malloc((N + 1) * sizeof(double));

	worst_interval = (struct interval *)malloc(sizeof(struct interval));
	worst_interval->left = left;
	worst_interval->right = right;
	if(CDF_OR_PDF == INTERP_CDF)
	{
		func_left = standard_astable_cdf(left, alpha, beta);
		func_right = standard_astable_cdf(right, alpha, beta);
		cdf_left = func_left;
		cdf_right = func_right;
	}
	else
	{
		func_left = standard_astable_pdf(left, alpha, beta);
		func_right = standard_astable_pdf(right, alpha, beta);
		cdf_left = standard_astable_cdf(left, alpha, beta);
		cdf_right = standard_astable_cdf(right, alpha, beta);
	}

	worst_interval->func_right = func_right;
	worst_interval->func_left = func_left;
	worst_interval->cdf_left = cdf_left;
	worst_interval->cdf_right = cdf_right;

	pqueue_enqueue(q, worst_interval);

	for(i = 0; i < N-1; i++)
	{
		worst_interval = (struct interval *)pqueue_dequeue(q);
		left = worst_interval->left;
		right = worst_interval->right;
		mid = .5 * (right + left);
		func_left = worst_interval->func_left;
		func_right = worst_interval->func_right;
		cdf_left = worst_interval->cdf_left;
		cdf_right = worst_interval->cdf_right;
		func_mid = (CDF_OR_PDF == INTERP_CDF ? standard_astable_cdf(mid, alpha, beta) : standard_astable_pdf(mid, alpha, beta));
		cdf_mid = (CDF_OR_PDF == INTERP_CDF ? func_mid : standard_astable_pdf(mid, alpha, beta));
		//printf("here: (%g, %g, %g), (%g, %g, %g), l-r: %g, fl-fr: %g\n", left, mid, right, func_left, func_mid, func_right, right-left, func_right-func_left);

		new_interval[0] = (struct interval *)malloc(sizeof(struct interval));
		new_interval[0]->left = left;
		new_interval[0]->right = mid;
		new_interval[0]->func_left = func_left;
		new_interval[0]->func_right = func_mid;
		new_interval[0]->cdf_left = cdf_left;
		new_interval[0]->cdf_right = cdf_mid;
		pqueue_enqueue(q, new_interval[0]);

		new_interval[1] = (struct interval *)malloc(sizeof(struct interval));
		new_interval[1]->left = mid;
		new_interval[1]->right = right;
		new_interval[1]->func_left = func_mid;
		new_interval[1]->func_right = func_right;
		new_interval[1]->cdf_left = cdf_mid;
		new_interval[1]->cdf_right = cdf_right;
		pqueue_enqueue(q, new_interval[1]);

		free(worst_interval);
	}

	/* Now sort them in increasing values of the left interval point. */
	q_x = pqueue_new(&interval_comp_x, N);
	i = 0;
	/* Take them off the old queue and add them to the new queue. */
	current_interval = NULL;
	while(q->size > 0)
	{
		i++;
		current_interval = (struct interval *)pqueue_dequeue(q);
		pqueue_enqueue(q_x, current_interval);
	}
	i = 0;
	/* Put them on the new queue. */
	while(q_x->size > 0)
	{
		if(i > 0) free(current_interval);
		current_interval = (struct interval *)pqueue_dequeue(q_x);
		x_sorted[i] = current_interval->left;
		y_sorted[i] = current_interval->func_left;
		//if(i > 0 && x_sorted[i] <= x_sorted[i-1]) printf("unsorted (%d): %g, %g -- %g, %g (a,b): (%g, %g)\n", i-1, x_sorted[i], x_sorted[i-1], y_sorted[i], y_sorted[i-1], alpha, beta);
		//else printf("  sorted (%d): %g, %g -- %g, %g (a,b): (%g, %g)\n", i-1, x_sorted[i], x_sorted[i-1], y_sorted[i], y_sorted[i-1], alpha, beta);
		i++;
	}
	x_sorted[i] = current_interval->right;
	y_sorted[i] = current_interval->func_right;
	free(current_interval);

	ip = (struct interpolant_params *)malloc(sizeof(struct interpolant_params));
	ip->alpha = alpha;
	ip->beta = beta;
	ip->ymax = y_sorted[N];
	ip->ymin = y_sorted[0];
	ip->size = N + 1;
	ip->acc = gsl_interp_accel_alloc();
	ip->spline = gsl_spline_alloc(gsl_interp_linear, N + 1);
	gsl_spline_init(ip->spline, x_sorted, y_sorted, N + 1);
	ip->limiting_coefficient = sin(M_PI * alpha / 2.) * gsl_sf_gamma(alpha) / M_PI;
	ip->cdf_or_pdf = CDF_OR_PDF;

	free(x_sorted);
	free(y_sorted);

	pqueue_delete(q);
	pqueue_delete(q_x);

	return ip;
}

void free_interpolant_params(struct interpolant_params *ip)
{
	gsl_spline_free(ip->spline);
	gsl_interp_accel_free(ip->acc);
	free(ip);
}
