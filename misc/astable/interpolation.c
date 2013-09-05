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
#include <gsl/gsl_sf.h>

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
	//TODO: find the expressions for when \beta=-1,x\to\infty and \beta=1,x\to-\infty
	if(x <= ip->xa[ip->size] && x >= ip->xa[0])
	{
		return gsl_spline_eval(ip->spline, (x - gamma) / lambda, ip->acc) / lambda;
	}
	else if(x > ip->xa[ip->size])
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
			return ip->ya[ip->size];
		}
	}
	else
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
			return ip->ya[0];
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
	double *x;
	double *y;
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

	if(left > 0) left = 0.;//We need the lower end point to be <= 0 in order for the limiting behavior code to work.

	q = pqueue_new(&interval_comp, N + 1);

	x = (double *)malloc((N + 2) * sizeof(double));
	y = (double *)malloc((N + 2) * sizeof(double));
	x_sorted = (double *)malloc((N + 2) * sizeof(double));
	y_sorted = (double *)malloc((N + 2) * sizeof(double));

	ip = (struct interpolant_params *)malloc(sizeof(struct interpolant_params));

	ip->alpha = alpha;
	ip->beta = beta;
	ip->acc = gsl_interp_accel_alloc();
	ip->spline = gsl_spline_alloc(gsl_interp_linear, N + 2);

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
	x[0] = left;
	x[1] = right;
	y[0] = func_left;
	y[1] = func_right;
	worst_interval->func_right = func_right;
	worst_interval->func_left = func_left;
	worst_interval->cdf_left = cdf_left;
	worst_interval->cdf_right = cdf_right;

	pqueue_enqueue(q, worst_interval);

	for(i = 0; i < N; i++)
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
		cdf_mid = (CDF_OR_PDF == INTERP_CDF ? func_mid : standard_astable_cdf(mid, alpha, beta));

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
	q_x = pqueue_new(&interval_comp_x, N + 1);
	i = 0;
	current_interval = NULL;//Tell the compiler to stop warning me.
	/* Take them off the old queue and add them to the new queue. */
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
		i++;
	}
	x_sorted[i] = current_interval->right;
	y_sorted[i] = current_interval->func_right;
	if(current_interval != NULL) free(current_interval);

	ip->xa = x_sorted;
	ip->ya = y_sorted;
	ip->size = N + 2;
	gsl_spline_init(ip->spline, ip->xa, ip->ya, ip->size);
	ip->limiting_coefficient = sin(M_PI * alpha / 2.) * gsl_sf_gamma(alpha) / M_PI;
	ip->cdf_or_pdf = CDF_OR_PDF;

	free(x);
	free(y);
	pqueue_delete(q);
	pqueue_delete(q_x);

	return ip;
}

void free_interpolant_params(struct interpolant_params *ip)
{
	gsl_spline_free(ip->spline);
	gsl_interp_accel_free(ip->acc);
	free(ip->xa);
	free(ip->ya);
	free(ip);
}
