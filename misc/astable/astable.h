/*
 * astable.h
 *
 *  Created on: Jul 11, 2013
 *      Author: Daniel O'Malley (omalled@lanl.gov)
 */

#ifndef ASTABLE_H_
#define ASTABLE_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>

#define INTERP_CDF 0
#define INTERP_PDF 1

struct integrand_params
{
	double alpha;
	double beta;
	double zeta;
	double theta0;
	double c1;
	double c2;
	double c3;
	double x;
};

struct interpolant_params
{
	double alpha;
	double beta;
	gsl_interp_accel *acc;
	gsl_spline *spline;
	double ymax;
	double ymin;
	int size;
	double limiting_coefficient;
	int cdf_or_pdf;
};

//Defined in astable.c
double astable_pdf(double x, double alpha, double beta, double gamma, double lambda);
double astable_cdf(double x, double alpha, double beta, double gamma, double lambda);
void setup_params(struct integrand_params *params, double x, double alpha, double beta);
double pdf_integrand(double theta, void *params);
double standard_astable_pdf(double x, double alpha, double beta);
double cdf_integrand(double theta, void *params);
double standard_astable_cdf(double x, double alpha, double beta);
inline double nolan_g(double theta, void *params);
inline double nolan_log_g(double theta, void *params);
inline double nolan_g_minus_one(double theta, void *params);
inline double exp_g_minus_half(double theta, void *params);
inline double nolan_log_V(double theta, void *params);
inline double nolan_V(double theta, void *params);
double bisection_solver(double (*f)(double, void *), double value, double a, double b, void *params, double tol);

//Defined in interpolation.c
void astable_cdf_interp(double x, double alpha, double beta, double gamma, double lambda, double *val);
struct interpolant_params *automate_interpolant(double alpha, double beta, double percentile, double abserr, int CDF_OR_PDF);
double interpolate(double x, double gamma, double lambda, struct interpolant_params *ip);
struct interpolant_params *setup_interpolant(double alpha, double beta, double left, double right, int N, int CDF_OR_PDF);
void free_interpolant_params(struct interpolant_params *ip);
int interval_comp(const void *d1, const void *d2);
int interval_comp_x(const void *d1, const void *d2);

#endif /* ASTABLE_H_ */
