/*
 * pdfs.c
 *
 *  Created on: Jul 11, 2013
 *      Author: Daniel O'Malley (omalled@lanl.gov)
 *
 *      All parameters are in Zolotarev's M parameterization.
 *      alpha - stability parameter
 *      beta - skewness parameter
 *      gamma - shift parameter
 *      lambda - scale parameter
 */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "astable.h"

#define MAXITER 1000
#define EPSREL 1e-10
#define EPSABS 1e-10

void setup_params(struct integrand_params *params, double x, double alpha, double beta)
{
	params->alpha = alpha;
	params->beta = beta;
	params->x = x;
	if(alpha == 1)
	{
		params->zeta = 0.;
		params->theta0 = M_PI / 2;
		params->c1 = 0.;
		params->c2 = .5 / fabs(beta);
		params->c3 = 1. / M_PI;
	}
	else
	{
		params->zeta = -beta * tan(alpha * M_PI / 2.);
		params->theta0 = atan(beta * tan(alpha * M_PI / 2)) / alpha;
		if(alpha < 1)
		{
			params->c1 = .5 - params->theta0 / M_PI;
			params->c3 = 1. / M_PI;
		}
		else
		{
			params->c1 = 1.;
			params->c3 = -1. / M_PI;
		}
		params->c2 = alpha / (M_PI * fabs(alpha - 1) * (x - params->zeta));
	}

	return;
}

double astable_pdf(double x, double alpha, double beta, double gamma, double lambda)
{
	return standard_astable_pdf((x - gamma) / lambda, alpha, beta) / lambda;
}

double astable_cdf(double x, double alpha, double beta, double gamma, double lambda)
{
	return standard_astable_cdf((x - gamma) / lambda, alpha, beta);
}

//A "standard" a-stable RV has
double standard_astable_pdf(double x, double alpha, double beta)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXITER);
	gsl_function f;
	int status;
	struct integrand_params params;
	double result1;
	double result2;
	double error;
	double high;
	double low;
	double theta2;

	status = 0;

	setup_params(&params, x, alpha, beta);

	if(alpha == 2.)
	{
		return exp(-x * x / 4.) / sqrt(4. * M_PI);
	}
	if(alpha == 1. && beta == 0.)
	{
		return 1. / (M_PI * (1. + x * x));
	}
	else if(alpha != 1 && fabs(x - params.zeta) < EPSABS)
	{
		return exp(lgamma(1 + 1 / params.alpha)) * cos(params.theta0) / (M_PI * pow(1 + params.zeta * params.zeta, .5 / params.alpha));
	}
	else if(alpha != 1 && x < params.zeta) return standard_astable_pdf(-x, alpha, -beta);
	else if(alpha == 1 && beta < 0) return standard_astable_pdf(-x, alpha, -beta);
	else
	{
		low = -params.theta0;
		high = low + .5 * (M_PI / 2. - low);
		theta2 = bisection_solver(&nolan_log_g, 0, low, high, (void *)&params, EPSABS);

		f.function = &pdf_integrand;
		f.params = &params;

		if(theta2 + params.theta0 < EPSABS) result1 = 0.;
		else status = gsl_integration_qag(&f, -params.theta0, theta2, EPSABS, EPSREL, MAXITER, GSL_INTEG_GAUSS61, w, &result1, &error);
		if(status != 0)
		{
			fprintf(stderr, "Error with the integration.\n");
		}
		if(M_PI / 2. - theta2 < EPSABS) result2 = 0.;
		else status = gsl_integration_qag(&f, theta2, M_PI / 2, EPSABS, EPSREL, MAXITER, GSL_INTEG_GAUSS61, w, &result2, &error);
		if(status != 0)
		{
			fprintf(stderr, "Error with the integration.\n");
		}
	}
	gsl_integration_workspace_free(w);

	return params.c2 * (result1 + result2);
}

void print_params(struct integrand_params *p)
{
	printf("a: %g, b %g, z: %g, t0 %g, x: %g\n", p->alpha, p->beta, p->zeta, p->theta0, p->x);
}

double standard_astable_cdf(double x, double alpha, double beta)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXITER);
	gsl_function f;
	int status;
	struct integrand_params params;
	double result1;
	double result2;
	double error;
	double high;
	double low;
	double theta2;

	status = 0;

	setup_params(&params, x, alpha, beta);

	if(alpha == 2.)
	{
		return .5 * (1. + erf(x / 2));
	}
	else if(alpha == 1. && beta == 0.)
	{
		return 0.5 + atan(x) / M_PI;
	}
	else if(fabs(x - params.zeta) < EPSABS) return 0.5 - params.theta0 / M_PI;
	else if(alpha != 1. && x < params.zeta) return 1. - standard_astable_cdf(-x, alpha, -beta);
	else if(alpha == 1 && beta < 0) return 1. - standard_astable_cdf(-x, alpha, -beta);
	else
	{
		low = -params.theta0;
		high = low + .5 * (M_PI / 2. - low);
		theta2 = bisection_solver(&nolan_log_g, log(log(2)), low, high, (void *)&params, EPSABS);

		f.function = &cdf_integrand;
		f.params = &params;

		if(theta2 + params.theta0 < EPSABS) result1 = 0.;
		else status = gsl_integration_qag(&f, -params.theta0, theta2, EPSABS, EPSREL, MAXITER, GSL_INTEG_GAUSS61, w, &result1, &error);
		if(status != 0)
		{
			fprintf(stderr, "Error with the integration.\n");
		}
		if(M_PI / 2. - theta2 < EPSABS) result2 = 0.;
		else status = gsl_integration_qag(&f, theta2, M_PI / 2, EPSABS, EPSREL, MAXITER, GSL_INTEG_GAUSS61, w, &result2, &error);
		if(status != 0)
		{
			fprintf(stderr, "Error with the integration.\n");
		}
	}
	gsl_integration_workspace_free(w);

	return params.c1 + params.c3 * (result1 + result2);
}

double pdf_integrand(double theta, void *params)
{
	double g;

	g = nolan_g(theta, params);

	if(isnan(exp(-g)) || isinf(exp(-g)))
	{
		print_params(params);
	}
	if(isinf(g) == 1) return 0.;
	else return g * exp(-g);
}

double cdf_integrand(double theta, void *params)
{
	double g;

	g = nolan_g(theta, params);

	if(isnan(exp(-g)) || isinf(exp(-g)))
	{
		print_params(params);
		printf("%g: %g %g -- %g %g %g\n", theta, g, exp(-g), nolan_V(theta, params), cos(theta), tan(theta));
	}

	return exp(-g);
}

//Solves for f(x)=value
double bisection_solver(double (*f)(double, void *), double value, double a, double b, void *params, double tol)
{
	double left;
	double mid;

	left = (*f)(a, params);
	while(fabs(a - b) > tol)
	{
		mid = (*f)(.5 * (a + b), params);
		if((left < value && mid < value) || (left > value && mid > value))
		{
			a = .5 * (a + b);
			left = mid;
		}
		else
		{
			b = .5 * (a + b);
		}
	}

	return .5 * (a + b);
}

inline double nolan_g(double theta, void *params)
{
	return exp(nolan_log_g(theta, params));
}

inline double nolan_log_g(double theta, void *params)
{
	struct integrand_params *p;
	double log_V;
	double log_g = -1;

	p = (struct integrand_params *)params;

	log_V = nolan_log_V(theta, params);
	if(p->alpha != 1.)
	{
		log_g = log_V + p->alpha / (p->alpha - 1) * log(p->x - p->zeta);
	}
	else if(p->beta != 0)
	{
		log_g = -M_PI * p->x / 2 / p->beta + log_V;
	}

	return log_g;
}

inline double nolan_g_minus_one(double theta, void *params)
{
	return nolan_g(theta, params) - 1.;
}

inline double nolan_log_V(double theta, void *params)
{
	struct integrand_params *p;
	double t1, t2, t3, t4;
	int inf;
	double log_V = -1.;

	p = (struct integrand_params *)params;

	if(p->alpha != 1)
	{
		t1 = log(fabs(cos(p->alpha * p->theta0)));
		t2 = log(fabs(cos(theta)));
		t3 = log(fabs(sin(p->alpha * (p->theta0 + theta))));
		t4 = log(fabs(cos(p->alpha * p->theta0 + (p->alpha - 1) * theta)));
		inf = 1 / (p->alpha - 1) * (isinf(t1) + isinf(t2) - p->alpha * isinf(t3)) + isinf(t4);
		if(inf < 0) log_V = -1. / 0.;
		else if(inf > 0) log_V = 1. / 0.;
		else log_V = 1 / (p->alpha - 1) * (t1 + t2 - p->alpha * t3) + t4;
	}
	else if(p->beta != 0.)
	{
		log_V = log(2 / M_PI * (M_PI / 2 + p->beta * theta)) - log(cos(theta)) + 1 / p->beta * (M_PI / 2 + p->beta * theta) * tan(theta);
		//printf("log_v: %g %g %g\n", log(2 / M_PI * (M_PI / 2 + p->beta * theta)), log(cos(theta)), 1 / p->beta * (M_PI / 2 + p->beta * theta) * tan(theta));
	}

	return log_V;
}

inline double nolan_V(double theta, void *params)
{
	struct integrand_params *p;
	double V = -1.;

	p = (struct integrand_params *)params;

	if(p->alpha != 1)
	{
		V = pow(fabs(cos(p->alpha * p->theta0)), 1 / (p->alpha - 1)) * pow(fabs(cos(theta) / sin(p->alpha * (p->theta0 - theta))), p->alpha / (p->alpha - 1)) * cos(p->alpha * p->theta0 + (p->alpha - 1) * theta) / cos(theta);
	}
	else if(p->beta != 0.)
	{
		V = 2 / M_PI * (M_PI / 2 + p->beta * theta) / cos(theta) * exp(1 / p->beta * (M_PI / 2 + p->beta * theta) * tan(theta));
	}

	return V;
}
