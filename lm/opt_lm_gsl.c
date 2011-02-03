#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "../mads.h"

#define FIT(i) gsl_vector_get(solver->x, i)
#define STDDEV(i) sqrt(gsl_matrix_get(covar,i,i))
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
int optimize_gsl( struct opt_data *op );
int func_gsl( const gsl_vector *x, void *data, gsl_vector *f );
int func_gsl_dx( const gsl_vector *x, void *data, gsl_matrix *jacobian );
int func_gsl_xdx( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J );
double func_gsl_deriv( double x, void *data );
int func_gsl_deriv_dx( const gsl_vector *x, void *data, gsl_matrix *J );
int func_gsl_deriv_xdx( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J );

/* Functions elsewhere */
double point_source( double x, double y, double z, double t, void *params );
double rectangle_source( double x, double y, double z, double t, void *params );
double rectangle_source_vz( double x, double y, double z, double t, void *params );
double box_source( double x, double y, double z, double t, void *params );
double int_point_source( double tau, void *params );
double int_rectangle_source( double tau, void *params );
double int_rectangle_source_vz( double tau, void *params );
double int_box_source( double tau, void *params );
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );

int lm_gsl( gsl_vector *opt_params, struct opt_data *op, gsl_matrix *jacobian, gsl_matrix *covar )
{
	int i, j, k, dof, iter = 0, gsl_deriv = 0;
	int status_i, status_d, status_g;
	const gsl_multifit_fdfsolver_type *solver_type;
	gsl_multifit_fdfsolver *solver;
	gsl_vector *opt_gradient;
	double *x_c, chi, eps, stddev_scale;
	gsl_multifit_function_fdf f;
	op->pd->var_current_gsl = gsl_vector_alloc( op->pd->nOptParam );
	op->od->obs_current_gsl = gsl_vector_alloc( op->od->nObs );
	opt_gradient = gsl_vector_alloc( op->pd->nOptParam );
	if(( x_c = ( double * ) malloc( op->pd->nParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); exit( 1 ); }
	f.f = &func_gsl; /* forward run only */
	if( strstr( op->cd->opt_method, "deriv" ) == NULL )
	{
		f.df = &func_gsl_dx; /* Forward numerical derivatives for GSL */
		f.fdf = &func_gsl_xdx; /* Forward numerical derivatives for GSL */
		printf( "Forward numerical derivatives\n" );
	}
	else
	{
		f.df = &func_gsl_deriv_dx; /* Forward numerical derivatives for GSL using GSL functions */
		f.fdf = &func_gsl_deriv_xdx; /* Forward numerical derivatives for GSL using GSL functions */
		printf( "Forward numerical derivatives for GSL using GSL functions\n" );
		gsl_deriv = 1;
	}
	f.n = op->od->nObs;
	f.p = op->pd->nOptParam;
	f.params = op;
	if( strstr( op->cd->opt_method, "uns" ) == NULL ) { printf( "Levenberg-Marquardt scalled GSL version (diagonal elements are not modified)\n" ); solver_type = gsl_multifit_fdfsolver_lmsder; } // LM method Scalled version
	else { printf( "Levenberg–Marquardt unscalled GSL version (diagonal elements set to 1)\n" ); solver_type = gsl_multifit_fdfsolver_lmder; } // LM method Unscalled version
	solver = gsl_multifit_fdfsolver_alloc( solver_type, op->od->nObs, op->pd->nOptParam );
	gsl_multifit_fdfsolver_set( solver, &f, opt_params );
	chi = gsl_blas_dnrm2( solver->f );
	printf( "Iteration %i: chi %g phi %g\n", iter, chi, chi * chi );
	do
	{
		iter++;
		do
		{
			status_i = gsl_multifit_fdfsolver_iterate( solver );
			if( op->cd->ldebug && status_i != GSL_SUCCESS ) printf( "Iteration status = %s\n", gsl_strerror( status_i ) );
		}
		while( status_i == GSL_CONTINUE );
		chi = gsl_blas_dnrm2( solver->f );
		if( op->cd->ldebug ) printf( "Iteration %i: chi %g phi %g\n", iter, chi, chi * chi );
		if( op->cd->ldebug >= 3 )
		{
			dof = op->od->nObs - op->pd->nOptParam;
			stddev_scale = chi / sqrt( dof );
			gsl_multifit_covar( solver->J, 0.0, covar );
			DeTransform( solver->x->data, op, x_c );
			for( i = 0; i < op->pd->nOptParam; i++ )
				printf( "%-40s : %g -> %g = %g +/- %g\n", op->pd->var_id[op->pd->var_index[i]], op->pd->var[op->pd->var_index[i]], x_c[i], FIT( i ), stddev_scale * STDDEV( i ) );
			for( k = 0, i = 0; i < op->wd->nW; i++ )
				for( j = 0; j < op->wd->nWellObs[i]; j++ )
					if( op->wd->obs_weight[i][j] > 0 )
					{
						printf( "%s(%g): %f\n", op->wd->id[i], op->wd->obs_time[i][j], solver->f->data[k] );
						k++;
					}
			printf( "Jacobian matrix\n" );
			printf( "%-25s :", "Observations" );
			for( k = 0; k < op->od->nObs; k++ )
			{
				i = op->od->well_index[k];
				j = op->od->time_index[k];
				printf( " %s(%g)", op->wd->id[i], op->wd->obs_time[i][j] );
			}
			printf( "\n" );
			for( k = i = 0; i < op->pd->nOptParam; i++ )
			{
				printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->od->nObs; j++, k++ )
				{
					eps = gsl_matrix_get( solver->J, j, i );
					if( fabs( eps ) > 1e3 ) printf( " %6.0e", eps );
					else printf( " %6.2f", eps );
				}
				printf( "\n" );
			}
		}
		status_d = gsl_multifit_test_delta( solver->dx, solver->x, 1e-3, 1e-2 );
		gsl_multifit_gradient( solver->J, solver->f, opt_gradient );
		status_g = gsl_multifit_test_gradient( opt_gradient, 1e-1 );
		if(( status_d == GSL_SUCCESS || status_g == GSL_SUCCESS ) )
		{
			f.df = &func_gsl_deriv_dx; // Numerical derivatives for GSL using GSL functions
			f.fdf = &func_gsl_deriv_xdx; // Numerical derivatives for GSL using GSL functions
			if( gsl_deriv == 0 ) { printf( "Forward numerical derivatives for GSL using GSL functions\n" ); gsl_deriv = 1; }
			else break;
		}
		if( op->cd->ldebug && iter >= op->cd->niter ) { printf( "Maximum number of iterations is exceeded (%d)\n", op->cd->niter ); break; }
		if( op->cd->ldebug && op->cd->eval >= op->cd->maxeval ) { printf( "Maximum number of evaluations is exceeded (%d)\n", op->cd->eval ); break; }
	}
	while( 1 );
	if( op->cd->ldebug ) printf( "Delta test status = %s\n", gsl_strerror( status_d ) );
	if( op->cd->ldebug ) printf( "Gradient test status = %s\n", gsl_strerror( status_g ) );
	gsl_matrix_memcpy( jacobian, solver->J );
	gsl_vector_memcpy( opt_params, solver->x );
	gsl_multifit_covar( solver->J, 0.0, covar );
	chi = gsl_blas_dnrm2( solver->f );
	gsl_multifit_fdfsolver_free( solver );
	gsl_vector_free( opt_gradient );
	gsl_vector_free( op->pd->var_current_gsl );
	gsl_vector_free( op->od->obs_current_gsl );
	free( x_c );
	op->phi = chi * chi;
	return( 1 );
}

int func_gsl( const gsl_vector *x, void *data, gsl_vector *f ) /* forward run for GSL */
{
	func( x->data, data, f->data );
	return GSL_SUCCESS;
}

double func_gsl_deriv( double x, void *data ) /* forward run for GSL deriative computation */
{
	struct opt_data *p = ( struct opt_data * )data;
	double r;
	gsl_vector_set( p->pd->var_current_gsl, p->cd->pderiv, x );
	func_gsl( p->pd->var_current_gsl, p, p->od->obs_current_gsl );
	r = gsl_vector_get( p->od->obs_current_gsl, p->cd->oderiv );
//	printf( "x %g result %g param %d obs %d\n", x, r, p->cd->pderiv, p->cd->oderiv );
	return( r );
}

int func_gsl_dx( const gsl_vector *x, void *data, gsl_matrix *jacobian ) /* Simple numerical derivatives for GSL */
{
	struct opt_data *p = ( struct opt_data * )data;
	gsl_vector *f_x, *f_xpdx, *x_c;
	double x_old, dx;
	int i, j;

	x_c = gsl_vector_alloc( p->pd->nOptParam );
	gsl_vector_memcpy( x_c, x );
	f_xpdx = gsl_vector_alloc( p->od->nObs );
	f_x = gsl_vector_alloc( p->od->nObs );
	func_gsl( x_c, data, f_x );
	for( j = 0; j < p->pd->nOptParam; j++ )
	{
		x_old = gsl_vector_get( x_c, j );
		if( p->cd->sintrans == 0 ) dx = p->pd->var_dx[j];
		else dx = p->cd->sindx;
		gsl_vector_set( x_c, j, x_old + dx );
		func_gsl( x_c, data, f_xpdx );
		gsl_vector_set( x_c, j, x_old );
		for( i = 0; i < p->od->nObs; i++ ) gsl_matrix_set( jacobian, i, j, ( gsl_vector_get( f_xpdx, i ) - gsl_vector_get( f_x, i ) ) / dx ); // first obs, second param
	}
	gsl_vector_free( f_x ); gsl_vector_free( f_xpdx ); gsl_vector_free( x_c );
	return GSL_SUCCESS;
}

int func_gsl_deriv_dx( const gsl_vector *x, void *data, gsl_matrix *J ) /* Numerical derivatives for GSL using GSL functions */
{
	gsl_function F;
	struct opt_data *p = ( struct opt_data * )data;
	double x_old, result, abserr, dx;
	int i, j;
//	printf("deriv xdx\n");
	F.function = &func_gsl_deriv;
	F.params = data;
	for( j = 0; j < p->pd->nOptParam; j++ )
	{
		gsl_vector_memcpy( p->pd->var_current_gsl, x );
//		printf( "Param %s\n", p->pd->var_id[p->pd->var_index[j]] );
		p->cd->pderiv = j;
		x_old = gsl_vector_get( x, j );
		if( p->cd->sintrans == 0 ) dx = p->pd->var_dx[j];
		else dx = p->cd->sindx;
		for( i = 0; i < p->od->nObs; i++ )
		{
			p->cd->oderiv = i;
			gsl_deriv_forward( &F, x_old, dx, &result, &abserr ); // Avoid using central when SIN transformation is applied
			gsl_matrix_set( J, i, j, result ); // first obs, second param
//			printf("grad param %s obs %d val %g %g grad %g err %g\n", p->pd->var_id[p->pd->var_index[j]], i, x_old, dx, result, abserr );
//			printf("grad param %s obs %d grad %g err %g\n", p->pd->var_id[p->pd->var_index[j]], i, result, abserr );
		}
	}
	p->cd->pderiv = p->cd->oderiv = -1;
	return GSL_SUCCESS;
}

int func_gsl_xdx( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J )
{
	func_gsl( x, data, f );
	func_gsl_dx( x, data, J ); /* Simple numerical derivatives for GSL */
	return GSL_SUCCESS;
}

int func_gsl_deriv_xdx( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J )
{
	func_gsl( x, data, f );
	func_gsl_deriv_dx( x, data, J ); /* Numerical derivatives for GSL using GSL functions */
	return GSL_SUCCESS;
}
