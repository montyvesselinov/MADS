// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Daniel O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
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

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../mads.h"
#include "astable/astable.h"

#define NUMITER 10000
#define EPSREL 1.E-7
#define EPSABS 1.E-9

double point_source( double x, double y, double z, double t, void *params );
double rectangle_source( double x, double y, double z, double t, void *params );
double point_source_triangle_time( double x, double y, double z, double t, void *params );
double int_point_source_triangle_time( double tau, void *params );
double rectangle_source_vz( double x, double y, double z, double t, void *params );
double box_source( double x, double y, double z, double t, void *params );
double gaussian_source_2d( double x, double y, double z, double t, void *params );
double gaussian_source_3d( double x, double y, double z, double t, void *params );
double box_source_levy_dispersion( double x, double y, double z, double t, void *params );
double int_point_source( double tau, void *params );
double int_rectangle_source( double tau, void *params );
double int_rectangle_source_vz( double tau, void *params );
double int_box_source( double tau, void *params );
double int_gaussian_source_2d( double tau, void *params );
double int_gaussian_source_3d( double tau, void *params );
double int_box_source_levy_dispersion( double tau, void *params );
// void handler(const char * reason, const char * file, int line, int gsl_errno);
// void handler(const char * reason, const char * file, int line, int gsl_errno) { }

double point_source( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	if( fabs( x - p->var[SOURCE_X] ) < DBL_EPSILON && fabs( y - p->var[SOURCE_Y] ) < DBL_EPSILON && fabs( z - p->var[SOURCE_Z] ) < DBL_EPSILON )
		return( p->var[FLUX] * 1e6 / ( 8 * pow( M_PI, 1.5 ) * p->var[POROSITY] * sqrt( p->var[AX] * p->var[AY] * p->var[AZ] * p->var[VX] * p->var[VX] * p->var[VX] ) ) );
	w = gsl_integration_workspace_alloc( NUMITER );
	time = t - p->var[TIME_INIT];
	F.function = &int_point_source;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
//	printf ("result			 = % .18f\n", result);
//	printf ("estimated error = % .18f\n", error);
//	printf ("intervals =  %d\n", (int) w->size);
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[FLUX] * 1e6 / ( 8 * pow( M_PI, 1.5 ) * p->var[POROSITY] * sqrt( p->var[AX] * p->var[AY] * p->var[AZ] * p->var[VX] * p->var[VX] * p->var[VX] ) ) * result );
}

double int_point_source( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_z = ( p->var[SOURCE_Z] );
	double rx, ry, rz, e1, ez, tau_d, tv, ts, tx, tz1, tz2;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - source_z );
	if( p->scaling_dispersion ) { tau_d = pow( tau, p->var[TSCALE_DISP] ); ts = pow( tau, -1.5 * p->var[TSCALE_DISP] ); }
	else { tau_d = tau; ts = pow( tau, -1.5 ); }
	tv = ( double ) 4 * tau_d * vx;
	rx = tv * ax;
	ry = tv * ay;
	rz = tv * az;
	tx = xe - tau * vx;
	e1 = exp( -tau * lambda - tx * tx / rx - ye * ye / ry );
	//tz1 = ze - source_z;
	//tz2 = ze + source_z;
	tz1 = ze;
	tz2 = ze + 2. * source_z;
	ez = exp( -tz1 * tz1 / rz ) + exp( -tz2 * tz2 / rz );
//	printf("tau %g %g %g %g\n",tau,ez,e1,ze);
	return( e1 * ez * ts );
}

double point_source_triangle_time( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	if( fabs( x - p->var[SOURCE_X] ) < DBL_EPSILON && fabs( y - p->var[SOURCE_Y] ) < DBL_EPSILON && fabs( z - p->var[SOURCE_Z] ) < DBL_EPSILON )
		return( p->var[FLUX] * 1e6 / ( 8 * pow( M_PI, 1.5 ) * p->var[POROSITY] * sqrt( p->var[AX] * p->var[AY] * p->var[AZ] * p->var[VX] * p->var[VX] * p->var[VX] ) ) );
	w = gsl_integration_workspace_alloc( NUMITER );
	time = t - p->var[TIME_INIT];
	F.function = &int_point_source_triangle_time;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
//	printf ("result			 = % .18f\n", result);
//	printf ("estimated error = % .18f\n", error);
//	printf ("intervals =  %d\n", (int) w->size);
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[FLUX] * 1e6 / ( 8 * pow( M_PI, 1.5 ) * p->var[POROSITY] * sqrt( p->var[AX] * p->var[AY] * p->var[AZ] * p->var[VX] * p->var[VX] * p->var[VX] ) ) * result );
}

double int_point_source_triangle_time( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_z = ( p->var[SOURCE_Z] );
	double rx, ry, rz, e1, ez, tau_d, tv, ts, tx, tz1, tz2;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - source_z );
	if( p->scaling_dispersion ) { tau_d = pow( tau, p->var[TSCALE_DISP] ); ts = pow( tau, -1.5 * p->var[TSCALE_DISP] ); }
	else { tau_d = tau; ts = pow( tau, -1.5 ); }
	tv = ( double ) 4 * tau_d * vx;
	rx = tv * ax;
	ry = tv * ay;
	rz = tv * az;
	tx = xe - tau * vx;
	e1 = exp( -tau * lambda - tx * tx / rx - ye * ye / ry );
	//tz1 = ze - source_z;
	//tz2 = ze + source_z;
	tz1 = ze;
	tz2 = ze + 2. * source_z;
	ez = exp( -tz1 * tz1 / rz ) + exp( -tz2 * tz2 / rz );
//	printf("tau %g %g %g %g\n",tau,ez,e1,ze);
	return( e1 * ez * ts );
}


double box_source( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	w = gsl_integration_workspace_alloc( NUMITER );
	F.function = &int_box_source;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 )
	{
		printf( "error: %s\n", gsl_strerror( status ) );
		//result = 0;
	}
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[FLUX] * 1e6 / ( 8. * p->var[POROSITY] * p->var[SOURCE_DX] * p->var[SOURCE_DY] * p->var[SOURCE_DZ] ) * result );
	// return( p->var[C0] * 1e6 / ( p->var[SOURCE_DX] * p->var[SOURCE_DY] * p->var[SOURCE_DZ] ) / ( 8. * p->var[POROSITY] ) * result );
}

double int_box_source( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double source_sizez = ( p->var[SOURCE_DZ] );
	double source_z = ( p->var[SOURCE_Z] );
	double rx, ry, rz, e1, ex, ey, ez, tau_d, tv;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - source_z );
	// if( p->debug >= 3 ) printf( "param %g %g %g %g %g %g %g %.12g %.12g %.12g %.12g\n", d, alpha, beta, xe, ye, x0, y0, p->xe, p->var[SOURCE_X], p->ye, p->var[SOURCE_Y] );
	// printf( "param %g %g %g %g %g %g\n", source_z, source_sizez, ze, p->ze, rz, az );
	tau_d = tau + p->var[NLC0] * sin( tau / p->var[NLC1] );
	if( p->scaling_dispersion ) tau_d = pow( tau_d, p->var[TSCALE_DISP] );
	else tau_d = tau_d;
	tv = ( double ) 4 * tau_d * vx;
	rx = sqrt( tv * ax );
	ry = sqrt( tv * ay );
	rz = sqrt( tv * az );
	e1 = exp( -tau * lambda );
	ex = erfc( ( xe - source_sizex / 2 - tau * vx ) / rx ) - erfc( ( xe + source_sizex / 2 - tau * vx ) / rx );
	ey = erfc( ( ye - source_sizey / 2 ) / ry ) - erfc( ( ye + source_sizey / 2 ) / ry );
	ez = erfc( ( ze - source_sizez ) / rz ) - erfc( ze / rz ) + erfc( ( ze + 2 * source_z ) / rz ) - erfc( ( ze + 2 * source_z + source_sizez ) / rz );
	// ez = erfc( ( ze - source_sizez ) / rz ) - erfc( ( ze + source_sizez ) / rz ) - erfc( ( ze - source_z ) / rz ) + erfc( ( ze + source_z ) / rz );
	// if( p->debug >= 3 ) printf( "int %g %g %g %g %g\n", tau, e1, ex, ey, ez );
	return( e1 * ex * ey * ez );
}

double rectangle_source( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	w = gsl_integration_workspace_alloc( NUMITER );
	F.function = &int_rectangle_source;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	//	printf("result %g ", result, var[C0], p );
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[FLUX] * 1e6 / ( 8. * sqrt( M_PI * p->var[AZ] * p->var[VX] )  * p->var[POROSITY] * p->var[SOURCE_DX] * p->var[SOURCE_DY] ) * result );
}

double int_rectangle_source( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double rx, ry, e1, ex, ey, ez, tau_d, tv;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - p->var[SOURCE_Z] );
	if( p->scaling_dispersion ) tau_d = pow( tau, p->var[TSCALE_DISP] );
	else tau_d = tau;
	tv = ( double ) 4 * tau_d * vx;
	rx = sqrt( tv * ax );
	ry = sqrt( tv * ay );
	e1 = exp( -tau * lambda );
	ez = exp( -ze * ze / ( tau * ( 4. * az * vx ) ) ) + exp( -( ze + 2 * p->var[SOURCE_Z] ) * ( ze + 2 * p->var[SOURCE_Z] ) );
	ex = erfc( ( xe - source_sizex / 2 - tau * vx ) / rx ) - erfc( ( xe + source_sizex / 2 - tau * vx ) / rx );
	ey = erfc( ( ye - source_sizey / 2 ) / ry ) - erfc( ( ye + source_sizey / 2 ) / ry );
	return( e1 * ex * ey * ez / sqrt( tau_d ) );
}

double rectangle_source_vz( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	w = gsl_integration_workspace_alloc( NUMITER );
	F.function = &int_rectangle_source_vz;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[FLUX] * 1e6 * p->var[VZ] / ( 4. * p->var[POROSITY] ) * result );
}

double int_rectangle_source_vz( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double vz = ( p->var[VZ] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double rx, ry, rz, e1, ex, ey, ez, tz, tv, tau_d, v;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - p->var[SOURCE_Z] );
	v = sqrt( vx * vx + vz * vz );
	if( p->scaling_dispersion ) tau_d = pow( tau, p->var[TSCALE_DISP] );
	else tau_d = tau;
	tv = ( double ) 4 * tau_d * v;
	rx = sqrt( tv * ax );
	ry = sqrt( tv * ay );
	rz = sqrt( tv * az );
	e1 = exp( -tau * lambda );
	ex = erfc( ( xe - source_sizex / 2 - tau * vx ) / rx ) - erfc( ( xe + source_sizex / 2 - tau * vx ) / rx );
	ey = erfc( ( ye - source_sizey / 2 ) / ry ) - erfc( ( ye + source_sizey / 2 ) / ry );
	tz = ze - tau * vz;
	ez = exp( -tz * tz / ( tau_d * ( 4. * az * v ) ) ) / sqrt( tau_d * ( M_PI * az * v ) ) -
		 vz / ( 2. * az * v ) * exp( vz * ze / ( az * v ) ) * erfc( ( ze + tau * vz ) / rz );
	return( e1 * ex * ey * ez );
}

double gaussian_source_2d( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	w  = gsl_integration_workspace_alloc( NUMITER );
	F.function = &int_gaussian_source_2d;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	//	printf("result %g ", result, var[C0], p );
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return p->var[FLUX] * 1e6 * result / ( p->var[POROSITY] * sqrt( M_PI * M_PI * M_PI * ( 2. * ( 2. * p->var[AX] * p->var[VX] + p->var[SOURCE_DX] * p->var[SOURCE_DX] ) ) * ( 2. * ( 2. * p->var[AY] * p->var[VX] + p->var[SOURCE_DY] * p->var[SOURCE_DY] ) ) * 4. * p->var[AZ] * p->var[VX] ) );
}

double int_gaussian_source_2d( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double d, alpha, beta, xe, ye, ze, x0, y0;
	double sx, sy, ez, tau_d, tv, ts;
	double varx, vary, varz;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - p->var[SOURCE_Z] );
	if( p->scaling_dispersion ) { tau_d = pow( tau, p->var[TSCALE_DISP] ); ts = pow( tau, -1.5 * p->var[TSCALE_DISP] ); }
	else { tau_d = tau; ts = pow( tau, -1.5 ); }
	tv = ( double ) 2 * tau_d * vx;
	varx = tv * ax + source_sizex * source_sizex;
	vary = tv * ay + source_sizey * source_sizey;
	varz = tv * az;
	sx =  -( xe - vx * tau ) * ( xe - vx * tau ) / ( 2 * varx );
	sy = -ye * ye / ( 2 * vary );
	ez = ( exp( -ze * ze / ( 2 * varz ) ) + exp( -( ze + 2 * p->var[SOURCE_Z] ) * ( ze + 2 * p->var[SOURCE_Z] ) / ( 2 * varz ) ) );
	return exp( -tau * lambda + sx + sy ) * ez * ts;
}

double gaussian_source_3d( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	w = gsl_integration_workspace_alloc( NUMITER );
	F.function = &int_gaussian_source_3d;
	F.params = p;
	gsl_set_error_handler_off();
	// gsl_set_error_handler(&handler);
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	//	printf("result %g ", result, var[C0], p );
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return p->var[FLUX] * 1e6 * result / ( p->var[POROSITY] * sqrt( M_PI * M_PI * M_PI * ( 2 * ( 2 * p->var[AX] * p->var[VX] + p->var[SOURCE_DX] * p->var[SOURCE_DX] ) ) * ( 2 * ( 2 * p->var[AY] * p->var[VX] + p->var[SOURCE_DY] * p->var[SOURCE_DY] ) ) * 2 * ( 2 * p->var[AZ] * p->var[VX] + p->var[SOURCE_DZ] * p->var[SOURCE_DZ] ) ) );
}

double int_gaussian_source_3d( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double source_sizez = ( p->var[SOURCE_DZ] );
	double d, alpha, beta, xe, ye, ze, x0, y0;
	double sx, sy, ez, tau_d, tv, ts;
	double varx, vary, varz;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - p->var[SOURCE_Z] );
	if( p->scaling_dispersion ) { tau_d = pow( tau, p->var[TSCALE_DISP] ); ts = pow( tau, -1.5 * p->var[TSCALE_DISP] ); }
	else { tau_d = tau; ts = pow( tau, -1.5 ); }
	tv = ( double ) 2 * tau_d * vx;
	varx = tv * ax + source_sizex * source_sizex;
	vary = tv * ay + source_sizey * source_sizey;
	varz = tv * az + source_sizez * source_sizez;
	sx = -( xe - vx * tau ) * ( xe - vx * tau ) / ( 2 * varx );
	sy = -ye * ye / ( 2 * vary );
	ez =  exp( sx + sy - tau * lambda - ze * ze / ( 2 * varz ) ) + exp( sx + sy - tau * lambda - ( ze + 2 * p->var[SOURCE_Z] ) * ( ze + 2 * p->var[SOURCE_Z] ) / ( 2 * varz ) ) ;
	return ez * ts;
}

double box_source_levy_dispersion( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct anal_data *p = ( struct anal_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	w = gsl_integration_workspace_alloc( NUMITER );
	F.function = &int_box_source_levy_dispersion;
	F.params = p;
	gsl_set_error_handler_off();
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 )
	{
		printf( "error: %s\n", gsl_strerror( status ) );
		//result = 0;
	}
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[FLUX] * 1e6 / ( 8. * p->var[POROSITY] * p->var[SOURCE_DX] * p->var[SOURCE_DY] * p->var[SOURCE_DZ] ) * result );
	// return( p->var[C0] * 1e6 / ( p->var[SOURCE_DX] * p->var[SOURCE_DY] * p->var[SOURCE_DZ] ) / ( 8. * p->var[POROSITY] ) * result );
}

double int_box_source_levy_dispersion( double tau, void *params )
{
	struct anal_data *p = ( struct anal_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double source_sizez = ( p->var[SOURCE_DZ] );
	double source_z = ( p->var[SOURCE_Z] );
	double alpha, beta;
	double lambda_x, lambda_y, lambda_z;
	double tau_d, tv;
	double angle, rot1, rot2, xe, ye, ze, x0, y0;
	double decay_factor, px, py, pz;
	alpha = p->var[ALPHA];
	beta = p->var[BETA];
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	angle = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	rot1 = cos( angle );
	rot2 = sin( angle );
	xe = x0 * rot1 - y0 * rot2;
	ye = x0 * rot2  + y0 * rot1;
	ze = ( p->ze - source_z );
	if( p->scaling_dispersion ) tau_d = pow( tau, p->var[TSCALE_DISP] );
	else tau_d = tau;
	tv = ( double ) 2 * tau_d * vx;
	lambda_x = pow( tv * ax, 1 / alpha );
	lambda_y = pow( tv * ay, 1 / alpha );
	lambda_z = pow( tv * az, 1 / alpha );
	decay_factor = exp( -tau * lambda );
	px = astable_cdf( xe + source_sizex / 2. - vx * tau, alpha, beta, 0., lambda_x ) - astable_cdf( xe - source_sizex / 2. - vx * tau, alpha, beta, 0., lambda_x );
	py = astable_cdf( ye + source_sizey / 2., alpha, beta, 0., lambda_y ) - astable_cdf( ye - source_sizey / 2., alpha, beta, 0., lambda_y );
	pz = astable_cdf( ze + source_sizez, alpha, beta, 0., lambda_z ) - astable_cdf( ze, alpha, beta, 0., lambda_z );
	// ez = erfc( ( ze - source_sizez ) / rz ) - erfc( ( ze + source_sizez ) / rz ) - erfc( ( ze - source_z ) / rz ) + erfc( ( ze + source_z ) / rz );
	// if( p->debug >= 3 ) printf( "int %g %g %g %g %g\n", tau, e1, ex, ey, ez );
	return( decay_factor * px * py * pz );
}
