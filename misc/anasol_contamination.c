#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../mads.h"

#define NUMITER 10000
#define EPSREL 1e-7
#define EPSABS 0

double point_source( double x, double y, double z, double t, void *params );
double rectangle_source( double x, double y, double z, double t, void *params );
double rectangle_source_vz( double x, double y, double z, double t, void *params );
double box_source( double x, double y, double z, double t, void *params );
double int_point_source( double tau, void *params );
double int_rectangle_source( double tau, void *params );
double int_rectangle_source_vz( double tau, void *params );
double int_box_source( double tau, void *params );

double point_source( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( NUMITER );
	gsl_function F;
	int status;
	struct calc_data *p = ( struct calc_data * )params;
	double result, error, time;
	if( t <= p->var[TIME_INIT] ) return( 0 );
	p->xe = x;
	p->ye = y;
	p->ze = z;
	time = t - p->var[TIME_INIT];
	F.function = &int_point_source;
	F.params = p;
	gsl_set_error_handler_off();
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
//	printf ("result			 = % .18f\n", result);
//	printf ("estimated error = % .18f\n", error);
//	printf ("intervals =  %d\n", w->size);
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[C0] * 1e6 / ( 44.546623974 * p->var[POROSITY] * sqrt( p->var[AX] * p->var[AY] * p->var[AZ] * p->var[VX] * p->var[VX] * p->var[VX] ) ) * result );
}

double int_point_source( double tau, void *params )
{
	struct calc_data *p = ( struct calc_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_z = ( p->var[SOURCE_Z] );
	double rx, ry, rz, e1, ez, tx, tz1, tz2;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - source_z );
	rx = 4.*tau * ax * vx;
	ry = 4.*tau * ay * vx;
	rz = 4.*tau * az * vx;
	tx = xe - tau * vx;
	e1 = exp( -tau * lambda - tx * tx / rx - ye * ye / ry );
	tz1 = ze - source_z;
	tz2 = ze + source_z;
	ez = exp( -tz1 * tz1 / rz ) + exp( -tz2 * tz2 / rz );
//	printf("tau %g %g %g %g\n",tau,ez,e1,ze);
	return( e1 * ez * pow( tau, -1.5 ) );
}

double box_source( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w;
	gsl_function F;
	int status;
	struct calc_data *p = ( struct calc_data * )params;
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
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[C0] * 1e6 / ( p->var[SOURCE_DX] * p->var[SOURCE_DY] * p->var[SOURCE_DZ] ) / ( 8. * p->var[POROSITY] ) * result );
}

double int_box_source( double tau, void *params )
{
	struct calc_data *p = ( struct calc_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double source_sizez = ( p->var[SOURCE_DZ] );
	double source_z = ( p->var[SOURCE_Z] );
	double rx, ry, rz, e1, ex, ey, ez;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - source_z );
//	if( p->debug >= 3 ) printf( "param %g %g %g %g %g %g %g %.12g %.12g %.12g %.12g\n", d, alpha, beta, xe, ye, x0, y0, p->xe, p->var[SOURCE_X], p->ye, p->var[SOURCE_Y] );
	rx = 2.*sqrt( tau * ax * vx );
	ry = 2.*sqrt( tau * ay * vx );
	rz = 2.*sqrt( tau * az * vx );
	e1 = exp( -tau * lambda );
	ex = erfc( ( xe - source_sizex / 2 - tau * vx ) / rx ) - erfc( ( xe + source_sizex / 2 - tau * vx ) / rx );
	ey = erfc( ( ye - source_sizey / 2 ) / ry ) - erfc( ( ye + source_sizey / 2 ) / ry );
	ez = erfc( ( ze - source_sizez ) / rz ) - erfc( ( ze + source_sizez ) / rz - erfc( ( ze - source_z ) / rz ) + erfc( ( ze + source_z ) / rz ) );
//	if( p->debug >= 3 ) printf( "int %g %g %g %g %g\n", tau, e1, ex, ey, ez );
	return( e1 * ex * ey * ez );
}

double rectangle_source( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( NUMITER );
	gsl_function F;
	int status;
	struct calc_data *p = ( struct calc_data * )params;
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
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	//	printf("result %g ", result, var[C0], p );
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[C0] * 1e6 * p->ze / ( 8. * sqrt( M_PI * p->var[AZ] * p->var[VX] ) ) * result );
}

double int_rectangle_source( double tau, void *params )
{
	struct calc_data *p = ( struct calc_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double rx, ry, rz, e1, ex, ey;
	double d, alpha, beta, xe, ye, ze, x0, y0;
	x0 = ( p->xe - p->var[SOURCE_X] );
	y0 = ( p->ye - p->var[SOURCE_Y] );
	d = ( -p->var[FLOW_ANGLE] * M_PI ) / 180;
	alpha = cos( d );
	beta = sin( d );
	xe = x0 * alpha - y0 * beta;
	ye = x0 * beta  + y0 * alpha;
	ze = ( p->ze - p->var[SOURCE_Z] );
	rx = 2.*sqrt( tau * ax * vx );
	ry = 2.*sqrt( tau * ay * vx );
	rz = 2.*sqrt( tau * az * vx );
	e1 = exp( -tau * lambda - ze * ze / ( tau * ( 4 * az * vx ) ) );
	ex = erfc( ( xe - source_sizex / 2 - tau * vx ) / rx ) - erfc( ( xe + source_sizex / 2 - tau * vx ) / rx );
	ey = erfc( ( ye - source_sizey / 2 ) / ry ) - erfc( ( ye + source_sizey / 2 ) / ry );
	return( e1 * ex * ey * pow( tau, -1.5 ) );
}

double rectangle_source_vz( double x, double y, double z, double t, void *params )
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( NUMITER );
	gsl_function F;
	int status;
	struct calc_data *p = ( struct calc_data * )params;
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
	if( t < p->var[TIME_END] )
		status = gsl_integration_qags( &F, 0, time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	else
		status = gsl_integration_qags( &F, time - ( p->var[TIME_END] - p->var[TIME_INIT] ), time, EPSABS, EPSREL, NUMITER, w, &result, &error );
	if( status != 0 ) result = 0;
	gsl_integration_workspace_free( w );
	// Concentrations are multiplied by 1e6 to convert in ppm!!!!!!!
	return( p->var[C0] * 1e6 * p->var[VZ] / ( 4. * p->var[POROSITY] ) * result );
}
double int_rectangle_source_vz( double tau, void *params )
{
	struct calc_data *p = ( struct calc_data * )params;
	double lambda = ( p->var[LAMBDA] );
	double vx = ( p->var[VX] );
	double vz = ( p->var[VZ] );
	double ax = ( p->var[AX] );
	double ay = ( p->var[AY] );
	double az = ( p->var[AZ] );
	double source_sizex = ( p->var[SOURCE_DX] );
	double source_sizey = ( p->var[SOURCE_DY] );
	double rx, ry, rz, e1, ex, ey, ez, tz, v;
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
	rx = 2.*sqrt( tau * ax * v );
	ry = 2.*sqrt( tau * ay * v );
	rz = 2.*sqrt( tau * az * v );
	e1 = exp( -tau * lambda );
	ex = erfc( ( xe - source_sizex / 2 - tau * vx ) / rx ) - erfc( ( xe + source_sizex / 2 - tau * vx ) / rx );
	ey = erfc( ( ye - source_sizey / 2 ) / ry ) - erfc( ( ye + source_sizey / 2 ) / ry );
	tz = ze - tau * vz;
	ez = 1. / sqrt( tau * ( M_PI * az * v ) ) * exp( -tz * tz / ( tau * ( 4 * az * v ) ) ) -
		 vz / ( 2 * az * v ) * exp( vz * ze / ( az * v ) ) * erfc( ( ze + tau * vz ) / rz );
	return( e1 * ex * ey * ez );
}
