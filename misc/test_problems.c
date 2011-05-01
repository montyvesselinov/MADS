#include <math.h>
#include <string.h>
#include "../mads.h"

int set_test_problems( struct opt_data *op );
double test_problems( int D, int function, double *x, int nObs, double *f );
char **char_matrix( int maxCols, int maxRows );

int set_test_problems( struct opt_data *op )
{
	int d;
	struct calc_data *cd;
	struct param_data *pd;
	struct obs_data *od;
	cd = op->cd;
	pd = op->pd;
	od = op->od;
	cd->sintrans = 0; // No sin transformations
	pd->nParam = pd->nOptParam = cd->test_func_dim;
	pd->nFlgParam = 0;
	pd->var_id = char_matrix( ( *pd ).nParam, 50 );
	pd->var = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	cd->var = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc( ( *pd ).nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc( ( *pd ).nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc( ( *pd ).nParam * sizeof( double ) );
	pd->var_index = ( int * ) malloc( ( *pd ).nOptParam * sizeof( int ) );
	pd->var_current = ( double * ) malloc( ( *pd ).nOptParam * sizeof( double ) );
	pd->var_best = ( double * ) malloc( ( *pd ).nOptParam * sizeof( double ) );
	od->nObs = 0; // assume no observations ...
	for( d = 0; d < cd->test_func_dim; d++ )
	{
		sprintf( pd->var_id[d], "Parameter #%d", d + 1 );
		pd->var[d] = cd->var[d] = pd->var_current[d] = pd->var_best[d] = 0;
		pd->var_max[d] = 100; pd->var_min[d] = -100; pd->var_log[d] = 0; pd->var_opt[d] = 1;
		if( cd->problem_type == ABAGUS ) pd->var_dx[d] = .1;
		else if( strcasestr( op->cd->opt_method, "lm" ) != NULL ) pd->var_dx[d] = 1;
		else pd->var_dx[d] = 0; // if not zero will force TRIBES to use discretized parameter space!
		pd->var_range[d] = pd->var_max[d] - pd->var_min[d];
		pd->var_index[d] = d;
	}
	od->nObs = 0; // modified for test problems with observations below
	switch( cd->test_func )
	{
		case 0: // Parabola (Sphere)
			printf( "Parabola (Sphere)" );
			od->nObs = cd->test_func_dim;
			break;
		case 1: // De Jong's Function 4
			printf( "De Jong's Function #4" );
			od->nObs = cd->test_func_dim;
			break;
		case 2: // Griewank
			printf( "Griewank" );
			od->nObs = cd->test_func_dim;
			if( cd->problem_type == ABAGUS )
				for( d = 0; d < cd->test_func_dim; d++ )
					pd->var[d] = cd->var[d] = pd->var_current[d] = pd->var_best[d] = 0; // global minimum at (0,0, ... )
			break;
		case 3: // Rosenbrock
			printf( "Rosenbrock" );
			od->nObs = cd->test_func_dim;
			if( cd->problem_type == ABAGUS )
				for( d = 0; d < cd->test_func_dim; d++ )
					pd->var[d] = cd->var[d] = pd->var_current[d] = pd->var_best[d] = 1; // global minimum at (1,1, ... )
			break;
		case 4: // Step
			printf( "Step" );
			od->nObs = cd->test_func_dim;
			break;
		case 6: //Foxholes 2D
			printf( "Foxholes 2D" );
			od->nObs = 25;
			if( cd->test_func_dim != 2 ) cd->test_func_dim = 2;
			break;
		case 7: // Polynomial fitting problem on [-100 100]
			printf( "Polynomial fitting" );
			break;
		case 8: // Clerc's f1, Alpine function, min 0
			printf( "Alpine function (Clerc's Function #1)" );
			od->nObs = cd->test_func_dim;
			break;
		case 9: // Rastrigin Minimum value 0. Solution (0,0 ...0)
			printf( "Rastrigin" );
			od->nObs = cd->test_func_dim;
			break;
		case 10: // Ackley
			printf( "Ackley" );
			break;
		case 13: // 2D Tripod function (Louis Gacogne) Search [-100, 100] min 0 on (0  -50)
			printf( "2D Tripod function" );
			if( cd->test_func_dim != 2 ) cd->test_func_dim = 2;
			break;
		case 17: // KrishnaKumar
			printf( "Krishna Kumar" );
			od->nObs = ( cd->test_func_dim - 1 ) * 2;
			break;
		case 18: // Eason 2D (usually on [-100,100] Minimum -1 on (pi,pi)
			printf( "Eason 2D " );
			if( cd->test_func_dim != 2 ) cd->test_func_dim = 2;
			break;
		case 33: // Rosenbrock
			printf( "Rosenbrock (with observations = (d-1)*2)" );
			od->nObs = ( cd->test_func_dim - 1 ) * 2;
			if( cd->problem_type == ABAGUS )
				for( d = 0; d < cd->test_func_dim; d++ )
					pd->var[d] = cd->var[d] = pd->var_current[d] = pd->var_best[d] = 1; // global minimum at (1,1, ... )
			break;
	}
	if( od->nObs > 1 )
	{
		od->obs_id = char_matrix( ( *od ).nObs, 50 );
		od->obs_target = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_weight = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_min = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_max = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_current = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->res = ( double * ) malloc( ( *od ).nObs * sizeof( double ) );
		od->obs_log = ( int * ) malloc( ( *od ).nObs * sizeof( int ) );
	}
	for( d = 0; d < od->nObs; d++ )
	{
		sprintf( od->obs_id[d], "Observation #%d", d + 1 );
		od->obs_target[d] = od->obs_max[d] = od->obs_min[d] = 0;
		od->obs_weight[d] = 1;
		od->obs_log[d] = 0;
	}
	return( 0 );
}

double test_problems( int D, int function, double *x, int nObs, double *o )
{
	// Evaluate the fitness value for the particle of rank s
	int d;
	int i, j, k;
	double f, p, xd, x1, x2;
	double sum1, sum2;
	double t0, tt, t1, E, pi;
	static int a[2][25] =
	{
		{ -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32 },
		{ -32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32 }
	}; // For Foxholes problem
	int const M = 60; // For polynomial fitting problem
	double py, y, dx;
	E = exp( 1 );
	pi = acos( -1 );
	switch( function )
	{
		case 0: // Parabola (Sphere)
			f = 0;
			p = 0; // Shift
			for( d = 0; d < D; d++ )
			{
				xd = x[d] - p;
				f += o[d] = xd * xd;
			}
			break;
		case 1: // De Jong's f4
			f = 0;
			p = 0; // Shift
			for( d = 0; d < D; d++ )
			{
				xd = x[d] - p;
				f += o[d] = ( d + 1 ) * xd * xd * xd * xd;
			}
			break;
		case 2: // Griewank
			f = 0;
			p = 1;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				f += o[d] = xd * xd / 4000;
				p *= cos( xd / sqrt( ( double ) d + 1 ) );
			}
			for( d = 0; d < D; d++ )
				o[d] += ( -p + 1 ) / D;
			f += -p + 1;
			break;
		case 3: // Rosenbrock
			f = 0;
			t0 = x[0];
			for( d = 1; d < D; d++ )
			{
				t1 = x[d];
				tt = ( double ) 1 - t0;
				f += o[d] = tt * tt;
				if( d == 1 ) { o[0] = o[1]; o[1] = 0; } // first element
				tt = t1 - t0 * t0;
				f += x1 = tt * tt * 100;
				o[d] += x1;
				t0 = t1;
			}
			break;
		case 4: // Step
			f = 0;
			for( d = 0; d < D; d++ )
				f += o[d] = x[d];
			break;
		case 6: //Foxholes 2D
			f = 0;
			for( j = 0; j < 25; j++ )
			{
				sum1 = 0;
				for( d = 0; d < 2; d++ )
					sum1 += pow( x[d] - a[d] [j], 6 );
				f += o[d] = ( double ) 1.0 / ( j + 1 + sum1 );
			}
			f = ( double ) 1.0 / ( 0.002 + f );
			break;
		case 7: // Polynomial fitting problem on [-100 100]
			f = 0;
			y = -1;
			dx = ( double ) 2.0 / M;
			for( i = 0; i <= M; i++ ) // M = 60
			{
				py = x[0];
				for( d = 1; d < D; d++ )
					py = y * py + x[d];
				if( py < -1 || py > 1 ) f += ( 1 - py ) * ( 1 - py );
				y += dx;
			}
			py = x[0];
			for( d = 1; d < D; d++ )
				py = -1.2 * py + x[d];
			py -= 72.661;
			if( py < 0 ) f += py * py;
			break;
		case 8: // Clerc's f1, Alpine function, min 0
			f = 0;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				f += o[d] = fabs( xd * sin( xd ) + 0.1 * xd );
			}
			break;
		case 9: // Rastrigin. Minimum value 0. Solution (0,0 ...0)
			k = 10;
			f = 0;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				f += o[d] = k + xd * xd - k * cos( 2 * pi * xd );
			}
			break;
		case 10: // Ackley
			sum1 = sum2 = 0;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				sum1 += xd * xd;
				sum2 += cos( 2 * pi * xd );
			}
			y = D;
			f = ( -20 * exp( -0.2 * sqrt( sum1 / y ) ) - exp( sum2 / y ) + 20 + E );
			break;
		case 13: // 2D Tripod function (Louis Gacogne) Search [-100, 100] min 0 on (0  -50)
			x1 = x[0];
			x2 = x[1];
			if( x2 < 0 )
				f = fabs( x1 ) + fabs( x2 + 50 );
			else
			{
				if( x1 < 0 )
					f = 1 + fabs( x1 + 50 ) + fabs( x2 - 50 );
				else
					f = 2 + fabs( x1 - 50 ) + fabs( x2 - 50 );
			}
			break;
		case 17: // KrishnaKumar
			f = 0;
			for( d = 0; d < D - 1; d++ )
			{
				o[i++] = x1 = sin( x[d] + x[d + 1] );
				o[i++] = x2 = sin( 2 * x[d] * x[d + 1] / 3 );
				f += x1 + x2;
			}
			break;
		case 18: // Eason 2D (usually on [-100,100] Minimum -1 on (pi,pi)
			x1 = x[0]; x2 = x[1];
			f = -cos( x1 ) * cos( x2 ) / exp( ( x1 - pi ) * ( x1 - pi ) + ( x2 - pi ) * ( x2 - pi ) );
			break;
		case 33: // Rosenbrock (with more observations)
			f = 0;
			t0 = x[0];
			for( i = 0, d = 1; d < D; d++ )
			{
				t1 = x[d];
				tt = ( double ) 1.0 - t0;
				f += o[i++] = tt * tt;
				tt = t1 - t0 * t0;
				f += o[i++] = tt * tt * 100;
				t0 = t1;
			}
			break;
	}
	for( d = 0; d < nObs; d++ )
		o[d] = sqrt( o[d] );
	return f;
}
