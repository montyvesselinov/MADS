#include <math.h>
#include "../mads.h"

int set_test_problems( struct calc_data *cd, struct param_data *pd );
double test_problems( int D, int function, double *x );
char **char_matrix( int maxCols, int maxRows );

int set_test_problems( struct calc_data *cd, struct param_data *pd )
{
	int d;
	cd->sintrans = 0; // No sin transformations
	pd->nParam = pd->nOptParam = cd->dim;
	pd->nFlgParam = 0;
	pd->var_id = char_matrix(( *pd ).nParam, 50 );
	pd->var = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	cd->var = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_opt = ( int * ) malloc(( *pd ).nParam * sizeof( int ) );
	pd->var_log = ( int * ) malloc(( *pd ).nParam * sizeof( int ) );
	pd->var_dx = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_min = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_max = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_range = ( double * ) malloc(( *pd ).nParam * sizeof( double ) );
	pd->var_index = ( int * ) malloc(( *pd ).nOptParam * sizeof( int ) );
	pd->var_current = ( double * ) malloc(( *pd ).nOptParam * sizeof( double ) );
	pd->var_best = ( double * ) malloc(( *pd ).nOptParam * sizeof( double ) );
	for( d = 0; d < cd->dim; d++ )
	{
		sprintf( pd->var_id[d], "Parameter #%d", d + 1 );
		pd->var[d] = cd->var[d] = pd->var_current[d] = pd->var_best[d] = 0;
		pd->var_max[d] = 100; pd->var_min[d] = -100; pd->var_log[d] = 0; pd->var_opt[d] = 1;
		if( cd->problem_type == ABAGUS ) pd->var_dx[d] = .1;
		else pd->var_dx[d] = 0; // if not zero will force TRIBES to use discretized parameter space!
		pd->var_range[d] = pd->var_max[d] - pd->var_min[d];
		pd->var_index[d] = d;
	}
	switch( cd->test )
	{
		case 0: // Parabola (Sphere)
			printf( "Parabola (Sphere)" );
			break;
		case 1: // De Jong's Function 4
			printf( "De Jong's Function #4" );
			break;
		case 2: // Griewank
			printf( "Griewank" );
			break;
		case 3: // Rosenbrock
			printf( "Rosenbrock" );
			if( cd->problem_type == ABAGUS )
				for( d = 0; d < cd->dim; d++ )
					pd->var[d] = cd->var[d] = pd->var_current[d] = pd->var_best[d] = 0;
			break;
		case 4: // Step
			printf( "Step" );
			break;
		case 6: //Foxholes 2D
			printf( "Foxholes 2D" );
			if( cd->dim != 2 ) cd->dim = 2;
			break;
		case 7: // Polynomial fitting problem on [-100 100]^9
			printf( "Polynomial fitting" );
			break;
		case 8: // Clerc's f1, Alpine function, min 0
			printf( "Alpine function (Clerc's Function #1)" );
			break;
		case 9: // Rastrigin Minimum value 0. Solution (0,0 ...0)
			printf( "Rastrigin" );
			break;
		case 10: // Ackley
			printf( "Ackley" );
			break;
		case 13: // 2D Tripod function (Louis Gacogne) Search [-100, 100]^2. min 0 on (0  -50)
			printf( "2D Tripod function" );
			if( cd->dim != 2 ) cd->dim = 2;
			break;
		case 17: // KrishnaKumar
			printf( "Krishna Kumar" );
			break;
		case 18: // Eason 2D (usually on [-100,100] Minimum -1 on (pi,pi)
			printf( "Eason 2D " );
			if( cd->dim != 2 ) cd->dim = 2;
			break;
	}
	return( 0 );
}

double test_problems( int D, int function, double *x )
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
	double py, y = -1, dx = ( double )M;
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
				f = f + xd * xd;
			}
			break;
		case 1: // De Jong's f4
			f = 0;
			p = 0; // Shift
			for( d = 0; d < D; d++ )
			{
				xd = x[d] - p;
				f = f + ( d + 1 ) * xd * xd * xd * xd;
			}
			break;
		case 2: // Griewank
			f = 0;
			p = 1;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				f = f + xd * xd;
				p = p * cos( xd / sqrt( d + 1 ) );
			}
			f = f / 4000 - p + 1;
			break;
		case 3: // Rosenbrock
			f = 0;
			t0 = x[0];
			for( d = 1; d < D; d++ )
			{
				t1 = x[d];
				tt = 1 - t0;
				f += tt * tt;
				tt = t1 - t0 * t0;
				f += 100 * tt * tt;
				t0 = t1;
			}
			break;
		case 4: // Step
			f = 0;
			for( d = 0; d < D; d++ ) f = f + ( int )x[d];
			break;
		case 6: //Foxholes 2D
			f = 0;
			for( j = 0; j < 25; j++ )
			{
				sum1 = 0;
				for( d = 0; d < 2; d++ )
				{
					sum1 = sum1 + pow( x[d] - a[d] [j], 6 );
				}
				f = f + 1 / ( j + 1 + sum1 );
			}
			f = 1 / ( 0.002 + f );
			break;
		case 7: // Polynomial fitting problem on [-100 100]^9
			f = 0;
			dx = 2 / dx;
			for( i = 0; i <= M; i++ )
			{
				py = x[0];
				for( d = 1; d < D; d++ )
				{
					py = y * py + x[d];
				}
				if( py < -1 || py > 1 ) f += ( 1 - py ) * ( 1 - py );
				y += dx;
			}
			py = x[0];
			for( d = 1; d < D; d++ ) py = 1.2 * py + x[d];
			py = py - 72.661;
			if( py < 0 ) f += py * py;
			py = x[0];
			for( d = 1; d < D; d++ ) py = -1.2 * py + x[d];
			py = py - 72.661;
			if( py < 0 ) f += py * py;
			break;
		case 8: // Clerc's f1, Alpine function, min 0
			f = 0;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				f += fabs( xd * sin( xd ) + 0.1 * xd );
			}
			break;
		case 9: // Rastrigin. Minimum value 0. Solution (0,0 ...0)
			k = 10;
			f = 0;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				f += xd * xd - k * cos( 2 * pi * xd );
			}
			f += D * k;
			break;
		case 10: // Ackley
			sum1 = 0;
			sum2 = 0;
			for( d = 0; d < D; d++ )
			{
				xd = x[d];
				sum1 += xd * xd;
				sum2 += cos( 2 * pi * xd );
			}
			y = D;
			f = ( -20 * exp( -0.2 * sqrt( sum1 / y ) ) - exp( sum2 / y ) + 20 + E );
			break;
		case 13: // 2D Tripod function (Louis Gacogne) Search [-100, 100]^2. min 0 on (0  -50)
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
				f += sin( x[d] + x[d + 1] ) + sin( 2 * x[d] * x[d + 1] / 3 );
			break;
		case 18: // Eason 2D (usually on [-100,100] Minimum -1 on (pi,pi)
			x1 = x[0]; x2 = x[1];
			f = -cos( x1 ) * cos( x2 ) / exp(( x1 - pi ) * ( x1 - pi ) + ( x2 - pi ) * ( x2 - pi ) );
			break;
	}
	return f;
}
