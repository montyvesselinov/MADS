// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
// http://www.ees.lanl.gov/staff/monty/codes/mads
// http://gitlab.com/monty/mads
//
// Licensing: GPLv3: http://www.gnu.org/licenses/gpl-3.0.html
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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

void lhs_imp_dist( int nvar, int npoint, int d, unsigned *seed, double x[] );
void lhs_center( int nvar, int npoint, unsigned *seed, double x[] );
void lhs_edge( int nvar, int npoint, unsigned *seed, double x[] );
void lhs_random( int nvar, int npoint, unsigned *seed, double x[] );
void smp_random( int nvar, int npoint, unsigned *seed, double x[] );
int get_seed();
int get_seed_old();
int nint( float x );
int int_max( int i1, int i2 );
int int_min( int i1, int i2 );
double float_uniform( unsigned *seed );
int *perm_uniform( int n, int base, unsigned *seed );
int int_uniform( int ilo, int ihi, unsigned *seed );
// Elsewhere ...
void tprintf( char const *fmt, ... );

void lhs_imp_dist( int nvar, int npoint, int d, unsigned *seed, double x_int[] )
{
	int *avail, best, count, i, j, k, *list, *point, point_index;
	double dist, min_all, min_can, opt;
	if( ( avail = ( int * ) malloc( nvar * npoint * sizeof( int ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return; }
	if( ( list = ( int * ) malloc( d * npoint * sizeof( int ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return; }
	if( ( point = ( int * ) malloc( nvar * d * npoint * sizeof( int ) ) ) == NULL )
	{ tprintf( "Not enough memory!\n" ); return; }
	opt = ( ( double ) npoint ) / pow( ( double ) npoint, ( double )( 1.0 / ( double ) nvar ) );
	for( i = 0; i < nvar; i++ )
		x_int[i + ( npoint - 1 )*nvar] = ( double ) int_uniform( 1, npoint, seed );
	for( j = 0; j < npoint; j++ )
		for( i = 0; i < nvar; i++ )
			avail[i + j * nvar] = j + 1;
	for( i = 0; i < nvar; i++ )
		avail[i + ( ( ( int ) x_int[i + ( npoint - 1 )*nvar] ) - 1 )*nvar] = npoint;
	for( count = npoint - 1; 2 <= count; count-- )
	{
		for( i = 0; i < nvar; i++ )
		{
			for( k = 0; k < d; k++ )
				for( j = 0; j < count; j++ )
					list[j + k * count] = avail[i + j * nvar];
			for( k = count * d - 1; 0 <= k; k-- )
			{
				point_index = int_uniform( 0, k, seed );
				point[i + k * nvar] = list[point_index];
				list[point_index] = list[k];
			}
		}
		min_all = 1e30;
		best = 0;
		for( k = 0; k < d * count; k++ )
		{
			min_can = 1e30;
			for( j = count; j < npoint; j++ )
			{
				dist = 0.0;
				for( i = 0; i < nvar; i++ )
					dist += ( point[i + k * nvar] - x_int[i + j * nvar] ) * ( point[i + k * nvar] - x_int[i + j * nvar] );
				dist = sqrt( dist );
				if( dist < min_can ) min_can = dist;
			}
			if( fabs( min_can - opt ) < min_all )
			{
				min_all = fabs( min_can - opt );
				best = k;
			}
		}
		for( i = 0; i < nvar; i++ )
			x_int[i + ( count - 1 )*nvar] = point[i + best * nvar];
		for( i = 0; i < nvar; i++ )
			for( j = 0; j < npoint; j++ )
				if( avail[i + j * nvar] == x_int[i + ( count - 1 )*nvar] ) avail[i + j * nvar] = avail[i + ( count - 1 ) * nvar];
	}
	for( i = 0; i < nvar; i++ )
		x_int[i + 0 * nvar] = avail[i + 0 * nvar];
	for( j = 0; j < npoint; j++ )
		for( i = 0; i < nvar; i++ )
			x_int[i + j * nvar] /= ( double ) npoint;
	free( avail ); free( list ); free( point );
}

void lhs_center( int nvar, int npoint, unsigned *seed, double x[] )
{
	int base = 0, i, j, k, *perm;
	for( k = i = 0; i < nvar; i++ )
	{
		perm = perm_uniform( npoint, base, seed );
		for( j = 0; j < npoint; j++ )
			x[k++] = ( ( ( double ) perm[j] ) + 0.5 ) / ( ( double ) npoint );
		free( perm );
	}
}

void smp_random( int nvar, int npoint, unsigned *seed, double x[] )
{
	int i, j, k;
	for( k = i = 0; i < nvar; i++ )
		for( j = 0; j < npoint; j++ )
			x[k++] = float_uniform( seed );
}

void lhs_random( int nvar, int npoint, unsigned *seed, double x[] )
{
	int base = 0, i, j, k, *perm;
	double r;
	for( k = i = 0; i < nvar; i++ )
	{
		perm = perm_uniform( npoint, base, seed );
		for( j = 0; j < npoint; j++ )
		{
			r = float_uniform( seed );
			x[k++] = ( ( ( double ) perm[j] ) + r ) / ( ( double ) npoint );
		}
		free( perm );
	}
}

void lhs_edge( int nvar, int npoint, unsigned *seed, double x[] )
{
	int base = 0, i, j, k, *perm;
	if( npoint == 1 )
		for( k = i = 0; i < nvar; i++ )
			x[k++] = 0.5;
	else
		for( k = i = 0; i < nvar; i++ )
		{
			perm = perm_uniform( npoint, base, seed );
			for( j = 0; j < npoint; j++ )
				x[k++] = ( ( double ) perm[j] ) / ( ( float )( npoint - 1 ) );
			free( perm );
		}
}

int *perm_uniform( int n, int base, unsigned *seed )
{
	int i, j, k, *p;
	if( ( p = ( int * ) malloc( n * sizeof( int ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( NULL ); }
	for( i = 0; i < n; i++ )
		p[i] = i + base;
	for( i = 0; i < n; i++ )
	{
		j = int_uniform( i, n - 1, seed );
		k    = p[i];
		p[i] = p[j];
		p[j] = k;
	}
	return p;
}

int nint( float x )
{
	int value;
	if( x < 0.0 ) value = - ( int )( fabs( x ) + 0.5 );
	else           value = ( int )( fabs( x ) + 0.5 );
	return value;
}

double float_uniform( unsigned *seed )
{
	int k;
	k = *seed / 127773;
	*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
	if( *seed < 0 ) *seed += 2147483647;
	return( ( double )( *seed ) * 4.656612875E-10 );
}

int int_uniform( int a, int b, unsigned *seed )
{
	int k;
	float r;
	int value;
	k = *seed / 127773;
	*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
	if( *seed < 0 ) *seed += 2147483647;
	r = ( float )( *seed ) * 4.656612875E-10;
	//  Scale R to lie between A-0.5 and B+0.5.
	r = ( 1.0 - r ) * ( ( float )( int_min( a, b ) ) - 0.5 ) +  r  * ( ( float )( int_max( a, b ) ) + 0.5 );
	//  Use rounding to convert R to an integer between A and B.
	value = nint( r );
	value = int_max( value, int_min( a, b ) );
	value = int_min( value, int_max( a, b ) );
	return( value );
}

int get_seed()
{
	struct timeval time;
	gettimeofday(&time, NULL);
	unsigned pid = getpid();
	unsigned int seed = ( time.tv_sec + pid ) * 1000 + time.tv_usec / 1000;
	return seed;
}

int get_seed_old()
{
	time_t clock;
	struct tm *lt;
	time_t tloc;
	int ihour, imin, isec, seed;
	clock = time( &tloc );
	lt = localtime( &clock );
	ihour = lt->tm_hour;
	if( ihour > 12 ) ihour -= 12;
	ihour--;
	imin = lt->tm_min;
	isec = lt->tm_sec;
	// ihour = 0; imin = 0; isec = 0;
	seed = isec + 60 * ( imin + 60 * ihour );
	// tprintf( "Seed: %d %d %d %d\n", ihour, imin, isec, seed );
	seed = ( int )( ( double )( seed / ( 60.0 * 60.0 * 12.0 ) ) * 2147483647 ); 	//  Remap *seed from [1,43200] to [1,2147483647].
	if( seed <= 0 ) seed = 1;
	tprintf( "Seed: %d\n", seed );
	return seed;
}

int int_max( int i1, int i2 )
{
	if( i2 < i1 ) return i1; else return i2;
}

int int_min( int i1, int i2 )
{
	if( i1 < i2 ) return i1; else return i2;
}
