// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
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

double epsilon( void );

double epsilon( void )
{
	double r = 1.0;
	while( 1.0 < ( double )( 1.0 + r ) ) r /= 2.0;
	return ( 2.0 * r );
}

int lu_decomposition( double a[], double ul[], int n )
{
	double x, rn;
	int ip, i, j, k, iq, ip1, ir;
	rn = ( double ) 1.0 / ( 16.0 * n );
	ip = ip1 = 0;
	for( i = 0; i < n; i++ )
	{
		iq = ip;
		for( ir = 0, j = 0; j <= i; j++, ir++ )
		{
			x = a[ip];
			if( j )
				for( k = iq; k <= ip1; k++ )
					x -= ul[k] * ul[ir++];
			if( j == i )
			{
				if( ( a[ip] + x * rn ) <= a[ip] )
					return( 1 );
				ul[ip] = ( double ) 1.0 / sqrt( x );
			}
			else
				ul[ip] = x * ul[ir];
			ip1 = ip++;
		}
	}
	return( 0 );
}

/* The Original LUDECP Function
int ludecp( double a[], double ul[], int n, double *d1, double *d2 )
{
	double x, rn;
	int ip, i, j, k, iq, ip1, ir;

	*d1 = 1.0;
	*d2 = 0.0;
	rn = 1.0 / ( (double) n * 16.0 );
	ip = 0;
	for( i = 0; i < n; i++ )
	{
		iq = ip;
		for( ir = 0, j = 0; j <= i; j++, ir++ )
		{
			x = a[ip];
			if( j )
				for( k = iq; k <= ip1; k++ )
					 x -= ul[k] * ul[ir++];
			if( j == i )
			{
				*d1 *= x;
				if( ( a[ip] + x * rn ) <= a[ip] )
					return( 129 );
				while( fabs( *d1 ) > 1.0 )
				{
					*d1 *= 0.0625;
					*d2 += 4.0;
				}
				while( fabs( *d1 ) < 0.0625 )
				{
					*d1 *= 16.0;
					*d2 -= 4.0;
				}
				ul[ip] = 1.0 / sqrt(x);
			}
			else
				ul[ip] = x * ul[ir];
			ip1 = ip++;
		}
	}
	return( 0 );
}
*/

void lu_elimination( double a[], double b[], int n, double x[] )
{
	double t;
	int i, ip, is, iw, k, kk;
	ip = 0;
	iw = -1;
	for( i = 0; i < n; i++ )
	{
		t = b[i];
		if( iw == -1 )
		{
			if( fabs( t ) > DBL_EPSILON ) iw = i;
			ip += i;
		}
		else
		{
			ip += iw;
			for( k = iw; k < i; )
				t -= a[ip++] * x[k++];
		}
		x[i] = t * a[ip++];
	}
	for( i = n - 1; i > -1; i-- )
	{
		is = --ip;
		t = x[i];
		for( kk = n, k = i + 1; k < n; k++ )
		{
			t -= a[is] * x[--kk];
			is -= kk;
		}
		x[i] = t * a[is];
	}
}
