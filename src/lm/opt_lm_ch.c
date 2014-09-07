// MADS: Model Analyses & Decision Support (v1.1) 2011
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MIN(X,Y) ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions elsewhere */
int lu_decomposition( double a[], double ul[], int n );
void lu_elimination( double a[], double b[], int n, double x[] );

int zxssqch( int func( double x[], void *, double f[] ), void *func_data,
			 int m, int n, int nsig, double eps, double delta, int MAXfn,
			 int iopt,  double parm[], double x[], double *ssq, double f[],
			 double xjac[], int ixjac, double xjtj[], int *infer )
{
	int     imjc, ieval, ibad, isw, iter, j, ijac, i, ier, k, l, is, js, icount, izero, label;
	double   al, cons2, dnorm, dsq, erl2, erl2x, f0, f0sq, g, hh, onesf0, ax,
			 prec, rel, rhh, sig, sqdif, ssqold, sum, xdif, xhold, up, xdabs,
			 relcon, delta2, temp;
//	double   f0sqs4;
	double  *work, *grad, *delx, *scall, *xnew, *xbad, *fplus, *fminus;
	f0 = f0sq = onesf0 = up = 0;
	ijac = icount = 0;
	sig = 6.3;
	ax = 0.1;
	ier = 0;
	if( m <= 0 || m > ixjac || n <= 0 || iopt < 0 || iopt > 2 ) { ier = 130; return( ier ); }
	if( iopt == 2 )
		if( parm[1] <= 1. || parm[0] <= 0 ) { ier = 130; return( ier ); }
	prec = pow( 10., -sig - 1. );
	rel = pow( 10., -sig * 0.5 );
	relcon = pow( 10., ( double ) - nsig );
	if( ( work = ( double * ) malloc( sizeof( double ) * ( n + 1 ) * n / 2 ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( grad = ( double * ) malloc( sizeof( double ) * n ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( delx = ( double * ) malloc( sizeof( double ) * n ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( scall = ( double * ) malloc( sizeof( double ) * n ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( xnew = ( double * ) malloc( sizeof( double ) * n ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( xbad = ( double * ) malloc( sizeof( double ) * n ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( fplus = ( double * ) malloc( sizeof( double ) * m ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	if( ( fminus = ( double * ) malloc( sizeof( double ) * m ) ) == NULL )
	{ printf( "Not enough memory!\n" ); return( 1 ); }
	imjc = ixjac - m;
	al = 1.;
	cons2 = 0.1;
	if( iopt != 0 )
	{
		if( iopt == 2 ) { al = parm[0]; f0 = parm[1]; up = parm[2]; cons2 = parm[3] * 0.5; }
		else { al = 0.01; f0 = 2.0; up = 1.2e2; }
		onesf0 = 1.0 / f0;
		f0sq = f0 * f0;
		// f0sqs4 = f0sq * f0sq * f0sq * f0sq;
	}
	ieval = 0;
	delta2 = delta * 0.5;
	erl2 = 1.0e10;
	ibad = -99;
	isw = 1;
	*infer = 0;
	ier = 0;
	iter = 0;
	for( j = 0; j < n; j++ ) delx[j] = 0;
	for( j = 0; j < n; j++ ) xnew[j] = x[j];
	func( xnew, func_data, fplus );
	ieval++;
	*ssq = 0;
	for( i = 0; i < m; i++ ) *ssq += fplus[i] * fplus[i];
	ssqold = *ssq;
	for( i = 0; i < m; i++ ) f[i] = fplus[i];
	printf( "Iteration: %5i phi: %g lambda: %g evaluations: %d\n", iter, *ssq, al, ieval );
	label = 55;
	/*
								main loop
	*/
	while( ier == 0 )
	{
		switch( label )
		{
			case 30:
				ssqold = *ssq;
				if( *infer > 0 || ijac >= n || iopt == 0 || icount > 0 ) { label = 55; continue; }
				else
				{
					//				printf( "ranking ..." );
					ijac++;
					dsq = 0;
					for( j = 0; j < n; j++ )
						dsq += delx[j] * delx[j];
					if( dsq > 0 )
					{
						for( i = 0; i < m; i++ )
						{
							g = f[i] - fminus[i];
							for( k = i, j = 0; j < n; j++, k += ixjac )
								g += xjac[k] * delx[j];
							g /= dsq;
							for( k = i, j = 0; j < n; j++, k += ixjac )
								xjac[k] -= g * delx[j];
						}
					}
					else { label = 55; printf( "delx magnitude is small\n" ); continue; }
				}
				break;
			case 55:
				ijac = 0;
				printf( "Iteration: %5i\n", ++iter );
				for( k = 0, j = 0; j < n; j++, k += imjc )
				{
					xdabs = fabs( x[j] );
					hh = rel * ( xdabs > ax ? xdabs : ax );
					xhold = x[j];
					x[j] += hh;
					func( x, func_data, fplus );
					ieval++;
					x[j] = xhold;
					if( isw == 1 )
					{
						/*  forward differences */
						rhh = 1.0 / hh;
						for( i = 0; i < m; i++, k++ )
							xjac[k] = ( fplus[i] - f[i] ) * rhh;
					}
					else
					{
						/*  central differences */
						x[j] = xhold - hh;
						func( x, func_data, fminus );
						ieval++;
						x[j] = xhold;
						rhh = 0.5 / hh;
						for( i = 0; i < m; i++, k++ )
							xjac[k] = ( fplus[i] - fminus[i] ) * rhh;
					}
				}
		}
		/*
		 *                              calculate gradient
		 */
		erl2x = erl2;
		erl2 = 0;
		for( k = 0, j = 0; j < n; j++, k += imjc )
		{
			sum = 0;
			for( i = 0; i < m; i++, k++ )
				sum += xjac[k] * f[i];
			grad[j] = sum;
			erl2 += sum * sum;
		}
		erl2 = sqrt( erl2 );
		/*
		 * 								convergence test for the norm of gradient
		 */
		if( ijac <= 0 )
		{
			if( erl2 <= delta2 ) *infer += 4;
			if( erl2 <= cons2 ) { isw = 2; printf( "central derivaives\n" ); }
		}
		/*
		 * 								calculate the lower triange of
		 * 								Jacobian(transposed) * Jacobian
		 */
		for( l = 0, is = 0, i = 0; i < n; i++, is += ixjac )
			for( j = 0, js = 0; j <= i; j++, js += ixjac )
			{
				sum = 0;
				for( k = 0; k < m; k++ )
					sum += xjac[is + k] * xjac[js + k];
				xjtj[l++] = sum;
			}
		/*
		 * 								convergence check
		 */
		if( *infer > 0 ) { ier = 1; break; }
		if( ieval >= MAXfn ) { ier = 133; break; }
		if( iopt == 0 )
		{
			dnorm = 0;
			for( k = 0, j = 0; j < n; j++, k += j + 1 )
			{
				scall[j] = sqrt( xjtj[k] );
				dnorm += xjtj[k] * xjtj[k];
			}
			dnorm = 1.0 / sqrt( dnorm );
			for( j = 0; j < n; j++ ) scall[j] *= dnorm * erl2;
		}
		else
		{
			/*  iopt !=0  */
			for( k = 0, j = 0; j < n; j++, k += j + 1 )
				scall[j] = xjtj[k];
		}
		/*
		 * 								add l-m factor to diagonal
		 */
		icount = 0;
		label = 140;
		while( 1 )
		{
			switch( label )
			{
				case 140:
					for( k = 0, i = 0; i < n; i++ )
					{
						for( j = 0; j <= i; j++, k++ )
							work[k] = xjtj[k];
						work[k - 1] += scall[i] * al;
						delx[i] = grad[i];
					}
					break;
				/*
				 * 								Cholesky decomposition
				 */
				case 155:
					ier = lu_decomposition( work, work, n );
					lu_elimination( work, delx, n, delx );
					if( ier != 0 )
					{
						ier = 0;
						if( ijac > 0 ) { label = 55; printf( "New iteration; singular hassian matrx\n" ); break; }
						if( ibad == 0 )
						{
							for( j = 0; j < n; j++ )
							{
								xhold = xbad[j];
								if( fabs( x[j] - xhold ) > relcon * MAX( ax, fabs( xhold ) ) )
									break;
							}
							if( j == n ) { ier = 132; break; }
						}
						if( ibad <= 0 )
						{
							for( j = 0; j < n; j++ ) xbad[j] = x[j];
							ibad = 1;
							if( iopt == 0 )
							{
								for( k = 0, i = 0; i < n; i++ )
								{
									for( j = 0; j <= i; j++, k++ )
										work[k] = xjtj[k];
									work[k - 1] = 1.5 * ( xjtj[k - 1] + al * erl2 * scall[i] ) + rel;
								}
								ibad = 2;
								label = 155;
								continue;
							}
							else
							{
								/* iopt != 0 */
								izero = 0;
								for( j = 0; j < n; j++ )
								{
									if( scall[j] > 0 ) continue;
									izero++;
									scall[j] = 1.0;
								}
								if( izero < n ) { label = 140; continue; }
								else { ier = 38; break; }
							}
						}
						else if( ibad >= 2 ) { ier = 129; break; }
						else { label = 190; continue; }
					}
					else
					{
						if( ibad != -99 ) ibad = 0;
						/*
						 * 								calculate sum of squares
						 */
						for( j = 0; j < n; j++ )
							xnew[j] = x[j] - delx[j];
						func( xnew, func_data, fplus );
						ieval++;
						*ssq = 0;
						for( i = 0; i < m; i++ )
							*ssq += fplus[i] * fplus[i];
						printf( "                 phi: %g lambda: %g evaluations: %d\n", *ssq, al, ieval );
						if( iopt != 0 )
						{
							if( *ssq > ssqold ) { label = 190; continue; }
							else if( icount == 0 ) al /= f0;
							if( erl2x > 0 )
							{
								g = erl2 / erl2x;
								if( erl2 < erl2x ) al *= ( onesf0 > g ? onesf0 : g );
								if( erl2 > erl2x ) al *= ( f0 < g ? f0 : g );
							}
							al = al > prec ? al : prec;
						}
						for( j = 0; j < n; j++ ) x[j] = xnew[j];
						for( i = 0; i < m; i++ )
						{
							fminus[i] = f[i];
							f[i] = fplus[i];
						}
						if( al > 5.0 ) { label = 30; break; }
						for( j = 0; j < n; j++ )
						{
							xdif = fabs( delx[j] ) / MAX( fabs( x[j] ), ax );
							if( xdif > relcon ) break;
						}
						if( j == n ) *infer = 1;
						sqdif = fabs( *ssq - ssqold ) / MAX( ssqold, ax );
						if( sqdif <= eps ) *infer += 2;
						label = 30;
						continue;
					}
					break;
				case 190:
					icount++;
					//				al*=f0; /* slow */
					al *= f0sq; /* fast */
					if( ijac == 0 || ( icount < 4 && al <= up ) )
					{
						if( al <= up ) { label = 140; continue; }
						else if( ibad == 1 ) { ier = 129; break; }
						else { ier = 39; break; }
					}
					//				else { al/=f0sq; printf("enough bad attempts (%d) or lambda is too large; (%d)\n",icount,ijac); ier = 39; break; }
					else { al /= f0sq; label = 55; printf( "enough bad attempts (%d) or lambda is too large; (%d)\n", icount, ijac ); continue; }
					break;
			}
			if( ier > 0 || label < 140 ) break;
		}
	}
	if( ier != 130 )
	{
		if( ier == 1 ) ier = 0;
		g = sig;
		for( j = 0; j < n; j++ )
		{
			xhold = fabs( delx[j] );
			if( xhold <= 0 ) continue;
			temp = ( ax > fabs( x[j] ) ) ? ax : fabs( x[j] );
			temp = -log10( xhold ) + log10( temp );
			g = g < temp ? g : temp;
		}
		if( n <= 2 ) { for( j = 0; j < n; j++ ) work[j + 5] = grad[j]; }
		work[0] = erl2 + erl2;
		work[1] = ieval;
		work[2] = g;
		work[3] = al;
		work[4] = iter;
	}
	free( work ); free( grad ); free( delx ); free( scall ); free( xnew ); free( xbad ); free( fplus ); free( fminus );
	return( ier );
}
