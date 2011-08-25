// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
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
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include "../mads.h"

#define MIN(X,Y) ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

#define TRUE  1
#define YES   1
#define NO    0
#define DAX 0.1
#define REL 1e-7
#define SIG 6.3

/* Functions elsewhere */
int lu_decomposition( double a[], double ul[], int n );
void lu_elimination( double a[], double b[], int n, double x[] );
double epsilon( void );

/* Purpose: Minimum of the sum of squares of m functions in n variables
 *          using a Finite Difference Levenberg-Marquardt algorithm
 * Arguments:
 * func  - A user supplied subroutine which calculates the residual vector
 *         f[1], f[2], ..., f[m] for given parameter values x[1], x[2], ...,
 *         x[n]. The calling sequence has the following form func(x, m, n, f)
 *         where x is a vector of length n and f is a vector of length m.
 *         Func must appear in an external statement in the calling program.
 *         Func must not alter the values of x[i], i=1, ..., n,  m,  or n.
 * m     - The number of residuals or observations (input)
 * n     - The number of unknown parameters (input).
 * nsig  - First convergence criterion. (input)
 *         Convergence condition satisfied if on two successive iterations, the
 *         parameter estimates agree, component by component, to nsig digits.
 * eps   - Second convergence criterion. (input)
 *         Convergence condition satisfied if, on two successive iterations the residual sum
 *         of squares estimates have relative difference less than or equal
 *         to eps. eps may be set to zero.
 * delta - Third convergence criterion. (input)
 *         Convergence condition satisfied if the (euclidean) norm of the
 *         approximate gradient is less than or equal to delta. Delta may be
 *         set to zero.
 *         Note, the iteration is terminated, and convergence is considered
 *         achieved, if any one of the three conditions is satisfied.
 * max_eval - Input maximum number of function evaluations
 *         (i.e., calls to subroutine func) allowed.
 *         The actual number of calls to func may exceed max_eval slightly.
 * iopt  - Input options parameter.
 *           iopt=0 implies Brown's algorithm without strict descent is desired.
 *           iopt=1 implies strict descent and default values for input vector parm are desired.
 *           iopt=2 implies strict descent is desired with user parameter choices in input vector parm.
 * parm  - Input vector of length 4 used only for iopt equal two.
 *         Parm[i] contains, when:
 *           i=0, the initial value of the Marquardt parameter used to scale the
 *                diagonal of the approximate jacTjac matrix, jacTjac, by the factor
 *                (1.0 + parm[0]). A small value gives a newton step, while a large
 *                value gives a steepest descent step.
 *                The default value for parm[0] is 0.01.
 *           i=1, The scaling factor used to modify the Marquardt parameter,
 *                which is decreased by parm[1] after an immediately successful
 *                descent direction, and increased by the square of parm[1] if not.
 *                parm[1] must be greater than one, and two is default.
 *           i=2, An upper bound for increasing the Marquardt parameter.
 *                The search for a descent point is abandoned if parm[2] is
 *                exceeded. parm[2] greater than 100.0 is recommended.
 *                Default is 120.0.
 *           i=3, Value for indicating when central rather than forward
 *                differencing is to be used for calculating the jacobian.
 *                The switch is made when the norm of the gradient of the sum
 *                of squares function becomes smaller than parm[3]. Central
 *                differencing is good in the vicinity of the solution, so
 *                parm[3] should be small.
 *                The default value is 0.10.
 * x     - Vector of length n containing parameter values.
 *         On input, x should contain the initial estimate of the location
 *         of the minimum.
 *         On output, x contains the final estimate of the location
 *         of the minimum.
 * phi   - Output scalar which is set to the residual sums of squares,
 *         f[1]^2 + ... + f[m]^2, for the final parameter estimates.
 * f     - Output vector of length m containing the *residuals* for the
 *         final parameter estimates.
 * jacobian  - Output m by n matrix containing the approximate jacobian at
 *         the output vector x.
 * ranking_countobian - Input row dimension of matrix jacobian exactly as specified in
 *         the dimension statement in the calling program.
 * jacTjac  - Output vector of length (n + 1) * n / 2 containing the n by n
 *         jacTjac matrix (transposed jacobian) * (jacobian) in symmetric storage mode.
 * Temp_Space  - Work vector of length 5 * n + 2 * m + (n + 1) * n / 2.
 *         On output, work[i] contains for
 *         i = 0, the norm of the gradient described under input parameters
 *                delta and parm[3].
 *         i = 1, the number of function evaluations required during the work[5]
 *                iterations.
 *         i = 2, the estimated number of significant digits in output vector x.
 *         i = 3, the final value of the Marquardt scaling parameter described
 *                under parm[0].
 *         i = 4, the number of iterations (i.e., changes to the x vector)
 *                performed.
 *         See programming notes for description of the latter elements of work.
 * infer - An integer that is set, on output, to indicate which convergence
 *         criterion was satisfied.
 *           infer = 0 indicates that convergence failed.
 *                     (ier gives further explanation.)
 *           infer = 1 indicates that the first criterion was satisfied.
 *           infer = 2 indicates that the second criterion was satisfied.
 *           infer = 4 indicates that the third criterion was satisfied.
 *         If more than one of the convergence criteria were satisfied on
 *         the final iteration, infer contains the corresponding sum.
 *         (e.g., infer = 3 implies first and second criteria satisfied simultaneously).
 * Return Parameter:
 */

int lm_opt( int func( double x[], void *data, double f[] ), int func_dx( double *x, double *f_x, void *data, double *jacobian ), void *data,
			int nObs, int nParam, int nsig, double eps, double delta, int max_eval, int max_iter,
			int iopt, double parm[], double x[], double *phi, double f[],
			double jacobian[], int nJacobian, double jacTjac[], int *infer )
{
	struct opt_data *func_data = ( struct opt_data * )data;
	int ieval, singular_count, central_derivatives, iter, iter_best, max_iter_best, max_phi_range, ranking_count, lambda_count, max_lambda_count, i, j, ier, k, l, is, js, iphi, bad_count, max_bad_count, lambda_up_count, izero, add, rank, debug, ofe_close;
	double lambda, limit_central_difference, dnorm, temp, temp2, grad_norm, grad_norm_old, factor, factor_squared, factor_pow8, dx, factor_recipical, phi_range, precision, rel, phi_diff, phi_old, phi_init, phi_best, phi_ibest, sum, x_diff, x_old, upper_bound, x_abs, relcon, rel_dx, *hessian, *grad, *x_update, *scale, *x_new, *x_best, *x_cur_best, *x_bad, *f_xpdx, *f_xmdx, *vphi, phi_range_min, phi_range_max, phi_cutoff;
	char filename[80];
	/* ERROR CHECKS FIRST EXECUTABLE STATEMENT */
	debug = func_data->cd->ldebug;
	phi_cutoff = func_data->cd->phi_cutoff;
	if( nObs <= 0 || nObs > nJacobian || nParam <= 0 || iopt < 0 || iopt > 2 || max_eval < 1 )
	{
		printf( "At least one of the input parameters (M, N, IOPT, max_eval) is specified incorrectly.\n" );
		return( 0 );
	}
	if( iopt == 2 )
		if( parm[1] <= 1.0 || parm[0] <= 0.0 )
		{
			printf( "At least one of the input parameters (PARM[1] or PARM[2]) is specified incorrectly.\n" );
			return( 0 );
		}
	precision = epsilon();
	/*
	printf( "p1 %g\n", precision );
	precision   = pow( 10.0, -sig - 1.0 ); // MACHINE DEPENDENT CONSTANTS
	printf( "p2 %g\n", precision );
	*/
	rel = sqrt( precision );
	phi_range = rel;
	rel_dx = sqrt( rel ) * 100;
	max_lambda_count = 5; // nax number of bad lambda's
	max_iter_best = max_iter / 2; // max number of interations since the current best one
	max_phi_range = max_iter / 4; // number of iterations to evaluate phi decline
	max_bad_count = MAX( 10, ( int ) nParam / 3 ); // maximum number of iterations to avoid Jacobian singuliarity
	if(( vphi = ( double * ) malloc( sizeof( double ) * max_eval ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( hessian = ( double * ) malloc( sizeof( double ) * ( nParam + 1 ) * nParam / 2 ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( grad = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( x_update = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( scale = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( x_new = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( x_cur_best = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( x_best = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( x_bad = ( double * ) malloc( sizeof( double ) * nParam ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( f_xpdx = ( double * ) malloc( sizeof( double ) * nObs ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
	if(( f_xmdx = ( double * ) malloc( sizeof( double ) * nObs ) ) == NULL )
		{ printf( "Not enough memory!\n" ); return( 1 ); }
//	imjc = nJacobian - nObs;
	if( iopt )                                      /* INITIALIZE VARIABLES */
	{
		switch( iopt )
		{
			case 1:
				lambda = 0.01;    factor = 2.0;
				upper_bound = 120.0;   limit_central_difference = 0.1; break;
			case 2:
				lambda = parm[0]; factor = parm[1];
				upper_bound = parm[2]; limit_central_difference = parm[3]; break;
		}
		factor_recipical = 1.0 / factor;
		factor_squared = factor * factor;
		factor_pow8 = factor_squared * factor_squared * factor_squared * factor_squared;
	}
	else { lambda = 1.0; limit_central_difference = 0.1; }
	grad_norm = 1.0e10;
	for( j = 0; j < nParam; j++ ) x_update[j] = 0.0;
	for( j = 0; j < nParam; j++ ) x_best[j] = x_new[j] = x[j];
	func( x_new, func_data, f_xpdx ); ieval = 1; iphi = 0;
	*phi = 0.0;
	for( i = 0; i < nObs; i++ ) { *phi += f_xpdx[i] * f_xpdx[i]; }
	if( debug ) printf( "Initial OF %g\n", *phi );
	vphi[iphi++] = phi_ibest = phi_best = phi_init = phi_old = *phi;
	for( i = 0; i < nObs; i++ ) f[i] = f_xpdx[i];
	lambda_count = bad_count = ranking_count = 0;
	singular_count = -99; central_derivatives = NO /*forwad derivatives */; rank = NO;
	iter = iter_best = 0;
	lambda_up_count = 0; //
	*infer = ier = 0;
	if( func_data->cd->odebug )
	{
		if( func_data->f_ofe == NULL )
		{
			if( func_data->counter > 0 ) sprintf( filename, "%s.%08d.ofe", func_data->root, func_data->counter );
			else sprintf( filename, "%s.ofe", func_data->root );
			if( func_data->cd->nretries > 1 && func_data->cd->retry_ind > 1 ) func_data->f_ofe = fopen( filename, "a" );
			else func_data->f_ofe = fopen( filename, "w" );
			ofe_close = 1;
		}
		else ofe_close = 0;
	}
	if( func_data->cd->check_success ) func_data->success = 0;
	while( !ier )                                 /* BEGIN OF THE MAIN LOOP */
	{
		if( rank ) /* BROYDEN'S RANK-ONE UPDATE OF JACOBIAN */
		{
			phi_old = *phi;
			if( debug > 2 ) printf( "Rank-one update of Jacobian matrix ...\n" );
			if( *infer > 0 || ranking_count >= nParam || iopt == 0 || bad_count > 0 )
			{
				if( debug > 2 ) printf( "Re-calculate Jacobian because ... %d %d %d\n", *infer, ranking_count, bad_count );
				rank = NO; continue; /* RE-CALCULATE JACOBIAN */
			}
			else if( debug > 2 ) printf( "Rank-one update of Jacobian matrix because ... %d %d %d\n", *infer, ranking_count, bad_count );
			ranking_count++;
			temp = 0.0;
			for( j = 0; j < nParam; j++ ) temp += x_update[j] * x_update[j];
			if( temp > DBL_EPSILON )
			{
				for( i = 0; i < nObs; i++ )
				{
					temp2 = f[i] - f_xmdx[i];
					for( k = i, j = 0; j < nParam; j++, k += nJacobian ) temp2 += jacobian[k] * x_update[j];
					temp2 /= temp;
					for( k = i, j = 0; j < nParam; j++, k += nJacobian ) jacobian[k] -= temp2 * x_update[j];
				}
			}
			else
			{
				if( debug ) printf( "Re-calculate Jacobian because because update vector norm is = 0 ...\n" );
				rank = NO; continue; /* RE-CALCULATE JACOBIAN */
			}
		}
		else /* CALCULATE JACOBIAN */
		{
			if( debug ) printf( "OF: last %g current best %g all-over best %g initial %g (evaluations %d)\n", *phi, phi_old, phi_best, phi_init, ieval );
			if( func_data->cd->odebug )
			{
				if( func_data->cd->standalone ) fprintf( func_data->f_ofe, "%d %g\n", func_data->cd->neval, phi_best ); // Print current best
				else fprintf( func_data->f_ofe, "%d %g\n", func_data->cd->neval, ( func_data->phi < phi_best ) ? func_data->phi : phi_best ); // Print overall best
				fflush( func_data->f_ofe );
			}
			if( fabs( phi_ibest - *phi ) / phi_ibest < phi_range ) iter_best++;
			else { phi_ibest = *phi; iter_best = 0; }
			if( ++iter > max_iter )
			{
				if( debug ) printf( "The number of maximum allowed iterations (%d) is exceeded.\n", max_iter );
				ier = 133;
				break; /* E X I T */
			}
			if( debug )
			{
				printf( "\nIteration %d", iter );
				if( iter > 1 && iter_best > 0 )  printf( " (%d following the current best result; max %d)\n", iter_best, max_iter_best );
				else printf( "\n" );
			}
			else if( debug > 2 )
			{
				if( central_derivatives == NO ) printf( "Calculate Jacobian using forward differences ...\n" );
				else printf( "Calculate Jacobian using central differences ...\n" ); // Do NOT use central difference with SIN transformation
			}
			lambda_count = ranking_count = 0;
			func_dx( x, f, func_data, jacobian );
			/*
			for( k = 0, j = 0; j < nParam; j++ ) // k += imjc
			{
				x_old = x[j];
				dx = REL * MAX( fabs( x_old) , DAX );
				x[j] += dx;
			//				printf( "Parameter %d: %g %g %g %g\n", j, x[j], ax, dx, rel );
				if( debug >= 4 ) printf( "Parameter %d: Upgradient\n", j + 1 );
				func( x, func_data, f_xpdx ); ieval++;
				x[j] = x_old;
				if( central_derivatives == NO )    // FORWARD DIFFERENCES
				{
					rdx = (double) 1.0 / dx;
					for( i = 0; i < nObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f[i] ) * rdx;
				}
				else                               // CENTRAL DIFFERENCES
				{
					x[j] -= dx;
					if( debug >= 4 ) printf( "Parameter %d: Downgradient\n", j + 1 );
					func( x, func_data, f_xmdx ); ieval++;
					x[j] = x_old;
					rdx = (double) 0.5 / dx;
					for( i = 0; i < nObs; i++, k++ ) jacobian[k] = ( f_xpdx[i] - f_xmdx[i] ) * rdx;
				}
			}
			*/
		}
		grad_norm_old = grad_norm;                 /* CALCULATE GRADIENT NORM */
		grad_norm = 0.0;
		for( k = 0, j = 0; j < nParam; j++ ) // k += imjc
		{
			sum = 0.0;
			for( i = 0; i < nObs; i++, k++ ) sum += jacobian[k] * f[i];
			grad[j] = sum;
			grad_norm += sum * sum;
		}
		grad_norm = sqrt( grad_norm ); /// Gradient norm equal to zero indicates that the current model predictions are equal to zero
		if( debug > 2 ) printf( "Gradient norm change %g -> %g\n", grad_norm_old, grad_norm );
		if( ranking_count == 0 )
		{
			if( grad_norm <= delta + DBL_EPSILON ) *infer += 4; /* CONVERGENCE TEST FOR NORM OF GRADIENT */
			if( grad_norm <= limit_central_difference ) { if( debug ) printf( "Central differences: gradient norm %g <= %g\n", grad_norm, limit_central_difference ); central_derivatives = YES; }
//			else { printf( "Forward differences: gradient norm %g > %g\n", grad_norm, limit_central_difference ); central_derivatives = NO; }
		}
		for( l = 0, is = 0, i = 0; i < nParam; i++, is += nJacobian ) /* CALCULATE THE LOWER/UPPER TRIANGLE OF */
			for( j = 0, js = 0; j <= i; j++, js += nJacobian )        /* JACOBIAN (TRANSPOSED) * JACOBIAN */
			{
				sum = 0.0;
				for( k = 0; k < nObs; k++ ) sum += jacobian[is + k] * jacobian[js + k];
				jacTjac[l++] = sum;
			}
		if( *infer > 0 ) /* CONVERGENCE CHECKS */
		{
			if( debug ) printf( "\nLevenberg-Marquardt optimization finished successfully!\n" );
			ier = 1;
			break; /* E X I T */
		}
		if( iter_best > max_iter_best )
		{
			if( debug ) printf( "The number of maximum allowed iterations (%d > %d) following the best achived result is exceeded.\n", iter_best, max_iter_best );
			ier = 133; *infer += 64;
			break; /* E X I T */
		}
		if( ieval >= max_eval )
		{
			if( debug ) printf( "The number of maximum allowed calls of the model function (%d) is exceeded.\n", max_eval );
			ier = 133; *infer += 32;
			break; /* E X I T */
		}
		if( iopt )                                /* COMPUTE SCALING VECTOR */
			for( k = 0, j = 0; j < nParam; k += ++j + 1 ) scale[j] = jacTjac[k];
		else
		{
			dnorm = 0.0;                 /* COMPUTE SCALING VECTOR AND NORM */
			for( k = 0, j = 0; j < nParam; k += ++j + 1 )
			{
				scale[j] = sqrt( jacTjac[k] );
				dnorm += jacTjac[k] * jacTjac[k];
			}
			dnorm = ( double ) 1.0 / sqrt( dnorm );
			for( j = 0; j < nParam; j++ )
				scale[j] *= dnorm * grad_norm;
		}
		bad_count = 0;
		add = YES;
		while( TRUE )                           /* BEGIN OF THE LOCAL WHILE */
		{
			if( add == YES )                         /* ADD L-M FACTOR TO DIAGONAL */
				for( k = 0, i = 0; i < nParam; i++, k++ )
				{
					l = k + i; // Upper limit
					for( ; k < l; k++ )
						hessian[k] = jacTjac[k];
					hessian[k] = jacTjac[k] + scale[i] * lambda;
					x_update[i] = grad[i];
				}
			if( lu_decomposition( hessian, hessian, nParam ) == 1 )         /* CHOLESKY DECOMPOSITION */
			{
				/* The matrix is algorithmically not positive definite */
				if( debug > 2 ) printf( "Hessian matrix is not positive definite\n" );
				if( ranking_count > 0 ) { if( debug ) printf( "New iteration; jacobian will be recomputed\n" ); rank = NO; break; /* REDO AND CALCULATE JACOBIAN */ }
				if( singular_count > 0 ) /* HAS CURRENT ITERATE CYCLED BACK TO THE LAST SINGULAR POINT */
				{
					for( j = 0; j < nParam; j++ )
					{
						x_old = x_bad[j];
						x_abs = fabs( x_old );
						if( fabs( x[j] - x_old ) > ( relcon * MAX( x_abs, DAX ) ) )
							break;
					}
					if( j == nParam )
					{
						if( debug ) printf( "After a successful recovery from a singular Jacobian\nthe solution has cycled back to the previous singularity.\n" );
						ier = 132; break; /* E X I T */
					}
				}
				if( singular_count <= 0 )
				{
					for( j = 0; j < nParam; j++ )     /* UPDATE THE BAD X VALUES */
						x_bad[j] = x[j];
					if( iopt )
					{
						/* REPLACE ZEROES ON jacTjac DIAGONAL */
						singular_count = 1;
						izero = 0;
						for( j = 0; j < nParam; j++ )
							if( scale[j] <= DBL_EPSILON )
							{
								izero++;
								scale[j] = 1.0;
							}
						if( izero < nParam )
						{
							if( debug ) printf( "Singularity Remedy: 0's replaced with 1's on diagonal of Hessian matrix\n" );
							add = YES; continue; /* ReDo the ADD loop */
						}
						else
						{
							if( debug ) printf( "Jacobian is ZERO. The solution X is a stationary point.\n" );
							ier = 38; break; /* E X I T */
						}
					}
					else
					{
						/* INCREASE DIAGONAL OF jacTjac */
						if( debug ) printf( "Singularity Remedy: Scale-up diagonal elements of Hessian matrix\n" );
						singular_count = 2;
						for( k = 0, i = 0; i < nParam; i++, k++ )
						{
							l = k + i; // Upper limit
							for( ; k < l; k++ )
								hessian[k] = jacTjac[k];
							hessian[k] = 1.5 * ( jacTjac[k] + lambda * grad_norm * scale[i] ) + rel;
						}
						add = NO; continue; /* ReDo the ADD loop */
					}
				}
				else if( singular_count == 2 ) /* iopt == 0 */
				{
					if( debug ) printf( "Singularity is detected in the Jacobian and recovery failed.\n" );
					ier = 129; break; /* E X I T */
				}
				else  /* singular_count == 1 && iopt != 1 */
				{
					bad_count++;
					if( ranking_count > 0 || bad_count > max_bad_count )
					{
						if( debug ) printf( "New iteration because singularity in Jacobian cannot be avoided\n" );
						rank = NO; break; /* New Iteration */
					}
//					lambda *= factor_squared; /* slow */
					lambda *= factor; /* fast */
					if( lambda <= upper_bound )
					{
						if( debug ) printf( "Increase lambda to %g to avoid singularity (attempt %d)\n", lambda, bad_count );
						add = YES; continue; /* INCREASE LAMBDA PARAMETER AND TRY AGAIN */
					}
					else
					{
						if( debug ) printf( "Singularity is detected in the Jacobian but lambda cannot be increased (%g > %g)\n", lambda, upper_bound );
						ier = 129; break; /* E X I T */
					}
				}
			}	/* END: The matrix is algorithmically NOT positive definite */
			else
			{
				/* The matrix is algorithmically positive definite */
				if( debug > 2 ) printf( "Hessian matrix is positive definite\n" );
				lu_elimination( hessian, x_update, nParam, x_update );
				if( singular_count != -99 ) singular_count = 0;
				for( j = 0; j < nParam; j++ ) x_new[j] = x[j] - x_update[j]; /* Test new parameters*/
//				for( j = 0; j < n; j++ ) printf( "lm %g ", x_new[j] ); printf( "\n" );
				if( func_data->cd->check_success ) func_data->success = 1;
				func( x_new, func_data, f_xpdx ); ieval++; /* CALCULATE SUM OF SQUARES */
				*phi = 0.0;
				for( i = 0; i < nObs; i++ ) *phi += f_xpdx[i] * f_xpdx[i];
				if( func_data->cd->check_success && func_data->success )
				{
					if( debug ) printf( "Predictions are within the predifined bounds!\n" );
					phi_best = *phi;
					for( j = 0; j < nParam; j++ ) x_best[j] = x_new[j];
					ier = 99; break; /* E X I T */
				}
				vphi[iphi++] = *phi;
				lambda_count++;
				if( debug ) printf( "phi %g lambda(%d) %g\n", *phi, lambda_count, lambda );
				if( iphi < 0 ) // AVOID THIS; not working properly at the moment
				{
					phi_range_max = phi_best; // initialize max search
					phi_range_min = HUGE_VAL; // initialize min search
					for( i = iphi - max_phi_range; i < iphi; i++ )
					{
//						printf( "p %g\n", vphi[i] );
						if( phi_range_max < vphi[i] ) phi_range_max = vphi[i];
						if( phi_range_min > vphi[i] ) phi_range_min = vphi[i];
					}
					temp = fabs( phi_range_min - phi_range_max ) / phi_range_min;
					if( temp < phi_range )
					{
						if( debug ) printf( "The last %d LM evaluations have OF in a close range (%g - %g; %g < %g)\n", max_phi_range, phi_range_min, phi_range_max, temp, phi_range );
						ier = 19; *infer += 16; break; /* E X I T */
					}
				}
				if( phi_best > *phi ) { phi_best = *phi; for( j = 0; j < nParam; j++ ) x_best[j] = x_new[j]; lambda_up_count = 0; }
				if( phi_old > *phi ) { for( j = 0; j < nParam; j++ ) x_cur_best[j] = x_new[j]; }
				if( *phi < phi_cutoff )
				{
					if( debug ) printf( "OF was successfully decreased below the predifined limit!\n" );
					ier = 29; *infer += 8; break; /* E X I T */
				}
				if( iopt )
				{
					if( *phi > phi_old )
					{
						bad_count++;
						if( bad_count > max_bad_count )
						{
							if( debug ) printf( "New iteration because OF cannot be reduced after %d attempts\n", max_bad_count );
							rank = NO; break; /* New Iteration */
						}
						lambda *= factor_squared;
						if( lambda <= upper_bound )
						{
							if( debug > 2 ) printf( "Increase lambda to %g to avoid increase in OF (%g > %g; attempt %d)\n", lambda, *phi, phi_old, bad_count );
//							rank = YES; break;
							add = YES; continue;
						}
						else
						{
							if( debug ) printf( "Marquardt lambda exceeded the maximally allowed value ..." );
							central_derivatives = YES;
							lambda_up_count++;
							if( lambda_up_count > 3 )
							{
								if( debug ) printf( "quitting after 3 attempts!\n" );
								ier = 39; break; /* E X I T */
							}
							if( debug ) printf( "try new iteration using central defivatives (attempt %d)\n", lambda_up_count );
							lambda /= factor;
							rank = NO; break; /* New Iteration */
						}
					}
					else if( bad_count == 0 )
					{
						lambda /= factor;
						if( debug > 2 ) printf( "Decrease lambda to %g because OF decreases (%g < %g)\n", lambda, *phi, phi_old );
						add = YES; continue; /* Same Iteration new lambda; NOT GOOD */
//						rank = YES; break; /* Same Iteration new lambda; EVEN WORST */
					}
					if( grad_norm_old > DBL_EPSILON )
					{
//						printf( "Lambda %g ", lambda );
						temp = grad_norm / grad_norm_old;
						if( grad_norm < grad_norm_old )      lambda *= MAX( factor_recipical, temp );
						else if( grad_norm > grad_norm_old ) lambda *= MIN( factor, temp );
						if( debug > 2 ) printf( "Lambda changed due to grad_norm %g (%g->%g)\n", lambda, grad_norm_old, grad_norm );
					}
					if( lambda < precision ) lambda = precision;
// 					add = YES; continue; // DOES NOT WORK
				}
				for( j = 0; j < nParam; j++ ) x[j] = x_cur_best[j];
				for( i = 0; i < nObs; i++ )
				{
					f_xmdx[i] = f[i];
					f[i] = f_xpdx[i];
				}
				for( j = 0; j < nParam; j++ )
				{
					x_abs = fabs( x[j] );
					x_diff = fabs( x_update[j] ) / MAX( x_abs, DAX );
					if( x_diff > relcon ) break;
				}
				if( j == nParam ) *infer += 1;
				phi_diff = fabs( *phi - phi_old ) / MAX( phi_old, DAX );
//				if( phi_diff <= eps ) *infer += 2;
				if( lambda > 1e2 && lambda_count > max_lambda_count ) { rank = NO; if( debug ) printf( "New interation because the lambda is large and there were enough attempts (%d < %d)\n", max_lambda_count, lambda_count ); break; }
				rank = YES; /* YES == standard LM; NO == may help to avoid local minimum */
				break;
			}	/* END: The matrix is algorithmically positive definite */
		}	/* END OF THE LOCAL WHILE */
	}	/* END OF THE MAIN LOOP */
	if( ier == 1 ) ier = 0;
	temp = HUGE_VAL;                     /* OUTPUT grad_norm, IEVAL, NSIG, AL, AND ITER */
	for( j = 0; j < nParam; j++ )
	{
		x_old = fabs( x_update[j] );
		if( x_old > 0.0 )
		{
			dx = log10( MAX( fabs( x[j] ), DAX ) ) - log10( x_old );
			if( temp < dx ) temp = dx;
		}
	}
	if( debug )
	{
		if( *infer == 0 ) printf( "No convergence criteria is satisfied!\n" );
		else
		{
			printf( "Convergence criteria is satisfied:\n" );
			if( *infer &  1 ) printf( "For two successive iterations, parameter estimates agree to %d leading digits [nsig]\n", nsig );
			if( *infer &  2 ) printf( "For two successive iterations, the relative OF difference is less than or equal to %g [eps]\n", eps );
			if( *infer &  4 ) printf( "Gradient norm is less than or equal to delta (%g < %g); gradient norm close to zero indicate that the current predictions are close to zero\n", grad_norm, delta );
			if( *infer &  8 ) printf( "OF (%g) is less the predifined cutoff limit (%g)\n", phi_best, phi_cutoff );
			if( *infer & 16 ) printf( "Last %d LM evaluations have OF in a close range (%g - %g; %g < %g)\n", max_phi_range, phi_range_min, phi_range_max, fabs( phi_range_min - phi_range_max ) / phi_range_min, phi_range );
		}
		printf( "\nInitial OF %g\n", phi_init );
		func( x_best, func_data, f_xpdx );
		*phi = 0.0;
		for( i = 0; i < nObs; i++ ) *phi += f_xpdx[i] * f_xpdx[i];
		printf( "Best OF %g\n", *phi );
	}
	*phi = func_data->phi = phi_best;
	for( j = 0; j < nParam; j++ ) x[j] = x_best[j];
	if( func_data->cd->odebug == 1 && ofe_close == 1 ) { fclose( func_data->f_ofe ); func_data->f_ofe = NULL; }
	free( vphi ); free( hessian ); free( grad ); free( x_update ); free( scale ); free( x_new ); free( x_cur_best ); free( x_best ); free( x_bad ); free( f_xpdx ); free( f_xmdx );
	return( ier );
}
