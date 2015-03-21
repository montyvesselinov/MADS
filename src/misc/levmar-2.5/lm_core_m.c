/////////////////////////////////////////////////////////////////////////////////
//
//  Levenberg - Marquardt non-linear minimization algorithm
//  Copyright (C) 2004  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////
//
// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
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

/* precision-specific definitions */
#include <string.h>
#include <math.h>
#include "../../mads.h"

#define LEVMAR_DER LM_ADD_PREFIX(levmar_der)
#define LEVMAR_DER2 LM_ADD_PREFIX(levmar_der2)
#define LEVMAR_DIF LM_ADD_PREFIX(levmar_dif)
#define LEVMAR_FDIF_FORW_JAC_APPROX LM_ADD_PREFIX(levmar_fdif_forw_jac_approx)
#define LEVMAR_FDIF_CENT_JAC_APPROX LM_ADD_PREFIX(levmar_fdif_cent_jac_approx)
#define LEVMAR_TRANS_MAT_MAT_MULT LM_ADD_PREFIX(levmar_trans_mat_mat_mult)
#define LEVMAR_L2NRMXMY LM_ADD_PREFIX(levmar_L2nrmxmy)
#define LEVMAR_COVAR LM_ADD_PREFIX(levmar_covar)

#ifdef HAVE_LAPACK
#define AX_EQ_B_LU LM_ADD_PREFIX(Ax_eq_b_LU)
#define AX_EQ_B_CHOL LM_ADD_PREFIX(Ax_eq_b_Chol)
#define AX_EQ_B_QR LM_ADD_PREFIX(Ax_eq_b_QR)
#define AX_EQ_B_QRLS LM_ADD_PREFIX(Ax_eq_b_QRLS)
#define AX_EQ_B_SVD LM_ADD_PREFIX(Ax_eq_b_SVD)
#define AX_EQ_B_BK LM_ADD_PREFIX(Ax_eq_b_BK)
#else
#define AX_EQ_B_LU LM_ADD_PREFIX(Ax_eq_b_LU_noLapack)
#endif /* HAVE_LAPACK */

int func_set( int n_sub, double *var_mat[], double *phi, double *f[], FILE *out, struct opt_data *op ); // parallel lambda search

/*
 * This function seeks the parameter vector p that best describes the measurements vector x.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized.
 *
 * This function requires an analytic Jacobian. In case the latter is unavailable,
 * use LEVMAR_DIF() bellow
 *
 * Returns the number of iterations (>=0) if successful, LM_ERROR if failed
 *
 * For more details, see K. Madsen, H.B. Nielsen and O. Tingleff's lecture notes on
 * non-linear least squares at http://www.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 */

// MADS calls LEVMAR_DER by DEFAULT

int LEVMAR_DER2(
	void ( *func )( LM_REAL *par_current, LM_REAL *obs_current, int nP, int nO, void *adata ), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
	void ( *jacf )( LM_REAL *par_current, LM_REAL *f, LM_REAL *j, int nP, int nO, void *adata ), /* function to evaluate the Jacobian \part x / \part p */
	LM_REAL *par_current,         /* I/O: initial parameter estimates. On output has the estimated solution */
	LM_REAL *obs_target,         /* I: measurement vector. NULL implies a zero vector */
	int nP,              /* I: parameter vector dimension (i.e. #unknowns) */
	int nO,              /* I: measurement vector dimension */
	int maxjac,          /* I: maximum number of iterations */
	LM_REAL opts[4],    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
		 * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used
		 */
	LM_REAL info[LM_INFO_SZ],
	/* O: information regarding the minimization. Set to NULL if don't care
	 * info[0]= ||e||_2 at initial p.
	 * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
	 * info[5]= # iterations,
	 * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
	 *                                 2 - stopped by small Dp
	 *                                 3 - stopped by itmax
	 *                                 4 - singular matrix. Restart from current p with increased mu
	 *                                 5 - no further error reduction is possible. Restart with increased mu
	 *                                 6 - stopped by small ||e||_2
	 *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
	 * info[7]= # function evaluations
	 * info[8]= # Jacobian evaluations
	 * info[9]= # linear systems solved, i.e. # attempts for reducing error
	 */
	LM_REAL *work,     /* working memory at least LM_DER_WORKSZ() reals large, allocated if NULL */
	LM_REAL *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */
	void *adata )       /* pointer to possibly additional data, passed uninterpreted to func & jacf.
		 * Set to NULL if not needed
		 */
{
	register int i, j, loop_count, l;
	int worksz, freework = 0, issolved, issolved1 = 0, success, odebug, change, computejac, changejac, maxnfev;
	struct opt_data *op = ( struct opt_data * ) adata;
	char filename[255];
	time_t time_start, time_end, time_elapsed, time_jacobian = 0, time_lambda = 0;
	int parallel_lambda_count = 0;
	/* temp work arrays */
	LM_REAL *obs_error,          /* nx1 */
			*obs_current,         /* \hat{x}_i, nx1 */
			*hx1,        /* used in acceleration comp (nx1) */
			*hx2,        /* used in acceleration comp (nx1) */
			*jac_min,
			*jac_max,
			*par_lam_last, *par_jac_last, *par_init,     /* old parameter set */
			*par_best, /* best parameter set */
			*obs_best, /* best observation set */
			*JTe,      /* J^T e_i mx1 */
			*jac,        /* nxm */
			*JTJ,    /* mxm */
			*par_change,         /* mx1 (=v) */
			*ephdp_plus,     /* residual used in acceleration computation (nx1) */
			*ephdp_minus,     /* residual used in acceleration computation (nx1) */
			*vvddr,      /* used to compute acceleration (nx1) */
			*jacTvv,     /* jacT*vvddr, mx1 */
			*acceleration,          /* acceleration (mx1) */
			*JTJ_diag,   /* diagonal of J^T J, mx1 */
			*par_update,        /* p + Dp, mx1 */
			*phDp_plus,       /* p + hDp, mx1 */
			*phDp_minus,      /* p - hDp, mx1 */
			*obs_update,        /* nx1 */
			*obs_error_update;       /* nx1, used only for holding a temporary e vector and when differentiating with central differences */
	int *jac_zero, *jac_zero_obs;
	int skipped, first, allzero;
	gsl_matrix *gsl_jacobian;
	int using_ffdif = 1;
	double acc_h = op->cd->lm_h;
	double lm_ratio = op->cd->lm_ratio; // alpha
	register LM_REAL lambda,  /* damping constant */
			 tmp = 0, /* mainly used in matrix & vector multiplications */
			 avRatio = 0.0; /* acceleration/velocity */
	LM_REAL phi_current, phi_jac_last, JTe_inf_norm, phi_update; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
	LM_REAL par_L2_norm, par_change_L2_norm = LM_REAL_MAX, acceleration_L2_norm, phi_change, par_change_L2_norm_sq = LM_REAL_MAX, dL;
	LM_REAL tau, eps1, eps2, eps2_sq, eps3, delta;
	LM_REAL phi_init, phi_best;
	int nu, nu2, stop = 0, nfev, njac = 0, nlss = 0, max_num_of_lambda_searches, jac_update_count, par_update_accepted = 1, newjac;
	max_num_of_lambda_searches = ( nP >= 10 ) ? nP : 10;
	if( op->cd->lm_num_lambda_searches > max_num_of_lambda_searches ) max_num_of_lambda_searches = op->cd->lm_num_lambda_searches;
	if( op->cd->lm_num_parallel_lambda > max_num_of_lambda_searches ) max_num_of_lambda_searches = op->cd->lm_num_parallel_lambda;
	if( op->cd->lm_nlamof > max_num_of_lambda_searches ) max_num_of_lambda_searches = op->cd->lm_nlamof;
	LM_REAL phi_jac_vector[maxjac], phi_lam_vector[max_num_of_lambda_searches];
	int phi_jac_count, phi_lam_count;
	int mu_big = 0, phi_decline = 0;
	int npl, ipar_max, ipar_min;
	double max_change, min_change, p_diff, fj;
	int imax, imin, omax, omin, ok;
	double max = 0, min = HUGE_VAL;
	double *phi_vector, **param_matrix, **obs_matrix;
	const int nm = nO * nP;
	int ( *linsolver )( LM_REAL * A, LM_REAL * B, LM_REAL * obs_target, int nP ) = NULL;
	lambda = JTe_inf_norm = par_L2_norm = 0.0; /* -Wall */
	jac_update_count = newjac = 0; /* -Wall */
	if( nO < nP )
	{
		// fprintf( stderr, LCAT( LEVMAR_DER, "(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n" ), n, m );
		// return LM_ERROR;
	}
	if( opts )
	{
		tau = opts[0];
		eps1 = opts[1];
		eps2 = opts[2];
		eps2_sq = opts[2] * opts[2];
		eps3 = opts[3];
		delta = opts[4];
		if( delta < 0.0 )
		{
			delta = -delta; /* make positive */
			using_ffdif = 0; /* use central differencing */
		}
	}
	else  // use default values
	{
		tau = LM_CNST( LM_INIT_MU );
		eps1 = LM_CNST( LM_STOP_THRESH );
		eps2 = LM_CNST( LM_STOP_THRESH );
		eps2_sq = LM_CNST( LM_STOP_THRESH ) * LM_CNST( LM_STOP_THRESH );
		eps3 = LM_CNST( LM_STOP_THRESH );
		delta = LM_CNST( LM_DIFF_DELTA );
	}
	if( !work )
	{
		worksz = LM_DIF_WORKSZ( nP, nO ); //4*n+4*nP + n*m + m*m;
		work = ( LM_REAL * )malloc( worksz * sizeof( LM_REAL ) ); /* allocate a big chunk in one step */
		if( !work )
		{
			tprintf( LCAT( LEVMAR_DIF, "ERROR (): memory allocation request failed\n" ) );
			return LM_ERROR;
		}
		freework = 1;
	}
	/* set up work arrays */
	obs_error = work;
	obs_current = obs_error + nO;
	JTe = obs_current + nO;
	jac = JTe + nP;
	JTJ = jac + nm;
	par_change = JTJ + nP * nP;
	JTJ_diag = par_change + nP;
	par_update = JTJ_diag + nP;
	obs_update = par_update + nP;
	obs_error_update = obs_update + nO;
	hx1 = ( LM_REAL * )malloc( nO * sizeof( LM_REAL ) );
	hx2 = ( LM_REAL * )malloc( nO * sizeof( LM_REAL ) );
	phDp_plus = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	phDp_minus = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	ephdp_plus = ( LM_REAL * )malloc( nO * sizeof( LM_REAL ) );
	ephdp_minus = ( LM_REAL * )malloc( nO * sizeof( LM_REAL ) );
	vvddr = ( LM_REAL * )malloc( nO * sizeof( LM_REAL ) );
	jacTvv = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	acceleration = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	jac_min = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	jac_max = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	jac_zero = ( int * )malloc( nP * sizeof( int ) );
	jac_zero_obs = ( int * )malloc( nO * sizeof( int ) );
	par_lam_last = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	par_jac_last = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	par_best = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	par_init = ( LM_REAL * )malloc( nP * sizeof( LM_REAL ) );
	obs_best = ( LM_REAL * )malloc( nO * sizeof( LM_REAL ) );
	if( op->cd->lm_eigen ) { gsl_jacobian = gsl_matrix_alloc( op->od->nTObs, op->pd->nOptParam ); }
	else gsl_jacobian = NULL;
	/* compute e=x - f(p) and its L2 norm */
	maxnfev = op->cd->maxeval - op->cd->neval;
	for( i = 0; i < nP; i++ )
		par_lam_last[i] = par_jac_last[i] = par_init[i] = par_best[i] = par_update[i] = par_current[i];
	( *func )( par_current, obs_current, nP, nO, adata ); nfev = 1;
	for( i = 0; i < nP; i++ )
		par_best[i] = par_current[i];
	for( i = 0; i < nO; i++ )
		obs_best[i] = obs_current[i];
	if( op->cd->lm_num_parallel_lambda > 0 )
	{
		if( ( phi_vector = ( double * ) malloc( op->cd->lm_num_parallel_lambda * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( param_matrix = double_matrix( op->cd->lm_num_parallel_lambda, nP ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( obs_matrix = double_matrix( op->cd->lm_num_parallel_lambda, nO ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		for( npl = 0; npl < op->cd->lm_num_parallel_lambda; npl++ )
		{
			phi_vector[npl] = 0;
			for( i = 0; i < nP; i++ )
				param_matrix[npl][i] = par_current[i];
			for( i = 0; i < nO; i++ )
				obs_matrix[npl][i] = obs_current[i];
		}
	}
	if( op->cd->check_success && op->success )
	{
		if( op->cd->ldebug ) tprintf( "SUCCESS: Model predictions are within predefined calibration ranges\n" );
		stop = 8;
	}
#ifdef HAVE_LAPACK
	/* 6 alternatives are available: LU, Cholesky, 2 variants of QR decomposition, SVD and LDLt.
	 * Cholesky is the fastest but might be inaccurate; QR is slower but more accurate;
	 * SVD is the slowest but most accurate; LU offers a tradeoff between accuracy and speed
	 */
	//linsolver = AX_EQ_B_BK; if( op->cd->ldebug ) tprintf( "BK decomposition\n" );
	//linsolver = AX_EQ_B_LU; if( op->cd->ldebug ) tprintf( "LU decomposition\n" );
	//linsolver = AX_EQ_B_CHOL; if( op->cd->ldebug ) tprintf( "Cholesky decomposition\n" );
	//linsolver = AX_EQ_B_QR; if( op->cd->ldebug ) tprintf( "QR decomposition\n" );
	//linsolver = (int (*)(LM_REAL *A, LM_REAL *B, LM_REAL *x, int m))AX_EQ_B_QRLS; if( op->cd->ldebug ) tprintf( "QRLS decomposition\n" );
	linsolver = AX_EQ_B_SVD; if( op->cd->ldebug ) tprintf( "LM using SVD decomposition\n" );
#else
	/* use the LU included with levmar */
	linsolver = AX_EQ_B_LU; if( op->cd->ldebug ) tprintf( "LM using LU decomposition\n" );
#endif /* HAVE_LAPACK */
	/* ### e=x-hx, p_eL2=||e|| */
#if 0
	phi_current = LEVMAR_L2NRMXMY( obs_error, obs_target, obs_current, nO );
#else
	for( i = 0, phi_current = 0.0; i < nO; i++ )
	{
		obs_error[i] = tmp = obs_target[i] - obs_current[i];
		phi_current += tmp * tmp;
	}
#endif
	if( op->cd->ldebug )
	{
		if( op->cd->lm_num_parallel_lambda )
		{
			tprintf( "LM with parallel lambda search\n" );
			tprintf( "LM lambda search will be performed for a series of %d lambdas in parallel.\n", op->cd->lm_num_parallel_lambda );
			tprintf( "LM lambda values will be increased/decrease by a factor of %g\n", op->cd->lm_mu );
			tprintf( "LM acceleration will not be applied even if selected!\n" );
			op->cd->lm_acc = 0;
		}
		else
		{
			tprintf( "LM without parallel lambda search\n" );
			if( op->cd->num_proc > 0 )
				tprintf( "WARNING: LM may perform better if parallel lambda search is evoked (nplambda>0)\n" );
		}
		if( op->cd->lm_acc ) tprintf( "LM with acceleration\n" );
		else tprintf( "LM without acceleration\n" );
		if( op->cd->lm_indir ) tprintf( "LM with indirect computation of lambda changes\n" );
		else tprintf( "LM with direct computation of lambda changes\n" );
	}
	phi_init = phi_best = phi_jac_last = phi_current;
	if( !LM_FINITE( phi_current ) ) stop = 7;
	nu = 32; /* force computation of J */
	if( op->cd->check_success ) success = 1; else success = 0;
	if( op->cd->odebug ) odebug = 1; else odebug = 0;
	computejac = 0;
	phi_jac_count = phi_lam_count = 0;
	if( op->cd->ldebug > 2 )
	{
		DeTransform( par_current, op, jac_min );
		tprintf( "Initial parameter estimates:\n" );
		for( i = 0; i < nP; i++ )
		{
			j = op->pd->var_index[i];
			tprintf( "%-20s:", op->pd->var_name[j] );
			if( op->pd->var_log[j] ) tprintf( " %g", pow( 10, jac_min[i] ) );
			else tprintf( " %g", jac_min[i] );
			tprintf( "\n" );
		}
	}
	if( op->cd->ldebug ) tprintf( "Initial OF %g\n", phi_current );
	else if( op->cd->lmstandalone ) tprintf( "OF %g -> ", phi_current );
	loop_count = -1;
	while( !stop )
	{
		sprintf( filename, "%s.quit", op->root );
		if( Ftest( filename ) == 0 ) { tprintf( "MADS quits! Termination request! %s\n", filename ); op->cd->quit = 1; stop = 99; break; }
		loop_count++;
		/* Note that p and e have been updated at a previous iteration */
		if( phi_current <= eps3 ) /* error is small */ // BELOW a cutoff value
		{
			if( op->cd->ldebug ) tprintf( "CONVERGED: OF below cutoff value (%g < %g)\n", phi_current, eps3 );
			stop = 6;
			break;
		}
		/* Compute the Jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
		 * The symmetry of J^T J is again exploited for speed
		 */
		if( ( par_update_accepted && nu > 16 ) || jac_update_count >= max_num_of_lambda_searches || mu_big || phi_decline || computejac || op->cd->lm_num_parallel_lambda + phi_lam_count >= max_num_of_lambda_searches ) /* compute difference approximation to J */
		{
			if( op->cd->ldebug && loop_count != 0 )
			{
				if( loop_count != 0 ) tprintf( "New Jacobian requested because:" );
				if( par_update_accepted && nu > 16 && loop_count != 0 ) tprintf( " > Lambda multiplication factor too large (nu = %d > 16); ", nu );
				if( jac_update_count >= max_num_of_lambda_searches ) tprintf( " > Maximum number of lambda iteration is reached (%d); ", max_num_of_lambda_searches );
				if( op->cd->lm_num_parallel_lambda > 0 )
				{
					if( op->cd->lm_num_parallel_lambda + phi_lam_count >= max_num_of_lambda_searches ) tprintf( " > Maximum number of lambda iteration will be reached if one more parallel search is performed (%d); ", max_num_of_lambda_searches );
				}
				else
				{
					if( phi_lam_count >= max_num_of_lambda_searches ) tprintf( " > Maximum number of lambda iteration is reached (%d); ", max_num_of_lambda_searches );
				}
				if( mu_big ) tprintf( " > Lambda is constrained (%g); ", 1e3 );
				if( phi_decline ) tprintf( " > OF estimate declined substantially (%g << %g)", phi_current, phi_jac_last );
				if( computejac ) tprintf( " > Linear solve OF estimates do not change substantially" );
				tprintf( "\n\n" );
			}
			if( njac >= maxjac ) { stop = 31; continue; }
			if( nfev >= maxnfev ) { stop = 32; continue; }
			computejac = 0;
			changejac = 1;
			mu_big = 0;
			phi_decline = 0;
			phi_lam_count = 0;
			phi_jac_last = phi_current;
			phi_jac_vector[phi_jac_count++] = phi_current;
			if( phi_jac_count > op->cd->lm_njacof )
			{
				tmp = HUGE_VAL;
				for( i = 0; i < phi_jac_count; i++ )
				{
					if( phi_jac_vector[i] < tmp ) tmp = phi_jac_vector[i];
					// tprintf( " %g", phi1[i] );
				}
				// tprintf( "d %g %d\n", tmp, phi1c );
				j = 0;
				for( i = 0; i < phi_jac_count; i++ )
				{
					// tprintf( "d %g %g %g\n", tmp, phi1[i], ( phi1[i] - tmp ) / tmp );
					if( ( phi_jac_vector[i] - tmp ) / tmp < 1 ) j++;
				}
				// tprintf( "d %g %d\n", tmp, j );
				if( j >= op->cd->lm_njacof )
				{
					if( op->cd->ldebug ) tprintf( "\nCONVERGED: %d Jacobian OF estimates are very close to the best current OF (%g)\n", j, tmp );
					stop = 9;
					break;
				}
			}
			if( success ) op->cd->check_success = 0;
			if( odebug ) op->cd->odebug = 0;
			time_start = time( NULL );
			if( using_ffdif ) /* use forward differences */
			{
				if( op->cd->ldebug > 5 ) tprintf( "Jacobian computed using forward differences\n" );
				jacf( par_current, obs_current, jac, nP, nO, adata );
				njac++; nfev += nP;
			}
			else /* use central differences */
			{
				if( op->cd->ldebug > 5 ) tprintf( "Jacobian computed using central differences\n" );
				LEVMAR_FDIF_CENT_JAC_APPROX( func, par_current, obs_update, obs_error_update, delta, jac, nP, nO, adata );
				njac++; nfev += 2 * nP;
			}
			time_end = time( NULL );
			time_elapsed = time_end - time_start;
			time_jacobian += time_elapsed;
			if( op->cd->tdebug )
			{
				if( time_elapsed > 86400 ) tprintf( "Jacobian PT = %g days (average %g days)\n", ( ( double ) time_elapsed / 86400 ), ( ( double ) time_jacobian / njac / 86400 ) );
				else if( time_elapsed > 3600 ) tprintf( "Jacobian PT = %g hours (average %g hours)\n", ( ( double ) time_elapsed / 3600 ), ( ( double ) time_jacobian / njac / 3600 ) );
				else if( time_elapsed > 60 ) tprintf( "Jacobian PT = %g minutes (average %g minutes)\n", ( ( double ) time_elapsed / 60 ), ( ( double ) time_jacobian / njac / 60 ) );
				else tprintf( "Jacobian PT = %ld seconds (average %g seconds)\n", time_elapsed, ( ( double ) time_jacobian / njac ) );
			}
			if( success ) op->cd->check_success = success;
			if( odebug ) op->cd->odebug = odebug;
			nu = op->cd->lm_nu; jac_update_count = 0; par_update_accepted = 0; newjac = 1;
			if( op->cd->ldebug )
			{
				if( op->cd->ldebug > 1 ) tprintf( "\n\n" );
				if( loop_count == 0 ) tprintf( "Jacobian #%d: Linear solves %d Evaluations %d OF %g lambda %g\n", njac, nlss, nfev, phi_current, tau );
				else tprintf( "Jacobian #%d: Linear solves %d Evaluations %d OF %g lambda %g\n", njac, nlss, nfev, phi_current, lambda );
				if( fabs( phi_best - phi_current ) > COMPARE_EPSILON ) tprintf( "Best OF %g\n", phi_best );
			}
			if( op->cd->ldebug >= 4 )
			{
				if( op->cd->ldebug >= 6 ) tprintf( "\nWeighted Jacobian matrix:\n" );
				else tprintf( "\nWeighted Jacobian matrix (observations with zero sensitivities are skipped):\n" );
				char *s;
				tprintf( "Parameters:" );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					s = op->pd->var_name[op->pd->var_index[i]];
					if( strlen( s ) < 7 ) tprintf( " %s", s );
					else tprintf( " p%d", i + 1 );
				}
				for( l = j = 0; j < op->od->nTObs; j++ )
				{
					if( op->cd->ldebug >= 5 || op->od->nTObs < 30 || ( j < 10 || j > op->od->nTObs - 10 ) )
					{
						int print = 0;
						if( op->od->obs_weight[j] <= DBL_EPSILON ) print = 0;
						else if( op->cd->ldebug >= 10 ) print = 1;
						else
						{
							for( i = 0; i < op->pd->nOptParam; i++, l++ )
								if( fabs( jac[l] ) > DBL_EPSILON ) print = 1;
							l -= op->pd->nOptParam; // reset jacobian counter
						}
						if( print )
						{
							tprintf( "\n%-12s: ", op->od->obs_id[j] );
							for( i = 0; i < op->pd->nOptParam; i++, l++ )
								tprintf( " %g", jac[l] );
						}
					}
					else l += op->pd->nOptParam;
					if( op->cd->ldebug == 4 && op->od->nTObs >= 30 && j == 11 ) tprintf( "\n..." );
				}
				tprintf( "\n" );
			}
			else if( op->cd->ldebug ) tprintf( "Jacobian matrix analysis ...\n" );
			// Jacobian matrix analyses: check for zero sensitivities & calculate absolute sensitivity ranges
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				jac_max[i] = jac_zero[i] = 0;
				jac_min[i] = HUGE_VAL;
			}
			for( j = 0; j < op->od->nTObs; j++ ) jac_zero_obs[j] = 0;
			max = 0; min = HUGE_VAL;
			skipped = 0;
			first = 1;
			imax = imin = omax = omin = 0;
			for( l = j = 0; j < op->od->nTObs; j++ )
			{
				if( op->od->obs_weight[j] > DBL_EPSILON )
				{
					allzero = 1;
					for( i = 0; i < op->pd->nOptParam; i++, l++ )
					{
						fj = fabs( jac[l] );
						if( fj < DBL_EPSILON ) { jac_zero[i]++; jac_zero_obs[j]++; }
						else
						{
							allzero = 0;
							if( jac_max[i] < fj ) jac_max[i] = fj;
							if( jac_min[i] > fj ) jac_min[i] = fj;
							if( max < fj ) { max = fj; imax = i; omax = j; }
							if( min > fj ) { min = fj; imin = i; omin = j; }
						}
					}
					if( allzero )
					{
						skipped++;
						if( op->cd->ldebug > 2 )
						{
							if( first ) { first = 0; tprintf( "WARNING: Observation(s) with zero sensitivities:" ); }
							tprintf( " %s", op->od->obs_id[j] );
						}
					}
				}
				else l += op->pd->nOptParam;
			}
			if( skipped && op->cd->ldebug > 2 ) tprintf( "\n" );
			if( skipped )
			{
				if( op->cd->ldebug )
				{
					if( skipped > 1 ) tprintf( "WARNING: %d observations have zero sensitivities!\n", skipped );
					else tprintf( "WARNING: %d observation has zero sensitivities!\n", skipped );
				}
				if( skipped >= op->od->nCObs )
				{
					tprintf( "\nERROR: All the calibration targets (%d) have zero sensitivities! LM quits! Potential model setup error or complex parameter space!\n", op->od->nCObs );
					stop = 10;
					break;
				}
			}
			// Jacobian matrix analyses: check for insensitive observations (redundant with above; not needed; for testing only)
			if( op->cd->ldebug >= 11 )
			{
				for( j = 0; j < op->od->nTObs; j++ )
					if( op->od->obs_weight[j] > DBL_EPSILON && jac_zero_obs[j] >= op->pd->nOptParam )
						tprintf( "WARNING: Model prediction \'%s\' is not impacted by model parameters!\n", op->od->obs_id[j] );
			}
			// Jacobian matrix analyses: check for insensitive parameters
			ok = 0;
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				if( jac_max[i] > DBL_EPSILON ) ok++;
				else { if( op->cd->ldebug ) tprintf( "WARNING: Model parameter \'%s\' is not impacting model predictions!\n", op->pd->var_name[op->pd->var_index[i]] ); }
			}
			if( ok == 0 )
			{
				if( op->cd->ldebug ) tprintf( "ERROR: None of the model parameters is impacting model predictions!\n" );
				stop = 10;
				break;
			}
			// Jacobian matrix analyses: print absolute sensitivity ranges
			if( op->cd->ldebug > 1 )
			{
				tprintf( "Highest absolute sensitivity (parameter #%d \'%s\' vs observation #%d \'%s\'): %g\n", imax + 1, op->pd->var_name[op->pd->var_index[imax]], omax + 1, op->od->obs_id[omax], max );
				tprintf( "Lowest  absolute sensitivity (parameter #%d \'%s\' vs observation #%d \'%s\'): %g\n", imin + 1, op->pd->var_name[op->pd->var_index[imin]], omin + 1, op->od->obs_id[omin], min );
				if( op->cd->ldebug > 2 )
				{
					char *s;
					tprintf( "\nTable of parameter sensitivities to observations:" );
					tprintf( "\nModel parameters          :" );
					for( i = 0; i < op->pd->nOptParam; i++ )
					{
						s = op->pd->var_name[op->pd->var_index[i]];
						if( strlen( s ) < 7 ) tprintf( " %s", s );
						else tprintf( " p%d", i );
					}
					tprintf( "\n" );
					tprintf( "Sensitivity zeros (out of %d):", op->od->nCObs ); for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( " %d", jac_zero[i] ); tprintf( "\n" );
					tprintf( "Min absolute sensitivities:" ); for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( " %g", jac_min[i] ); tprintf( "\n" );
					tprintf( "Max absolute sensitivities:" ); for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( " %g", jac_max[i] ); tprintf( "\n" );
				}
			}
			// Jacobian matrix analyses: calculate & print overall sensitivity ranges
			if( op->cd->ldebug > 3 )
			{
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					jac_max[i] = -HUGE_VAL;
					jac_min[i] = HUGE_VAL;
				}
				for( l = j = 0; j < op->od->nTObs; j++ )
				{
					if( op->od->obs_weight[j] <= DBL_EPSILON ) { l += op->pd->nOptParam; continue; }
					for( i = 0; i < op->pd->nOptParam; i++, l++ )
					{
						if( jac_max[i] < jac[l] ) jac_max[i] = jac[l];
						if( jac_min[i] > jac[l] ) jac_min[i] = jac[l];
					}
				}
				tprintf( "Min overall  sensitivities:" ); for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( " %g", jac_min[i] ); tprintf( "\n" );
				tprintf( "Max overall  sensitivities:" ); for( i = 0; i < op->pd->nOptParam; i++ ) tprintf( " %g", jac_max[i] ); tprintf( "\n\n" );
			}
			if( op->cd->lm_eigen )
			{
				for( l = j = 0; j < op->od->nTObs; j++ )
					for( i = 0; i < op->pd->nOptParam; i++ )
						gsl_matrix_set( gsl_jacobian, j, i, jac[l++] ); // LEVMAR is using different jacobian order
				DeTransform( par_current, op, jac_min );
				for( i = 0; i < op->pd->nOptParam; i++ )
					op->pd->var[op->pd->var_index[i]] = jac_min[i];
				op->cd->lm_eigen--;
				op->phi = phi_current;
				eigen( op, obs_current, gsl_jacobian, NULL );
				save_results( 0, "", op, op->gd );
				op->cd->lm_eigen++;
			}
		}
		if( newjac ) /* Jacobian has changed, recompute J^T J, J^t e, etc */
		{
			newjac = 0;
			/* J^T J, J^T e */
			if( nm <= __BLOCKSZ__SQ ) // this is a small problem
			{
				/* J^T*J_ij = \sum_l J^T_il * J_lj = \sum_l J_li * J_lj.
				 * Thus, the product J^T J can be computed using an outer loop for
				 * l that adds J_li*J_lj to each element ij of the result. Note that
				 * with this scheme, the accesses to J and JtJ are always along rows,
				 * therefore induces less cache misses compared to the straightforward
				 * algorithm for computing the product (i.e., l loop is innermost one).
				 * A similar scheme applies to the computation of J^T e.
				 * However, for large minimization problems (i.e., involving a large number
				 * of unknowns and measurements) for which J/J^T J rows are too large to
				 * fit in the L1 cache, even this scheme incures many cache misses. In
				 * such cases, a cache-efficient blocking scheme is preferable.
				 *
				 * Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
				 * performance problem.
				 *
				 * Note that the non-blocking algorithm is faster on small
				 * problems since in this case it avoids the overheads of blocking.
				 */
				register int l, im;
				register LM_REAL alpha, *jaclm;
				for( i = nP * nP; i-- > 0; ) /* looping downwards saves a few computations */
					JTJ[i] = 0.0;
				for( i = nP; i-- > 0; )
					JTe[i] = 0.0;
				for( l = nO; l-- > 0; )
				{
					jaclm = jac + l * nP;
					for( i = nP; i-- > 0; )
					{
						im = i * nP;
						alpha = jaclm[i]; //jac[l*nP+i];
						for( j = i + 1; j-- > 0; ) /* j<=i computes lower triangular part only */
							JTJ[im + j] += jaclm[j] * alpha; //jac[l*nP+j]
						JTe[i] += alpha * obs_error[l]; /* J^T e */
					}
				}
				for( i = nP; i-- > 0; ) /* copy to upper part */
					for( j = i + 1; j < nP; ++j )
						JTJ[i * nP + j] = JTJ[j * nP + i];
			}
			else // this is a large problem
			{
				LEVMAR_TRANS_MAT_MAT_MULT( jac, JTJ, nO, nP ); /* Cache efficient computation of J^T J based on blocking */
				for( i = 0; i < nP; i++ ) /* cache efficient computation of J^T e */
					JTe[i] = 0.0;
				for( i = 0; i < nO; i++ )
				{
					register LM_REAL *jacrow;
					for( l = 0, jacrow = jac + i * nP, tmp = obs_error[i]; l < nP; ++l )
						JTe[l] += jacrow[l] * tmp;
				}
			}
			for( i = 0, par_L2_norm = JTe_inf_norm = 0.0; i < nP; i++ ) /* Compute ||J^T e||_inf and ||p||^2 */
			{
				if( JTe_inf_norm < ( tmp = FABS( JTe[i] ) ) ) JTe_inf_norm = tmp;
				JTJ_diag[i] = JTJ[i * nP + i]; /* save diagonal entries so that augmentation can be later canceled */
				par_L2_norm += par_current[i] * par_current[i];
			}
			if( op->cd->ldebug && ( changejac && ( loop_count > 0 ) ) )
			{
				DeTransform( par_current, op, jac_min );
				DeTransform( par_jac_last, op, jac_max );
				max_change = 0; min_change = HUGE_VAL;
				ipar_min = ipar_max = 0;
				for( i = 0 ; i < nP; i++ )
				{
					j = op->pd->var_index[i];
					if( op->pd->var_log[j] ) p_diff = fabs( pow( 10, jac_max[i] ) - pow( 10, jac_min[i] ) );
					else p_diff = fabs( jac_max[i] - jac_min[i] );
					if( max_change < p_diff ) { max_change = p_diff; ipar_max = i; }
					if( min_change > p_diff ) { min_change = p_diff; ipar_min = i; }
				}
				if( max_change < DBL_EPSILON )
					tprintf( "Parameters did not change between last two jacobian iterations.\n" );
				else
				{
					tprintf( "Parameter with maximum absolute estimate change between last two jacobian iterations: %-30s (%g)\n", op->pd->var_name[op->pd->var_index[ipar_max]], max_change );
					tprintf( "Parameter with minimum absolute estimate change between last two jacobian iterations: %-30s (%g)\n", op->pd->var_name[op->pd->var_index[ipar_min]], min_change );
				}
				for( i = 0 ; i < nP; i++ ) /* update parameter estimates */
					par_jac_last[i] = par_current[i];
			}
			/* check for convergence */
			if( op->cd->ldebug > 6 ) tprintf( "\nTest for convergence: Magnitude of the largest J^T e component (||J^T e||_inf = %g < %g to converge)\n", JTe_inf_norm, eps1 );
			if( ( JTe_inf_norm <= eps1 ) )
			{
				par_change_L2_norm = 0.0; /* no increment for p in this case */
				if( op->cd->ldebug ) tprintf( "CONVERGED: Largest J^T e component ||J^T e||_inf is too small (%g < %g)\n", JTe_inf_norm, eps1 );
				stop = 1;
				break;
			}
			changejac = 0;
		}
		if( loop_count == 0 ) // compute initial damping factor
		{
			for( i = 0, tmp = LM_REAL_MIN; i < nP; i++ )
			{
				if( JTJ_diag[i] > tmp ) tmp = JTJ_diag[i]; /* find max diagonal element */
				if( JTJ_diag[i] < DBL_EPSILON && op->cd->ldebug && op->cd->paranoid ) tprintf( "WARNING: Parameter %s is NOT impacting model predictions (JTJ diagonal %g)\n", op->pd->var_name[op->pd->var_index[i]], JTJ_diag[i] );
				else if( op->cd->ldebug > 2 && JTJ_diag[i] < 1 ) tprintf( "WARNING: Parameter %s is not very sensitive (JTJ diagonal %g)\n", op->pd->var_name[op->pd->var_index[i]], JTJ_diag[i] );
			}
			if( tmp > 1e4 ) tmp = 1e4;
			lambda = tau * tmp;
			if( op->cd->ldebug ) tprintf( "Computed initial lambda %g (%g, %g)\n", lambda, tau, tmp );
		}
		if( op->cd->lm_num_parallel_lambda > 0 ) // Parallel lambda search
		{
			tprintf( "Parallel lambda search ...\n" );
			double lambda_current, lambda_up = lambda, lambda_down = lambda;
			for( npl = 0; npl < op->cd->lm_num_parallel_lambda; npl++ )
			{
				phi_vector[npl] = HUGE_VAL;
				if( npl == 0 )
					lambda_current = lambda;
				else if( npl % 2 ) // even number
					lambda_current = lambda_down *= op->cd->lm_mu;
				else // odd number
					lambda_current = lambda_up /= op->cd->lm_mu;
				for( i = 0; i < nP; i++ )
					JTJ[i * nP + i] += lambda_current;  // Add lambda to the matrix diagonal
				/* solve augmented equations */
				issolved = linsolver( JTJ, JTe, par_change, nP ); ++nlss;
				if( issolved )
					for( i = 0; i < nP; i++ )
						param_matrix[npl][i] = par_current[i] + par_change[i];
				else
					tprintf( "WARNING: Linear solver failed!\n" );
				for( i = 0; i < nP; i++ )
					JTJ[i * nP + i] -= lambda_current;  // Subtract lambda from the matrix diagonal
				tprintf( "Parallel lambda #%d = %g ...\n", npl + 1, lambda_current );
			}
			tprintf( "Parallel execution of %d lambda searches ...\n", op->cd->lm_num_parallel_lambda );
			time_start = time( NULL );
			func_set( op->cd->lm_num_parallel_lambda, param_matrix, phi_vector, obs_matrix, ( FILE * ) mads_output, adata );
			time_end = time( NULL );
			tprintf( "Done.\n" );
			time_elapsed = time_end - time_start;
			time_lambda += time_elapsed;
			parallel_lambda_count++;
			if( op->cd->tdebug )
			{
				if( time_elapsed > 86400 ) tprintf( "Parallel lambda total PT = %g days (average %g days)\n", ( ( double ) time_elapsed / 86400 ), ( ( double ) time_lambda / parallel_lambda_count / 86400 ) );
				else if( time_elapsed > 3600 ) tprintf( "Parallel lambda total PT = %g hours (average %g hours)\n", ( ( double ) time_elapsed / 3600 ), ( ( double ) time_lambda / parallel_lambda_count / 3600 ) );
				else if( time_elapsed > 60 ) tprintf( "Parallel lambda total PT = %g minutes (average %g minutes)\n", ( ( double ) time_elapsed / 60 ), ( ( double ) time_lambda / parallel_lambda_count / 60 ) );
				else tprintf( "Parallel lambda total PT = %ld seconds (average %g seconds)\n", time_elapsed, ( ( double ) time_lambda / parallel_lambda_count ) );
			}
			for( npl = 0; npl < op->cd->lm_num_parallel_lambda; npl++ )
				tprintf( "Parallel lambda #%d => OF %g ...\n", npl + 1, phi_vector[npl] );
			int npl_min = 0;
			double phi_alpha_min = HUGE_VAL;
			for( npl = 0; npl < op->cd->lm_num_parallel_lambda; npl++ )
			{
				phi_lam_vector[phi_lam_count++] = phi_vector[npl];
				if( phi_vector[npl] < phi_alpha_min ) { phi_alpha_min = phi_vector[npl]; npl_min = npl; }
			}
			tprintf( "Parallel OF (total): %d\n", phi_lam_count );
			phi_lam_vector[phi_lam_count + npl_min] = phi_lam_vector[phi_lam_count];
			phi_lam_vector[phi_lam_count] = phi_alpha_min;
			for( i = 0 ; i < nP; i++ )
				par_update[i] = param_matrix[npl_min][i];
			for( i = 0, par_change_L2_norm = 0.0; i < nP; i++ )
			{
				par_change[i] = tmp = param_matrix[npl_min][i] - par_current[i];
				par_change_L2_norm += tmp * tmp;
			}
			par_change_L2_norm_sq = sqrt( par_change_L2_norm );
			for( i = 0 ; i < nO; i++ )
				obs_update[i] = obs_matrix[npl_min][i];
			if( npl_min != 0 )
			{
				if( npl_min % 2 ) // even number
					lambda *= pow( op->cd->lm_mu, npl_min / 2 + 1 );
				else // odd number
					lambda *= pow( op->cd->lm_mu, -( npl_min + 1 ) / 2 );
			}
			tprintf( "Best Parallel lambda #%d = %g OF = %g\n", npl_min + 1, lambda, phi_alpha_min );
			phi_update = phi_alpha_min;
		}
		else
		{
			/* determine increment using adaptive damping augment normal equations */
			for( i = 0; i < nP; i++ )
				JTJ[i * nP + i] += lambda;  // Add lambda to the matrix diagonal
			/* solve augmented equations */
			issolved = linsolver( JTJ, JTe, par_change, nP ); ++nlss;
			if( issolved )
			{
				if( op->cd->lm_acc ) // Acceleration
				{
					register LM_REAL beta, *jacT;
					for( i = 0; i < nP; i++ )
						jacTvv[i] = 0.0;
					for( i = 0; i < nP; i++ )
					{
						phDp_plus[i] = par_current[i] + acc_h * par_change[i];
						phDp_minus[i] = par_current[i] - acc_h * par_change[i];
					}
					change = 0;
					if( op->cd->compute_phi ) { op->cd->compute_phi = 0; change = 1; }
					( *func )( phDp_plus, hx1, nP, nO, adata ); nfev++;
					for( i = 0; i < nO; i++ )
						ephdp_plus[i] = obs_target[i] - hx1[i];
					( *func )( phDp_minus, hx2, nP, nO, adata ); nfev++;
					if( change ) op->cd->compute_phi = 1;
					for( i = 0; i < nO; i++ )
						ephdp_minus[i] = obs_target[i] - hx2[i];
					for( i = 0; i < nO; i++ )
						vvddr[i] = ( ephdp_plus[i] - 2.0 * obs_error[i] + ephdp_minus[i] ) / ( acc_h * acc_h );
					for( j = nO; j-- > 0; )
					{
						jacT = jac + j * nP;
						for( i = nP; i-- > 0; )
						{
							beta = jacT[i]; //jac[l*nP+i];
							/* J^T*vvddr */
							jacTvv[i] += beta * vvddr[j];
						}
					}
					issolved1 = linsolver( JTJ, jacTvv, acceleration, nP ); ++nlss;
				}
				/* compute p's new estimate and ||Dp||^2 */
				for( i = 0, par_change_L2_norm = 0.0; i < nP; i++ )
				{
					par_update[i] = par_current[i] + ( tmp = par_change[i] );
					par_change_L2_norm += tmp * tmp;
				}
				if( op->cd->lm_acc && issolved1 )
				{
					for( i = 0, acceleration_L2_norm = 0.0; i < nP; i++ )
					{
						par_update[i] += 0.5 * ( tmp = acceleration[i] );
						acceleration_L2_norm += tmp * tmp;
					}
					avRatio = acceleration_L2_norm / par_change_L2_norm;
					double par_change_L2_norm_new = 0;
					for( i = 0; i < nP; i++ )
					{
						tmp = par_change[i] + 0.5 * acceleration[i];
						par_change_L2_norm_new += tmp * tmp;
					}
					if( op->cd->ldebug > 1 ) tprintf( "LM acceleration performed (with acceleration %g vs without acceleration %g)\n", par_change_L2_norm_new, par_change_L2_norm );
					par_change_L2_norm = par_change_L2_norm_new;
				}
				par_change_L2_norm_sq = sqrt( par_change_L2_norm );
			}
			if( op->cd->ldebug > 6 ) tprintf( "Test for convergence: Magnitude of the relative change in parameter space (%g < %g = %g * %g  to converge)\n", par_change_L2_norm_sq, eps2_sq * par_L2_norm, eps2_sq, par_L2_norm );
			if( par_change_L2_norm_sq <= eps2_sq * par_L2_norm ) /* relative change in p is small, stop */
			{
				if( op->cd->ldebug ) tprintf( "CONVERGED: Relative change in parameter space is small (%g < %g)\n", par_change_L2_norm_sq, eps2_sq * par_L2_norm );
				stop = 2;
				break;
			}
			if( op->cd->ldebug > 6 ) tprintf( "Test for convergence: Solution singularity (%g > %g to converge)\n", par_change_L2_norm_sq, ( par_L2_norm + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) );
			if( par_change_L2_norm_sq >= ( par_L2_norm + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) ) /* almost singular */
			{
				if( op->cd->ldebug ) tprintf( "CONVERGED: almost singular solution (%g > %g)\n", par_change_L2_norm_sq, ( par_L2_norm + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) );
				stop = 4;
				break;
			}
		}
		if( op->cd->lm_num_parallel_lambda == 0 ) // if parallel this is already executed
		{
			( *func )( par_update, obs_update, nP, nO, adata ); ++nfev; /* evaluate function at p + Dp */
		}
#if 0
		if( op->cd->solution_type[0] == TEST ) // this is a test; not needed in general
			for( i = 0; i < nO; i++ )
				obs_update[i] = sqrt( obs_update[i] );
#endif
		/* compute ||e(pDp)||_2 */
		/* ### wrk2=x-wrk, pDp_eL2=||wrk2|| */
#if 0
		phi_update = LEVMAR_L2NRMXMY( obs_error_update, obs_target, obs_update, nO );
#else
		for( i = 0, phi_update = 0.0; i < nO; i++ )
		{
			obs_error_update[i] = tmp = obs_target[i] - obs_update[i];
#if 1
			phi_update += tmp * tmp;
#else
			if( op->cd->solution_type[0] == TEST ) // this is a test; not needed in general
				phi_update += obs_update[i];
			else
				phi_update += tmp * tmp;
#endif
		}
#endif
		if( op->cd->ldebug == 1 ) tprintf( "OF %g lambda %g\n", phi_update, lambda );
		else if( op->cd->ldebug > 1 )
		{
			tprintf( "\nLambda search #%d: OF %g lambda %g \n", phi_lam_count + 1, phi_update, lambda );
			DeTransform( par_update, op, jac_min );
			DeTransform( par_lam_last, op, jac_max );
			max_change = 0; min_change = HUGE_VAL;
			ipar_max = ipar_min = 0;
			for( i = 0; i < nP; i++ )
			{
				j = op->pd->var_index[i];
				if( op->pd->var_log[j] ) p_diff = fabs( pow( 10, jac_min[i] ) - pow( 10, jac_max[i] ) );
				else p_diff = fabs( jac_max[i] - jac_min[i] );
				if( max_change < p_diff ) { max_change = p_diff; ipar_max = i; }
				if( min_change > p_diff ) { min_change = p_diff; ipar_min = i; }
			}
			if( op->cd->ldebug > 2 )
			{
				tprintf( "Current parameter estimates:\n" );
				for( i = 0; i < nP; i++ )
				{
					j = op->pd->var_index[i];
					tprintf( "%-30s:", op->pd->var_name[j] );
					if( op->cd->ldebug > 3 )
					{
						if( op->pd->var_log[j] ) tprintf( " %12g", pow( 10, jac_max[i] ) );
						else tprintf( " %12g", jac_max[i] );
						tprintf( " =>" );
					}
					if( op->pd->var_log[j] ) tprintf( " %12g", pow( 10, jac_min[i] ) );
					else tprintf( " %12g", jac_min[i] );
					if( op->cd->ldebug > 3 )
					{
						tprintf( " change" );
						if( op->pd->var_log[j] ) tprintf( " %12g",  pow( 10, jac_min[i] ) - pow( 10, jac_max[i] ) );
						else tprintf( " %12g", jac_min[i] - jac_max[i] );
					}
					tprintf( "  (JTJ diagonal term %g)", JTJ_diag[i] );
					if( JTJ_diag[i] < DBL_EPSILON ) tprintf( " WARNING: not impacting model predictions" );
					else if( JTJ_diag[i] < 1 ) tprintf( " not very sensitive" );
					tprintf( "\n" );
				}
			}
			if( max_change < DBL_EPSILON )
				tprintf( "Parameters did not change between last two lambda searches.\n" );
			else
			{
				tprintf( "Parameter with maximum absolute estimate change between last two lambda searches: %-30s (%g)\n", op->pd->var_name[op->pd->var_index[ipar_max]], max_change );
				tprintf( "Parameter with minimum absolute estimate change between last two lambda searches: %-30s (%g)\n", op->pd->var_name[op->pd->var_index[ipar_min]], min_change );
			}
			for( i = 0 ; i < nP; i++ ) /* update parameter estimates */
				par_lam_last[i] = par_update[i];
		}
		if( phi_update < phi_best )
		{
			phi_best = phi_update;
			for( i = 0; i < nP; i++ )
				par_best[i] = par_update[i];
			for( i = 0; i < nO; i++ )
				obs_best[i] = obs_update[i];
			if( op->cd->ldebug >= 1 ) tprintf( "New Best OF %g\n", phi_best );
		}
		DeTransform( par_update, op, jac_min );
		for( i = 0; i < op->pd->nOptParam; i++ )
			op->pd->var[op->pd->var_index[i]] = jac_min[i];
		op->phi = phi_update;
		save_results( 0, "", op, op->gd );
		if( op->cd->lm_num_parallel_lambda == 0 )
			phi_lam_vector[phi_lam_count++] = phi_update;
		if( phi_lam_count > op->cd->lm_nlamof )
		{
			tmp = HUGE_VAL;
			j = 0;
			for( i = 0; i < phi_lam_count; i++ )
			{
				if( phi_lam_vector[i] < tmp )
				{
					tmp = phi_lam_vector[i];
					if( i < phi_lam_count - op->cd->lm_nlamof ) j = 1;
					else j = 0;
				}
				// tprintf( " %g", phi1[i] );
			}
			// for( i = 0; i < phi2c; i++ )
			// tprintf( " %g", phi2[i] );
			// tprintf( "a %g %d\n", tmp, phi2c );
			for( i = phi_lam_count - op->cd->lm_nlamof; i < phi_lam_count; i++ )
			{
				if( ( phi_lam_vector[i] - tmp ) / tmp < 1 ) j++;
				// tprintf( "a %g %g %g\n", tmp, phi2[i], ( phi2[i] - tmp ) / tmp );
			}
			if( j >= op->cd->lm_nlamof )
				computejac = 1; // New jacobian: Lambda search OF are very similar
		}
		if( phi_update <= eps3 ) /* error is small */ // BELOW a cutoff value
		{
			if( op->cd->ldebug ) tprintf( "CONVERGED: OF below cutoff value (%g < %g)\n", phi_update, eps3 );
			for( i = 0; i < nP; i++ ) par_current[i] = par_update[i];
			phi_current = phi_update;
			stop = 6;
			break;
		}
		if( op->cd->check_success && op->success )
		{
			if( op->cd->ldebug ) tprintf( "CONVERGED: Predictions are within predefined calibration ranges\n" );
			for( i = 0; i < nP; i++ ) par_current[i] = par_update[i];
			phi_current = phi_update;
			stop = 8;
			break;
		}
		if( !LM_FINITE( phi_update ) )
		{
			/* sum of squares is not finite, most probably due to a user error.
			 * This check makes sure that the loop terminates early in the case
			 * of invalid input. Thanks to Steve Danauskas for suggesting it */
			if( op->cd->ldebug ) tprintf( "CONVERGED: sum of squares is not finite, most probably due to a user error\n" );
			stop = 7;
			break;
		}
		if( op->cd->quit == 1 ) { stop = 99; break; } // Terminate
		if( op->cd->lm_indir ) // original
		{
			tmp = phi_jac_last / phi_update; // original code
			// tprintf( "Original tmp = %g\n", tmp );
			if( tmp > op->cd->lm_ofdecline ) phi_decline = 1; /* recompute jacobian because OF decreased */
		}
		else
		{
			tmp = phi_current / phi_update;
			// tprintf( "Leif tmp = %g\n", tmp );
			if( tmp > op->cd->lm_ofdecline && avRatio < lm_ratio ) phi_decline = 1; /* recompute jacobian because OF decreased */
		}
		if( op->cd->ldebug > 10 ) tprintf( "OF (%g ? %g) ...\n", phi_current, phi_update );
		phi_change = phi_current - phi_update; // Difference between current and previous OF
		if( par_update_accepted || phi_change > 0.0 ) /* update jacobian because OF decreases */
		{
			if( op->cd->ldebug > 3 ) tprintf( "Update jacobian because OF decreases (%g > %g) ...\n", phi_current, phi_update );
			for( i = 0; i < nO; i++ )
			{
				for( l = 0, tmp = 0.0; l < nP; ++l )
					tmp += jac[i * nP + l] * par_change[l]; /* (J * Dp)[i] */
				tmp = ( obs_update[i] - obs_current[i] - tmp ) / par_change_L2_norm; /* (f(p+dp)[i] - f(p)[i] - (J * Dp)[i])/(dp^T*dp) */
				for( j = 0; j < nP; ++j )
					jac[i * nP + j] += tmp * par_change[j];
			}
			newjac = 1;
			jac_update_count++;
		}
		if( op->cd->lm_indir )
		{
			for( i = 0, dL = 0.0; i < nP; i++ )
				dL += par_change[i] * ( lambda * par_change[i] + JTe[i] );
			if( op->cd->ldebug > 3 ) tprintf( "dL %g\n", dL );
		}
		else dL = ( double ) 1.0;
		if( dL > 0.0 && phi_change > 0.0 && avRatio < lm_ratio ) /* reduction in error, increment is accepted */
		{
			if( op->cd->lm_indir )
			{
				tmp = ( LM_CNST( 2.0 ) * phi_change / dL - LM_CNST( 1.0 ) );
				tmp = LM_CNST( 1.0 ) - tmp * tmp * tmp;
				tmp = ( ( tmp >= LM_CNST( ONE_THIRD ) ) ? tmp : LM_CNST( ONE_THIRD ) );
			}
			else tmp = op->cd->lm_mu;
			lambda = lambda * tmp; // change lambda
			/*
								if( mu > 1e3 )
								{
									if( mu_constrained > 500000 )
									{
										if( op->cd->ldebug ) tprintf( "lambda has been already constrained; new iteration" );
										mu_big = 1;
										mu_constrained = 0;
									}
									else
									{
										// mu = 1e3;
										mu_constrained++;
										if( op->cd->ldebug ) tprintf( "lambda is constrained to be less than %g", mu );
									}
								}
								else mu_constrained = 0;
			 */
			if( op->cd->ldebug > 3 ) tprintf( "Reduction in OF (%g > %g) ... Lambda multiplied by factor %g\n", phi_current, phi_update, tmp );
			nu = op->cd->lm_nu;
			for( i = 0 ; i < nP; i++ ) /* update parameter estimates */
				par_current[i] = par_update[i];
			for( i = 0; i < nO; i++ ) /* update e, hx and ||e||_2 */
			{
				obs_error[i] = obs_error_update[i]; //x[i]-wrk[i];
				obs_current[i] = obs_update[i];
			}
			phi_current = phi_update; // Update OF
			par_update_accepted = 1;
			continue; // Solve for a new lambda
		}
		/* if this point is reached, either the linear system could not be solved or
		 * the error did not reduce; in any case, the increment must be rejected */
		lambda *= nu; // increase lambda
		/*
				if( mu > 1e3 )
				{
					if( mu_constrained > 5000000 )
					{
						if( op->cd->ldebug ) tprintf( "lambda has been already constrained; new iteration" );
						mu_big = 1;
						mu_constrained = 0;
					}
					else
					{
						// mu = 1e3;
						mu_constrained++;
						if( op->cd->ldebug ) tprintf( "lambda is constrained to be less than %g", mu );
					}
				}
				else mu_constrained = 0;
		 */
		if( op->cd->ldebug > 3 ) tprintf( "Increase in OF (%g < %g) ... Lambda increased by factor %d\n", phi_current, phi_update, nu );
		if( op->cd->lm_indir )
		{
			nu2 = nu << 1; // 2 * nu;
			// if( op->cd->ldebug > 3 ) tprintf( "NU %d %d\n", nu2, nu );
			if( nu2 <= nu )
			{
				if( op->cd->ldebug ) tprintf( "CONVERGED: lambda multiplication factor has wrapped around (overflown)\n" );
				stop = 5;
				break;
			}
			nu = nu2;
		}
		else { computejac = 0; }
		for( i = 0; i < nP; i++ ) /* restore diagonal J^T J entries */
			JTJ[i * nP + i] = JTJ_diag[i];
		if( op->cd->lm_num_parallel_lambda > 0 && op->cd->lm_num_parallel_lambda + phi_lam_count >= max_num_of_lambda_searches )
		{
			for( i = 0 ; i < nP; i++ ) /* update parameter estimates */
				par_current[i] = par_update[i];
			for( i = 0; i < nO; i++ ) /* update e, hx and ||e||_2 */
			{
				obs_error[i] = obs_error_update[i]; //x[i]-wrk[i];
				obs_current[i] = obs_update[i];
			}
			phi_current = phi_update; // Update OF
		}
	} // END OF OPTIMIZATION LOOP
	if( fabs( phi_best - phi_current ) > COMPARE_EPSILON ) tprintf( "Best OF %g\n", phi_best );
	op->phi = phi_current = phi_best;
	for( i = 0; i < nP; i++ )
		par_current[i] = par_best[i];
	for( i = 0; i < nO; i++ )
		op->od->obs_current[i] = op->od->res[i] = obs_best[i];
	if( op->cd->ldebug > 3 )
	{
		DeTransform( par_best, op, jac_min );
		DeTransform( par_init, op, jac_max );
		max_change = 0; min_change = HUGE_VAL;
		ipar_max = ipar_min = 0;
		for( i = 0 ; i < nP; i++ )
		{
			j = op->pd->var_index[i];
			if( op->pd->var_log[j] ) p_diff = fabs( pow( 10, jac_min[i] ) - pow( 10, jac_max[i] ) );
			else p_diff = fabs( jac_max[i] - jac_min[i] );
			if( max_change < p_diff ) { max_change = p_diff; ipar_max = i; }
			if( min_change > p_diff ) { min_change = p_diff; ipar_min = i; }
		}
		if( max_change < DBL_EPSILON )
			tprintf( "Parameters did not change.\n" );
		else
		{
			tprintf( "Parameter with maximum absolute estimate change: %-30s (%g)\n", op->pd->var_name[op->pd->var_index[ipar_max]], max_change );
			tprintf( "Parameter with minimum absolute estimate change: %-30s (%g)\n", op->pd->var_name[op->pd->var_index[ipar_min]], min_change );
		}
	}
	for( i = 0; i < nP; i++ ) /* restore diagonal J^T J entries */
		JTJ[i * nP + i] = JTJ_diag[i];
	if( info || op->cd->ldebug )
		for( i = 0, tmp = LM_REAL_MIN; i < nP; i++ )
			if( tmp < JTJ[i * nP + i] ) tmp = JTJ[i * nP + i];
	if( info )
	{
		info[0] = phi_init;
		info[1] = phi_current;
		info[2] = JTe_inf_norm;
		info[3] = par_change_L2_norm;
		info[4] = lambda / tmp;
		info[5] = ( LM_REAL )loop_count;
		info[6] = ( LM_REAL )stop;
		info[7] = ( LM_REAL )nfev;
		info[8] = ( LM_REAL )njac;
		info[9] = ( LM_REAL )nlss;
	}
	if( op->cd->ldebug )
	{
		tprintf( "Final OF %g\n", phi_current );
		tprintf( "LM optimization is completed. Reason: " );
		switch( stop )
		{
			case 1: tprintf( "small gradient J^T e (%g)\n", JTe_inf_norm ); break;
			case 2: tprintf( "small Dp (%g)\n", par_change_L2_norm ); break;
			case 31: tprintf( "maximum number of jacobian iterations is reached (lmiter=%d)\n", maxjac ); break;
			case 32: tprintf( "maximum number of functional evaluations is exceeded (eval=%d; %d > %d)\n", op->cd->maxeval, nfev, maxnfev ); break;
			case 4: tprintf( "singular matrix. Restart from current p with increased mu (current mu=%g; mu/max[J^T J]_ii=%g)\n", lambda, lambda / tmp ); break;
			case 5: tprintf( "no further error reduction is possible. Restart with increased mu (current mu=%g; mu/max[J^T J]_ii=%g)\n", lambda, lambda / tmp ); break;
			case 6: tprintf( "small ||e||_2; OF below cutoff value (%g < %g)\n", phi_current, eps3 ); break;
			case 7: tprintf( "invalid (i.e. NaN or Inf) values returned by the solver (func). This is a user error\n" ); break;
			case 8: tprintf( "model predictions are within predefined calibration ranges\n" ); break;
			case 9: tprintf( "small OF changes\n" ); break;
			case 10: tprintf( "all jacobian matrix elements are equal to zero\n" ); break;
			case 99: tprintf( "MADS forced to quit\n" ); break;
			default: tprintf( "UNKNOWN flag: %d\n", stop ); break;
		}
		if( op->cd->ldebug > 14 )
		{
			tprintf( "Levenberg-Marquardt Optimization completed after %g iteration (reason %g) (returned value %d)\n", info[5], info[6], ( stop != 4 && stop != 7 ) ?  loop_count : LM_ERROR );
			tprintf( "initial phi %g final phi %g ||J^T e||_inf %g ||Dp||_2 %g mu/max[J^T J]_ii %g\n", info[0], info[1], info[2], info[3], info[4] );
			tprintf( "function evaluation %g jacobian evaluations %g linear systems solved %g\n", info[7], info[8], info[9] );
		}
	}
	else if( op->cd->lmstandalone ) { if( op->cd->lmstandalone > 1 ) tprintf( "%g\n", phi_current ); else tprintf( "%g ", phi_current ); }
	/* covariance matrix */
	if( covar )
	{
		LEVMAR_COVAR( JTJ, covar, phi_current, nP, nO );
	}
	if( freework ) free( work );
	free( hx1 ); free( hx2 ); free( phDp_plus ); free( phDp_minus ); free( ephdp_plus ); free( ephdp_minus );
	free( vvddr ); free( jacTvv ); free( acceleration ); free( jac_min ); free( jac_max ); free( jac_zero ); free( jac_zero_obs );
	free( par_lam_last ); free( par_jac_last ); free( par_best ); free( par_init );
	if( op->cd->lm_eigen ) gsl_matrix_free( gsl_jacobian );
	if( op->cd->lm_num_parallel_lambda > 0 )
	{
		free( phi_vector );
		free_matrix( ( void ** ) param_matrix, op->cd->lm_num_parallel_lambda );
		free_matrix( ( void ** ) obs_matrix, op->cd->lm_num_parallel_lambda );
	}
#ifdef LINSOLVERS_RETAIN_MEMORY
	if( linsolver )( *linsolver )( NULL, NULL, NULL, 0 );
#endif
	return ( stop != 4 && stop != 7 ) ?  loop_count : LM_ERROR;
}

int LEVMAR_DER(
	void ( *func )( LM_REAL *p, LM_REAL *hx, int m, int n, void *adata ), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
	void ( *jacf )( LM_REAL *p, LM_REAL *j, int m, int n, void *adata ), /* function to evaluate the Jacobian \part x / \part p */
	LM_REAL *p,         /* I/O: initial parameter estimates. On output has the estimated solution */
	LM_REAL *x,         /* I: measurement vector. NULL implies a zero vector */
	int m,              /* I: parameter vector dimension (i.e. #unknowns) */
	int n,              /* I: measurement vector dimension */
	int itmax,          /* I: maximum number of iterations */
	LM_REAL opts[4],    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
		 * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used
		 */
	LM_REAL info[LM_INFO_SZ],
	/* O: information regarding the minimization. Set to NULL if don't care
	 * info[0]= ||e||_2 at initial p.
	 * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
	 * info[5]= # iterations,
	 * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
	 *                                 2 - stopped by small Dp
	 *                                 3 - stopped by itmax
	 *                                 4 - singular matrix. Restart from current p with increased mu
	 *                                 5 - no further error reduction is possible. Restart with increased mu
	 *                                 6 - stopped by small ||e||_2
	 *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
	 * info[7]= # function evaluations
	 * info[8]= # Jacobian evaluations
	 * info[9]= # linear systems solved, i.e. # attempts for reducing error
	 */
	LM_REAL *work,     /* working memory at least LM_DER_WORKSZ() reals large, allocated if NULL */
	LM_REAL *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */
	void *adata )       /* pointer to possibly additional data, passed uninterpreted to func & jacf.
		 * Set to NULL if not needed
		 */
{
	return 1;
}


/* Secant version of the LEVMAR_DER() function above: the Jacobian is approximated with
 * the aid of finite differences (forward or central, see the comment for the opts argument)
 */
int LEVMAR_DIF(
	void ( *func )( LM_REAL *p, LM_REAL *hx, int m, int n, void *adata ), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
	LM_REAL *p,         /* I/O: initial parameter estimates. On output has the estimated solution */
	LM_REAL *x,         /* I: measurement vector. NULL implies a zero vector */
	int m,              /* I: parameter vector dimension (i.e. #unknowns) */
	int n,              /* I: measurement vector dimension */
	int itmax,          /* I: maximum number of iterations */
	LM_REAL opts[5],    /* I: opts[0-4] = minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \delta]. Respectively the
		 * scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2 and
		 * the step used in difference approximation to the Jacobian. Set to NULL for defaults to be used.
		 * If \delta<0, the Jacobian is approximated with central differences which are more accurate
		 * (but slower!) compared to the forward differences employed by default.
		 */
	LM_REAL info[LM_INFO_SZ],
	/* O: information regarding the minimization. Set to NULL if don't care
	 * info[0]= ||e||_2 at initial p.
	 * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
	 * info[5]= # iterations,
	 * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
	 *                                 2 - stopped by small Dp
	 *                                 3 - stopped by itmax
	 *                                 4 - singular matrix. Restart from current p with increased mu
	 *                                 5 - no further error reduction is possible. Restart with increased mu
	 *                                 6 - stopped by small ||e||_2
	 *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
	 * info[7]= # function evaluations
	 * info[8]= # Jacobian evaluations
	 * info[9]= # linear systems solved, i.e. # attempts for reducing error
	 */
	LM_REAL *work,     /* working memory at least LM_DIF_WORKSZ() reals large, allocated if NULL */
	LM_REAL *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */
	void *adata )       /* pointer to possibly additional data, passed uninterpreted to func.
		 * Set to NULL if not needed
		 */
{
	return 1;
}

/* undefine everything. THIS MUST REMAIN AT THE END OF THE FILE */
#undef LEVMAR_DER
#undef LEVMAR_DIF
#undef LEVMAR_FDIF_FORW_JAC_APPROX
#undef LEVMAR_FDIF_CENT_JAC_APPROX
#undef LEVMAR_COVAR
#undef LEVMAR_TRANS_MAT_MAT_MULT
#undef LEVMAR_L2NRMXMY
#undef AX_EQ_B_LU
#undef AX_EQ_B_CHOL
#undef AX_EQ_B_QR
#undef AX_EQ_B_QRLS
#undef AX_EQ_B_SVD
#undef AX_EQ_B_BK
