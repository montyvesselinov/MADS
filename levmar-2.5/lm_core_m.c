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

/* precision-specific definitions */
#include <string.h>

#define LEVMAR_DER LM_ADD_PREFIX(levmar_der)
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

int LEVMAR_DER(
	void ( *func )( LM_REAL *p, LM_REAL *hx, int m, int n, void *adata ), /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
	void ( *jacf )( LM_REAL *p, LM_REAL *f, LM_REAL *j, int m, int n, void *adata ), /* function to evaluate the Jacobian \part x / \part p */
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
	register int i, j, k, l;
	int worksz, freework = 0, issolved, issolved1, success, odebug, change, computejac, changejac, kmax, maxnfev;
	struct opt_data *op = ( struct opt_data * ) adata;
	/* temp work arrays */
	LM_REAL *e,          /* nx1 */
			*hx,         /* \hat{x}_i, nx1 */
			*hx1,        /* used in acceleration comp (nx1) */
			*hx2,        /* used in acceleration comp (nx1) */
			*jac_min,
			*jac_max,
			*p_old, *p_old2,      /* old parameter set */
			*jacTe,      /* J^T e_i mx1 */
			*jac,        /* nxm */
			*jacTjac,    /* mxm */
			*Dp,         /* mx1 (=v) */
			*ephdp_plus,     /* residual used in acceleration computation (nx1) */
			*ephdp_minus,     /* residual used in acceleration computation (nx1) */
			*vvddr,      /* used to compute acceleration (nx1) */
			*jacTvv,     /* jacT*vvddr, mx1 */
			*a,          /* acceleration (mx1) */
			*diag_jacTjac,   /* diagonal of J^T J, mx1 */
			*pDp,        /* p + Dp, mx1 */
			*phDp_plus,       /* p + hDp, mx1 */
			*phDp_minus,      /* p - hDp, mx1 */
			*wrk,        /* nx1 */
			*wrk2;       /* nx1, used only for holding a temporary e vector and when differentiating with central differences */
	int using_ffdif = 1;
	double acc_h = op->cd->lm_h;
	double lm_ratio = op->cd->lm_ratio; // alpha
	register LM_REAL mu,  /* damping constant */
			 tmp, /* mainly used in matrix & vector multiplications */
			 avRatio = 0.0; /* acceleration/velocity */
	LM_REAL p_eL2, p_eL2_old, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
	LM_REAL p_L2, Dp_L2 = LM_REAL_MAX, a_L2, dF, Dpa_L2, dL;
	LM_REAL tau, eps1, eps2, eps2_sq, eps3, delta;
	LM_REAL init_p_eL2;
	int nu, nu2, stop = 0, nfev, njap = 0, nlss = 0, K = ( m >= 10 ) ? m : 10, updjac, updp = 1, newjac;
	LM_REAL phi1[itmax], phi2[K];
	int phi1c, phi2c;
	int mu_big = 0, phi_decline = 0;
	const int nm = n * m;
	int ( *linsolver )( LM_REAL * A, LM_REAL * B, LM_REAL * x, int m ) = NULL;
	mu = jacTe_inf = p_L2 = 0.0; /* -Wall */
	updjac = newjac = 0; /* -Wall */
	if( n < m )
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
		worksz = LM_DIF_WORKSZ( m, n ); //4*n+4*m + n*m + m*m;
		work = ( LM_REAL * )malloc( worksz * sizeof( LM_REAL ) ); /* allocate a big chunk in one step */
		if( !work )
		{
			fprintf( stderr, LCAT( LEVMAR_DIF, "(): memory allocation request failed\n" ) );
			return LM_ERROR;
		}
		freework = 1;
	}
	/* set up work arrays */
	e = work;
	hx = e + n;
	jacTe = hx + n;
	jac = jacTe + m;
	jacTjac = jac + nm;
	Dp = jacTjac + m * m;
	diag_jacTjac = Dp + m;
	pDp = diag_jacTjac + m;
	wrk = pDp + m;
	wrk2 = wrk + n;
	hx1 = ( LM_REAL * )malloc( n * sizeof( LM_REAL ) );
	hx2 = ( LM_REAL * )malloc( n * sizeof( LM_REAL ) );
	phDp_plus = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	phDp_minus = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	ephdp_plus = ( LM_REAL * )malloc( n * sizeof( LM_REAL ) );
	ephdp_minus = ( LM_REAL * )malloc( n * sizeof( LM_REAL ) );
	vvddr = ( LM_REAL * )malloc( n * sizeof( LM_REAL ) );
	jacTvv = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	a = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	jac_min = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	jac_max = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	p_old = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	p_old2 = ( LM_REAL * )malloc( m * sizeof( LM_REAL ) );
	/* compute e=x - f(p) and its L2 norm */
	maxnfev = op->cd->maxeval - op->cd->neval;
	for( i = 0; i < m; ++i )
		p_old[i] = p_old2[i] = pDp[i] = p[i];
	( *func )( p, hx, m, n, adata ); nfev = 1;
	if( op->cd->check_success && op->success )
	{
		if( op->cd->ldebug ) printf( "SUCCESS: Model predictions are within predefined calibration ranges\n" );
		stop = 8;
	}
#ifdef HAVE_LAPACK
	/* 6 alternatives are available: LU, Cholesky, 2 variants of QR decomposition, SVD and LDLt.
	 * Cholesky is the fastest but might be inaccurate; QR is slower but more accurate;
	 * SVD is the slowest but most accurate; LU offers a tradeoff between accuracy and speed
	 */
	//linsolver = AX_EQ_B_BK; if( op->cd->ldebug ) printf( "BK decomposition\n" );
	//linsolver = AX_EQ_B_LU; if( op->cd->ldebug ) printf( "LU decomposition\n" );
	//linsolver = AX_EQ_B_CHOL; if( op->cd->ldebug ) printf( "Cholesky decomposition\n" );
	//linsolver = AX_EQ_B_QR; if( op->cd->ldebug ) printf( "QR decomposition\n" );
	//linsolver = (int (*)(LM_REAL *A, LM_REAL *B, LM_REAL *x, int m))AX_EQ_B_QRLS; if( op->cd->ldebug ) printf( "QRLS decomposition\n" );
	linsolver = AX_EQ_B_SVD; if( op->cd->ldebug ) printf( "LM using SVD decomposition\n" );
#else
	/* use the LU included with levmar */
	linsolver = AX_EQ_B_LU; if( op->cd->ldebug ) printf( "LM using LU decomposition\n" );
#endif /* HAVE_LAPACK */
	/* ### e=x-hx, p_eL2=||e|| */
#if 0
	p_eL2 = LEVMAR_L2NRMXMY( e, x, hx, n );
#else
	for( i = 0, p_eL2 = 0.0; i < n; ++i )
	{
		e[i] = tmp = x[i] - hx[i];
		p_eL2 += tmp * tmp;
	}
#endif
	if( op->cd->ldebug )
	{
		if( op->cd->lm_acc ) printf( "LM with acceleration\n" );
		else printf( "LM without LM acceleration\n" );
		if( op->cd->lm_indir ) printf( "LM with indirect computation of lambda changes\n" );
		else printf( "LM with direct computation of lambda changes\n" );
	}
	init_p_eL2 = p_eL2_old = p_eL2;
	if( !LM_FINITE( p_eL2 ) ) stop = 7;
	nu = 20; /* force computation of J */
	if( op->cd->check_success ) success = 1; else success = 0;
	if( op->cd->odebug ) odebug = 1; else odebug = 0;
	computejac = 0;
	kmax = itmax * 100;
	phi1c = phi2c = 0;
	if( op->cd->ldebug > 2 )
	{
		DeTransform( p, op, jac_min );
		printf( "Initial parameter estimates:\n" );
		for( i = 0; i < m; i++ )
		{
			j = op->pd->var_index[i];
			printf( "%-20s:", op->pd->var_id[j] );
			if( op->pd->var_log[j] ) printf( " %g", pow( 10, jac_min[i] ) );
			else printf( " %g", jac_min[i] );
			printf( "\n" );
		}
	}
	if( op->cd->ldebug ) printf( "Initial evaluation: OF %g\n", p_eL2 );
	else if( op->cd->standalone ) { printf( "OF %g -> ", p_eL2 ); fflush( stdout ); }
	for( k = 0; k < kmax && !stop; ++k )
	{
		/* Note that p and e have been updated at a previous iteration */
		if( p_eL2 <= eps3 ) /* error is small */ // BELOW a cutoff value
		{
			if( op->cd->ldebug ) printf( "CONVERGED: OF below cutoff value (%g < %g)\n", p_eL2, eps3 );
			stop = 6;
			break;
		}
		/* Compute the Jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
		 * The symmetry of J^T J is again exploited for speed
		 */
		if( ( updp && nu > 16 ) || updjac >= K || mu_big || phi_decline || computejac ) /* compute difference approximation to J */
		{
			if( op->cd->ldebug && k != 0 )
			{
				if( k != 0 ) printf( "New Jacobian requested because: " );
				if( updp && nu > 16 && k != 0 ) printf( "Lambda multiplication factor too large (nu = %d > 16); ", nu );
				if( updjac >= K ) printf( "Maximum number of lambda iteration is reached (%d); ", K );
				if( mu_big ) printf( "Lambda is constrained (%g); ", 1e3 );
				if( phi_decline ) printf( "OF estimate declined substantially (%g << %g)", p_eL2, p_eL2_old );
				if( computejac ) printf( "Linear solve OF estimates do not change substantially" );
				printf( "\n\n" );
			}
			if( njap >= itmax ) { stop = 31; continue; }
			if( nfev >= maxnfev ) { stop = 32; continue; }
			computejac = 0;
			changejac = 1;
			mu_big = 0;
			phi_decline = 0;
			phi2c = 0;
			p_eL2_old = p_eL2;
			phi1[phi1c++] = p_eL2;
			if( phi1c > op->cd->lm_njacof )
			{
				tmp = HUGE_VAL;
				for( i = 0; i < phi1c; i++ )
				{
					if( phi1[i] < tmp ) tmp = phi1[i];
					// printf( " %g", phi1[i] );
				}
				// printf( "d %g %d\n", tmp, phi1c );
				j = 0;
				for( i = 0; i < phi1c; i++ )
				{
					// printf( "d %g %g %g\n", tmp, phi1[i], ( phi1[i] - tmp ) / tmp );
					if( ( phi1[i] - tmp ) / tmp < 1 ) j++;
				}
				// printf( "d %g %d\n", tmp, j );
				if( j >= op->cd->lm_njacof )
				{
					if( op->cd->ldebug ) printf( "\nCONVERGED: %d Jacobian OF estimates are very close to the best current OF (%g)\n", j, tmp );
					stop = 9;
					break;
				}
			}
			if( success ) op->cd->check_success = 0;
			if( odebug ) op->cd->odebug = 0;
			if( using_ffdif ) /* use forward differences */
			{
				jacf( p, hx, jac, m, n, adata );
				++njap; nfev += m;
			}
			else  /* use central differences */
			{
				if( op->cd->ldebug ) printf( "Central Differences\n" );
				LEVMAR_FDIF_CENT_JAC_APPROX( func, p, wrk, wrk2, delta, jac, m, n, adata );
				++njap; nfev += 2 * m;
			}
			if( success ) op->cd->check_success = success;
			if( odebug ) op->cd->odebug = odebug;
			nu = op->cd->lm_nu; updjac = 0; updp = 0; newjac = 1;
			if( op->cd->ldebug )
			{
				if( k == 0 ) printf( "Jacobians %d Linear solves %d Evaluations %d OF %g lambda %g\n", njap, nlss, nfev, p_eL2, tau );
				else printf( "Jacobians %d Linear solves %d Evaluations %d OF %g lambda %g\n", njap, nlss, nfev, p_eL2, mu );
			}
			if( op->cd->ldebug >= 4 )
			{
				printf( "\nJacobian matrix:\n" );
				char *s;
				printf( "Parameter:" ); for( i = 0; i < op->pd->nOptParam; i++ )
				{
					s = op->pd->var_id[op->pd->var_index[i]];
					if( strlen( s ) < 7 ) printf( " %s", s );
					else printf( " p%d", i );
				}
				for( l = j = 0; j < op->od->nObs; j++ )
				{
					if( op->cd->ldebug >= 5 || op->od->nObs < 30 || ( j < 10 || j > op->od->nObs - 10 ) )
					{
						int print = 0;
						for( i = 0; i < op->pd->nOptParam; i++, l++ )
							if( fabs( jac[l] ) > DBL_EPSILON ) print = 1;
						if( print )
						{
							l -= op->pd->nOptParam; // reset jacobian counter
							printf( "\n%-12s: ", op->od->obs_id[j] );
							for( i = 0; i < op->pd->nOptParam; i++, l++ )
								printf( " %g", jac[l] );
						}
					}
					else l += op->pd->nOptParam;
					if( op->cd->ldebug == 4 && op->od->nObs >= 30 && j == 11 ) printf( "\n..." );
				}
				printf( "\n" );
			}
			if( op->cd->ldebug > 1 )
			{
				int skipped = 0;
				int first = 1;
				int allzero = 1;
				for( l = j = 0; j < op->od->nObs; j++ )
				{
					if( op->od->obs_weight[j] > DBL_EPSILON )
					{
						for( i = 0; i < op->pd->nOptParam; i++, l++ )
							if( fabs( jac[l] ) > DBL_EPSILON ) allzero = 0;
						if( allzero )
						{
							skipped++;
							if( op->cd->ldebug > 2 )
							{
								if( first ) { first = 0; printf( "WARNING: Observation(s) with zero sensitivities:" ); }
								printf( " %s", op->od->obs_id[j] );
							}
						}
					}
					else l += op->pd->nOptParam;
				}
				if( skipped && op->cd->ldebug > 2 ) printf( "\n" );
				if( skipped )
				{
					if( skipped > 1 ) printf( "\nWARNING: %d observations have zero sensitivities!\n", skipped );
					else printf( "\nWARNING: %d observation has zero sensitivities!\n", skipped );
				}
			}
			if( op->cd->ldebug > 1 )
			{
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					jac_max[i] = 0;
					jac_min[i] = HUGE_VAL;
				}
				int imax, imin, omax, omin;
				double max = 0, min = HUGE_VAL;
				for( l = j = 0; j < op->od->nObs; j++ )
				{
					if( !( op->od->obs_weight > 0 ) ) continue;
					for( i = 0; i < op->pd->nOptParam; i++, l++ )
					{
						double fj = fabs( jac[l] );
						if( jac_max[i] < fj ) jac_max[i] = fj;
						if( jac_min[i] > fj ) jac_min[i] = fj;
						if( fj > max ) { max = fj; imax = i; omax = j; }
						if( fj < min && fj > DBL_EPSILON ) { min = fj; imin = i; omin = j; }
					}
				}
				printf( "Highest absolute sensitivity - parameter: %s (p%d) - observation: %s (o%d): %g\n", op->pd->var_id[op->pd->var_index[imax]], imax + 1, op->od->obs_id[omax], omax + 1, max );
				printf( "Lowest  absolute sensitivity - parameter: %s (p%d) - observation: %s (o%d): %g\n", op->pd->var_id[op->pd->var_index[imin]], imin + 1, op->od->obs_id[omin], omin + 1, min );
				if( op->cd->ldebug > 2 )
				{
					char *s;
					printf( "\nParameter:" ); for( i = 0; i < op->pd->nOptParam; i++ )
					{
						s = op->pd->var_id[op->pd->var_index[i]];
						if( strlen( s ) < 7 ) printf( " %s", s );
						else printf( " p%d", i );
					}
					printf( "\n" );
					printf( "Min absolute observation sensitivity:" ); for( i = 0; i < op->pd->nOptParam; i++ ) printf( " %g", jac_min[i] ); printf( "\n" );
					printf( "Max absolute observation sensitivity:" ); for( i = 0; i < op->pd->nOptParam; i++ ) printf( " %g", jac_max[i] ); printf( "\n" );
				}
			}
			if( op->cd->ldebug > 3 )
			{
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					jac_max[i] = -HUGE_VAL;
					jac_min[i] = HUGE_VAL;
				}
				for( l = j = 0; j < op->od->nObs; j++ )
				{
					for( i = 0; i < op->pd->nOptParam; i++, l++ )
					{
						if( jac_max[i] < jac[l] ) jac_max[i] = jac[l];
						if( jac_min[i] > jac[l] ) jac_min[i] = jac[l];
					}
				}
				printf( "Min observation sensitivity:" ); for( i = 0; i < op->pd->nOptParam; i++ ) printf( " %g", jac_min[i] ); printf( "\n" );
				printf( "Max observation sensitivity:" ); for( i = 0; i < op->pd->nOptParam; i++ ) printf( " %g", jac_max[i] ); printf( "\n\n" );
				int ok = 0;
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					if( fabs( jac_max[i] ) > DBL_EPSILON || fabs( jac_min[i] ) > DBL_EPSILON ) ok++;
					else { if( op->cd->ldebug ) printf( "WARNING: Model parameter \'%s\' is not impacting model predictions!\n", op->pd->var_id[op->pd->var_index[i]] ); }
				}
				if( ok == 0 && !( op->cd->ldebug && op->cd->standalone ) )
					printf( "WARNING: None of the model parameters is impacting model predictions!\n" );
			}
			else
			{
				int ok = 0;
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					if( jac_max[i] > DBL_EPSILON ) ok++;
					else { if( op->cd->ldebug ) printf( "WARNING: Model parameter \'%s\' is not impacting model predictions\n", op->pd->var_id[op->pd->var_index[i]] ); }
				}
				if( ok == 0 && !( op->cd->ldebug && op->cd->standalone ) )
					printf( "WARNING: None of the model parameters is impacting model predictions!\n" );
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
				/* looping downwards saves a few computations */
				for( i = m * m; i-- > 0; )
					jacTjac[i] = 0.0;
				for( i = m; i-- > 0; )
					jacTe[i] = 0.0;
				for( l = n; l-- > 0; )
				{
					jaclm = jac + l * m;
					for( i = m; i-- > 0; )
					{
						im = i * m;
						alpha = jaclm[i]; //jac[l*m+i];
						for( j = i + 1; j-- > 0; ) /* j<=i computes lower triangular part only */
							jacTjac[im + j] += jaclm[j] * alpha; //jac[l*m+j]
						/* J^T e */
						jacTe[i] += alpha * e[l];
					}
				}
				for( i = m; i-- > 0; ) /* copy to upper part */
					for( j = i + 1; j < m; ++j )
						jacTjac[i * m + j] = jacTjac[j * m + i];
			}
			else  // this is a large problem
			{
				/* Cache efficient computation of J^T J based on blocking
				 */
				LEVMAR_TRANS_MAT_MAT_MULT( jac, jacTjac, n, m );
				/* cache efficient computation of J^T e */
				for( i = 0; i < m; ++i )
					jacTe[i] = 0.0;
				for( i = 0; i < n; ++i )
				{
					register LM_REAL *jacrow;
					for( l = 0, jacrow = jac + i * m, tmp = e[i]; l < m; ++l )
						jacTe[l] += jacrow[l] * tmp;
				}
			}
			/* Compute ||J^T e||_inf and ||p||^2 */
			for( i = 0, p_L2 = jacTe_inf = 0.0; i < m; ++i )
			{
				if( jacTe_inf < ( tmp = FABS( jacTe[i] ) ) ) jacTe_inf = tmp;
				diag_jacTjac[i] = jacTjac[i * m + i]; /* save diagonal entries so that augmentation can be later canceled */
				p_L2 += p[i] * p[i];
			}
			if( op->cd->ldebug && ( changejac && ( k > 0 ) ) )
			{
				DeTransform( p, op, jac_min );
				DeTransform( p_old2, op, jac_max );
				int ipar_max, ipar_min;
				double max_change = 0, min_change = HUGE_VAL, p_diff;
				for( i = 0 ; i < m; ++i )
				{
					p_diff = fabs( jac_max[i] - jac_min[i] );
					if( max_change < p_diff ) { max_change = p_diff; ipar_max = i; }
					if( min_change > p_diff ) { min_change = p_diff; ipar_min = i; }
				}
				printf( "Parameter with maximum estimate change between jacobian iterations: %s (%g)\n", op->pd->var_id[op->pd->var_index[ipar_max]], jac_max[ipar_max] - jac_min[ipar_max] );
				printf( "Parameter with minimum estimate change between jacobian iterations: %s (%g)\n", op->pd->var_id[op->pd->var_index[ipar_min]], jac_max[ipar_min] - jac_min[ipar_min] );
				for( i = 0 ; i < m; ++i ) /* update p's estimate */
					p_old2[i] = p[i];
			}
			/* check for convergence */
			if( op->cd->ldebug > 6 ) printf( "\nTest for convergence: Magnitude of the largest J^T e component (||J^T e||_inf = %g < %g to converge)\n", jacTe_inf, eps1 );
			if( ( jacTe_inf <= eps1 ) )
			{
				Dp_L2 = 0.0; /* no increment for p in this case */
				if( op->cd->ldebug ) printf( "CONVERGED: Largest J^T e component ||J^T e||_inf is too small (%g < %g)\n", jacTe_inf, eps1 );
				stop = 1;
				break;
			}
			changejac = 0;
		}
		else printf( "\n" );
		/* compute initial damping factor */
		if( k == 0 )
		{
			for( i = 0, tmp = LM_REAL_MIN; i < m; ++i )
			{
				if( diag_jacTjac[i] > tmp ) tmp = diag_jacTjac[i]; /* find max diagonal element */
				if( op->cd->ldebug > 2 && diag_jacTjac[i] < 1 ) printf( "Parameter %s is not very sensitive (JTJ diagonal %g)\n", op->pd->var_id[op->pd->var_index[i]], diag_jacTjac[i] );
				if( diag_jacTjac[i] < DBL_EPSILON ) printf( "WARNING: Parameter %s is not impacting model predictions (JTJ diagonal %g)\n", op->pd->var_id[op->pd->var_index[i]], diag_jacTjac[i] );
			}
			if( tmp > 1e4 ) tmp = 1e4;
			mu = tau * tmp;
			if( op->cd->ldebug ) printf( "Computed initial lambda %g\n", mu );
		}
		/* determine increment using adaptive damping */
		/* augment normal equations */
		for( i = 0; i < m; ++i )
			jacTjac[i * m + i] += mu;  // Add lambda to the matrix diagonal
		/* solve augmented equations */
		issolved = linsolver( jacTjac, jacTe, Dp, m ); ++nlss;
		if( issolved )
		{
			if( op->cd->lm_acc ) // Acceleration
			{
				register LM_REAL beta, *jacT;
				for( i = 0; i < m; ++i )
					jacTvv[i] = 0.0;
				for( i = 0; i < m; ++i )
				{
					phDp_plus[i] = p[i] + acc_h * Dp[i];
					phDp_minus[i] = p[i] - acc_h * Dp[i];
				}
				change = 0;
				if( op->cd->compute_phi ) { op->cd->compute_phi = 0; change = 1; }
				( *func )( phDp_plus, hx1, m, n, adata ); nfev++;
				for( i = 0; i < n; ++i )
					ephdp_plus[i] = x[i] - hx1[i];
				( *func )( phDp_minus, hx2, m, n, adata ); nfev++;
				if( change ) op->cd->compute_phi = 1;
				for( i = 0; i < n; ++i )
					ephdp_minus[i] = x[i] - hx2[i];
				for( i = 0; i < n; ++i )
					vvddr[i] = ( ephdp_plus[i] - 2.0 * e[i] + ephdp_minus[i] ) / ( acc_h * acc_h );
				for( j = n; j-- > 0; )
				{
					jacT = jac + j * m;
					for( i = m; i-- > 0; )
					{
						beta = jacT[i]; //jac[l*m+i];
						/* J^T*vvddr */
						jacTvv[i] += beta * vvddr[j];
					}
				}
				issolved1 = linsolver( jacTjac, jacTvv, a, m ); ++nlss;
			}
			/* compute p's new estimate and ||Dp||^2 */
			for( i = 0, Dp_L2 = 0.0; i < m; ++i )
			{
				pDp[i] = p[i] + ( tmp = Dp[i] );
				Dp_L2 += tmp * tmp;
			}
			Dpa_L2 = Dp_L2;
			if( op->cd->lm_acc && issolved1 )
			{
				for( i = 0, a_L2 = 0.0; i < m; ++i )
				{
					pDp[i] += 0.5 * ( tmp = a[i] );
					a_L2 += tmp * tmp;
				}
				avRatio = a_L2 / Dp_L2;
				for( i = 0, Dpa_L2 = 0.0; i < m; ++i )
				{
					tmp = Dp[i] + 0.5 * a[i];
					Dpa_L2 += tmp * tmp;
				}
				if( op->cd->ldebug > 1 ) printf( "LM acceleration performed (with acceleration %g vs without acceleration %g)\n", Dpa_L2, Dp_L2 );
			}
			Dpa_L2 = sqrt( Dpa_L2 );
			if( op->cd->ldebug > 6 ) printf( "Test for convergence: Magnitude of the relative change in parameter space (%g < %g = %g * %g  to converge)\n", Dpa_L2, eps2_sq * p_L2, eps2_sq, p_L2 );
			if( Dpa_L2 <= eps2_sq * p_L2 ) /* relative change in p is small, stop */
			{
				if( op->cd->ldebug ) printf( "CONVERGED: Relative change in parameter space is small (%g < %g)\n", Dpa_L2, eps2_sq * p_L2 );
				stop = 2;
				break;
			}
			if( op->cd->ldebug > 6 ) printf( "Test for convergence: solution singularity (%g > %g to converge)\n", Dpa_L2, ( p_L2 + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) );
			if( Dpa_L2 >= ( p_L2 + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) ) /* almost singular */
			{
				if( op->cd->ldebug ) printf( "CONVERGED: almost singular solution (%g > %g)\n", Dpa_L2, ( p_L2 + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) );
				stop = 4;
				break;
			}
			( *func )( pDp, wrk, m, n, adata ); ++nfev; /* evaluate function at p + Dp */
#if 0
			if( op->cd->solution_type[0] == TEST ) // this is a test; not needed in general
				for( i = 0; i < n; i++ )
					wrk[i] = sqrt( wrk[i] );
#endif
			/* compute ||e(pDp)||_2 */
			/* ### wrk2=x-wrk, pDp_eL2=||wrk2|| */
#if 0
			pDp_eL2 = LEVMAR_L2NRMXMY( wrk2, x, wrk, n );
#else
			for( i = 0, pDp_eL2 = 0.0; i < n; ++i )
			{
				wrk2[i] = tmp = x[i] - wrk[i];
#if 1
				pDp_eL2 += tmp * tmp;
#else
				if( op->cd->solution_type[0] == TEST ) // this is a test; not needed in general
					pDp_eL2 += wrk[i];
				else
					pDp_eL2 += tmp * tmp;
#endif
			}
#endif
			if( op->cd->ldebug == 1 ) printf( "OF %g lambda %g\n", pDp_eL2, mu );
			else if( op->cd->ldebug > 1 )
			{
				printf( "\nLinear solve (lambda search) #%d: OF %g lambda %g \n", phi2c + 1, pDp_eL2, mu );
				DeTransform( pDp, op, jac_min );
				DeTransform( p_old, op, jac_max );
				int ipar_max, ipar_min;
				double max_change = 0, min_change = HUGE_VAL, p_diff;
				for( i = 0 ; i < m; i++ )
				{
					p_diff = fabs( jac_max[i] - jac_min[i] );
					if( max_change < p_diff ) { max_change = p_diff; ipar_max = i; }
					if( min_change > p_diff ) { min_change = p_diff; ipar_min = i; }
				}
				if( op->cd->ldebug > 2 )
				{
					printf( "Current parameter estimates:\n" );
					for( i = 0; i < m; i++ )
					{
						j = op->pd->var_index[i];
						printf( "%-20s:", op->pd->var_id[j] );
						if( op->cd->ldebug > 3 )
						{
							if( op->pd->var_log[j] ) printf( " %12g", pow( 10, jac_max[i] ) );
							else printf( " %12g", jac_max[i] );
							printf( " =>" );
						}
						if( op->pd->var_log[j] ) printf( " %12g", pow( 10, jac_min[i] ) );
						else printf( " %12g", jac_min[i] );
						if( op->cd->ldebug > 3 )
						{
							printf( " change" );
							if( op->pd->var_log[j] ) printf( " %12g",  pow( 10, jac_min[i] ) - pow( 10, jac_max[i] ) );
							else printf( " %12g", jac_min[i] - jac_max[i] );
						}
						printf( "  (JTJ diagonal term %g)", diag_jacTjac[i] );
						if( diag_jacTjac[i] < 1 ) printf( " not very sensitive" );
						else if( diag_jacTjac[i] < DBL_EPSILON ) printf( " WARNING: not impacting model predictions" );
						printf( "\n" );
					}
				}
				printf( "Parameter with maximum estimate change: %-20s (%g)\n", op->pd->var_id[op->pd->var_index[ipar_max]], jac_max[ipar_max] - jac_min[ipar_max] );
				printf( "Parameter with minimum estimate change: %-20s (%g)\n", op->pd->var_id[op->pd->var_index[ipar_min]], jac_max[ipar_min] - jac_min[ipar_min] );
				for( i = 0 ; i < m; i++ ) /* update p's estimate */
					p_old[i] = pDp[i];
			}
			phi2[phi2c++] = pDp_eL2;
			if( phi2c > op->cd->lm_nlamof )
			{
				tmp = HUGE_VAL;
				for( i = 0; i < phi2c; i++ )
				{
					if( phi2[i] < tmp )
					{
						tmp = phi2[i];
						if( i < phi2c - op->cd->lm_nlamof ) j = 1;
						else j = 0;
					}
					// printf( " %g", phi1[i] );
				}
				// for( i = 0; i < phi2c; i++ )
				// printf( " %g", phi2[i] );
				// printf( "a %g %d\n", tmp, phi2c );
				for( i = phi2c - op->cd->lm_nlamof; i < phi2c; i++ )
				{
					if( ( phi2[i] - tmp ) / tmp < 1 ) j++;
					// printf( "a %g %g %g\n", tmp, phi2[i], ( phi2[i] - tmp ) / tmp );
				}
				if( j >= op->cd->lm_nlamof )
					computejac = 1; // New jacobian: Lambda search OF are very similar
			}
			if( pDp_eL2 <= eps3 ) /* error is small */ // BELOW a cutoff value
			{
				if( op->cd->ldebug ) printf( "CONVERGED: OF below cutoff value (%g < %g)\n", pDp_eL2, eps3 );
				for( i = 0; i < m; i++ ) p[i] = pDp[i];
				p_eL2 = pDp_eL2;
				stop = 6;
				break;
			}
			if( op->cd->check_success && op->success )
			{
				if( op->cd->ldebug ) printf( "CONVERGED: Predictions are within predefined calibration ranges\n" );
				for( i = 0; i < m; i++ ) p[i] = pDp[i];
				p_eL2 = pDp_eL2;
				stop = 8;
				break;
			}
			if( !LM_FINITE( pDp_eL2 ) )
			{
				/* sum of squares is not finite, most probably due to a user error.
				 * This check makes sure that the loop terminates early in the case
				 * of invalid input. Thanks to Steve Danauskas for suggesting it */
				if( op->cd->ldebug ) printf( "CONVERGED: sum of squares is not finite, most probably due to a user error\n" );
				stop = 7;
				break;
			}
			if( op->cd->lm_indir ) // original
			{
				tmp = p_eL2_old / pDp_eL2; // original code
				//				printf( "Original tmp = %g\n", tmp );
				if( tmp > op->cd->lm_ofdecline ) phi_decline = 1; /* recompute jacobian because OF decreased */
			}
			else
			{
				tmp = p_eL2 / pDp_eL2;
				//				printf( "Leif tmp = %g\n", tmp );
				if( tmp > op->cd->lm_ofdecline && avRatio < lm_ratio ) phi_decline = 1; /* recompute jacobian because OF decreased */
			}
			dF = p_eL2 - pDp_eL2; // Difference between current and previous OF
			if( updp || dF > 0.0 ) /* update jacobian because OF increases */
			{
				for( i = 0; i < n; ++i )
				{
					for( l = 0, tmp = 0.0; l < m; ++l )
						tmp += jac[i * m + l] * Dp[l]; /* (J * Dp)[i] */
					tmp = ( wrk[i] - hx[i] - tmp ) / Dp_L2; /* (f(p+dp)[i] - f(p)[i] - (J * Dp)[i])/(dp^T*dp) */
					for( j = 0; j < m; ++j )
						jac[i * m + j] += tmp * Dp[j];
				}
				newjac = 1;
				++updjac;
			}
			if( op->cd->lm_indir )
			{
				for( i = 0, dL = 0.0; i < m; ++i )
					dL += Dp[i] * ( mu * Dp[i] + jacTe[i] );
			}
			else dL = ( double ) 1.0;
			if( dL > 0.0 && dF > 0.0 && avRatio < lm_ratio ) /* reduction in error, increment is accepted */
			{
				if( op->cd->lm_indir )
				{
					tmp = ( LM_CNST( 2.0 ) * dF / dL - LM_CNST( 1.0 ) );
					tmp = LM_CNST( 1.0 ) - tmp * tmp * tmp;
					tmp = ( ( tmp >= LM_CNST( ONE_THIRD ) ) ? tmp : LM_CNST( ONE_THIRD ) );
				}
				else tmp = op->cd->lm_mu;
				mu = mu * tmp; // change lambda
				/*
								if( mu > 1e3 )
								{
									if( mu_constrained > 500000 )
									{
										if( op->cd->ldebug ) printf( "lambda has been already constrained; new iteration" );
										mu_big = 1;
										mu_constrained = 0;
									}
									else
									{
										// mu = 1e3;
										mu_constrained++;
										if( op->cd->ldebug ) printf( "lambda is constrained to be less than %g", mu );
									}
								}
								else mu_constrained = 0;
				 */
				if( op->cd->ldebug > 3 ) printf( "OF change factor tmp (%g)\n", tmp );
				nu = op->cd->lm_nu;
				for( i = 0 ; i < m; ++i ) /* update p's estimate */
					p[i] = pDp[i];
				for( i = 0; i < n; ++i ) /* update e, hx and ||e||_2 */
				{
					e[i] = wrk2[i]; //x[i]-wrk[i];
					hx[i] = wrk[i];
				}
				p_eL2 = pDp_eL2; // Update OF
				updp = 1;
				continue; // Solve for a new lambda
			}
		}
		/* if this point is reached, either the linear system could not be solved or
		 * the error did not reduce; in any case, the increment must be rejected
		 */
		mu *= nu; // increase lambda
		/*
				if( mu > 1e3 )
				{
					if( mu_constrained > 5000000 )
					{
						if( op->cd->ldebug ) printf( "lambda has been already constrained; new iteration" );
						mu_big = 1;
						mu_constrained = 0;
					}
					else
					{
						// mu = 1e3;
						mu_constrained++;
						if( op->cd->ldebug ) printf( "lambda is constrained to be less than %g", mu );
					}
				}
				else mu_constrained = 0;
		 */
		if( op->cd->ldebug > 3 ) printf( "OF change factor (nu) %d\n", nu );
		if( op->cd->lm_indir )
		{
			nu2 = nu << 1; // 2 * nu;
			if( nu2 <= nu ) /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
			{
				if( op->cd->ldebug ) printf( "CONVERGED: lambda multiplication factor has wrapped around (overflown)\n" );
				stop = 5;
				break;
			}
			nu = nu2;
		}
		else { computejac = 0; }
		for( i = 0; i < m; ++i ) /* restore diagonal J^T J entries */
			jacTjac[i * m + i] = diag_jacTjac[i];
	}
	if( op->cd->ldebug > 3 )
	{
		DeTransform( pDp, op, jac_min );
		DeTransform( p_old, op, jac_max );
		int ipar_max, ipar_min;
		double max_change = 0, min_change = HUGE_VAL, p_diff;
		for( i = 0 ; i < m; ++i )
		{
			p_diff = fabs( jac_max[i] - jac_min[i] );
			if( max_change < p_diff ) { max_change = p_diff; ipar_max = i; }
			if( min_change > p_diff ) { min_change = p_diff; ipar_min = i; }
		}
		printf( "Parameter with maximum estimate change: %s (%g)\n", op->pd->var_id[op->pd->var_index[ipar_max]], jac_max[ipar_max] - jac_min[ipar_max] );
		printf( "Parameter with minimum estimate change: %s (%g)\n", op->pd->var_id[op->pd->var_index[ipar_min]], jac_max[ipar_min] - jac_min[ipar_min] );
		for( i = 0 ; i < m; ++i ) /* update p's estimate */
			p_old[i] = pDp[i];
	}
	if( k >= kmax && stop == 0 ) stop = 3;
	for( i = 0; i < m; ++i ) /* restore diagonal J^T J entries */
		jacTjac[i * m + i] = diag_jacTjac[i];
	if( info )
	{
		info[0] = init_p_eL2;
		info[1] = p_eL2;
		info[2] = jacTe_inf;
		info[3] = Dp_L2;
		for( i = 0, tmp = LM_REAL_MIN; i < m; ++i )
			if( tmp < jacTjac[i * m + i] ) tmp = jacTjac[i * m + i];
		info[4] = mu / tmp;
		info[5] = ( LM_REAL )k;
		info[6] = ( LM_REAL )stop;
		info[7] = ( LM_REAL )nfev;
		info[8] = ( LM_REAL )njap;
		info[9] = ( LM_REAL )nlss;
	}
	if( op->cd->ldebug )
	{
		printf( "LM optimization is completed. Reason: " );
		switch( stop )
		{
			case 1: printf( "small gradient J^T e (%g)\n", jacTe_inf ); break;
			case 2: printf( "small Dp (%g)\n", Dp_L2 ); break;
			case 3: printf( "maximum number of LevMar iterations is exceeded (kmax=%d)\n", kmax ); break;
			case 31: printf( "maximum number of jacobian iterations is exceeded (lmiter=%d)\n", itmax ); break;
			case 32: printf( "maximum number of functional evaluations is exceeded (eval=%d; %d > %d)\n", op->cd->maxeval, nfev, maxnfev ); break;
			case 4: printf( "singular matrix. Restart from current p with increased mu (current mu=%g; mu/max[J^T J]_ii=%g)\n", mu, mu / tmp ); break;
			case 5: printf( "no further error reduction is possible. Restart with increased mu (current mu=%g; mu/max[J^T J]_ii=%g)\n", mu, mu / tmp ); break;
			case 6: printf( "small ||e||_2; OF below cutoff value (%g < %g)\n", p_eL2, eps3 ); break;
			case 7: printf( "invalid (i.e. NaN or Inf) values returned by the solver (func). This is a user error\n" ); break;
			case 8: printf( "model predictions are within predefined calibration ranges\n" ); break;
			case 9: printf( "small OF changes\n" ); break;
			default: printf( "UNKNOWN flag: %d\n", stop ); break;
		}
		if( op->cd->ldebug > 2 )
		{
			printf( "Levenberg-Marquardt Optimization completed after %g iteration (reason %g) (returned value %d)\n", info[5], info[6], ( stop != 4 && stop != 7 ) ?  k : LM_ERROR );
			printf( "initial phi %g final phi %g ||J^T e||_inf %g ||Dp||_2 %g mu/max[J^T J]_ii %g\n", info[0], info[1], info[2], info[3], info[4] );
			printf( "function evaluation %g jacobian evaluations %g linear systems solved %g\n", info[7], info[8], info[9] );
		}
	}
	else if( op->cd->standalone ) { printf( "%g", p_eL2 ); fflush( stdout ); }
	/* covariance matrix */
	if( covar )
	{
		LEVMAR_COVAR( jacTjac, covar, p_eL2, m, n );
	}
	if( freework ) free( work );
#ifdef LINSOLVERS_RETAIN_MEMORY
	if( linsolver )( *linsolver )( NULL, NULL, NULL, 0 );
#endif
	return ( stop != 4 && stop != 7 ) ?  k : LM_ERROR;
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
	register int i, j, k, l;
	int worksz, freework = 0, issolved, success, odebug, kmax, maxnfev;
	struct opt_data *op = ( struct opt_data * ) adata;
	/* temp work arrays */
	LM_REAL *e,          /* nx1 */
			*hx,         /* \hat{x}_i, nx1 */
			*jacTe,      /* J^T e_i mx1 */
			*jac,        /* nxm */
			*jacTjac,    /* mxm */
			*Dp,         /* mx1 */
			*diag_jacTjac,   /* diagonal of J^T J, mx1 */
			*pDp,        /* p + Dp, mx1 */
			*wrk,        /* nx1 */
			*wrk2;       /* nx1, used only for holding a temporary e vector and when differentiating with central differences */
	int using_ffdif = 1;
	register LM_REAL mu,  /* damping constant */
			 tmp; /* mainly used in matrix & vector multiplications */
	LM_REAL p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
	LM_REAL p_L2, Dp_L2 = LM_REAL_MAX, dF, dL;
	LM_REAL tau, eps1, eps2, eps2_sq, eps3, delta;
	LM_REAL init_p_eL2;
	int nu, nu2, stop = 0, nfev, njap = 0, nlss = 0, K = ( m >= 10 ) ? m : 10, updjac, updp = 1, newjac;
	const int nm = n * m;
	int ( *linsolver )( LM_REAL * A, LM_REAL * B, LM_REAL * x, int m ) = NULL;
	mu = jacTe_inf = p_L2 = 0.0;
	updjac = newjac = 0;
	if( n < m )
	{
		fprintf( stderr, LCAT( LEVMAR_DIF, "(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n" ), n, m );
		return LM_ERROR;
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
		worksz = LM_DIF_WORKSZ( m, n ); //4*n+4*m + n*m + m*m;
		work = ( LM_REAL * )malloc( worksz * sizeof( LM_REAL ) ); /* allocate a big chunk in one step */
		if( !work )
		{
			fprintf( stderr, LCAT( LEVMAR_DIF, "(): memory allocation request failed\n" ) );
			return LM_ERROR;
		}
		freework = 1;
	}
	/* set up work arrays */
	e = work;
	hx = e + n;
	jacTe = hx + n;
	jac = jacTe + m;
	jacTjac = jac + nm;
	Dp = jacTjac + m * m;
	diag_jacTjac = Dp + m;
	pDp = diag_jacTjac + m;
	wrk = pDp + m;
	wrk2 = wrk + n;
	/* compute e=x - f(p) and its L2 norm */
	maxnfev = op->cd->maxeval - op->cd->neval;
	( *func )( p, hx, m, n, adata ); nfev = 1;
	if( op->cd->check_success && op->success )
	{
		if( op->cd->ldebug ) printf( "SUCCESS: Model predictions are within predefined calibration ranges\n" );
		stop = 8;
	}
#ifdef HAVE_LAPACK
	/* 6 alternatives are available: LU, Cholesky, 2 variants of QR decomposition, SVD and LDLt.
	 * Cholesky is the fastest but might be inaccurate; QR is slower but more accurate;
	 * SVD is the slowest but most accurate; LU offers a tradeoff between accuracy and speed
	 */
	//linsolver = AX_EQ_B_BK; if( op->cd->ldebug ) printf( "BK decomposition\n" );
	//linsolver = AX_EQ_B_LU; if( op->cd->ldebug ) printf( "LU decomposition\n" );
	//linsolver = AX_EQ_B_CHOL; if( op->cd->ldebug ) printf( "Cholesky decomposition\n" );
	//linsolver = AX_EQ_B_QR; if( op->cd->ldebug ) printf( "QR decomposition\n" );
	//linsolver = (int (*)(LM_REAL *A, LM_REAL *B, LM_REAL *x, int m))AX_EQ_B_QRLS; if( op->cd->ldebug ) printf( "QRLS decomposition\n" );
	linsolver = AX_EQ_B_SVD; if( op->cd->ldebug ) printf( "LM using: SVD decomposition\n" );
#else
	/* use the LU included with levmar */
	linsolver = AX_EQ_B_LU; if( op->cd->ldebug ) printf( "LM using: LU decomposition\n" );
#endif
	/* ### e=x-hx, p_eL2=||e|| */
#if 0
	p_eL2 = LEVMAR_L2NRMXMY( e, x, hx, n );
#else
	for( i = 0, p_eL2 = 0.0; i < n; ++i )
	{
		e[i] = tmp = x[i] - hx[i];
		p_eL2 += tmp * tmp;
	}
#endif
	init_p_eL2 = p_eL2;
	if( !LM_FINITE( p_eL2 ) ) stop = 7;
	nu = 20; /* force computation of J */
	if( op->cd->check_success ) success = 1; else success = 0;
	if( op->cd->odebug ) odebug = 1; else odebug = 0;
	kmax = itmax * 100;
	if( op->cd->ldebug ) printf( "Initial evaluation: OF %g\n", p_eL2 );
	for( k = 0; k < kmax && !stop; ++k )
	{
		/* Compute the Jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
		 * The symmetry of J^T J is again exploited for speed
		 */
		if( ( updp && nu > 16 ) || updjac == K ) /* compute difference approximation to J */
		{
			if( op->cd->ldebug )
			{
				printf( "New Jacobian requested because: " );
				if( nu > 16 ) printf( "Lambda multiplication factor too large (nu = %d > 16); ", nu );
				if( updjac == K ) printf( "Maximum number of lambda iteration is reached (%d); ", K );
				printf( "\n\n" );
			}
			if( njap >= itmax ) { stop = 31; continue; }
			if( nfev >= maxnfev ) { stop = 32; continue; }
			if( success ) op->cd->check_success = 0;
			if( odebug ) op->cd->odebug = 0;
			if( using_ffdif ) /* use forward differences */
			{
				LEVMAR_FDIF_FORW_JAC_APPROX( func, p, hx, wrk, delta, jac, m, n, adata );
				++njap; nfev += m;
			}
			else  /* use central differences */
			{
				if( op->cd->ldebug ) printf( "Central Differences\n" );
				LEVMAR_FDIF_CENT_JAC_APPROX( func, p, wrk, wrk2, delta, jac, m, n, adata );
				++njap; nfev += 2 * m;
			}
			if( success ) op->cd->check_success = 1;
			if( odebug ) op->cd->odebug = 1;
			nu = 10; updjac = 0; updp = 0; newjac = 1;
			if( op->cd->ldebug )
			{
				if( k == 0 ) printf( "Jacobians %d Linear solves %d Evaluations %d OF %g lambda %g\n", njap, nlss, nfev, p_eL2, tau );
				else printf( "Jacobians %d Linear solves %d Evaluations %d OF %g lambda %g\n", njap, nlss, nfev, p_eL2, mu );
			}
			if( op->cd->ldebug > 13 )
			{
				printf( "Jacobian matrix:\n" );
				for( l = j = 0; j < op->od->nObs; j++ )
				{
					if( op->cd->ldebug > 14 || op->od->nObs < 50 || ( j < 20 || j > op->od->nObs - 20 ) )
					{
						printf( "Observation %i: ", j + 1 );
						for( i = 0; i < op->pd->nOptParam; i++, l++ )
							printf( " %g", jac[l] );
					}
					else
						for( i = 0; i < op->pd->nOptParam; i++, l++ );
					if( ( !( op->cd->ldebug > 14 ) || op->od->nObs > 50 ) && j == 21 ) printf( "...\n" );
					printf( "\n" );
				}
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
				/* looping downwards saves a few computations */
				for( i = m * m; i-- > 0; )
					jacTjac[i] = 0.0;
				for( i = m; i-- > 0; )
					jacTe[i] = 0.0;
				for( l = n; l-- > 0; )
				{
					jaclm = jac + l * m;
					for( i = m; i-- > 0; )
					{
						im = i * m;
						alpha = jaclm[i]; //jac[l*m+i];
						for( j = i + 1; j-- > 0; ) /* j<=i computes lower triangular part only */
							jacTjac[im + j] += jaclm[j] * alpha; //jac[l*m+j]
						/* J^T e */
						jacTe[i] += alpha * e[l];
					}
				}
				for( i = m; i-- > 0; ) /* copy to upper part */
					for( j = i + 1; j < m; ++j )
						jacTjac[i * m + j] = jacTjac[j * m + i];
			}
			else  // this is a large problem
			{
				/* Cache efficient computation of J^T J based on blocking
				 */
				LEVMAR_TRANS_MAT_MAT_MULT( jac, jacTjac, n, m );
				/* cache efficient computation of J^T e */
				for( i = 0; i < m; ++i )
					jacTe[i] = 0.0;
				for( i = 0; i < n; ++i )
				{
					register LM_REAL *jacrow;
					for( l = 0, jacrow = jac + i * m, tmp = e[i]; l < m; ++l )
						jacTe[l] += jacrow[l] * tmp;
				}
			}
			/* Compute ||J^T e||_inf and ||p||^2 */
			for( i = 0, p_L2 = jacTe_inf = 0.0; i < m; ++i )
			{
				if( jacTe_inf < ( tmp = FABS( jacTe[i] ) ) ) jacTe_inf = tmp;
				diag_jacTjac[i] = jacTjac[i * m + i]; /* save diagonal entries so that augmentation can be later canceled */
				p_L2 += p[i] * p[i];
			}
			//p_L2=sqrt(p_L2);
			if( op->cd->ldebug > 4 ) printf( "||J^T e||_inf %g (%g < %g to converge)\n", jacTe_inf, jacTe_inf, eps1 );
			/* check for convergence */
			if( ( jacTe_inf <= eps1 ) )
			{
				Dp_L2 = 0.0; /* no increment for p in this case */
				if( op->cd->ldebug ) printf( "CONVERGED: ||J^T e||_inf is too small (%g < %g)\n", jacTe_inf, eps1 );
				stop = 1;
				break;
			}
		}
		/* compute initial damping factor */
		if( k == 0 )
		{
			for( i = 0, tmp = LM_REAL_MIN; i < m; ++i )
				if( diag_jacTjac[i] > tmp ) tmp = diag_jacTjac[i]; /* find max diagonal element */
			mu = tau * tmp;
			if( op->cd->ldebug ) printf( "Computed initial lambda %g\n", mu );
		}
		/* determine increment using adaptive damping */
		/* augment normal equations */
		for( i = 0; i < m; ++i )
			jacTjac[i * m + i] += mu; // Add lambda to the matrix diaganol
		/* solve augmented equations */
		issolved = linsolver( jacTjac, jacTe, Dp, m ); ++nlss;
		if( issolved )
		{
			/* compute p's new estimate and ||Dp||^2 */
			for( i = 0, Dp_L2 = 0.0; i < m; ++i )
			{
				pDp[i] = p[i] + ( tmp = Dp[i] );
				Dp_L2 += tmp * tmp;
			}
			// Dp_L2=sqrt(Dp_L2);
			if( op->cd->ldebug > 4 ) printf( "Relative change in the OF %g (%g < %g to converge)\n", Dp_L2, Dp_L2, eps2_sq * p_L2 );
			if( Dp_L2 <= eps2_sq * p_L2 ) /* relative change in p is small, stop */
			{
				//if(Dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
				if( op->cd->ldebug ) printf( "CONVERGED: Relative change in the OF is small (%g < %g)\n", Dp_L2, eps2_sq * p_L2 );
				stop = 2;
				break;
			}
			if( Dp_L2 >= ( p_L2 + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) ) /* almost singular */
			{
				//if(Dp_L2>=(p_L2+eps2)/LM_CNST(EPSILON)){ /* almost singular */
				if( op->cd->ldebug ) printf( "CONVERGED: almost singular solution (%g > %g)\n", Dp_L2, ( p_L2 + eps2 ) / ( LM_CNST( EPSILON )*LM_CNST( EPSILON ) ) );
				stop = 4;
				break;
			}
			( *func )( pDp, wrk, m, n, adata ); ++nfev; /* evaluate function at p + Dp */
			/* compute ||e(pDp)||_2 */
			/* ### wrk2=x-wrk, pDp_eL2=||wrk2|| */
#if 1
			pDp_eL2 = LEVMAR_L2NRMXMY( wrk2, x, wrk, n );
#else
			for( i = 0, pDp_eL2 = 0.0; i < n; ++i )
			{
				wrk2[i] = tmp = x[i] - wrk[i];
				pDp_eL2 += tmp * tmp;
			}
#endif
			if( op->cd->ldebug ) printf( "OF %g lambda %g ", pDp_eL2, mu );
			if( pDp_eL2 <= eps3 ) /* error is small */ // BELOW a cutoff value
			{
				if( op->cd->ldebug ) printf( "CONVERGED: OF below cutoff value (%g < %g)\n", p_eL2, eps3 );
				for( i = 0; i < m; i++ ) p[i] = pDp[i];
				p_eL2 = pDp_eL2;
				stop = 6;
				break;
			}
			if( op->cd->check_success && op->success )
			{
				if( op->cd->ldebug ) printf( "CONVERGED: Predictions are within predefined calibration ranges\n" );
				for( i = 0; i < m; i++ ) p[i] = pDp[i];
				p_eL2 = pDp_eL2;
				stop = 8;
				break;
			}
			if( !LM_FINITE( pDp_eL2 ) )
			{
				/* sum of squares is not finite, most probably due to a user error.
				 * This check makes sure that the loop terminates early in the case
				 * of invalid input. Thanks to Steve Danauskas for suggesting it */
				if( op->cd->ldebug ) printf( "CONVERGED: sum of squares is not finite, most probably due to a user error\n" );
				stop = 7;
				break;
			}
			dF = p_eL2 - pDp_eL2; // Difference between current and previous OF
			// printf( "dF = %g (%g - %g)\n", dF, p_eL2, pDp_eL2 );
			if( updp || dF > 0.0 ) /* update jac because OF increases */
			{
				for( i = 0; i < n; ++i )
				{
					for( l = 0, tmp = 0.0; l < m; ++l )
						tmp += jac[i * m + l] * Dp[l]; /* (J * Dp)[i] */
					tmp = ( wrk[i] - hx[i] - tmp ) / Dp_L2; /* (f(p+dp)[i] - f(p)[i] - (J * Dp)[i])/(dp^T*dp) */
					for( j = 0; j < m; ++j )
						jac[i * m + j] += tmp * Dp[j];
				}
				++updjac;
				newjac = 1;
			}
			for( i = 0, dL = 0.0; i < m; ++i )
				dL += Dp[i] * ( mu * Dp[i] + jacTe[i] );
			if( dL > 0.0 && dF > 0.0 ) /* reduction in error, increment is accepted */
			{
				tmp = ( LM_CNST( 2.0 ) * dF / dL - LM_CNST( 1.0 ) );
				tmp = LM_CNST( 1.0 ) - tmp * tmp * tmp;
				tmp = ( ( tmp >= LM_CNST( ONE_THIRD ) ) ? tmp : LM_CNST( ONE_THIRD ) );
				mu = mu * tmp; // change lambda
				if( op->cd->ldebug > 1 ) printf( "change factor (tmp) %g\n", tmp );
				else if( op->cd->ldebug ) printf( "\n" );
				nu = 2;
				for( i = 0 ; i < m; ++i ) /* update p's estimate */
					p[i] = pDp[i];
				for( i = 0; i < n; ++i ) /* update e, hx and ||e||_2 */
				{
					e[i] = wrk2[i]; //x[i]-wrk[i];
					hx[i] = wrk[i];
				}
				p_eL2 = pDp_eL2; // Update OF
				updp = 1;
				continue; // Solve for a new lambda
			}
		}
		/* if this point is reached, either the linear system could not be solved or
		 * the error did not reduce; in any case, the increment must be rejected
		 */
		mu *= nu; // increase lambda
		if( op->cd->ldebug > 1 ) printf( "change factor (nu) %d\n", nu );
		else if( op->cd->ldebug ) printf( "\n" );
		nu2 = nu << 1; // 2*nu;
		if( nu2 <= nu ) /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
		{
			if( op->cd->ldebug ) printf( "CONVERGED: lambda multiplication factor has wrapped around (overflown)\n" );
			stop = 5;
			break;
		}
		nu = nu2;
		for( i = 0; i < m; ++i ) /* restore diagonal J^T J entries */
			jacTjac[i * m + i] = diag_jacTjac[i];
	}
	if( k >= kmax && stop == 0 ) stop = 3;
	for( i = 0; i < m; ++i ) /* restore diagonal J^T J entries */
		jacTjac[i * m + i] = diag_jacTjac[i];
	if( info )
	{
		info[0] = init_p_eL2;
		info[1] = p_eL2;
		info[2] = jacTe_inf;
		info[3] = Dp_L2;
		for( i = 0, tmp = LM_REAL_MIN; i < m; ++i )
			if( tmp < jacTjac[i * m + i] ) tmp = jacTjac[i * m + i];
		info[4] = mu / tmp;
		info[5] = ( LM_REAL )k;
		info[6] = ( LM_REAL )stop;
		info[7] = ( LM_REAL )nfev;
		info[8] = ( LM_REAL )njap;
		info[9] = ( LM_REAL )nlss;
	}
	if( op->cd->ldebug )
	{
		printf( "LM optimization is completed. Reason: " );
		switch( stop )
		{
			case 1: printf( "small gradient J^T e (%g)\n", jacTe_inf ); break;
			case 2: printf( "small Dp (%g)\n", Dp_L2 ); break;
			case 3: printf( "maximum number of LevMar iterations is exceeded (kmax=%d)\n", kmax ); break;
			case 31: printf( "maximum number of jacobian iterations is exceeded (lmiter=%d)\n", itmax ); break;
			case 32: printf( "maximum number of functional evaluations is exceeded (eval=%d; %d > %d)\n", op->cd->maxeval, nfev, maxnfev ); break;
			case 4: printf( "singular matrix. Restart from current p with increased mu (current mu=%g; mu/max[J^T J]_ii=%g)\n", mu, mu / tmp ); break;
			case 5: printf( "no further error reduction is possible. Restart with increased mu (current mu=%g; mu/max[J^T J]_ii=%g)\n", mu, mu / tmp ); break;
			case 6: printf( "small ||e||_2; OF below cutoff value (%g < %g)\n", p_eL2, eps3 ); break;
			case 7: printf( "invalid (i.e. NaN or Inf) values returned by the solver (func). This is a user error\n" ); break;
			case 8: printf( "model predictions are within predefined calibration ranges\n" ); break;
			default: printf( "UNKNOWN flag: %d\n", stop ); break;
		}
		if( op->cd->ldebug > 2 )
		{
			printf( "Levenberg-Marquardt Optimization completed after %g iteration (reason %g) (returned value %d)\n", info[5], info[6], ( stop != 4 && stop != 7 ) ?  k : LM_ERROR );
			printf( "initial phi %g final phi %g ||J^T e||_inf %g ||Dp||_2 %g mu/max[J^T J]_ii %g\n", info[0], info[1], info[2], info[3], info[4] );
			printf( "function evaluation %g jacobian evaluations %g linear systems solved %g\n", info[7], info[8], info[9] );
		}
	}
	/* covariance matrix */
	if( covar )
	{
		LEVMAR_COVAR( jacTjac, covar, p_eL2, m, n );
	}
	if( freework ) free( work );
#ifdef LINSOLVERS_RETAIN_MEMORY
	if( linsolver )( *linsolver )( NULL, NULL, NULL, 0 );
#endif
	return ( stop != 4 && stop != 7 ) ?  k : LM_ERROR;
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
