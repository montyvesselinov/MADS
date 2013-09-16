// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
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


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "mads.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#ifdef MATHEVAL
#include <matheval.h>
#endif

/* Functions here */
int check_mads_problem( char *filename );
int set_param_id( struct opt_data *op );
int set_param_names( struct opt_data *op );
void init_params( struct opt_data *op );
int parse_cmd_debug( char *buf );
int parse_cmd( char *buf, struct calc_data *cd );
int load_problem( char *filename, int argn, char *argv[], struct opt_data *op );
int save_problem( char *filename, struct opt_data *op );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc2( char *filename, char *filename2, struct opt_data *op );
void compute_btc( char *filename, struct opt_data *op );
static char *strsave( const char *s, const char *lim );
char **shellpath( void );
void freeshellpath( char *shellpath[] );
unsigned maxpathlen( char *path[], const char *base );
void execvepath( char *path[], const char *base, char *const argv[], char *const envp[] );
int count_lines( char *filename );
int count_cols( char *filename, int row );
char *timestamp(); // create time stamp
char *datestamp(); // create date stamp
char *str_replace( char *orig, char *rep, char *with ); // replace all string occurrences
int set_optimized_params( struct opt_data *op );
int map_obs( struct opt_data *op );
int map_well_obs( struct opt_data *op );

/* Functions elsewhere */
char **char_matrix( int maxCols, int maxRows );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
int set_test_problems( struct opt_data *op );
void *malloc_check( const char *what, size_t n );
int Ftest( char *filename );
FILE *Fread( char *filename );
void removeChars( char *str, char *garbage );
char *white_trim( char *x );
void white_skip( char **s );

int check_mads_problem( char *filename )
{
	FILE *infile;
	char buf[5000], *word;
	char *separator = " \t\n";
	int c;
	if( ( infile = fopen( filename, "r" ) ) == NULL ) return( 0 ); // No file
	fscanf( infile, "%1000[^\n]s", buf );
	fclose( infile );
	// check the first line in the file
	for( c = 0, word = strtok( buf, separator ); word; c++, word = strtok( NULL, separator ) )
	{
		if( !strncasecmp( word, "yaml", 4 ) ) return( 1 ); // YAML file
		if( !strncmp( word, "---", 3 ) && c == 0 ) return( 1 ); // YAML file
		if( !strncmp( word, "{", 1 ) ) return( 1 ); // YAML file
	}
	return( 2 ); // Not YAML file
}

int set_param_id( struct opt_data *op )
{
	op->cd->num_aquifer_params = NUM_ANAL_PARAMS_AQUIFER;
	op->cd->num_source_params = NUM_ANAL_PARAMS_SOURCE;
	if( ( op->ad->var = ( double * ) malloc( ( op->cd->num_source_params + op->cd->num_aquifer_params ) * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	op->sd->param_id = char_matrix( op->cd->num_source_params, 6 );
	op->qd->param_id = char_matrix( op->cd->num_aquifer_params, 6 );
	strcpy( op->sd->param_id[0], "x" ); strcpy( op->sd->param_id[1], "y" ); strcpy( op->sd->param_id[2], "z" );
	strcpy( op->sd->param_id[3], "dx" ); strcpy( op->sd->param_id[4], "dy" ); strcpy( op->sd->param_id[5], "dz" );
	strcpy( op->sd->param_id[6], "f" ); strcpy( op->sd->param_id[7], "t0" ); strcpy( op->sd->param_id[8], "t1" );
	strcpy( op->qd->param_id[0], "n" ); strcpy( op->qd->param_id[1], "rf" ); strcpy( op->qd->param_id[2], "lambda" );
	strcpy( op->qd->param_id[3], "tetha" ); strcpy( op->qd->param_id[4], "vx" ); strcpy( op->qd->param_id[5], "vy" ); strcpy( op->qd->param_id[6], "vz" );
	strcpy( op->qd->param_id[7], "ax" ); strcpy( op->qd->param_id[8], "ay" ); strcpy( op->qd->param_id[9], "az" );
	strcpy( op->qd->param_id[10], "ts_dsp" ); strcpy( op->qd->param_id[11], "ts_adv" ); strcpy( op->qd->param_id[12], "ts_rct" );
	strcpy( op->qd->param_id[13], "alpha" ); strcpy( op->qd->param_id[14], "beta" ); strcpy( op->qd->param_id[15], "nlc0" );
	op->sd->param_name = char_matrix( op->cd->num_source_params, 50 );
	op->qd->param_name = char_matrix( op->cd->num_aquifer_params, 50 );
	strcpy( op->sd->param_name[0], "Source x coordinate [L]" ); strcpy( op->sd->param_name[1], "Source y coordinate [L]" ); strcpy( op->sd->param_name[2], "Source z coordinate [L]" );
	strcpy( op->sd->param_name[3], "Source x dimension [L]" ); strcpy( op->sd->param_name[4], "Source y dimension [L]" ); strcpy( op->sd->param_name[5], "Source z dimension [L]" );
	strcpy( op->sd->param_name[6], "Contaminant flux [M/T]" ); strcpy( op->sd->param_name[7], "Start Time [T]" ); strcpy( op->sd->param_name[8], "End Time [T]" );
	strcpy( op->qd->param_name[0], "Porosity [L3/L3]" ); strcpy( op->qd->param_name[1], "Retardation Factor [-]" ); strcpy( op->qd->param_name[2], "Half-life decay [1/T]" );
	strcpy( op->qd->param_name[3], "Flow Angle [degrees]" ); strcpy( op->qd->param_name[4], "Pore x velocity [L/T]" ); strcpy( op->qd->param_name[5], "Pore y velocity [L/T]" ); strcpy( op->qd->param_name[6], "Pore z velocity [L/T]" );
	strcpy( op->qd->param_name[7], "Dispersivity x [L]" ); strcpy( op->qd->param_name[8], "Dispersivity y [L]" ); strcpy( op->qd->param_name[9], "Dispersivity z [L]" );
	strcpy( op->qd->param_name[10], "Time Scale Dispersivity [-]" ); strcpy( op->qd->param_name[11], "Time Scale Advection [-]" ); strcpy( op->qd->param_name[12], "Time Scale Reaction [-]" );
	strcpy( op->qd->param_name[13], "Levy alpha [-]" ); strcpy( op->qd->param_name[14], "Levy beta [-]" ); strcpy( op->qd->param_name[15], "NLC 0 [-]" );
	return( 1 );
}

int set_param_names( struct opt_data *op )
{
	int i, j, k;
	if( op->cd->num_sources == 1 )
	{
		for( j = 0; j < op->cd->num_source_params; j++ )
		{
			strcpy( op->pd->var_name[j], op->sd->param_name[j] );
			strcpy( op->pd->var_id[j], op->qd->param_id[j] );
		}
		k = op->cd->num_source_params;
	}
	else
		for( k = 0, i = 0; i < op->cd->num_sources; i++ )
		{
			for( k = 0, j = 0; j < op->cd->num_source_params; j++, k++ )
			{
				sprintf( op->pd->var_name[k], "%s #%d", op->sd->param_name[j], i + 1 );
				strcpy( op->pd->var_id[k], op->sd->param_id[j] );
			}
		}
	for( i = 0; i < op->cd->num_aquifer_params; i++, k++ )
	{
		strcpy( op->pd->var_name[k], op->qd->param_name[i] );
		strcpy( op->pd->var_id[k], op->qd->param_id[i] );
	}
	return( 1 );
}

void init_params( struct opt_data *op )
{
	struct opt_data *p = ( struct opt_data * ) op;
	struct param_data *pd;
	struct calc_data *cd;
	int k;
	cd = p->cd;
	pd = p->pd;
	if( cd->num_sources > 1 ) k = cd->num_source_params * ( cd->num_sources - 1 );
	else k = 0;
	cd->var[k + TSCALE_DISP] = pd->var[k + TSCALE_DISP] = 1;
	cd->var[k + TSCALE_ADV] = pd->var[k + TSCALE_ADV] = 1;
	cd->var[k + TSCALE_REACT] = pd->var[k + TSCALE_REACT] = 1;
	pd->var_opt[k + TSCALE_DISP] = 0; pd->var_opt[k + TSCALE_ADV] = 0; pd->var_opt[k + TSCALE_REACT] = 0;
	pd->var_log[k + TSCALE_DISP] = 0; pd->var_log[k + TSCALE_ADV] = 0; pd->var_log[k + TSCALE_REACT] = 0;
	pd->var_dx[k + TSCALE_DISP] = 0.1; pd->var_dx[k + TSCALE_ADV] = 0.1; pd->var_dx[k + TSCALE_REACT] = 0.1;
	pd->var_min[k + TSCALE_DISP] = 0.1; pd->var_min[k + TSCALE_ADV] = 0.1; pd->var_min[k + TSCALE_REACT] = 0.1;
	pd->var_max[k + TSCALE_DISP] = 10; pd->var_max[k + TSCALE_ADV] = 10; pd->var_max[k + TSCALE_REACT] = 10;

	cd->var[k + ALPHA] = pd->var[k + ALPHA] = 2.;
	cd->var[k + BETA] = pd->var[k + BETA] = 0.;
	cd->var[k + NLC0] = pd->var[k + NLC0] = 1.;
	pd->var_opt[k + ALPHA] = 0; pd->var_opt[k + BETA] = 0; pd->var_opt[k + NLC0] = 0;
	pd->var_log[k + ALPHA] = 0; pd->var_log[k + BETA] = 0; pd->var_log[k + NLC0] = 0;
	pd->var_dx[k + ALPHA] = 0.1; pd->var_dx[k + BETA] = 0.1; pd->var_dx[k + NLC0] = 0.1;
	pd->var_min[k + ALPHA] = 0.5; pd->var_min[k + BETA] = -1.; pd->var_min[k + NLC0] = 0.1;
	pd->var_max[k + ALPHA] = 2; pd->var_max[k + BETA] = 1.; pd->var_max[k + NLC0] = 2.;
}

int parse_cmd_debug( char *buf )
{
	int debug;
	char *sep = " \t\n", *word;
	for( word = strtok( buf, sep ); word; word = strtok( NULL, sep ) )
	{
		if( !strncasecmp( word, "debug", 5 ) ) { if( sscanf( word, "debug=%d", &debug ) == 0 || debug == 0 ) debug = 1; } // Global debug
	}
	return( debug );
}

int parse_cmd( char *buf, struct calc_data *cd )
{
	int w;
	char *sep = " \t\n", *word;
	cd->opt_method = ( char * ) malloc( 50 * sizeof( char ) ); cd->opt_method[0] = 0;
	cd->smp_method = ( char * ) malloc( 50 * sizeof( char ) ); cd->smp_method[0] = 0;
	cd->paran_method = ( char * ) malloc( 50 * sizeof( char ) ); cd->paran_method[0] = 0;
	cd->infile = ( char * ) malloc( 255 * sizeof( char ) ); cd->infile[0] = 0;
	cd->resultsfile = ( char * ) malloc( 255 * sizeof( char ) ); cd->resultsfile[0] = 0;
	cd->resultscase = 0;
	cd->restart_zip_file = ( char * ) malloc( 255 * sizeof( char ) ); cd->restart_zip_file[0] = 0;
	strcpy( cd->opt_method, "lm" );
	cd->problem_type = UNKNOWN;
	cd->calib_type = SIMPLE;
	cd->solution_type[0] = EXTERNAL;
	cd->levy = 0;
	cd->objfunc_type = SSR;
	cd->check_success = 0;
	cd->c_background = 0;
	cd->debug = 0;
	cd->fdebug = 0;
	cd->ldebug = 0;
	cd->lm_eigen = 0;
	cd->pdebug = 0;
	cd->mdebug = 0;
	cd->odebug = 0;
	cd->insdebug = 0;
	cd->tpldebug = 0;
	cd->pardebug = 0;
	cd->sintrans = 1;
	cd->plogtrans = -1;
	cd->ologtrans = -1;
	cd->oweight = -1;
	cd->num_proc = -1;
	cd->restart = 1;
	cd->nreal = 0;
	cd->niter = 0;
	cd->disp_tied = 0;
	cd->disp_scaled = 0;
	cd->save = 0;
	cd->pargen = 0;
	cd->init_particles = -1;
	cd->nretries = 0;
	cd->seed = cd->seed_init = 0;
	cd->maxeval = 5000;
	cd->paranoid = 0;
	cd->phi_cutoff = 0;
	cd->parerror = -1;
	cd->obserror = -1;
	cd->obsdomain = -1;
	cd->obsstep = 0;
	cd->obsrange = 0;
	cd->test_func = -1;
	cd->test_func_dim = 0;
	cd->energy = 0;
	cd->ireal = 0;
	cd->sindx = 0;
	cd->lindx = 0;
	cd->pardx = 0;
	cd->pardomain = 100;
	cd->lm_factor = 10.0;
	cd->lm_acc = 0;
	cd->lm_mu = 0.1;
	cd->lm_nu = 2; // default value
	cd->lm_h = 0.1;
	cd->lm_ratio = 0.75 * 0.75;
	cd->lm_ofdecline = 2; // TODO Leif prefers to be 1 and works better for Cr problem; Original code is 2 and works better for test cases
	cd->lm_error = 1e-5;
	cd->lm_indir = 1;
	cd->lm_njacof = 5;
	cd->lm_nlamof = 3;
	cd->squads = 0;
	cd->test_func_npar = cd->test_func_nobs = 0;
	cd->obs_int = 2;
	quiet = 0;
	for( word = strtok( buf, sep ); word; word = strtok( NULL, sep ) )
	{
		w = 0;
		if( !strncasecmp( word, "yaml", 5 ) ) { w = 1; cd->yaml = 1; }
		if( !strncasecmp( word, "quiet", 5 ) ) { w = 1; quiet = 1; }
		if( !strncasecmp( word, "check", 5 ) ) { w = 1; cd->problem_type = CHECK; }
		if( !strncasecmp( word, "create", 6 ) ) { w = 1; cd->problem_type = CREATE; }
		if( !strncasecmp( word, "forward", 7 ) ) { w = 1; cd->problem_type = FORWARD; }
		if( !strncasecmp( word, "calib", 5 ) ) { w = 1; cd->problem_type = CALIBRATE; }
		if( !strncasecmp( word, "lsens", 5 ) ) { w = 1; if( cd->problem_type == CALIBRATE ) cd->lm_eigen = 1; else cd->problem_type = LOCALSENS; }
		if( !strncasecmp( word, "eigen", 5 ) ) { w = 1; if( cd->problem_type == CALIBRATE ) cd->lm_eigen = 1; else cd->problem_type = EIGEN; }
		if( !strncasecmp( word, "lmeigen", 7 ) ) { w = 1; cd->problem_type = CALIBRATE; sscanf( word, "lmeigen=%d", &cd->lm_eigen ); if( cd->lm_eigen == 0 ) cd->lm_eigen = 1; }
		if( !strncasecmp( word, "monte", 5 ) ) { w = 1; cd->problem_type = MONTECARLO; }
		if( !strncasecmp( word, "gsens", 5 ) ) { w = 1; cd->problem_type = GLOBALSENS; cd->gsa_type = SOBOL; }
		if( !strncasecmp( word, "sobol", 5 ) ) { w = 1; cd->problem_type = GLOBALSENS; cd->gsa_type = SOBOL; }
		if( !strncasecmp( word, "salt", 4 ) ) { w = 1; cd->problem_type = GLOBALSENS; cd->gsa_type = SALTELLI; }
		if( !strncasecmp( word, "moat", 4 ) ) { w = 1; cd->problem_type = GLOBALSENS; cd->gsa_type = MOAT; }
		if( !strncasecmp( word, "glue", 4 ) ) { w = 1; cd->problem_type = GLUE; }
		if( !strncasecmp( word, "abagus", 6 ) ) { w = 1; cd->problem_type = ABAGUS; }
		if( !strncasecmp( word, "infogap", 7 ) ) { w = 1; cd->problem_type = INFOGAP; }
		if( !strncasecmp( word, "postpua", 7 ) ) { w = 1; cd->problem_type = POSTPUA; }
		if( !strncasecmp( word, "single", 6 ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = SIMPLE; }
		if( !strncasecmp( word, "simple", 6 ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = SIMPLE; }
		if( !strncasecmp( word, "igpd", 4 ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = IGPD; }
		if( !strncasecmp( word, "ppsd", 4 ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = PPSD; }
		if( !strncasecmp( word, "igrnd", 5 ) ) { w = 1; cd->problem_type = CALIBRATE; cd->calib_type = IGRND; }
		if( !strncasecmp( word, "energy=", 7 ) ) { w = 1; sscanf( word, "energy=%d", &cd->energy ); }
		if( !strncasecmp( word, "background=", 11 ) ) { w = 1; sscanf( word, "background=%lf", &cd->c_background ); }
		if( !strncasecmp( word, "lmfactor=", 9 ) ) { w = 1; sscanf( word, "lmfactor=%lf", &cd->lm_factor ); }
		if( !strncasecmp( word, "lmacc", 5 ) ) { w = 1; cd->lm_acc = 1; }
		if( !strncasecmp( word, "lmratio=", 8 ) ) { w = 1; sscanf( word, "lmratio=%lf", &cd->lm_ratio ); }
		if( !strncasecmp( word, "lmofdecline=", 12 ) ) { w = 1; sscanf( word, "lmofdecline=%lf", &cd->lm_ofdecline ); }
		if( !strncasecmp( word, "lmerror=", 8 ) ) { w = 1; sscanf( word, "lmerror=%lf", &cd->lm_error ); }
		if( !strncasecmp( word, "lmh=", 4 ) ) { w = 1; sscanf( word, "lmh=%lf", &cd->lm_h ); cd->lm_acc = 1; }
		if( !strncasecmp( word, "lmind", 5 ) ) { w = 1; cd->lm_indir = 1; cd->lm_ofdecline = 2; }
		if( !strncasecmp( word, "lmdir", 5 ) ) { w = 1; cd->lm_indir = 0; cd->lm_ofdecline = 1; }
		if( !strncasecmp( word, "lmmu=", 5 ) ) { w = 1; sscanf( word, "lmmu=%lf", &cd->lm_mu ); }
		if( !strncasecmp( word, "lmnu=", 5 ) ) { w = 1; sscanf( word, "lmnu=%d", &cd->lm_nu ); }
		if( !strncasecmp( word, "lmiter=", 7 ) ) { w = 1; sscanf( word, "lmiter=%d", &cd->niter ); }
		if( !strncasecmp( word, "lmnlamof=", 9 ) ) { w = 1; sscanf( word, "lmnlamof=%d", &cd->lm_nlamof ); }
		if( !strncasecmp( word, "lmnjacof=", 9 ) ) { w = 1; sscanf( word, "lmnjacof=%d", &cd->lm_njacof ); }
		if( !strncasecmp( word, "infile=", 7 ) ) { w = 1; sscanf( word, "infile=%s", cd->infile ); }
		if( !strncasecmp( word, "tied", 4 ) ) { w = 1; cd->disp_tied = 1; } // Tied dispersivities
		if( !strncasecmp( word, "disp_tied", 9 ) ) { w = 1; cd->disp_tied = 1; } // Tied dispersivities
		if( !strncasecmp( word, "disp_scaled", 10 ) ) { w = 1; if( sscanf( word, "disp_scaled=%d", &cd->disp_scaled ) != 1 ) cd->disp_scaled = 1; } // Scaled dispersivities
		if( !strncasecmp( word, "save", 4 ) ) { w = 1; cd->save = 1; }
		if( !strncasecmp( word, "pargen", 6 ) ) { w = 1; cd->pargen = 1; }
		if( !strncasecmp( word, "real=", 5 ) ) { w = 1; sscanf( word, "real=%d", &cd->nreal ); }
		if( !strncasecmp( word, "eval=", 5 ) ) { w = 1; sscanf( word, "eval=%d", &cd->maxeval ); }
		if( !strncasecmp( word, "case=", 5 ) ) { w = 1; sscanf( word, "case=%d", &cd->ireal ); }
		if( !strncasecmp( word, "retry", 5 ) ) { w = 1; sscanf( word, "retry=%d", &cd->nretries ); }
		if( !strncasecmp( word, "particles", 9 ) ) { w = 1; if( sscanf( word, "particles=%d", &cd->init_particles ) != 1 ) cd->init_particles = -1; }
		if( !strncasecmp( word, "opt=", 4 ) ) { w = 1; if( cd->problem_type == UNKNOWN ) cd->problem_type = CALIBRATE; sscanf( word, "opt=%s", cd->opt_method ); }
		if( !strncasecmp( word, "smp=", 4 ) ) { w = 1; sscanf( word, "smp=%s", cd->smp_method ); }
		if( !strncasecmp( word, "mslm", 4 ) ) { w = 1; cd->paranoid = 1; sscanf( word, "mslm=%s", cd->paran_method ); }
		if( !strncasecmp( word, "paran=", 6 ) ) { w = 1; cd->paranoid = 1; sscanf( word, "paran=%s", cd->paran_method ); } // legacy
		if( !strncasecmp( word, "nosin", 5 ) ) { w = 1; cd->sintrans = 0; }
		if( !strncasecmp( word, "plog", 4 ) ) { w = 1; if( sscanf( word, "plog=%d", &cd->plogtrans ) != 1 ) cd->plogtrans = 1; }
		if( !strncasecmp( word, "olog", 4 ) ) { w = 1; if( sscanf( word, "olog=%d", &cd->ologtrans ) != 1 ) cd->ologtrans = 1; }
		if( !strncasecmp( word, "oweight", 7 ) ) { w = 1; if( sscanf( word, "oweight=%d", &cd->oweight ) != 1 ) cd->oweight = 1; }
		if( !strncasecmp( word, "cutoff=", 7 ) ) { w = 1; sscanf( word, "cutoff=%lf", &cd->phi_cutoff ); cd->check_success = cd->obsrange = 0; cd->obserror = cd->parerror = ( double ) 0; }
		if( !strncasecmp( word, "obsrange", 8 ) ) { w = 1; cd->check_success = 1; cd->obsrange = 1; cd->phi_cutoff = cd->obserror = cd->parerror = ( double ) 0; }
		if( !strncasecmp( word, "success", 7 ) ) { w = 1; cd->check_success = 1; cd->obsrange = 1; cd->phi_cutoff = cd->obserror = cd->parerror = ( double ) 0; } // legacy
		if( !strncasecmp( word, "truth", 5 ) ) { w = 1; sscanf( word, "truth=%lf", &cd->parerror ); cd->check_success = 1; cd->obsrange = 0; cd->phi_cutoff = cd->obserror = ( double ) 0; if( cd->parerror < DBL_EPSILON ) cd->parerror = 0.1; } // legacy
		if( !strncasecmp( word, "sindx=", 5 ) ) { w = 1; cd->sintrans = 1; sscanf( word, "sindx=%lf", &cd->sindx ); if( cd->sindx < DBL_EPSILON ) cd->sindx = 0.0000001; }
		if( !strncasecmp( word, "lindx=", 5 ) ) { w = 1; cd->sintrans = 0; sscanf( word, "lindx=%lf", &cd->lindx ); if( cd->lindx < DBL_EPSILON ) cd->lindx = 0.001; }
		if( !strncasecmp( word, "pardx", 5 ) ) { w = 1; cd->sintrans = 0; sscanf( word, "pardx=%lf", &cd->pardx ); if( cd->pardx < DBL_EPSILON ) cd->pardx = 0.1; }
		if( !strncasecmp( word, "parerror", 8 ) ) { w = 1; sscanf( word, "parerror=%lf", &cd->parerror ); cd->check_success = 1; cd->obsrange = 0; cd->phi_cutoff = cd->obserror = ( double ) 0; if( cd->parerror < DBL_EPSILON ) cd->parerror = 0.1; }
		if( !strncasecmp( word, "obserror", 8 ) ) { w = 1; sscanf( word, "obserror=%lf", &cd->obserror ); cd->check_success = 1; cd->obsrange = 0; cd->phi_cutoff =  cd->parerror = ( double ) 0; if( cd->obserror < DBL_EPSILON ) cd->obserror = 0.1; }
		if( !strncasecmp( word, "pardomain=", 10 ) ) { w = 1; cd->sintrans = 0; sscanf( word, "pardomain=%lf", &cd->pardomain ); if( cd->pardomain < DBL_EPSILON ) cd->pardomain = ( double ) 100; }
		if( !strncasecmp( word, "obsdomain=", 10 ) ) { w = 1; sscanf( word, "obsdomain=%lf", &cd->obsdomain ); if( cd->obsdomain < DBL_EPSILON ) cd->pardomain = ( double ) 0; }
		if( !strncasecmp( word, "obsstep=", 8 ) ) { w = 1; sscanf( word, "obsstep=%lf", &cd->obsstep ); if( fabs( cd->obsstep ) < DBL_EPSILON ) cd->obsstep = ( double ) 0; else cd->problem_type = INFOGAP; }
		if( !strncasecmp( word, "seed=", 5 ) ) { w = 1; sscanf( word, "seed=%d", &cd->seed ); cd->seed_init = cd->seed; }
		if( !strncasecmp( word, "np", 2 ) ) { w = 1; cd->num_proc = 0; sscanf( word, "np=%d", &cd->num_proc ); if( cd->num_proc <= 0 ) cd->num_proc = 0; }
		if( !strncasecmp( word, "restart", 7 ) ) { w = 1; sscanf( word, "restart=%d", &cd->restart ); if( cd->restart < 0 || cd->restart > 1 ) cd->restart = -1; }
		if( !strncasecmp( word, "rstfile=", 8 ) ) { w = 1; sscanf( word, "rstfile=%s", cd->restart_zip_file ); cd->restart = -1; }
		if( !strncasecmp( word, "resultsfile=", 12 ) ) { w = 1; sscanf( word, "resultsfile=%s", cd->resultsfile ); cd->problem_type = FORWARD; }
		if( !strncasecmp( word, "resultscase=", 12 ) ) { w = 1; sscanf( word, "resultscase=%d", &cd->resultscase ); }
		if( !strncasecmp( word, "debug", 5 ) ) { w = 1; if( sscanf( word, "debug=%d", &cd->debug ) == 0 || cd->debug == 0 ) cd->debug = 1; } // Global debug
		if( !strncasecmp( word, "fdebug", 6 ) ) { w = 1; sscanf( word, "fdebug=%d", &cd->fdebug ); if( cd->fdebug == 0 ) cd->fdebug = 1; }
		if( !strncasecmp( word, "ldebug", 6 ) ) { w = 1; sscanf( word, "ldebug=%d", &cd->ldebug ); if( cd->ldebug == 0 ) cd->ldebug = 1; }
		if( !strncasecmp( word, "pdebug", 6 ) ) { w = 1; sscanf( word, "pdebug=%d", &cd->pdebug ); if( cd->pdebug == 0 ) cd->pdebug = 1; }
		if( !strncasecmp( word, "mdebug", 6 ) ) { w = 1; sscanf( word, "mdebug=%d", &cd->mdebug ); if( cd->mdebug == 0 ) cd->mdebug = 1; }
		if( !strncasecmp( word, "odebug", 6 ) ) { w = 1; sscanf( word, "odebug=%d", &cd->odebug ); if( cd->odebug != 1 ) cd->odebug = 1; }
		if( !strncasecmp( word, "insdebug", 8 ) ) { w = 1; sscanf( word, "insdebug=%d", &cd->insdebug ); if( cd->insdebug == 0 ) cd->insdebug = 1; }
		if( !strncasecmp( word, "tpldebug", 8 ) ) { w = 1; sscanf( word, "tpldebug=%d", &cd->tpldebug ); if( cd->tpldebug == 0 ) cd->tpldebug = 1; }
		if( !strncasecmp( word, "pardebug", 8 ) ) { w = 1; sscanf( word, "pardebug=%d", &cd->pardebug ); if( cd->pardebug == 0 ) cd->pardebug = 1; }
		if( !strncasecmp( word, "ssr", 3 ) ) { w = 1; cd->objfunc_type = SSR; }
		if( !strncasecmp( word, "ssd0", 4 ) ) { w = 1; cd->objfunc_type = SSD0; }
		if( !strncasecmp( word, "ssdx", 4 ) ) { w = 1; cd->objfunc_type = SSDX; }
		if( !strncasecmp( word, "ssda", 4 ) ) { w = 1; cd->objfunc_type = SSDA; }
		if( !strncasecmp( word, "ssdr", 4 ) ) { w = 1; cd->objfunc_type = SSDR; }
		if( !strncasecmp( word, "test", 4 ) ) { w = 1; cd->test_func = 1; cd->test_func_dim = 2; sscanf( word, "test=%d", &cd->test_func ); cd->solution_type[0] = TEST; }
		if( !strncasecmp( word, "dim=", 4 ) ) { w = 1; sscanf( word, "dim=%d", &cd->test_func_dim ); if( cd->test_func_dim < 2 ) cd->test_func_dim = 2; }
		if( !strncasecmp( word, "npar=", 5 ) ) { w = 1; sscanf( word, "npar=%d", &cd->test_func_npar ); }
		if( !strncasecmp( word, "nobs=", 5 ) ) { w = 1; sscanf( word, "nobs=%d", &cd->test_func_nobs ); }
		if( !strncasecmp( word, "poi", 3 ) ) { w = 1; cd->solution_type[0] = POINT; }
		if( !strncasecmp( word, "gau", 3 ) ) { w = 1; if( strcasestr( word, "2" ) ) cd->solution_type[0] = GAUSSIAN2D; else cd->solution_type[0] = GAUSSIAN3D; }
		if( !strncasecmp( word, "rec", 3 ) ) { w = 1; if( strcasestr( word, "ver" ) ) cd->solution_type[0] = PLANE3D; else cd->solution_type[0] = PLANE; }
		if( !strncasecmp( word, "box", 3 ) ) { w = 1; cd->solution_type[0] = BOX; }
		if( !strncasecmp( word, "levy", 4 ) ) { w = 1; cd->levy = 1; }
		if( !strncasecmp( word, "point_tri", 9 ) ) { w = 1; cd->solution_type[0] = POINT_TRIANGLE_TIME; }
		if( !strncasecmp( word, "obs_int=", 8 ) ) { w = 1; sscanf( word, "obs_int=%d", &cd->obs_int ); if( cd->obs_int > 2 || cd->obs_int < 1 ) cd->obs_int = 2; }
		if( !strncasecmp( word, "paran", 5 ) ) { w = 1; cd->paranoid = 1; } // legacy
		if( strcasestr( word, "_ms" ) ) { w = 1; cd->paranoid = 1; } // legacy
		if( w == 0 ) { tprintf( "\nERROR: Unknown keyword \'%s\'!\n", word ); return( -1 ); }
	}
	if( cd->seed != 0 ) cd->seed *= -1; // Modify the seed to show that is imported
	if( cd->seed_init != 0 ) cd->seed_init *= -1; // Modify the seed to show that is imported
	if( cd->problem_type == UNKNOWN ) { cd->problem_type = CALIBRATE; cd->calib_type = SIMPLE; }
	if( cd->nreal == 0 && ( cd->problem_type == MONTECARLO || cd->calib_type == IGRND || cd->problem_type == GLOBALSENS || cd->problem_type == ABAGUS ) ) cd->nreal = 100;
	if( cd->nretries > 0 && strncasecmp( cd->opt_method, "lm", 2 ) == 0 ) cd->paranoid = 1;
	if( cd->test_func > 0 )
	{
		tprintf( "\nTest Function %d ", cd->test_func );
		if( cd->test_func < 40 )
		{
			tprintf( "Dimensionality %d ", cd->test_func );
			if( cd->test_func_nobs > 0 ) tprintf( "Observations %d", cd->test_func_nobs );
		}
		else
		{
			if( cd->test_func_npar > 0 ) tprintf( "Parameters %d ", cd->test_func_npar );
			if( cd->test_func_nobs > 0 ) tprintf( "Observations %d", cd->test_func_nobs );
		}
		tprintf( "\n" );
	}
	tprintf( "\nProblem type: " );
	switch( cd->problem_type )
	{
		case CHECK: tprintf( "check model setup and input/output files (no model execution)" ); break;
		case CREATE: tprintf( "create a calibration input file based on a forward run (no calibration)" ); break;
		case FORWARD: tprintf( "forward run (no calibration)" ); break;
		case CALIBRATE: tprintf( "calibration" ); break;
		case LOCALSENS: tprintf( "sensitivity analysis" ); break;
		case EIGEN: tprintf( "eigen analysis" ); break;
		case MONTECARLO: tprintf( "monte-carlo analysis (realizations = %d)", cd->nreal ); break;
		case GLOBALSENS: tprintf( "global sensitivity analysis (realizations = %d)", cd->nreal ); break;
		case ABAGUS: tprintf( "abagus: agent-based global uncertainty and sensitivity analysis" ); break;
		case GLUE: tprintf( "glue: Generalized Likelihood Uncertainty Estimation: GLUE runs currently postprocess ABAGUS results" ); break;
		case INFOGAP: tprintf( "Info-gap decision analysis" ); break;
		case POSTPUA: tprintf( "predictive uncertainty analysis of sampling results" ); break;
		default: tprintf( "WARNING: unknown problem type; calibration assumed" ); cd->problem_type = CALIBRATE; break;
	}
	tprintf( "\n" );
	if( cd->resultsfile[0] != 0 )
	{
		if( Ftest( cd->resultsfile ) == 0 )
		{
			tprintf( "\nModel analyses based on previously saved results in file %s\n", cd->resultsfile );
			int implemented = 0;
			if( cd->phi_cutoff > DBL_EPSILON ) { tprintf( "Model analyses for cases with phi < %g.", cd->phi_cutoff ); implemented = 1; }
			if( cd->obsrange ) { tprintf( "Model analyses for successful cases." ); implemented = 1; }
			if( !implemented )
			{
				if( cd->resultscase == 0 ) cd->resultscase = 1;
				if( cd->resultscase > 0 ) tprintf( "Model analyses for case #%d", cd->resultscase );
				else tprintf( "Model analyses for first %d cases", cd->resultscase );
			}
			tprintf( "\n" );
		}
		else { tprintf( "\nERROR Results file %s cannot be opened \n", cd->resultsfile ); cd->resultscase = 0; cd->resultsfile[0] = 0; }
	}
	if( cd->problem_type == CALIBRATE )
	{
		tprintf( "\nCalibration technique: " );
		switch( cd->calib_type )
		{
			case IGRND: tprintf( "sequential calibration using a set of random initial values (realizations = %d)", cd->nreal ); break;
			case IGPD: tprintf( "sequential calibration using a set discretized initial values" ); break;
			case PPSD: tprintf( "sequential calibration using partial parameter parameter discretization" ); break;
			case SIMPLE: tprintf( "single calibration using initial guesses provided in the input file" ); break;
			default: tprintf( "WARNING: unknown calibration type!\nASSUMED: single calibration using initial guesses provided in the input file" ); cd->calib_type = SIMPLE; break;
		}
		tprintf( "\n" );
		if( cd->lm_eigen > 0 ) tprintf( "Eigen analysis will be performed for the final optimization results\n" );
	}
	tprintf( "\nOptimization method: opt=%s --- ", cd->opt_method );
	if( strncasecmp( cd->opt_method, "squad", 5 ) == 0 || ( strcasestr( cd->opt_method, "pso" ) && strcasestr( cd->opt_method, "lm" ) ) )
	{
		tprintf( "SQUADS: Coupled Particle-Swarm and Levenberg-Marquardt optimization (Vesselinov & Harp, 2012)\n" );
		cd->squads = 1;
		cd->lmstandalone = 0;
	}
	else if( strncasecmp( cd->opt_method, "lm", 2 ) == 0 )
	{
		if( cd->paranoid ) tprintf( "Multi-Start Levenberg-Marquardt optimization\n" );
		else tprintf( "Levenberg-Marquardt optimization\n" );
		if( cd->calib_type == SIMPLE && cd->nretries <= 1 && !( fabs( cd->obsstep ) > DBL_EPSILON ) ) cd->ldebug = cd->lm_eigen = 1;
	}
	else if( strcasestr( cd->opt_method, "pso" ) || strncasecmp( cd->opt_method, "swarm", 5 ) == 0 || strncasecmp( cd->opt_method, "tribe", 5 ) == 0 )
	{
		tprintf( "Particle-Swarm optimization" );
		if( strcasestr( cd->opt_method, "apso" ) || strncasecmp( cd->opt_method, "tribe", 5 ) == 0 ) tprintf( " TRIBES-D (Clerc 2004; http://clerc.maurice.free.fr/pso)\n" );
		else tprintf( " Standard2006 (http://clerc.maurice.free.fr/pso)\n" );
	}
	else { tprintf( "WARNING: Unknown method (opt=%s)! Levenberg-Marquardt optimization assumed\n", cd->opt_method ); strcpy( cd->opt_method, "lm" ); }
	if( cd->nretries > 0 ) tprintf( "Number of calibration retries = %d\n", cd->nretries );
	if( cd->niter < 0 ) cd->niter = 0;
	if( cd->niter > 0 ) tprintf( "Number of Levenberg-Marquardt iterations = %d\n", cd->niter );
	else tprintf( "Number of Levenberg-Marquardt iterations = will be computed internally\n" );
	if( strcasestr( cd->opt_method, "apso" ) || strcasestr( cd->opt_method, "tribe" ) || strcasestr( cd->opt_method, "squad" ) )
	{
		if( cd->init_particles > 1 ) tprintf( "Number of particles = %d\n", cd->init_particles );
		if( cd->init_particles == -1 ) tprintf( "Number of particles = will be computed internally\n" );
	}
	tprintf( "\nGlobal termination criteria:\n" );
	tprintf( "1: Maximum number of evaluations = %d\n", cd->maxeval );
	tprintf( "2: Objective function cutoff value: " );
	if( cd->phi_cutoff <= DBL_EPSILON ) tprintf( "NOT implemented (ADD keyword cutoff=[value] to implement)\n" );
	else tprintf( "%g\n", cd->phi_cutoff );
	tprintf( "3: Observations within predefined calibration ranges or an absolute observation error: " );
	if( cd->obsrange ) tprintf( "implemented using calibration ranges (keyword 'obsrange')\n" );
	else if( cd->obserror > 0 ) tprintf( "implemented using a predefined absolute error (keyword 'obserror=%g')\n", cd->obserror );
	else tprintf( "NOT implemented (ADD keyword 'obsrange' or 'obserror' to implement)\n" );
	tprintf( "4: Parameters within a predefined absolute error from known 'true' values: " );
	if( cd->parerror > 0 ) tprintf( "implemented (keyword 'parerror=%g')\n", cd->parerror );
	else tprintf( "NOT implemented (ADD keyword 'parerror' to implement)\n" );
	tprintf( "Objective function: " );
	if( cd->test_func > 0 )
		tprintf( "test function %d", cd->test_func );
	else
	{
		switch( cd->objfunc_type )
		{
			case SSR: tprintf( "sum of squared residuals" ); break;
			case SSDR: tprintf( "sum of squared discrepancies and squared residuals" ); break;
			case SSDA: tprintf( "sum of squared discrepancies and residuals" ); break;
			case SSD0: tprintf( "sum of squared discrepancies" ); break;
			case SSDX: tprintf( "sum of squared discrepancies increased to get within the bounds" ); break;
			default: tprintf( "unknown value; sum of squared residuals assumed" ); cd->objfunc_type = SSR; break;
		}
	}
	tprintf( "\n" );
	if( cd->sintrans == 0 ) tprintf( "\nSin transformation of the model parameters: NOT applied (keyword 'nosin')!\n" );
	else tprintf( "\nSin transformation of the model parameters: applied (ADD keyword 'nosin' to remove)\n" );
	if( cd->plogtrans == 1 ) tprintf( "\nLog transformation enforced on all parameters!\n" );
	else if( cd->plogtrans == 0 ) tprintf( "\nLog transformation is not applied to any parameters!\n" );
	if( cd->ologtrans == 1 ) tprintf( "\nLog transformation enforced on all calibration targets!\n" );
	else if( cd->ologtrans == 0 ) tprintf( "\nLog transformation is not applied to any calibration targets!\n" );
	if( cd->oweight == 1 ) tprintf( "\nUnit residual weights for all calibration targets!\n" );
	else if( cd->oweight == 0 ) tprintf( "\nWARNING: Zero residual weights for all calibration targets (for testing and manipulation purposes)!\n" );
	else if( cd->oweight == 2 ) tprintf( "\nResidual weights reversely proportional to observations for all calibration targets!\n" );
	if( cd->problem_type == ABAGUS )
	{
		if( cd->energy == 0 ) { cd->energy = 10000; tprintf( "\nInitial particle energy set to default value: %d\n", cd->energy );}
		else tprintf( "\nInitial particle energy set to: %d\n", cd->energy );
		cd->sintrans = 0; tprintf( "\nsine transformation disabled for ABAGUS runs" );
	}
	if( cd->problem_type == ABAGUS && cd->infile[0] != 0 ) { tprintf( "\nResults in %s to be read into kdtree\n", cd->infile );}
	if( cd->problem_type == INFOGAP && ( fabs( cd->obsstep ) < DBL_EPSILON ) && cd->infile[0] == 0 ) { tprintf( "\nInfile must be specified for infogap run\n" ); return( 0 ); }
	if( cd->problem_type == POSTPUA && cd->infile[0] == 0 ) { tprintf( "\nInfile must be specified for postpua run\n" ); return( 0 ); }
	if( cd->smp_method[0] != 0 )
	{
		tprintf( "\nSampling method: " );
		if( strncasecmp( cd->smp_method, "olhs", 4 ) == 0 ) tprintf( "Optimal Latin Hyper Cube (LHS) (if real <= 500 IDLHS; if real > 500 LHS)\n" );
		else if( strncasecmp( cd->smp_method, "lhs", 3 ) == 0 ) tprintf( "Standard Latin Hyper Cube (LHS)\n" );
		else if( strncasecmp( cd->smp_method, "idlhs", 5 ) == 0 ) tprintf( "Improved Distance Latin Hyper Cube (IDLHS)\n" );
		else if( strncasecmp( cd->smp_method, "random", 5 ) == 0 ) tprintf( "Pure random\n" );
		else { tprintf( "WARNING: Unknown (rnd=%s); Optimal Latin Hyper Cube selected (if real <= 500 IDLHS; if real > 500 LHS)\n", cd->smp_method ); strcpy( cd->smp_method, "olhs" ); }
	}
	if( cd->paran_method[0] != 0 )
	{
		tprintf( "\nParanoid Sampling method: " );
		if( strncasecmp( cd->paran_method, "olhs", 4 ) == 0 ) tprintf( "Optimal Latin Hyper Cube (LHS) (if real <= 500 IDLHS; if real > 500 LHS)\n" );
		else if( strncasecmp( cd->paran_method, "lhs", 3 ) == 0 ) tprintf( "Standard Latin Hyper Cube (LHS)\n" );
		else if( strncasecmp( cd->paran_method, "idlhs", 5 ) == 0 ) tprintf( "Improved Distance Latin Hyper Cube (IDLHS)\n" );
		else if( strncasecmp( cd->paran_method, "random", 5 ) == 0 ) tprintf( "Pure random\n" );
		else { tprintf( "WARNING: Unknown (rnd=%s); Optimal Latin Hyper Cube selected (if real <= 500 IDLHS; if real > 500 LHS)\n", cd->paran_method ); strcpy( cd->paran_method, "olhs" ); }
	}
	tprintf( "\nGlobal debug (verbosity) level: debug=%d\n", cd->debug );
	if( cd->debug || cd->fdebug ) tprintf( "Debug (verbosity) level for the analytical model evaluations: fdebug=%d\n", cd->fdebug );
	if( cd->debug )
	{
		tprintf( "Debug (verbosity) level for Levenberg-Marquardt optimization progress: ldebug=%d\n", cd->ldebug );
		tprintf( "Debug (verbosity) level for Particle-Swarm optimization progress: pdebug=%d\n", cd->pdebug );
		tprintf( "Debug (verbosity) level for objective function progress: odebug=%d\n", cd->odebug );
	}
	if( ( cd->debug || cd->mdebug ) && cd->problem_type != CREATE && cd->problem_type != EIGEN )
		tprintf( "Debug (verbosity) level for random sets: mdebug=%d\n", cd->mdebug );
	if( ( cd->debug || cd->pardebug ) && cd->num_proc > 1 )
		tprintf( "Debug (verbosity) level for parallel execution: pardebug=%d\n", cd->pardebug );
	if( ( cd->debug || cd->tpldebug || cd->insdebug ) && cd->solution_type[0] == EXTERNAL )
	{
		tprintf( "Debug (verbosity) level for template file: tpldebug=%d\n", cd->tpldebug );
		tprintf( "Debug (verbosity) level for instruction file: insdebug=%d\n", cd->insdebug );
	}
	tprintf( "\n" );
	return( 1 );
}

int load_problem( char *filename, int argn, char *argv[], struct opt_data *op )
{
	FILE *infile, *infile2;
	//	FILE *infileb;
	double d;
	char buf[5000], *file, **path, exec[1000], *word, *start, charecter;
	char *separator = " \t\n";
	int  i, j, k, c, l1, l2, bad_data, status, nofile = 0, skip = 0, short_names_printed;
	struct calc_data *cd;
	struct param_data *pd;
	struct regul_data *rd;
	struct obs_data *od;
	struct obs_data *preds;
	struct well_data *wd;
	struct grid_data *gd;
	struct extrn_data *ed;
	char **expvar_names;
	int expvar_count;
	cd = op->cd;
	pd = op->pd;
	rd = op->rd;
	od = op->od;
	preds = op->preds;
	wd = op->wd;
	gd = op->gd;
	ed = op->ed;
	pd->nParam = pd->nFixParam = pd->nFlgParam = pd->nOptParam = pd->nExpParam = 0;
	od->nObs = od->nTObs = od->nCObs = 0;
	gd->min_t = gd->time = 0;
	// IMPORTANT
	// internal problem: nCObs = nObs
	// external problem: nTObs = nObs
	// internal test problem: nTObs = nCObs = nObs
	wd->nW = 0;
	ed->ntpl = ed->nins = 0;
	gd->min_t = 0;
	bad_data = 0;
	if( ( infile = fopen( filename, "r" ) ) == NULL )
	{
		sprintf( filename, "%s.in", op->root );
		if( ( infile = fopen( filename, "r" ) ) == NULL )
			nofile = 1;
		else
			tprintf( "WARNING: File \'%s\' is opened to read problem information!\n", filename );
	}
	if( nofile == 0 )
	{
		// Read commands in the file
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ":" ); fscanf( infile, "%1000[^\n]s\n", buf ); fscanf( infile, "\n" );
		if( sscanf( buf, "%i", &pd->nParam ) == 1 ) { buf[0] = 0; skip = 1; }
	}
	else buf[0] = 0;
	// Add commands provided as arguments
	for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
	cd->solution_type = ( int * ) malloc( sizeof( int ) );
	if( parse_cmd( buf, cd ) == -1 ) return( -1 );
	od->include_predictions = 1;
	if( cd->problem_type == INFOGAP ) od->include_predictions = 0;
	if( fabs( cd->obsstep ) > DBL_EPSILON ) od->include_predictions = 1;
	// Read Solution Type
	cd->solution_id = ( char * ) malloc( 150 * sizeof( char ) ); // Needed only to save text MADS files
	cd->solution_id[0] = 0;
	if( nofile == 0 && skip == 0 ) { fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ":" ); fgets( cd->solution_id, 150, infile ); /*fscanf( infile, "%s\n", cd->solution_id );*/ }
	strcpy( buf, cd->solution_id );
	for( c = 0, word = strtok( buf, separator ); word; c++, word = strtok( NULL, separator ) )
		if( cd->debug > 1 ) tprintf( "Model #%d %s\n", c + 1, word );
	cd->num_sources = c;
	if( cd->num_sources > 1 )
	{
		tprintf( "Number of analytical solutions: %d\n", cd->num_sources );
		free( cd->solution_type );
		cd->solution_type = ( int * ) malloc( cd->num_sources * sizeof( int ) );
	}
	strcpy( buf, cd->solution_id );
	for( c = 0, word = strtok( buf, separator ); word; c++, word = strtok( NULL, separator ) )
	{
		sscanf( word, "%d", &cd->solution_type[c] );
		if( !strncasecmp( word, "ext", 3 ) ) { cd->solution_type[c] = EXTERNAL; if( cd->num_sources > 1 ) { tprintf( "ERROR: Multiple solutions can be only internal; no external!\n" ); bad_data = 1; } }
		if( !strncasecmp( word, "poi", 3 ) ) cd->solution_type[c] = POINT;
		if( !strncasecmp( word, "gau", 3 ) ) { if( strcasestr( word, "2" ) ) cd->solution_type[c] = GAUSSIAN2D; else cd->solution_type[c] = GAUSSIAN3D; }
		if( !strncasecmp( word, "rec", 3 ) ) { if( strcasestr( word, "ver" ) ) cd->solution_type[c] = PLANE3D; else cd->solution_type[c] = PLANE; }
		if( !strncasecmp( word, "box", 3 ) ) cd->solution_type[c] = BOX;
		if( !strncasecmp( word, "point_tri", 9 ) ) cd->solution_type[c] = POINT_TRIANGLE_TIME;
		if( !strncasecmp( word, "test", 4 ) || cd->test_func >= 0 ) { cd->solution_type[c] = TEST; od->nTObs = 0; if( cd->num_sources > 1 ) { tprintf( "ERROR: Multiple solutions can be only internal; no test functions!\n" ); bad_data = 1; } }
	}
	if( cd->num_sources == 0 && cd->test_func >= 0 ) { cd->num_sources = 1; cd->solution_type[0] = TEST; od->nTObs = 0; if( cd->num_sources > 1 ) { tprintf( "ERROR: Multiple solutions can be only internal; no test functions!\n" ); bad_data = 1; } }
	if( bad_data ) return( -1 );
	if( nofile )
	{
		if( cd->solution_type[0] != TEST || ( cd->problem_type == INFOGAP && ( fabs( cd->obsstep ) < DBL_EPSILON ) ) )
		{
			tprintf( "File \'%s\' cannot be opened to read problem information!\n", filename );
			tprintf( "ERROR: Input file is needed!\n\n" );
			bad_data = 1;
			return( -1 );
		}
	}
	cd->solution_id[0] = 0; // Regenerate the solution ID to save in the text MADS output file
	if( cd->num_sources > 1 ) tprintf( "\nModels:" );
	else tprintf( "Model: " );
	for( c = 0; c < cd->num_sources; c++ )
	{
		if( cd->num_sources > 1 ) tprintf( " (%d) ", c + 1 );
		switch( cd->solution_type[c] )
		{
			case EXTERNAL: { tprintf( "external" ); strcat( cd->solution_id, "external" ); break; }
			case POINT: { tprintf( "internal point contaminant source" ); strcat( cd->solution_id, "point" ); break; }
			case PLANE: { tprintf( "internal rectangular contaminant source" ); strcat( cd->solution_id, "rect" ); break; }
			case GAUSSIAN2D: { tprintf( "internal planar (2d) gaussian contaminant source" ); strcat( cd->solution_id, "gaussian_2d" ); break; }
			case GAUSSIAN3D: { tprintf( "internal spatial (3d) gaussian contaminant source" ); strcat( cd->solution_id, "gaussian_3d" ); break; }
			case PLANE3D: { tprintf( "internal rectangular contaminant source with vertical flow component" ); strcat( cd->solution_id, "rect_vert" ); break; }
			case BOX: { tprintf( "internal box contaminant source" ); strcat( cd->solution_id, "box" ); break; }
			case POINT_TRIANGLE_TIME: { tprintf( "internal point contaminant source with triangle shape in time" ); strcat( cd->solution_id, "point_tri" ); break; }
			case TEST: { tprintf( "internal test optimization problem #%d: ", cd->test_func ); set_test_problems( op ); sprintf( cd->solution_id, "test=%d", cd->test_func ); break; }
			default: tprintf( "WARNING! UNDEFINED model type!" ); break;
		}
		if( cd->num_sources > 1 ) { strcat( cd->solution_id, " " ); tprintf( ";" ); }
	}
	if( cd->c_background ) tprintf( " | background concentration = %g", cd->c_background );
	tprintf( "\n" );
	//
	if( cd->solution_type[0] == TEST ) return( 1 ); // DONE IF INTERNAL TEST PROBLEM
	//
	if( cd->solution_type[0] != EXTERNAL )
	{
		if( cd->disp_scaled ) tprintf( "Longitudinal dispersivity is scaled!\n" );
		if( cd->disp_scaled > 1 && !cd->disp_tied ) tprintf( "Transverse dispersivities are scaled!\n" );
		else if( cd->disp_tied ) tprintf( "Transverse dispersivities are tied!\n" );
		else tprintf( "Transverse dispersivities are neither tied or scaled!\n" );
	}
	rd->nRegul = 0;
	// ------------------------------------------------------------ Reading parameters ----------------------------------------------------------------
	if( skip == 0 ) fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %i\n", &pd->nParam );
	tprintf( "\nNumber of model parameters: %d\n", pd->nParam );
	if( cd->solution_type[0] != TEST && cd->solution_type[0] != EXTERNAL )
	{
		set_param_id( op ); // set analytical parameter id's
		pd->nAnalParam = cd->num_source_params * cd->num_sources + cd->num_aquifer_params;
		if( pd->nAnalParam != pd->nParam )
		{
			tprintf( "WARNING: Internal analytical solver expects %d parameters (%d != %d)!\n", pd->nAnalParam , pd->nAnalParam , pd->nParam );
			// bad_data = 1; TODO revisit this; currently the code does not check for consistency
			// return( -1 );
		}
		pd->var_id = char_matrix( pd->nAnalParam , 10 );
		pd->var_name = char_matrix( pd->nAnalParam , 50 );
		for( i = 0; i < pd->nAnalParam ; i++ )
			pd->var_id[i][0] = pd->var_name[i][0] = 0;
		pd->var = ( double * ) malloc( pd->nAnalParam  * sizeof( double ) );
		cd->var = ( double * ) malloc( pd->nAnalParam  * sizeof( double ) );
		pd->var_opt = ( int * ) malloc( pd->nAnalParam  * sizeof( int ) );
		pd->var_log = ( int * ) malloc( pd->nAnalParam  * sizeof( int ) );
		pd->var_dx = ( double * ) malloc( pd->nAnalParam  * sizeof( double ) );
		pd->var_min = ( double * ) malloc( pd->nAnalParam  * sizeof( double ) );
		pd->var_max = ( double * ) malloc( pd->nAnalParam  * sizeof( double ) );
		pd->var_range = ( double * ) malloc( pd->nAnalParam  * sizeof( double ) );
		pd->param_expressions_index = ( int * ) malloc( pd->nAnalParam  * sizeof( int ) );
		pd->param_expression = ( void ** ) malloc( pd->nAnalParam  * sizeof( void * ) );
		init_params( op );
	}
	else
	{
		pd->nAnalParam = pd->nParam; // just in case
		pd->var_name = char_matrix( pd->nParam, 50 );
		pd->var = ( double * ) malloc( pd->nParam * sizeof( double ) );
		cd->var = ( double * ) malloc( pd->nParam * sizeof( double ) );
		pd->var_opt = ( int * ) malloc( pd->nParam * sizeof( int ) );
		pd->var_log = ( int * ) malloc( pd->nParam * sizeof( int ) );
		pd->var_dx = ( double * ) malloc( pd->nParam * sizeof( double ) );
		pd->var_min = ( double * ) malloc( pd->nParam * sizeof( double ) );
		pd->var_max = ( double * ) malloc( pd->nParam * sizeof( double ) );
		pd->var_range = ( double * ) malloc( pd->nParam * sizeof( double ) );
		pd->param_expressions_index = ( int * ) malloc( pd->nParam * sizeof( int ) );
		pd->param_expression = ( void ** ) malloc( pd->nParam * sizeof( void * ) );
	}
	pd->nOptParam = pd->nFlgParam = 0;
	for( i = 0; i < pd->nParam; i++ )
	{
		pd->var[i] = 0;
		fscanf( infile, "%50[^:=]s", pd->var_name[i] );
		fscanf( infile, "%c", &charecter );
		white_skip( &pd->var_name[i] );
		white_trim( pd->var_name[i] );
		if( charecter == ':' ) // regular parameter
		{
			if( cd->debug ) tprintf( "%-27s: ", pd->var_name[i] );
			if( fscanf( infile, "%lf %d %d %lf %lf %lf\n", &pd->var[i], &pd->var_opt[i], &pd->var_log[i], &pd->var_dx[i], &pd->var_min[i], &pd->var_max[i] ) != 6 )
			{
				tprintf( "ERROR: Specific parameter values expected for parameter \"%s\":\n", pd->var_name[i] );
				tprintf( "       initial value (float), optimization flag (int), log-transformation flag (int), dx (float), min (float), max (float)\n" );
				bad_data = 1;
				return( -1 );
			}
			cd->var[i] = pd->var[i];
			if( cd->debug ) tprintf( "init %9g opt %1d log %1d step %7g min %9g max %9g\n", pd->var[i], pd->var_opt[i], pd->var_log[i], pd->var_dx[i], pd->var_min[i], pd->var_max[i] );
			if( pd->var_opt[i] == 1 ) pd->nOptParam++;
			else if( pd->var_opt[i] == 2 )
			{
				pd->nFlgParam++;
				if( cd->calib_type != PPSD ) pd->nOptParam++;
			}
			if( pd->var_opt[i] >= 1 )
			{
				if( pd->var_max[i] < pd->var[i] || pd->var_min[i] > pd->var[i] )
				{
					tprintf( "ERROR: Parameter initial value is outside the specified min/max range! " );
					tprintf( "Parameter %s: %g min %g max %g\n", pd->var_name[i], pd->var[i], pd->var_min[i], pd->var_max[i] );
					bad_data = 1;
				}
				if( pd->var_max[i] < pd->var_min[i] )
				{
					tprintf( "ERROR: Parameter min/max range is not correctly specified! " );
					tprintf( "Parameter %s: min %g max %g\n", pd->var_name[i], pd->var_min[i], pd->var_max[i] );
					bad_data = 1;
				}
				if( cd->plogtrans == 1 ) pd->var_log[i] = 1;
				else if( cd->plogtrans == 0 ) pd->var_log[i] = 0;
				if( pd->var_log[i] == 1 )
				{
					if( pd->var_min[i] < 0 || pd->var[i] < 0 )
					{
						tprintf( "ERROR: Parameter cannot be log transformed (negative values)!\n" );
						tprintf( "Parameter %s: min %g max %g\n", pd->var_name[i], pd->var_min[i], pd->var_max[i] );
						if( cd->plogtrans ) { pd->var_log[i] = 0; pd->var_range[i] = pd->var_max[i] - pd->var_min[i]; continue; }
						else bad_data = 1;
					}
					if( pd->var_dx[i] < 2 ) d = ( pd->var_max[i] - pd->var_min[i] ) / pd->var_dx[i];
					if( pd->var[i] < DBL_EPSILON ) pd->var[i] = DBL_EPSILON;
					if( pd->var_min[i] < DBL_EPSILON ) pd->var_min[i] = DBL_EPSILON;
					pd->var[i] = log10( pd->var[i] );
					pd->var_min[i] = log10( pd->var_min[i] );
					pd->var_max[i] = log10( pd->var_max[i] );
					if( pd->var_dx[i] < 2 ) pd->var_dx[i] = ( pd->var_max[i] - pd->var_min[i] ) / d;
					else pd->var_dx[i] = log10( pd->var_dx[i] );
				}
				pd->var_range[i] = pd->var_max[i] - pd->var_min[i];
				if( pd->var_dx[i] > DBL_EPSILON ) cd->pardx = 1; // discretization is ON
			}
		}
		else if( charecter == '=' )
		{
			if( cd->debug ) tprintf( "%-27s =", pd->var_name[i] );
			pd->var_opt[i] = pd->var_log[i] = 0;
			fscanf( infile, "%1000[^\n]s", buf );
			fscanf( infile, "\n" );
			if( cd->debug ) tprintf( " %s", buf );
#ifdef MATHEVAL
			pd->param_expressions_index[pd->nExpParam] = i;
			pd->param_expression[pd->nExpParam] = evaluator_create( buf );
			assert( pd->param_expression[pd->nExpParam] );
			evaluator_get_variables( pd->param_expression[pd->nExpParam], &expvar_names, &expvar_count );
#else
			expvar_count = 0;
#endif
			if( expvar_count > 0 )
			{
				if( cd->debug )
				{
					tprintf( " -> variables:" );
					for( j = 0; j < expvar_count; j++ )
						tprintf( " %s", expvar_names[j] );
					tprintf( "\n" );
				}
				pd->var_opt[i] = -1;
				pd->nExpParam++;
			}
			else
			{
#ifdef MATHEVAL
				pd->var[i] = cd->var[i] = evaluator_evaluate_x( pd->param_expression[pd->nExpParam], 0 );
				if( cd->debug ) tprintf( " = %g (NO variables; fixed parameter)\n", pd->var[i] );
#else
				tprintf( " MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
			}
		}
	}
	// ------------------------------------------------------------ Reading parameters from previously saved results in a file ----------------------------------------------------------------
	if( cd->resultscase > 0 )
	{
		bad_data = 0;
		tprintf( "\nModel parameters initiated based on previously saved results in file %s (case %d)\n", cd->resultsfile, cd->resultscase );
		infile2 = Fread( cd->resultsfile );
		for( i = 0; i < cd->resultscase; i++ )
		{
			if( fgets( buf, sizeof buf, infile2 ) == NULL )
			{
				bad_data = 1;
				tprintf( "ERROR reading model parameters initiated based on previously saved results in file %s (case %d)\n", cd->resultsfile, cd->resultscase );
				return( 0 );
			}
			// else tprintf( "%s\n", buf );
		}
		fclose( infile2 );
		int caseid;
		sscanf( buf, "%d", &caseid );
		tprintf( "Case ID %d in %s (case %d)\n", caseid, cd->resultsfile, cd->resultscase );
		start = strstr( buf, ": OF" );
		if( start == NULL ) tprintf( "WARNING Objective function value is missing in %s (case %d)\n", cd->resultsfile, cd->resultscase );
		else
		{
			double phi;
			sscanf( start, ": OF %lg", &phi );
			tprintf( "Objective function = %g\n", phi );
		}
		start = strstr( buf, "success" );
		if( start == NULL ) tprintf( "WARNING Success value is missing in %s (case %d)\n", cd->resultsfile, cd->resultscase );
		else
		{
			int success;
			sscanf( start, "success %d", &success );
			tprintf( "Success = %d\n", success );
		}
		start = strstr( buf, "final var" );
		if( start == NULL )
		{
			bad_data = 1;
			tprintf( "ERROR Final model parameters cannot be located in %s (case %d)\n", cd->resultsfile, cd->resultscase );
			return( 0 );
		}
		// tprintf( "%s\n", start );
		strcpy( buf, start );
		tprintf( "\nInitialized model parameters:\n" );
		for( k = 0, i = 0, c = -2, word = strtok( buf, separator ); word; c++, word = strtok( NULL, separator ) )
		{
			if( c > -1 )
			{
				// tprintf( "Par #%d %s\n", c + 1, word );
				while( pd->var_opt[i] == 0 ) i++;
				sscanf( word, "%lf", &pd->var[i] );
				k++;
				if( pd->var_log[i] ) pd->var[i] = log10( pd->var[i] );
				tprintf( "%s %g\n", pd->var_name[i], pd->var[i] );
				i++;
			}
		}
		tprintf( "Number of initialized parameters = %d\n\n", k );
		if( pd->nOptParam != k )
		{
			bad_data = 1;
			tprintf( "ERROR Number of optimized (%d) and initialized (%d) parameters in %s (case %d) do not match\n", pd->nOptParam, k, cd->resultsfile, cd->resultscase );
			return( 0 );
		}
	}
	/*
		if( (*cd).problem_type != 0 )
		{
			// Read binary data if available
			sprintf( buf, "%s.bin", filename );
			if ( ( infileb = fopen( buf, "r" ) ) == NULL )
				tprintf( "Binary file %s cannot be opened to read problem information!\n", buf );
			else
			{
				tprintf( "Binary file %s is found and parameter values are read\n", buf );
				fread( (*pd).var, sizeof((*pd).var[i]), (*pd).nParam, infileb);
				for( i = 0; i < (*pd).nParam; i++ )
						tprintf( "%-27s: binary init %15.12g\n", pd->var_id[i], (*pd).var[i] );
				fclose( infileb );
			}
		}
	 */
	if( cd->solution_type[0] == EXTERNAL ) // check for consistent parameter names
	{
		pd->var_id = ( char ** ) malloc( pd->nParam * sizeof( char * ) );
		for( i = 0; i < pd->nParam; i++ )
		{
			pd->var_id[i] = pd->var_name[i];
			if( strchr( pd->var_name[i], ' ' ) || strchr( pd->var_name[i], '\t' ) ) { tprintf( "ERROR: \'%s\' - invalid parameter name (contains empty space)\n", pd->var_name[i] ); bad_data = 1; }
			l1 = strlen( pd->var_name[i] );
			if( l1 == 0 ) { tprintf( "ERROR: \'%s\' empty parameter name\n", pd->var_name[i] ); bad_data = 1; }
			for( j = i + 1; j < pd->nParam; j++ )
			{
				l2 = strlen( pd->var_name[j] );
				if( l1 == l2 && strcmp( pd->var_name[i], pd->var_name[j] ) == 0 ) { tprintf( "ERROR: %d (\'%s\') and %d (\'%s\') parameter names are the same\n", i + 1, pd->var_name[i], j + 1, pd->var_name[j] ); bad_data = 1; }
			}
		}
	}
	else // set short names for the internal model parameters
	{
		set_param_id( op );
		if( cd->num_sources == 1 )
		{
			for( i = 0; i < cd->num_source_params; i++ )
			{
				sprintf( pd->var_id[i], "%s", op->sd->param_id[i] );
//				tprintf( "%s %s\n", pd->var_id[i],op->sd->param_id[i] );
//				tprintf( "%s %s\n", pd->var_name[i],op->sd->param_name[i] );
				if( !strcasestr( pd->var_name[i], op->sd->param_name[i] ) )
				{
					tprintf( "WARNING: Parameter name \"%s\" did not match expected \"%s\"! Potential input error!\n", pd->var_name[i], op->sd->param_name[i] );
					sprintf( pd->var_name[i], "%s", op->sd->param_name[i] );
				}
			}
		}
		else
		{
			for( j = c = 0; c < cd->num_sources; c++ )
				for( i = 0; i < cd->num_source_params; i++, j++ )
				{
//					tprintf( "%s\n", pd->var_id[j] );
					sprintf( pd->var_id[j], "%s_%d", op->sd->param_id[i], c + 1 );
					if( !strcasestr( pd->var_name[j], op->sd->param_name[i] ) )
					{
						tprintf( "WARNING: Parameter name \"%s\" did not match expected \"%s\"! Potential input error!\n", pd->var_name[j], op->sd->param_name[i] );
						sprintf( pd->var_name[j], "%s", op->sd->param_name[i] );
					}
				}
		}
		j = cd->num_sources * cd->num_source_params;
		for( i = 0; i < cd->num_aquifer_params; i++, j++ )
		{
			sprintf( pd->var_id[j], "%s", op->qd->param_id[i] );
			// tprintf( "%s\n", pd->var_id[j] );
			// tprintf( "%s %s\n", pd->var_name[j],op->qd->param_name[i] );
			if( !strcasestr( pd->var_name[j], op->qd->param_name[i] ) )
			{
				tprintf( "WARNING: Parameter name \"%s\" did not match expected \"%s\"! Potential input error!\n", pd->var_name[j], op->qd->param_name[i] );
				sprintf( pd->var_name[j], "%s", op->qd->param_name[i] );
			}
		}
		if( cd->num_sources > 1 ) k = cd->num_source_params * ( cd->num_sources - 1 );
		else k = 0;
		if( op->cd->debug )
		{
			tprintf( "\nParameter ID's:\n" );
			for( i = 0; i < pd->nParam; i++ )
				tprintf( "%d %s\n", i + 1, pd->var_id[i] );
		}
	}
	if( cd->debug ) tprintf( "\n" );
	set_optimized_params( op );
	// ------------------------------------------------------------ Set parameters with computational expressions (coupled or tied parameters) ----------------------------------------------------------------
	short_names_printed = 0;
	if( pd->nExpParam > 0 )
	{
#ifndef MATHEVAL
		tprintf( "WARNING: MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
		for( i = 0; i < pd->nExpParam; i++ )
		{
#ifdef MATHEVAL
			evaluator_get_variables( pd->param_expression[i], &expvar_names, &expvar_count );
#else
			expvar_count = 0;
#endif
			for( j = 0; j < expvar_count; j++ )
			{
				l1 = strlen( expvar_names[j] );
				status = 0;
				for( k = 0; k < pd->nParam; k++ )
				{
					word = pd->var_id[k];
					l2 = strlen( word );
					if( l1 == l2 && strcmp( expvar_names[j], word ) == 0 ) { status = 1; break; }
				}
#ifdef MATHEVAL
				if( status == 0 ) { tprintf( "ERROR: parameter name \'%s\' in expression \'%s\' for parameter \'%s\' is not defined!\n", expvar_names[j], evaluator_get_string( pd->param_expression[i] ), pd->var_name[pd->param_expressions_index[i]] ); bad_data = 1; }
#endif
			}
		}
		if( bad_data ) return( 0 );
		for( i = 0; i < pd->nParam; i++ )
		{
			if( pd->var_opt[i] > 0 && pd->var_log[i] == 1 ) cd->var[i] = pow( 10, pd->var[i] );
			else cd->var[i] = pd->var[i];
			if( cd->debug )
			{
				if( cd->solution_type[0] == EXTERNAL ) tprintf( "%-27s: %.12g\n", pd->var_id[i], cd->var[i] );
				else tprintf( "%-27s: %-12s: %.12g\n", pd->var_name[i], pd->var_id[i], cd->var[i] );
			}
		}
#ifdef MATHEVAL
		for( i = 0; i < pd->nExpParam; i++ )
		{
			k = pd->param_expressions_index[i];
			pd->var[k] = cd->var[k] = evaluator_evaluate( pd->param_expression[i], pd->nParam, pd->var_id, cd->var );
		}
#endif
		if( cd->debug )
		{
			short_names_printed = 1; // short names printed in the loop above
			for( i = 0; i < pd->nExpParam; i++ )
			{
				k = pd->param_expressions_index[i];
				tprintf( "%-27s= ", pd->var_name[k] );
#ifdef MATHEVAL
				tprintf( "%s", evaluator_get_string( pd->param_expression[i] ) );
				pd->var[k] = cd->var[k] = evaluator_evaluate( pd->param_expression[i], pd->nParam, pd->var_id, cd->var );
				tprintf( " = %g\n", pd->var[k] );
#else
				tprintf( "MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
			}
		}
	}
	if( pd->nParam == 0 || ( pd->nOptParam == 0 && pd->nFlgParam == 0 ) ) { tprintf( "\nERROR: Number of model parameters is zero!\n\n" ); bad_data = 1; }
	if( bad_data ) return( 0 );
	fscanf( infile, "%1000[^:]s", buf );
	// ------------------------------------------------------------ Reading Observations ----------------------------------------------------------------
	// ------------------------------------------------------------ Reading external problem ----------------------------------------------------------------
	if( cd->solution_type[0] == EXTERNAL )
	{
		fscanf( infile, ": %i\n", &od->nObs );
		tprintf( "Number of total observations = %d\n", od->nObs );
		od->nTObs = od->nObs;
		od->obs_id = char_matrix( od->nTObs, 50 );
		od->obs_target = ( double * ) malloc( od->nTObs * sizeof( double ) );
		od->obs_weight = ( double * ) malloc( od->nTObs * sizeof( double ) );
		od->obs_min = ( double * ) malloc( od->nTObs * sizeof( double ) );
		od->obs_max = ( double * ) malloc( od->nTObs * sizeof( double ) );
		od->obs_current = ( double * ) malloc( od->nTObs * sizeof( double ) );
		od->obs_best = ( double * ) malloc( od->nTObs * sizeof( double ) );
		// if( ( od->obs_best = ( double * ) malloc( od->nObs * sizeof( double ) ) ) == NULL ) tprintf( "***\nNO MEMORY!!!!\n***\n" );
		od->res = ( double * ) malloc( od->nTObs * sizeof( double ) );
		od->obs_log = ( int * ) malloc( od->nTObs * sizeof( int ) );
		od->nCObs = 0;
		preds->nTObs = 0; // TODO INFOGAP and GLUE analysis for external problems
		for( i = 0; i < od->nObs; i++ )
		{
			od->obs_min[i] = -1e6; od->obs_max[i] = 1e6; od->obs_weight[i] = 1; od->obs_log[i] = 0;
			fscanf( infile, "%s %lf %lf %d %lf %lf\n", od->obs_id[i], &od->obs_target[i], &od->obs_weight[i], &od->obs_log[i], &od->obs_min[i], &od->obs_max[i] );
			if( cd->obsdomain > DBL_EPSILON && &od->obs_weight[i] > 0 ) { od->obs_min[i] = od->obs_target[i] - cd->obsdomain; od->obs_max[i] = od->obs_target[i] + cd->obsdomain; }
			if( od->obs_max[i] < od->obs_target[i] || od->obs_min[i] > od->obs_target[i] )
			{
				tprintf( "ERROR: Observation target is outside the specified min/max range! " );
				tprintf( "Observation %s: %g min %g max %g\n", od->obs_id[i], od->obs_target[i], od->obs_min[i], od->obs_max[i] );
				bad_data = 1;
			}
			if( od->obs_max[i] <= od->obs_min[i] )
			{
				tprintf( "ERROR: Calibration range is not correctly specified! " );
				tprintf( "Observation %s: min %g max %g\n", od->obs_id[i], od->obs_min[i], od->obs_max[i] );
				bad_data = 1;
			}
			if( cd->ologtrans == 1 ) od->obs_log[i] = 1;
			else if( cd->ologtrans == 0 ) od->obs_log[i] = 0;
			if( cd->oweight == 1 ) od->obs_weight[i] = 1;
			else if( cd->oweight == 0 ) od->obs_weight[i] = 0;
			else if( cd->oweight == 2 ) { if( fabs( od->obs_target[i] ) > DBL_EPSILON ) od->obs_weight[i] = ( double ) 1.0 / od->obs_target[i]; else od->obs_weight[i] = HUGE_VAL; }
			if( od->obs_weight[i] > DBL_EPSILON ) od->nCObs++;
			if( od->obs_weight[i] < -DBL_EPSILON ) { preds->nTObs++; if( od->include_predictions ) od->nCObs++; } // Predictions have negative weights
		}
		tprintf( "Number of calibration targets = %d\n", od->nCObs );
		if( bad_data ) return( 0 );
		if( cd->debug )
		{
			tprintf( "\n" );
			for( i = 0; i < od->nObs; i++ )
			{
				if( cd->debug > 10 || od->nObs <= 50 || ( i < 20 || i > od->nObs - 20 ) )
					tprintf( "%-20s: %15g weight %7g log %1d acceptable range: min %15g max %15g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
				if( ( !( cd->debug > 10 ) || od->nObs > 50 ) && i == 21 ) tprintf( "...\n" );
			}
			for( i = od->nObs; i < od->nTObs; i++ )
				tprintf( "%-20s: %15g weight %7g log %1d acceptable range: min %15g max %15g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
		}
		for( i = 0; i < od->nObs; i++ )
			for( j = i + 1; j < od->nObs; j++ )
				if( strcmp( od->obs_id[i], od->obs_id[j] ) == 0 )
				{
					tprintf( "ERROR: Observation names #%i (%s) and #%i (%s) are identical!\n", i + 1, od->obs_id[i], j + 1, od->obs_id[j] );
					bad_data = 1;
				}
	}
	else
	{
		// ------------------------------------------------------------ Reading internal problem ----------------------------------------------------------------
		fscanf( infile, ": %i\n", &wd->nW );
		tprintf( "Number of wells: %d\n", wd->nW );
		wd->id = char_matrix( wd->nW, 40 );
		if( ( wd->x = ( double * ) malloc( wd->nW * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->y = ( double * ) malloc( wd->nW * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->z1 = ( double * ) malloc( wd->nW * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->z2 = ( double * ) malloc( wd->nW * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->nWellObs = ( int * ) malloc( wd->nW * sizeof( int ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->obs_time = ( double ** ) malloc( wd->nW * sizeof( double * ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->obs_target = ( double ** ) malloc( wd->nW * sizeof( double * ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->obs_weight = ( double ** ) malloc( wd->nW * sizeof( double * ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->obs_log = ( int ** ) malloc( wd->nW * sizeof( int * ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->obs_min = ( double ** ) malloc( wd->nW * sizeof( double * ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		if( ( wd->obs_max = ( double ** ) malloc( wd->nW * sizeof( double * ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
		od->nObs = preds->nTObs = 0;
		if( cd->debug ) tprintf( "\nObservation data:\n" );
		for( i = 0; i < wd->nW; i++ )
		{
			status = fscanf( infile, "%s %lf %lf %lf %lf %i ", wd->id[i], &wd->x[i], &wd->y[i], &wd->z1[i], &wd->z2[i], &wd->nWellObs[i] );
			if( status != 6 ) { tprintf( "ERROR: Well %s data provided in the input file %s is incomplete; input file error!\n", wd->id[i], filename ); bad_data = 1; }
			if( cd->debug ) tprintf( "Well %-6s x %8g y %8g z0 %6g z1 %6g nObs %2i ", wd->id[i], wd->x[i], wd->y[i], wd->z1[i], wd->z2[i], wd->nWellObs[i] );
			if( wd->nWellObs[i] <= 0 ) { if( cd->debug ) tprintf( "WARNING: no observations!\n" ); fscanf( infile, "%lf %lf %lf %i %lf %lf\n", &d, &d, &d, &j, &d, &d ); continue; }
			if( ( wd->obs_time[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( wd->obs_target[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( wd->obs_weight[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( wd->obs_log[i] = ( int * ) malloc( wd->nWellObs[i] * sizeof( int ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( wd->obs_min[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			if( ( wd->obs_max[i] = ( double * ) malloc( wd->nWellObs[i] * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
			for( j = 0; j < wd->nWellObs[i]; j++ )
			{
				wd->obs_min[i][j] = -1e6; wd->obs_max[i][j] = 1e6; wd->obs_weight[i][j] = 1; wd->obs_log[i][j] = 0;
				status = fscanf( infile, "%lf %lf %lf %i %lf %lf\n", &wd->obs_time[i][j], &wd->obs_target[i][j], &wd->obs_weight[i][j], &wd->obs_log[i][j], &wd->obs_min[i][j], &wd->obs_max[i][j] );
				if( status != 6 )
				{
					tprintf( "ERROR:\tObservation data provided for well %s in the input file %s is incomplete; input file error!\n", wd->id[i], filename );
					tprintf( "\tWell %-6s x %8g y %8g z0 %6g z1 %6g nObs %2i\n", wd->id[i], wd->x[i], wd->y[i], wd->z1[i], wd->z2[i], wd->nWellObs[i] );
					tprintf( "\tObservation #%d: time %5g concentration %5g weight %7g log %1d acceptable range: min %5g max %5g\n\n", j + 1, wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_weight[i][j], wd->obs_log[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
					bad_data = 1;
				}
				if( cd->obsdomain > DBL_EPSILON && wd->obs_weight[i][j] > DBL_EPSILON ) { wd->obs_min[i][j] = wd->obs_target[i][j] - cd->obsdomain; wd->obs_max[i][j] = wd->obs_target[i][j] + cd->obsdomain; }
				if( cd->ologtrans == 1 ) wd->obs_log[i][j] = 1;
				else if( cd->ologtrans == 0 ) wd->obs_log[i][j] = 0;
				if( cd->oweight == 1 ) wd->obs_weight[i][j] = 1;
				else if( cd->oweight == 0 ) wd->obs_weight[i][j] = 0;
				else if( cd->oweight == 2 ) { if( fabs( wd->obs_target[i][j] ) > DBL_EPSILON ) wd->obs_weight[i][j] = ( double ) 1.0 / wd->obs_target[i][j]; else wd->obs_weight[i][j] = HUGE_VAL; }
				if( cd->debug )
					tprintf( "t %5g c %5g weight %7g log %1d acceptable range: min %5g max %5g\n", wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_weight[i][j], wd->obs_log[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
				if( wd->obs_max[i][j] < wd->obs_target[i][j] || wd->obs_min[i][j] > wd->obs_target[i][j] )
				{
					tprintf( "ERROR: Observation target is outside the specified min/max range! " );
					tprintf( "Observation %s(%g): %g min %g max %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
					bad_data = 1;
				}
				if( wd->obs_max[i][j] <= wd->obs_min[i][j] )
				{
					tprintf( "ERROR: Calibration range is not correctly specified! " );
					tprintf( "Observation %s(%g): min %g max %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
					bad_data = 1;
				}
				if( wd->obs_weight[i][j] > DBL_EPSILON ) od->nObs++;
				if( wd->obs_weight[i][j] < -DBL_EPSILON ) { preds->nTObs++; if( od->include_predictions ) od->nObs++; } // Predictions have negative weights
				if( j + 1 < wd->nWellObs[i] ) { fscanf( infile, "\t\t" ); if( cd->debug ) tprintf( "\t\t\t\t\t\t\t      " ); }
			}
		}
		od->nCObs = od->nObs;
		map_well_obs( op );
	}
	// ------------------------------------------------------------ Reading regularization terms ----------------------------------------------------------------
	fscanf( infile, "%1000[^:]s", buf );
	if( !strncasecmp( buf, "Number of regul", 15 ) )
	{
		if( cd->debug && short_names_printed == 0 ) // if param names and values are not already printed (above for param expressions)
		{
			tprintf( "\nParameter values for computation of the regularization terms:\n" );
			for( i = 0; i < pd->nParam; i++ )
			{
				if( pd->var_opt[i] > 0 && pd->var_log[i] == 1 ) cd->var[i] = pow( 10, pd->var[i] );
				else cd->var[i] = pd->var[i];
				if( cd->debug )
				{
					if( cd->solution_type[0] == EXTERNAL ) tprintf( "%-27s: %.12g\n", pd->var_name[i], cd->var[i] );
					else tprintf( "%-27s: %-12s: %.12g\n", pd->var_name[i], pd->var_id[i], cd->var[i] );
				}
			}
		}
		fscanf( infile, ": %i\n", &rd->nRegul );
		if( cd->debug ) tprintf( "\n" );
		tprintf( "Number of regularization terms = %d\n", rd->nRegul );
		rd->regul_expression = ( void ** ) malloc( rd->nRegul * sizeof( void * ) );
		rd->regul_id = char_matrix( rd->nRegul, 10 );
		rd->regul_target = ( double * ) malloc( rd->nRegul * sizeof( double ) );
		rd->regul_weight = ( double * ) malloc( rd->nRegul * sizeof( double ) );
		rd->regul_min = ( double * ) malloc( rd->nRegul * sizeof( double ) );
		rd->regul_max = ( double * ) malloc( rd->nRegul * sizeof( double ) );
		rd->regul_log = ( int * ) malloc( rd->nRegul * sizeof( int ) );
#ifdef MATHEVAL
		rd->regul_nMap = pd->nAnalParam + od->nObs; // rd->nRegul is not needed
		rd->regul_map_id = ( char ** ) malloc( ( rd->regul_nMap ) * sizeof( char * ) );
		rd->regul_map_val = ( double * ) malloc( ( rd->regul_nMap + rd->nRegul ) * sizeof( double ) ); // rd->nRegul added to accommodate cd->obs_current
#endif
#ifndef MATHEVAL
		tprintf( "WARNING: MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
		for( i = 0; i < rd->nRegul; i++ )
		{
			sprintf( rd->regul_id[i], "reg%d", i + 1 );
			fscanf( infile, "%1000[^=]s", buf );
			fscanf( infile, "= %lf %lf %i %lf %lf\n", &rd->regul_target[i], &rd->regul_weight[i], &rd->regul_log[i], &rd->regul_min[i], &rd->regul_max[i] );
			if( cd->debug ) tprintf( "%-12s: target %g weight %g log %i min %g max %g : equation %s", rd->regul_id[i], rd->regul_target[i], rd->regul_weight[i], rd->regul_log[i], rd->regul_min[i], rd->regul_max[i], buf );
			if( !( rd->regul_weight[i] > DBL_EPSILON ) ) tprintf( " WARNING Weight <= 0 " );
#ifdef MATHEVAL
			rd->regul_expression[i] = evaluator_create( buf );
			assert( rd->regul_expression[i] );
			evaluator_get_variables( rd->regul_expression[i], &expvar_names, &expvar_count );
#else
			expvar_count = 0;
#endif
			if( expvar_count > 0 )
			{
				if( cd->debug )
				{
					tprintf( " -> variables:" );
					for( j = 0; j < expvar_count; j++ )
						tprintf( " %s", expvar_names[j] );
					tprintf( "\n" );
				}
				for( j = 0; j < expvar_count; j++ )
				{
					l1 = strlen( expvar_names[j] );
					status = 0;
					for( k = 0; k < pd->nParam; k++ )
						if( !strncasecmp( expvar_names[j], pd->var_id[k], l1 ) ) { status = 1; break; }
					for( k = 0; k < od->nObs; k++ )
						if( !strncasecmp( expvar_names[j], od->obs_id[k], l1 ) ) { status = 1; break; }
#ifdef MATHEVAL
					if( status == 0 ) { tprintf( "ERROR: parameter name \'%s\' in regularization term \'%s\' is not defined!\n", expvar_names[j], evaluator_get_string( rd->regul_expression[i] ) ); bad_data = 1; }
#endif
				}
			}
#ifdef MATHEVAL
			else { tprintf( "ERROR: no variables\n" ); bad_data = 1; }
#endif
		}
#ifdef MATHEVAL
		for( k = 0; k < pd->nAnalParam; k++ ) { rd->regul_map_id[k] = pd->var_id[k]; rd->regul_map_val[k] = cd->var[k]; }
		for( i = pd->nAnalParam, k = 0; k < od->nObs; k++, i++ ) { rd->regul_map_id[i] = od->obs_id[k]; rd->regul_map_val[i] = od->obs_current[k] = od->obs_target[k]; }
		// free( cd->var ); // TODO Linux fails to free cd->var; needs debugging
		free( od->obs_current );
		cd->var = &rd->regul_map_val[0];
		od->obs_current = &rd->regul_map_val[pd->nAnalParam];
		// for( k = 0; k < rd->regul_nMap; k++ ) { tprintf( "%s %g\n", rd->regul_map_id[k], rd->regul_map_val[k] ); }
#endif
		if( cd->debug )
		{
#ifdef MATHEVAL
			tprintf( "Regularization expressions evaluated (initial values applied for parameters).\n" );
#endif
			for( i = 0; i < rd->nRegul; i++ )
			{
				tprintf( "%-12s= ", rd->regul_id[i] );
#ifdef MATHEVAL
				tprintf( "%s", evaluator_get_string( rd->regul_expression[i] ) );
				d = evaluator_evaluate( rd->regul_expression[i], rd->regul_nMap, rd->regul_map_id, rd->regul_map_val );
				tprintf( " = %g\n", d );
#else
				tprintf( "MathEval is not installed; expressions cannot be evaluated.\n" );
#endif
			}
		}
		fscanf( infile, "%1000[^:]s", buf );
	}
	else tprintf( "Number of regularization terms = %d\n", rd->nRegul );
	if( bad_data ) return ( 0 );
	// add regularization terms to observations
	if( rd->nRegul > 0 )
	{
		tprintf( "Number of total observations & regularizations = %d\n", od->nTObs );
		map_obs( op ); // add regularizations to the observations
		if( cd->debug )
			for( i = 0; i < rd->nRegul; i++ )
				tprintf( "%s: %g weight %g", rd->regul_id[i], rd->regul_target[i], rd->regul_weight[i] );
	}
	// ------------------------------------------------------------ Set predictions ----------------------------------------------------------------
	if( preds->nTObs > 0 ) // TODO add regularization in INFOGAP and GLUE analysis
	{
		if( cd->problem_type == INFOGAP ) tprintf( "Number of performance criterion predictions for info-gap analysis = %d\n", preds->nTObs );
		else tprintf( "Number of predictions = %d\n", preds->nTObs );
		preds->obs_index = ( int * ) malloc( preds->nTObs * sizeof( int ) );
		preds->obs_target = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		preds->obs_current = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		preds->obs_best = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		preds->obs_well_index = ( int * ) malloc( preds->nTObs * sizeof( int ) );
		preds->obs_time_index = ( int * ) malloc( preds->nTObs * sizeof( int ) );
		preds->obs_id = char_matrix( preds->nTObs, 50 );
		preds->obs_weight = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		preds->obs_min = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		preds->obs_max = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		preds->obs_log = ( int * ) malloc( preds->nTObs * sizeof( int ) );
		preds->res = ( double * ) malloc( preds->nTObs * sizeof( double ) );
		for( c = k = i = 0; i < wd->nW; i++ )
			for( j = 0; j < wd->nWellObs[i]; j++ )
			{
				if( fabs( wd->obs_weight[i][j] ) > DBL_EPSILON ) c++;
				if( wd->obs_weight[i][j] < -DBL_EPSILON )
				{
					preds->obs_index[k] = c - 1;
					preds->obs_target[k] = wd->obs_target[i][j];
					preds->obs_weight[k] = 1.0;
					preds->obs_min[k] = wd->obs_target[i][j];
					preds->obs_max[k] = wd->obs_target[i][j];
					preds->obs_log[k] = wd->obs_log[i][j];
					preds->obs_well_index[k] = i;
					preds->obs_time_index[k] = j;
					sprintf( preds->obs_id[k], "%s(%g)", wd->id[i], wd->obs_time[i][j] );
					if( cd->debug ) tprintf( "%s(%g): %g weight %g\n", wd->id[i], wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_weight[i][j] );
					k++;
				}
			}
	}
	else
	{
		tprintf( "Number of predictions = %d\n", preds->nTObs );
		if( cd->problem_type == INFOGAP ) // INFOGAP problem
		{
			tprintf( "\nERROR: Weight of at least one observation must be set as performance criterion prediction\nby setting weight to -1 for Info-gap analysis\n\n" );
			bad_data = 1;
		}
		if( cd->problem_type == GLUE ) // GLUE problem
		{
			tprintf( "\nERROR: Weight of at least one observation must be set as a prediction\nby setting weight to -1 for GLUE analysis\n\n" );
			bad_data = 1;
		}
	}
	if( bad_data ) return( 0 );
	// ------------------------------------------------------------ Reading external problem ----------------------------------------------------------------
	if( cd->solution_type[0] == EXTERNAL )
	{
		// Executable Command Line
		ed->cmdline = ( char * ) malloc( 80 * sizeof( char ) );
		fscanf( infile, ": " ); fgets( ed->cmdline, 80, infile );
		ed->cmdline[strlen( ed->cmdline ) - 1] = 0;
		tprintf( "Execution command: %s\n", ed->cmdline );
		if( sscanf( ed->cmdline, "%i", &i ) == -1 )
		{
			tprintf( "ERROR: Execution command is not valid!\n" );
			bad_data = 1;
		}
		strcpy( buf, ed->cmdline );
		file = &buf[0];
		file = strsep( &file, " \t" );
		k = 1;
		if( access( file, X_OK ) == -1 )
		{
			k = 0;
			if( file[0] == '~' )
			{
				sprintf( exec, "%s/%s", getenv( "HOME" ), &file[1] );
				if( access( exec, X_OK ) == 0 )
					k = 1;
			}
			if( k == 0 )
			{
				path = shellpath();
				for( i = 0; k == 0 && path[i]; i++ )
				{
					if( cd->debug > 2 ) tprintf( "%s\n", path[i] );
					sprintf( exec, "%s/%s", path[i], file );
					if( access( exec, X_OK ) == 0 )
						k = 1;
				}
			}
		}
		if( k == 0 )
		{
			tprintf( "ERROR: Program \'%s\' does not exist or cannot be executed!\n", file );
			bad_data = 1;
		}
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %d\n", &ed->ntpl );
		ed->fn_tpl = char_matrix( ed->ntpl, 80 );
		ed->fn_out = char_matrix( ed->ntpl, 80 );
		for( i = 0; i < ed->ntpl; i++ )
		{
			fscanf( infile, "%s %s\n", ed->fn_tpl[i], ed->fn_out[i] );
			if( access( ed->fn_tpl[i], R_OK ) == -1 )
			{
				tprintf( "ERROR: File \'%s\' does not exist!\n", ed->fn_tpl[i] );
				bad_data = 1;
			}
		}
		tprintf( "External files:\n" );
		tprintf( "- to provide current model parameters:\n" );
		for( i = 0; i < ed->ntpl; i++ )
			tprintf( "%s -> %s\n", ed->fn_tpl[i], ed->fn_out[i] );
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %d\n", &ed->nins );
		ed->fn_ins = char_matrix( ed->nins, 80 );
		ed->fn_obs = char_matrix( ed->nins, 80 );
		for( i = 0; i < ed->nins; i++ )
		{
			fscanf( infile, "%s %s\n", ed->fn_ins[i], ed->fn_obs[i] );
			if( access( ed->fn_ins[i], R_OK ) == -1 )
			{
				tprintf( "ERROR: File \'%s\' does not exist!\n", ed->fn_ins[i] );
				bad_data = 1;
			}
		}
		tprintf( "- to read current model predictions:\n" );
		for( i = 0; i < ed->nins; i++ )
			tprintf( "%s <- %s\n", ed->fn_ins[i], ed->fn_obs[i] );
		fclose( infile );
		gd->min_t = gd->time = 0;
		tprintf( "\n" );
		if( bad_data ) return ( 0 );
	}
	// ------------------------------------------------------------ Reading internal problem ----------------------------------------------------------------
	else
	{
		// ------------------------------------------------------------ Read grid and breakthrough computational data ----------------------------------------------------------------
		fscanf( infile, ": %lf\n", &gd->time );
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %i %i %i\n", &gd->nx, &gd->ny, &gd->nz );
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &gd->min_x, &gd->min_y, &gd->min_z );
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &gd->max_x, &gd->max_y, &gd->max_z );
		fscanf( infile, "%1000[^:]s", buf ); fscanf( infile, ": %lf %lf %lf\n", &gd->min_t, &gd->max_t, &gd->dt );
		fclose( infile );
		if( cd->debug )
		{
			tprintf( "\nGrid Time: %g\n", gd->time );
			tprintf( "Grid lines: %i %i %i\n", gd->nx, gd->ny, gd->nz );
			tprintf( "Grid Minimums: %g %g %g\n", gd->min_x, gd->min_y, gd->min_z );
			tprintf( "Grid Maximums: %g %g %g\n", gd->max_x, gd->max_y, gd->max_z );
		}
		if( gd->nx == 1 ) gd->dx = 0;
		else gd->dx = ( gd->max_x - gd->min_x ) / ( gd->nx - 1 );
		if( gd->ny == 1 ) gd->dy = 0;
		else gd->dy = ( gd->max_y - gd->min_y ) / ( gd->ny - 1 );
		// if(gd->nz == 1 ) gd->dz = gd->max_z - gd->min_z ); // In this way compute_grid computed for min_z
		if( gd->nz == 1 ) gd->dz = 0;
		else gd->dz = ( gd->max_z - gd->min_z ) / ( gd->nz - 1 );
		gd->nt = 1 + ( int )( ( double )( gd->max_t - gd->min_t ) / gd->dt );
		if( gd->nt < 0 ) gd->nt = 0;
		if( cd->debug ) tprintf( "Breakthrough-curve time window: %g %g %g number of time steps: %g\n", gd->min_t, gd->max_t, gd->dt, gd->nt );
	}
	return( 1 );
}

int save_problem( char *filename, struct opt_data *op )
{
	struct calc_data *cd;
	struct param_data *pd;
	struct regul_data *rd;
	struct obs_data *od;
	struct well_data *wd;
	struct grid_data *gd;
	struct extrn_data *ed;
	cd = op->cd;
	pd = op->pd;
	rd = op->rd;
	od = op->od;
	wd = op->wd;
	gd = op->gd;
	ed = op->ed;
	FILE *outfile;
	//	FILE *outfileb;
	//	char buf[100];
	int  i, j;
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		tprintf( "File \'%s\' cannot be opened to save the problem information!\n", filename );
		return( 0 );
	}
	/*
		sprintf( buf, "%s.bin", filename );
		if ( ( outfileb = fopen( buf, "wb" ) ) == NULL )
		 tprintf( "Binary file %s cannot be opened to save problem information!\n", buf );
		else
		{
			fwrite( (void *) (*pd).var, sizeof((*pd).var[i]), (*pd).nParam, outfileb );
			fclose( outfileb );
		}
	 */
	i = 0;
	fprintf( outfile, "Problem: " );
	switch( cd->problem_type )
	{
		case CREATE: fprintf( outfile, "create" ); break;
		case FORWARD: fprintf( outfile, "forward" ); break;
		case CALIBRATE: fprintf( outfile, "calibration" ); break;
		case LOCALSENS: fprintf( outfile, "lsens" ); break;
		case GLOBALSENS: fprintf( outfile, "gsens" ); break;
		case EIGEN: fprintf( outfile, "eigen" ); break;
		case MONTECARLO: fprintf( outfile, "montecarlo real=%d", cd->nreal ); break;
		case ABAGUS: fprintf( outfile, " abagus energy=%d", cd->energy ); break;
		case POSTPUA: fprintf( outfile, " postpua" ); break;
	}
	if( cd->debug > 0 ) fprintf( outfile, " debug=%d", cd->debug );
	if( cd->fdebug > 0 ) fprintf( outfile, " fdebug=%d", cd->fdebug );
	if( cd->ldebug > 0 ) fprintf( outfile, " ldebug=%d", cd->ldebug );
	if( cd->pdebug > 0 ) fprintf( outfile, " pdebug=%d", cd->pdebug );
	if( cd->mdebug > 0 ) fprintf( outfile, " mdebug=%d", cd->mdebug );
	if( cd->odebug > 0 ) fprintf( outfile, " odebug=%d", cd->odebug );
	if( cd->insdebug > 0 ) fprintf( outfile, " insdebug=%d", cd->insdebug );
	if( cd->tpldebug > 0 ) fprintf( outfile, " tpldebug=%d", cd->tpldebug );
	if( cd->pardebug > 0 ) fprintf( outfile, " pardebug=%d", cd->pardebug );
	if( cd->test_func >= 0 ) fprintf( outfile, " test=%d", cd->test_func );
	if( cd->test_func_dim > 2 ) fprintf( outfile, " dim=%d", cd->test_func_dim );
	if( cd->phi_cutoff > 0 ) fprintf( outfile, " cutoff=%g", cd->phi_cutoff );
	if( cd->sintrans ) { if( cd->sindx > DBL_EPSILON ) fprintf( outfile, " sindx=%g", cd->sindx ); }
	else { if( cd->lindx > DBL_EPSILON ) fprintf( outfile, " lindx=%g", cd->lindx ); }
	// if( cd->pardx > DBL_EPSILON ) fprintf( outfile, " pardx=%g", cd->pardx ); TODO when to print pardx?
	if( cd->check_success ) fprintf( outfile, " success" );
	fprintf( outfile, " " );
	switch( cd->calib_type )
	{
		case SIMPLE: fprintf( outfile, "single" ); break;
		case PPSD: fprintf( outfile, "ppsd" ); break;
		case IGRND: fprintf( outfile, "igrnd real=%d", cd->nreal ); break;
		case IGPD: fprintf( outfile, "igpd" ); break;
	}
	fprintf( outfile, " eval=%d", cd->maxeval );
	if( cd->opt_method[0] != 0 ) fprintf( outfile, " opt=%s", cd->opt_method );
	if( cd->c_background > 0 ) fprintf( outfile, " background=%g", cd->c_background );
	if( cd->disp_tied ) fprintf( outfile, " disp_tied" );
	if( cd->disp_scaled ) fprintf( outfile, " disp_scaled" );
	if( cd->save ) fprintf( outfile, " save" );
	if( cd->seed_init < 0 ) fprintf( outfile, " seed=%d", cd->seed_init * -1 );
	if( cd->nretries > 0 ) fprintf( outfile, " retry=%d", cd->nretries );
	if( cd->init_particles > 1 ) fprintf( outfile, " particles=%d", cd->init_particles );
	else if( cd->init_particles < 0 ) fprintf( outfile, " particles" );
	if( cd->lm_eigen ) fprintf( outfile, " lmeigen=%d", cd->lm_eigen );
	if( cd->niter > 0 ) fprintf( outfile, " iter=%d", cd->niter );
	if( cd->smp_method[0] != 0 ) fprintf( outfile, " rnd=%s", cd->smp_method );
	if( cd->paran_method[0] != 0 ) fprintf( outfile, " paran=%s", cd->paran_method );
	fprintf( outfile, " " );
	switch( cd->objfunc_type )
	{
		case SSR: fprintf( outfile, "ssr" ); break;
		case SSDR: fprintf( outfile, "ssdr" ); break;
		case SSD0: fprintf( outfile, "ssd0" ); break;
		case SSDX: fprintf( outfile, "ssdx" ); break;
		case SSDA: fprintf( outfile, "ssda" ); break;
	}
	fprintf( outfile, "\n" );
	fprintf( outfile, "Solution: %s\n", cd->solution_id );
	fprintf( outfile, "Number of parameters: %i\n", pd->nParam );
	for( i = j = 0; i < pd->nParam; i++ )
	{
		if( pd->var_opt[i] == -1 ) // tied parameter
#ifdef MATHEVAL
			fprintf( outfile, "%s= %s\n", pd->var_name[i], evaluator_get_string( pd->param_expression[j++] ) );
#else
			fprintf( outfile, "%s= MathEval is not installed; expressions cannot be evaluated\n", pd->var_name[i] );
#endif
		else if( pd->var_opt[i] >= 1 && pd->var_log[i] == 1 ) // optimized log transformed parameter
			fprintf( outfile, "%s: %.15g %d %d %g %g %g\n", pd->var_name[i], pow( 10, pd->var[i] ), pd->var_opt[i], pd->var_log[i], pow( 10, pd->var_dx[i] ), pow( 10, pd->var_min[i] ), pow( 10, pd->var_max[i] ) );
		else // fixed or not log-transformed parameter
			fprintf( outfile, "%s: %.15g %d %d %g %g %g\n", pd->var_name[i], pd->var[i], pd->var_opt[i], pd->var_log[i], pd->var_dx[i], pd->var_min[i], pd->var_max[i] );
	}
	if( cd->solution_type[0] != EXTERNAL )
	{
		fprintf( outfile, "Number of wells: %i\n", wd->nW );
		for( i = 0; i < wd->nW; i++ )
		{
			fprintf( outfile, "%s %.15g %.15g %g %g %i ", wd->id[i], wd->x[i], wd->y[i], wd->z1[i], wd->z2[i], wd->nWellObs[i] );
			if( wd->nWellObs[i] > 1 ) { fprintf( outfile, "\n" ); }
			for( j = 0; j < wd->nWellObs[i]; j++ )
				fprintf( outfile, "%g %g %g %i %g %g\n", wd->obs_time[i][j], wd->obs_target[i][j], wd->obs_weight[i][j], wd->obs_log[i][j], wd->obs_min[i][j], wd->obs_max[i][j] );
		}
	}
	else
	{
		fprintf( outfile, "Number of observations: %i\n", od->nTObs );
		for( i = 0; i < od->nTObs; i++ )
			fprintf( outfile, "%s %g %g %d %g %g\n", od->obs_id[i], od->obs_target[i], od->obs_weight[i], od->obs_log[i], od->obs_min[i], od->obs_max[i] );
		fprintf( outfile, "Execution command: %s\n", ed->cmdline );
		fprintf( outfile, "Number of execution templates: %d\n", ed->ntpl );
		for( i = 0; i < ed->ntpl; i++ )
			fprintf( outfile, "%s %s\n", ed->fn_tpl[i], ed->fn_out[i] );
		fprintf( outfile, "Number of execution instructions: %d\n", ed->nins );
		for( i = 0; i < ed->nins; i++ )
			fprintf( outfile, "%s %s\n", ed->fn_ins[i], ed->fn_obs[i] );
	}
	if( rd->nRegul > 0 )
	{
		fprintf( outfile, "Number of regularization terms: %d\n", rd->nRegul );
		for( i = 0; i < rd->nRegul; i++ )
#ifdef MATHEVAL
			fprintf( outfile, "%s = %g %g %i %g %g\n", evaluator_get_string( rd->regul_expression[i] ), rd->regul_target[i], rd->regul_weight[i], rd->regul_log[i], rd->regul_min[i], rd->regul_max[i] );
#else
			fprintf( outfile, "Regularization term #%d = %g %g %i %g %g\n", i, rd->regul_target[i], rd->regul_weight[i], rd->regul_log[i], rd->regul_min[i], rd->regul_max[i] );
#endif
	}
	if( cd->solution_type[0] != EXTERNAL )
	{
		fprintf( outfile, "Grid Time: %g\n", gd->time );
		fprintf( outfile, "Grid lines: %i %i %i\n", gd->nx, gd->ny, gd->nz );
		fprintf( outfile, "Grid Minimums: %g %g %g\n", gd->min_x, gd->min_y, gd->min_z );
		fprintf( outfile, "Grid Maximums: %g %g %g\n", gd->max_x, gd->max_y, gd->max_z );
		fprintf( outfile, "Time Window: %g %g %g\n", gd->min_t, gd->max_t, gd->dt );
	}
	fclose( outfile );
	return( 1 );
}

void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd )
{
	FILE *outfile;
	double x, y, z, t, c, dx, dy, max_conc, max_x, max_y, max_z;
	int  i, j, k, nx, ny;
	if( gd->time <= 0 ) return;
	t = gd->time;
	nx = ny = 15;
	dx = ( gd->max_x - gd->min_x ) / nx;
	dy = ( gd->max_y - gd->min_y ) / ny;
	z = gd->min_z;
	for( j = 0; j < ny; j++ )
	{
		y = gd->max_y - j * dy;
		for( i = 0; i < nx; i++ )
		{
			x = i * dx + gd->min_x;
			c = func_solver1( x, y, z, t, ( void * ) cd );
			tprintf( "%6.0g ", c );
		}
		tprintf( "\n" );
	}
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		tprintf( "Output file %s cannot be opened!\n", filename );
		return;
	}
	fprintf( outfile, "# vtk DataFile Version 2.0\n" );
	fprintf( outfile, "MADS output\n" );
	fprintf( outfile, "ASCII\n" );
	fprintf( outfile, "DATASET STRUCTURED_POINTS\n" );
	fprintf( outfile, "DIMENSIONS %d %d %d\n", gd->nx, gd->ny, gd->nz );
	fprintf( outfile, "ORIGIN %g %g %g\n", gd->min_x, gd->min_y, gd->min_z ); /* negative z */
	//	fprintf( outfile, "ORIGIN %g %g %g\n", gd->min_x, gd->min_y, -gd->max_z ); /* negative z */
	fprintf( outfile, "SPACING %g %g %g\n", gd->dx, gd->dy, gd->dz );
	fprintf( outfile, "POINT_DATA %d\n", gd->nx * gd->ny * gd->nz );
	fprintf( outfile, "SCALARS Cr_concentration_ppb float 1\n" );
	fprintf( outfile, "LOOKUP_TABLE default\n" );
	max_conc = max_x = max_y = max_z = 0;
	for( k = 0; k < gd->nz; k++ )
	{
		//		z = gd->max_z - k * gd->dz; /* reverse z */
		z = k * gd->dz + gd->min_z;
		for( j = 0; j < gd->ny; j++ )
		{
			y = j * gd->dy + gd->min_y;
			for( i = 0; i < gd->nx; i++ )
			{
				x = i * gd->dx + gd->min_x;
				c = func_solver1( x, y, z, t, ( void * ) cd );
				fprintf( outfile, "%g ", c );
				if( max_conc < c ) { max_conc = c; max_x = x; max_y = y; max_z = z; }
			}
			fprintf( outfile, "\n" );
		}
		fprintf( outfile, "\n" );
	}
	tprintf( "Max Conc = %g @ (%g,%g,%g)\n", max_conc, max_x, max_y, max_z );
	fclose( outfile );
	tprintf( "Spatial concentration data saved in %s.\n", filename );
}

void compute_btc2( char *filename, char *filename2, struct opt_data *op )
{
	FILE *outfile;
	double time, time_expected, c, max_source_conc, max_source_time, max_source_x, max_source_y, *max_conc, *max_time, d, x0, y0, alpha, beta, xe, ye, v, v_apparent;
	int  i, k, s, p;
	struct grid_data *gd;
	gd = op->gd;
	if( gd->min_t <= 0 ) return;
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		tprintf( "Output file %s cannot be opened!\n", filename );
		return;
	}
	tprintf( "\n" );
	max_time = ( double * ) malloc( op->wd->nW * sizeof( double ) );
	max_conc = ( double * ) malloc( op->wd->nW * sizeof( double ) );
	max_source_x = max_source_y = max_source_conc = 0;
	max_source_time = op->pd->var[TIME_INIT];
	fprintf( outfile, "variables = \"Time [a]\"" );
	for( i = 0; i < op->wd->nW; i++ )
	{
		fprintf( outfile, " \"%s\"", op->wd->id[i] );
		max_conc[i] = 0;
		max_time[i] = -1;
	}
	for( s = 0; s < op->cd->num_sources; s++ )
		fprintf( outfile, " \"S%d\"", s + 1 );
	fprintf( outfile, "\n" );
	for( k = 0; k < gd->nt; k++ )
	{
		time = gd->min_t + gd->dt * k;
		fprintf( outfile, "%g", time );
		for( i = 0; i < op->wd->nW; i++ )
		{
			c = func_solver( op->wd->x[i], op->wd->y[i], op->wd->z1[i], op->wd->z2[i], time, ( void * ) op->cd );
			fprintf( outfile, " %g", c );
			if( max_conc[i] < c ) { max_conc[i] = c; max_time[i] = time; }
		}
		for( s = 0; s < op->cd->num_sources; s++ )
		{
			p = s * NUM_ANAL_PARAMS_SOURCE;
			c = func_solver( op->pd->var[p + SOURCE_X], op->pd->var[p + SOURCE_Y], op->pd->var[p + SOURCE_Z], op->pd->var[p + SOURCE_Z], time, ( void * ) op->cd );
			fprintf( outfile, " %g", c );
			if( max_source_conc < c ) { max_source_conc = c; max_source_x = op->pd->var[p + SOURCE_X]; max_source_y = op->pd->var[p + SOURCE_Y]; max_source_time = time; }
		}
		fprintf( outfile, "\n" );
	}
	fclose( outfile );
	tprintf( "Concentration breakthrough data saved in %s\n", filename );
	tprintf( "\nPeak source concentration (x = %g, y = %g, t = %g) = %g\n", max_source_x, max_source_y, max_source_time, max_source_conc );
	if( ( outfile = fopen( filename2, "w" ) ) == NULL )
	{
		tprintf( "Output file %s cannot be opened!\n", filename );
		return;
	}
	for( i = 0; i < op->pd->nParam; i++ ) // IMPORTANT: Take care of log transformed variable
		if( op->pd->var_opt[i] >= 1 && op->pd->var_log[i] == 1 )
			op->pd->var[i] = pow( 10, op->pd->var[i] );
	i = ( op->cd->num_sources - 1 ) * NUM_ANAL_PARAMS_SOURCE;
	v = sqrt( op->pd->var[i + VX] * op->pd->var[i + VX] + op->pd->var[i + VY] * op->pd->var[i + VY] + op->pd->var[i + VZ] * op->pd->var[i + VZ] ); // Flow velocity
	// tprintf( "Flow velocity pointers = (%d %d %d)\n", i+VX, i+VY, i+VZ );
	tprintf( "Flow velocity = %g (%g %g %g)\n", v, op->pd->var[i + VX], op->pd->var[i + VY], op->pd->var[i + VZ] );
	for( i = 0; i < op->wd->nW; i++ )
	{
		x0 = ( op->wd->x[i] - max_source_x );
		y0 = ( op->wd->y[i] - max_source_y );
		d = ( -op->pd->var[( op->cd->num_sources - 1 ) * NUM_ANAL_PARAMS_SOURCE + FLOW_ANGLE] * M_PI ) / 180;
		alpha = cos( d );
		beta = sin( d );
		xe = x0 * alpha - y0 * beta;
		ye = x0 * beta  + y0 * alpha;
		d = sqrt( xe * xe + ye * ye );
		time = max_time[i] - max_source_time;
		if( time > DBL_EPSILON ) v_apparent = d / time; else { v_apparent = -1; if( time < 0 ) time = -1; };
		c = max_conc[i] / max_source_conc; // Normalized concentration
		if( v > DBL_EPSILON ) time_expected = d / v; else time_expected = -1;
		tprintf( "%s\tPeak Concentration  = %12g (%12g) @ time %12g (%12g expected %12g) velocity = %12g (%12g) distance = %12g\n", op->wd->id[i], max_conc[i], c, max_time[i], time, time_expected, v_apparent, v, d );
		fprintf( outfile, "%s\tPeak Conc = %12g (%12g) @ time %12g (%12g expected %12g) velocity = %12g (%12g) distance = %12g\n", op->wd->id[i], max_conc[i], c, max_time[i], time, time_expected, v_apparent, v, d );
	}
	fclose( outfile );
	tprintf( "Concentration peak data saved in %s\n", filename2 );
	for( i = 0; i < op->pd->nParam; i++ ) // IMPORTANT: Take care of log transformed variable
		if( op->pd->var_opt[i] >= 1 && op->pd->var_log[i] == 1 )
			op->pd->var[i] = log10( op->pd->var[i] );
	free( max_conc );
	free( max_time );
}

void compute_btc( char *filename, struct opt_data *op )
{
	FILE *outfile;
	double time, c, max_conc, max_time;
	int  i, k;
	struct grid_data *gd;
	gd = op->gd;
	if( gd->min_t <= 0 ) return;
	if( ( outfile = fopen( filename, "w" ) ) == NULL )
	{
		tprintf( " Output file %s cannot be opened!\n", filename );
		return;
	}
	tprintf( "\n" );
	for( i = 0; i < op->wd->nW; i++ )
	{
		c = time = max_conc = max_time = 0;
		for( k = 0; k < gd->nt; k++ )
		{
			time = gd->min_t + gd->dt * k;
			c = func_solver( op->wd->x[i], op->wd->y[i], op->wd->z1[i], op->wd->z2[i], time, ( void * ) op->cd );
			fprintf( outfile, "%s %g %g\n", op->wd->id[i], time, c );
			if( max_conc < c ) { max_conc = c; max_time = time; }
		}
		tprintf( "%s\tPeak Conc = %12.4g @ time %12g\tConc = %12.4f @ time %12g\n", op->wd->id[i], max_conc, max_time, c, time );
	}
	fclose( outfile );
	tprintf( "Concentration breakthrough data saved in %s\n", filename );
}

static char *strsave( const char *s, const char *lim )
{
	if( lim == NULL )
		lim = s + strlen( s );
	char *p = ( char * ) malloc_check( "save string", lim - s + 1 );
	strncpy( p, s, lim - s );
	p[lim - s] = '\0';
	return p;
}

char **shellpath( void )
{
	const char *path = getenv( "PATH" );
	if( !path )
		path = "/bin:/usr/bin:/usr/local/bin";
	char **vector = // size is overkill
		( char ** ) malloc_check( "hold path elements", strlen( path ) * sizeof( *vector ) );
	const char *p = path;
	int next = 0;
	while( p )
	{
		char *q = strchr( p, ':' );
		vector[next++] = strsave( p, q );
		p = q ? q + 1 : NULL;
	}
	vector[next] = NULL;
	return vector;
}

void freeshellpath( char *shellpath[] )
{
	int i;
	for( i = 0; shellpath[i]; i++ )
		free( shellpath[i] );
	free( shellpath );
}

unsigned maxpathlen( char *path[], const char *base )
{
	unsigned blen = strlen( base );
	unsigned n = 0;
	int i;
	for( i = 0; path[i]; i++ )
	{
		unsigned pn = strlen( path[i] );
		if( pn > n ) n = pn;
	}
	return blen + n + 1;
}

void execvepath( char *path[], const char *base, char *const argv[], char *const envp[] )
{
	int i;
	if( strchr( base, '/' ) )
		execve( base, argv, envp );
	else
	{
		size_t maxlen = maxpathlen( path, base ) + 1;
		char *buf = ( char * ) malloc_check( "hold path", maxlen );
		for( i = 0; path[i]; i++ )
		{
			snprintf( buf, maxlen, "%s/%s", path[i], base );
			execve( buf, argv, envp );
		}
	}
}

int count_lines( char *filename )
{
	int nol = 0;
	FILE *fl;
	char buf[1000];
	fl = fopen( filename, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", filename ); exit( 0 ); }
	while( ( fgets( buf, sizeof buf, fl ) ) != NULL ) nol++;
	fclose( fl );
	return nol;
}

int count_cols( char *filename, int row )
{
	int ncol = 0, i, n = 0;
	FILE *fl;
	char buf[1000], entry[16], *ln;
	if( row < 1 ) { printf( "\nInput Error: Row number %d < 1\n", row ); return( 0 ); }
	fl = fopen( filename, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", filename ); return( 0 ); }
	for( i = 0; i < row; i++ ) ln = fgets( buf, sizeof buf, fl );
	while( sscanf( ln, "%10s%n", entry, &n ) == 1 )
	{
		ncol++;
		ln += n;
	}
	fclose( fl );
	return ncol;
}

char *timestamp()
{
	time_t raw_time;
	struct tm *ptr_ts;
	char *datetime;
	datetime = ( char * ) malloc( 10 * sizeof( char ) );
	time( &raw_time );
	ptr_ts = localtime( &raw_time );
	// printf( "%s\n", asctime( ptr_ts ) );
	sprintf( datetime, "%02d:%02d:%02d", ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
	return( datetime );
}

char *datestamp()
{
	time_t raw_time;
	struct tm *ptr_ts;
	char *datetime;
	datetime = ( char * ) malloc( 16 * sizeof( char ) );
	time( &raw_time );
	ptr_ts = localtime( &raw_time );
	// printf( "%s\n", asctime( ptr_ts ) );
	sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
	return( datetime );
}

// You must free the result if result is non-NULL.
char *str_replace( char *orig, char *rep, char *with )
{
	char *result; // the return string
	char *ins;    // the next insert point
	char *tmp;    // varies
	int len_rep;  // length of rep
	int len_with; // length of with
	int len_front; // distance between rep and end of last rep
	int count;    // number of replacements
	if( !orig ) return NULL;
	if( !rep || !( len_rep = strlen( rep ) ) ) return NULL;
	if( !( ins = strstr( orig, rep ) ) ) return NULL;
	if( !with ) with = "";
	len_with = strlen( with );
	for( count = 0; ( tmp = strstr( ins, rep ) ); ++count )
		ins = tmp + len_rep;
	// first time through the loop, all the variable are set correctly
	// from here on,
	//    tmp points to the end of the result string
	//    ins points to the next occurrence of rep in orig
	//    orig points to the remainder of orig after "end of rep"
	tmp = result = ( char * ) malloc( strlen( orig ) + ( len_with - len_rep ) * count + 1 );
	if( !result ) return NULL;
	while( count-- )
	{
		ins = strstr( orig, rep );
		len_front = ins - orig;
		tmp = strncpy( tmp, orig, len_front ) + len_front;
		tmp = strcpy( tmp, with ) + len_with;
		orig += len_front + len_rep; // move to next "end of rep"
	}
	strcpy( tmp, orig );
	return result;
}

void tprintf( char const *fmt, ... )
{
	va_list ap;
	if( !quiet )
	{
		va_start( ap, fmt );
		vprintf( fmt, ap );
		va_end( ap );
		fflush( stdout );
	}
	va_start( ap, fmt );
	vfprintf( mads_output, fmt, ap );
	va_end( ap );
	fflush( mads_output );
}

int set_optimized_params( struct opt_data *op )
{
	struct calc_data *cd;
	struct param_data *pd;
	int i, k, bad_data = 0;
	cd = op->cd;
	pd = op->pd;
	if( cd->problem_type == CALIBRATE && pd->nFlgParam == 0 )
	{
		if( cd->calib_type == PPSD )
		{
			tprintf( "\nERROR: Partial parameter-space discretization (PPSD) is selected.\nHowever no parameters are flagged! Use optimization code value = 2 to flag model parameters.\n\n" );
			bad_data = 1;
		}
		if( cd->calib_type == IGPD )
		{
			tprintf( "\nERROR: Partial parameter-space discretization of initial guesses (IGPD) is selected.\nHowever no parameters are flagged! Use optimization code value = 2 to flag model parameters.\n\n" );
			bad_data = 1;
		}
	}
	if( cd->debug ) tprintf( "\n" );
	tprintf( "Number of optimized parameters = %d\n", pd->nOptParam );
	pd->var_index = ( int * ) malloc( pd->nOptParam * sizeof( int ) );
	pd->var_current = ( double * ) malloc( pd->nOptParam * sizeof( double ) );
	pd->var_best = ( double * ) malloc( pd->nOptParam * sizeof( double ) );
	for( k = i = 0; i < pd->nParam; i++ )
		if( pd->var_opt[i] == 1 || ( pd->var_opt[i] > 1 && cd->calib_type != PPSD ) )
		{
			if( cd->debug ) tprintf( "%-27s:%-6s: init %9g step %8.3g min %9g max %9g\n", pd->var_name[i], pd->var_id[i], pd->var[i], pd->var_dx[i], pd->var_min[i], pd->var_max[i] );
			pd->var_index[k++] = i;
		}
	if( cd->debug ) tprintf( "\n" );
	tprintf( "Number of flagged parameters = %d\n", pd->nFlgParam );
	if( cd->debug )
	{
		for( i = 0; i < pd->nParam; i++ )
			if( pd->var_opt[i] == 2 )
				tprintf( "%-27s:%-6s: init %9g step %6g min %9g max %9g\n", pd->var_name[i], pd->var[i], pd->var_id[i], pd->var_dx[i], pd->var_min[i], pd->var_max[i] );
	}
	pd->nIgnParam = 0;
	for( i = 0; i < pd->nParam; i++ )
		if( pd->var_name[i][0] == 0 ) pd->nIgnParam++;
	pd->nFixParam = pd->nParam - pd->nOptParam - pd->nFlgParam - pd->nExpParam - pd->nIgnParam;
	if( pd->nFixParam == 0 && cd->debug ) tprintf( "\nNO fixed parameters\n" );
	else
	{
		if( cd->debug ) tprintf( "\n" );
		tprintf( "Number of fixed parameters = %d\n", pd->nFixParam );
		if( cd->debug )
		{
			for( i = 0; i < pd->nParam; i++ )
				if( pd->var_opt[i] == 0 && pd->var_name[i][0] != 0 )
					tprintf( "%-27s:%-6s: %g\n", pd->var_name[i], pd->var_id[i],  pd->var[i] );
		}
	}
	tprintf( "Number of parameters with computational expressions (coupled or tied parameters) = %d\n", pd->nExpParam );
	if( bad_data ) return( -1 );
	return( 1 );
}

int map_obs( struct opt_data *op )
{
	struct obs_data od2;
	struct param_data *pd;
	struct obs_data *od;
	struct regul_data *rd;
	int i, k;
	od = op->od;
	rd = op->rd;
	pd = op->pd;
	if( rd->nRegul < 1 ) return( 1 );
	od->nTObs = od->nObs + rd->nRegul;
	od2.obs_id = char_matrix( od->nTObs, 50 );
	od2.obs_target = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od2.obs_weight = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od2.obs_min = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od2.obs_max = ( double * ) malloc( od->nTObs * sizeof( double ) );
	// od2.obs_current = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od2.obs_best = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od2.res = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od2.obs_log = ( int * ) malloc( od->nTObs * sizeof( int ) );
	for( i = 0; i < od->nObs; i++ )
	{
		strcpy( od2.obs_id[i], od->obs_id[i] );
		od2.obs_target[i] = od->obs_target[i];
		od2.obs_weight[i] = od->obs_weight[i];
		od2.obs_min[i] = od->obs_min[i];
		od2.obs_max[i] = od->obs_max[i];
		// od2.obs_current[i] = od->obs_current[i];
		od2.obs_best[i] = od->obs_best[i];
		od2.res[i] = od->res[i];
		od2.obs_log[i] = od->obs_log[i];
	}
	free_matrix( ( void ** ) od->obs_id, od->nObs );
	free( od->obs_target );
	free( od->obs_weight );
	free( od->obs_min );
	free( od->obs_max );
	// if( rd->nRegul == 0 ) free( od->obs_current ); // Already freed if there are regularization terms ...
	free( od->obs_best );
	free( od->res );
	free( od->obs_log );
	for( k = od->nObs, i = 0; i < rd->nRegul; i++, k++ ) // add regularization terms
	{
		strcpy( od2.obs_id[k], rd->regul_id[i] );
		od2.obs_target[k] = rd->regul_target[i];
		od2.obs_weight[k] = rd->regul_weight[i];
		od2.obs_min[k] = rd->regul_min[i];
		od2.obs_max[k] = rd->regul_max[i];
		od2.obs_log[k] = rd->regul_log[i];
	}
	od->obs_id = od2.obs_id;
	od->obs_target = od2.obs_target;
	od->obs_weight = od2.obs_weight;
	od->obs_min = od2.obs_min;
	od->obs_max = od2.obs_max;
	// od->obs_current = od2.obs_current;
	od->obs_best = od2.obs_best;
	od->res = od2.res;
	od->obs_log = od2.obs_log;
	for( i = pd->nAnalParam, k = 0; k < od->nObs; k++, i++ ) rd->regul_map_id[i] = od->obs_id[k]; // remap the pointer to the new observation id array
	// for( k = 0; k < rd->regul_nMap; k++ ) tprintf( "%s %g\n", rd->regul_map_id[k], rd->regul_map_val[k] );
	return( 1 );
}

int map_well_obs( struct opt_data *op )
{
	struct calc_data *cd;
	struct obs_data *od;
	struct obs_data *preds;
	struct well_data *wd;
	struct param_data *pd;
	struct regul_data *rd;
	int i, j, k, bad_data = 0;
	cd = op->cd;
	pd = op->pd;
	od = op->od;
	rd = op->rd;
	preds = op->preds;
	wd = op->wd;
	od->nCObs = od->nObs;
	for( i = 0; i < wd->nW; i++ )
	{
		if( wd->nWellObs[i] <= 0 )
			tprintf( "WARNING: Well %s has no observations!\n", wd->id[i] );
		for( j = 0; j < wd->nWellObs[i]; j++ )
			if( wd->obs_time[i][j] < DBL_EPSILON )
				tprintf( "WARNING: Observation #%d time for well %s is too small (%g); potential error in the input file %s!\n", j + 1, wd->id[i], wd->obs_time[i][j], op->filename );
		for( j = i + 1; j < wd->nW; j++ )
			if( strcmp( wd->id[i], wd->id[j] ) == 0 )
				tprintf( "WARNING: Well names #%i (%s) and #%i (%s) are identical!\n", i + 1, wd->id[i], j + 1, wd->id[j] );
	}
	if( od->nObs == 0 )
	{
		if( cd->problem_type != FORWARD && cd->problem_type != MONTECARLO )
		{ tprintf( "\nERROR: Number of calibration targets is equal to zero!\n\n" ); bad_data = 1; }
		else tprintf( "\nWARNING: Number of calibration targets is equal to zero!\n\n" );
	}
	if( bad_data ) return( 0 );
	if( cd->debug > 2 )
	{
		double d = ( pd->var[FLOW_ANGLE] * M_PI ) / 180;
		double alpha = cos( d );
		double beta = sin( d );
		tprintf( "\nCoordinate transformation of the observation points relative to the source:\n" );
		for( i = 0; i < wd->nW; i++ )
		{
			double x0 = wd->x[i] - pd->var[SOURCE_X];
			double y0 = wd->y[i] - pd->var[SOURCE_Y];
			double x = x0 * alpha - y0 * beta;
			double y = x0 * beta  + y0 * alpha;
			tprintf( "Well %10s %.15g %.15g : %.15g %.15g\n", wd->id[i], wd->x[i], wd->y[i], x, y );
		}
	}
	if( cd->debug ) tprintf( "\n" );
	tprintf( "Number of calibration targets = %d", od->nObs );
	if( preds->nTObs )
	{
		if( od->include_predictions ) tprintf( " (including predictions; observations with weight < 0)\n" );
		else tprintf( " (excluding predictions; observations with weight < 0)\n" );
	}
	else tprintf( "\n" );
	od->nTObs = od->nObs + rd->nRegul;
	if( rd->nRegul > 0 ) tprintf( "Number of total calibration targets (including regularization terms) = %d", od->nTObs );
	if( op->pd->nOptParam > op->od->nTObs ) { tprintf( "WARNING: Number of optimized model parameters is greater than number of observations and regularizations (%d>%d)\n", op->pd->nOptParam, op->od->nTObs ); }
	od->obs_target = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_current = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_best = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_id = char_matrix( od->nTObs, 50 );
	od->pred_id = char_matrix( od->nTObs, 50 );
	od->obs_weight = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_min = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_max = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->res = ( double * ) malloc( od->nTObs * sizeof( double ) );
	od->obs_log = ( int * ) malloc( od->nTObs * sizeof( int ) );
	od->obs_well_index = ( int * ) malloc( od->nTObs * sizeof( int ) );
	od->obs_time_index = ( int * ) malloc( od->nTObs * sizeof( int ) );
	for( k = i = 0; i < wd->nW; i++ )
		for( j = 0; j < wd->nWellObs[i]; j++ )
			if( ( wd->obs_weight[i][j] > DBL_EPSILON ) || ( od->include_predictions && wd->obs_weight[i][j] < -DBL_EPSILON ) )
			{
				od->obs_target[k] = wd->obs_target[i][j];
				od->obs_weight[k] = wd->obs_weight[i][j];
				od->obs_min[k] = wd->obs_min[i][j];
				od->obs_max[k] = wd->obs_max[i][j];
				od->obs_log[k] = wd->obs_log[i][j];
				od->obs_well_index[k] = i;
				od->obs_time_index[k] = j;
				sprintf( od->obs_id[k], "%s_%g", wd->id[i], wd->obs_time[i][j] );
				removeChars( od->obs_id[k], "()-*/" );
				if( cd->debug ) tprintf( "%s: %g weight %g\n", od->obs_id[k], od->obs_target[k], od->obs_weight[k] );
				k++;
			}
	for( i = 0; i < rd->nRegul; i++, k++ ) // add regularization targets
	{
		strcpy( od->obs_id[k], rd->regul_id[i] );
		od->obs_target[k] = rd->regul_target[i];
		od->obs_weight[k] = rd->regul_weight[i];
		od->obs_min[k] = rd->regul_min[i];
		od->obs_max[k] = rd->regul_max[i];
		od->obs_log[k] = rd->regul_log[i];
		od->obs_well_index[k] = -1;
		od->obs_time_index[k] = -1;
		if( cd->debug ) tprintf( "%s: %g weight %g", rd->regul_id[i], rd->regul_target[i], rd->regul_weight[i] );
	}
	if( cd->debug ) tprintf( "\n" );
	if( bad_data ) return( -1 );
	return( 1 );
}
