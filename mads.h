// MADS: Model Analyses & Decision Support (v1.1) 2011
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
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

#ifndef MADS_H
#define MADS_H

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

enum PROBLEM_TYPE {UNKNOWN = -3, CHECK, CREATE, FORWARD, CALIBRATE, LOCALSENS, EIGEN, MONTECARLO, GLOBALSENS, ABAGUS, INFOGAP, POSTPUA, GLUE };
enum CALIBRATION_TYPE {SIMPLE, PPSD, IGPD, IGRND};
enum OBJFUNC_TYPE {SSR = 0, SSDR, SSD0, SSDX, SSDA, SCR };
enum SOLUTION_TYPE {TEST = -2, EXTERNAL = -1, POINT = 0, PLANE = 1, PLANE3D = 2, BOX = 3 };
#define NUM_ANAL_PARAMS 19
#define NUM_ANAL_PARAMS_SOURCE 9
enum PARAM_TAGS {SOURCE_X = 0, SOURCE_Y, SOURCE_Z, SOURCE_DX, SOURCE_DY, SOURCE_DZ, C0, TIME_INIT, TIME_END, POROSITY, KD, LAMBDA, FLOW_ANGLE, VX, VY, VZ, AX, AY, AZ };

int (*func_global)( double *x, void *data, double *f ); // global pointer to the model evaluation func (external or internal)
void tprintf( char const *fmt, ... );
FILE *mads_output;
int quiet;

#define COMPARE_EPSILON pow( FLT_EPSILON, (double) 1/2 ) // EPSILON FOR BOUND COMPARISON

struct opt_data // TODO class MADS (in C++)
{
	double phi; // current objective function (can be applied as termination criteria)
	int success; // success (can be applied as termination criteria)
	int global_success; // global success for TRIBES and SQUADS (can be applied as termination criteria)
	int counter; // current run counter (Monte Carlo / paranoid / retry / igrd / igpd ... )
	char *root; // problem name (filename root)
	char *label; // problem label (for output file generation)
	char *datetime_stamp; // date & time of the simulation
	FILE *f_ofe; // runtime output file with current best objective function
	struct param_data *pd; // parameters subclass
	struct regul_data *rd; // regularization data
	struct obs_data *od; // observations subclass
	struct obs_data *preds; // predictions subclass (special observations)
	struct well_data *wd; // well-data subclass
	struct calc_data *cd; // calculation parameters subclass
	struct grid_data *gd; // grid subclass to compute model predictions
	struct extrn_data *ed; // parameter subclass for external simulations
	struct anal_data *ad;
	// TODO model of the MADS functions should be part of this class
};

struct calc_data // calculation parameters; TODO some of the flags can be boolean type
{
	int problem_type; // problem type: forward, calibration, ...
	int calib_type; // calibration type: simple, igpd, ...
	int paranoid; // paranoid calibration
	int num_solutions; // number of internal solutions
	int *solution_type; // external / internal (box, ... )
	int nlmo; // number of LM calls
	int nreal; // number of realizations 
	int ireal; // execution of specific realization (case)
	int nretries; // number of paranoid retries
	int retry_ind; // retry evaluation counter
	int init_particles; // number of tribes (squads, tribes)
	int niter; // number of iterations
	int neval; // current number of evaluations (can be applied as termination criteria)
	int njac; // current number of jacobian evaluations
	int maxeval; // maximum number of evaluations (termination criteria)
	int lmstandalone; // flag standalone LM run (yes/no)
	int squads;
	int seed; //random seed
	int seed_init; //random seed
	int test_func; // test function id
	int test_func_dim; // test function dimensionality
	int test_func_npar; // test function number of parameters
	int test_func_nobs; // test function number of observations
	int num_proc; // number of processors for parallel execution
	int restart; // flag restart for parallel jobs
	int njob; // number of parallel jobs
	int energy; // starting energy for pssa particles
	int disp_tied;
	int disp_scaled;
	int save;
	int pargen;
	double lm_factor;
	int lm_acc; // to accelerate or npt LM
	int lm_indir;
	double lm_mu;
	int lm_nu;
	int lm_nlamof;
	int lm_njacof;
	double lm_h;
	double lm_ratio;
	double lm_ofdecline;
	double lm_error;
	double c_background;
	char *solution_id; // solution identifier (name)
	char *opt_method; // optimization method identifier
	char *smp_method; // sampling method identifier
	char *paran_method; // sampling method identifier for paranoid runs (same options as for smp_method)
	char *mydir; // working directory
	char *mydir_hosts; // directories for the parallel jobs
	char **paral_hosts; // parallel host identifier (name)
	time_t time_infile; // time of the input problem file (*.mads)
	char *datetime_infile; // date & time of the input problem file (*.mads) (equivalent to time_infile but in different format)
	char *restart_zip_file; // filename of the zip restart file
	char *infile; // old results file from pssa to be read in to initialize kdtree
	int resultscase; // read specific case
	char *resultsfile; // read existing results file
	double phi_cutoff; // objective function cutoff value (termination criteria)
	int obsrange; // flag; observations are within predefined ranges (termination criteria)
	double obserror; // absolute error from the known 'true' observation (termination criteria)
	double parerror; // absolute error from the known 'true' model parameters (if known; e.g. in the case of test functions; termination criteria)
	double sindx; // increments to compute model parameter gradients in the sin transformed parameter space
	double lindx; // increments to compute model parameter gradients in the linear (not sin) transformed parameter space
	double pardx; // parameter space discretization step for pso, tribes, squads and abagus
	double pardomain; // parameter space domain size for test functions
	double obsdomain; // observation space domain size for info-gap analyses
	double obsstep; // observation space domain step for info-gap analyses
	int debug; // various debug / verbosity levels
	int fdebug;
	int ldebug;
	int pdebug;
	int mdebug;
	int odebug;
	int tpldebug;
	int insdebug;
	int pardebug;
	int lm_eigen; // flag for eigen analysis
	int sintrans; // flag for sin transformation
	int plogtrans; // flag for log transformation of all parameters
	int ologtrans; // flag for log transformation of all observations
	int oweight; // flag for weights of all observations
	int pderiv; // internal flag for computation of parameter derivatives
	int oderiv; // internal flag for computation of observation derivatives
	int objfunc_type; // objective function type
	int check_success; // flag to check success
	int compute_phi; // flag to compute objective function
	double xe; // x coordinate; needed only for the functions during integration (can be a subclass)
	double ye; // y coordinate; needed only for the functions during integration (can be a subclass)
	double ze; // z coordinate; needed only for the functions during integration (can be a subclass)
	double te; // t coordinate; needed only for the functions during integration (can be a subclass)
	double *var; // optimized model parameters; needed only for the functions during integration (can be a subclass)
};

struct param_data // data structure for model parameters
{
	int nParam; // number of parameters
	int nOptParam; // number of optimized parameters
	int nFlgParam; // number of special (flagged) parameters
	int nFixParam; // number of fixed parameters
	int nExpParam; // number of parameters with computational expressions (i.e. tied parameters)
	int *var_index; // parameter index array
	char **var_id; // parameter identifier (name)
	char **var_id_short; // short parameter identifier (name) for the case of internal models only
	double *var; // parameter value (initial/final)
	int *var_opt; // parameter flag
	int *var_log; // flag for log transformation
	double *var_dx; // parameter increase/decrease step (discretization/derivatives)
	double *var_min; // parameter min value
	double *var_max; // parameter max value
	double *var_range; // parameter range: range = max - min
	double *var_current; // parameter value (current)
	double *var_best; // parameter value (current best)
	double *var_truth; // true parameter value
	int *param_expressions_index; // index array for parameters with computational expressions (i.e. tied parameters)
	void **param_expressions; // math expressions for "tied" parameters
	gsl_vector *var_current_gsl; // current model parameters as GSL vector
};

struct regul_data // data structure for model parameters
{
	int nRegul; // number of regularization terms
	void **regul_expressions; // math expressions for the regularization terms
	char **regul_id; // regularization identifier (name)
	double *regul_target; // regularization value (target)
	double *regul_weight; // regularization weight
	int *regul_log; // flag for log transformation
	double *regul_min; // regularization min
	double *regul_max; // regularization max
};

struct obs_data // data structure for observation data (EXTERNAL PROBLEM)
{
	int nTObs; // total number of observations and regularization terms; nGObs = nObs + nRegul
	int nObs; // number of observations (internal problem: nObs = nCObs; external problem: nObs > nCObs (all observations) )
	int nCObs; // total number of calibration targets observations with weight greater than zero
	int nPred; // number of performance criterion prediction
	// observations
	char **obs_id; // observation identifier (name) 
	double *obs_target; // observation value (target)
	double *obs_current; // current model predicted observation; NOTE: redundant
	double *obs_best; // current model predicted observation; NOTE: redundant
	int *obs_index;
	int *obs_log; // flag for log transformation
	int *obs_well_index; // well index for observations
	int *obs_time_index; // time index for observations
	double *obs_weight; // observation weight
	double *obs_min; // observation min
	double *obs_max; // observation max
	double *res; // current residual: res = obs_target - obs_current
	gsl_vector *obs_current_gsl; // current model predicted observation as GSL vector; NOTE: redundant
	// predictions
	char **pred_id; // performance criterion identifier (name)
	double *pred_crit; // value of prediction criterion
	int *pred_well_index; // well index for prediction
	int *pred_time_index; // time index for prediction
};

struct well_data // data structure for well data (INTERNAL PROBLEM)
{
	int    nW; // number of wells
	char   **id; // well ID's
	double *x; // well coordinates
	double *y;
	double *z1;
	double *z2;
	double *xa; // adjusted well coordinates
	double *ya;
	double *za1;
	double *za2;
	int    *nWellObs; // number of observations at the wells
	double **obs_time; // observation time
	double **obs_target; // observation value (target)
	int **obs_log; // log transformation
	double **obs_weight; // observation weight
	double **obs_min; // observation min
	double **obs_max; // observation max
};

struct grid_data // data structure for model predictions along a grid
{
	int nx;
	int ny;
	int nz;
	int nt;
	double time;
	double min_x;
	double min_y;
	double min_z;
	double min_t;
	double max_x;
	double max_y;
	double max_z;
	double max_t;
	double dx;
	double dy;
	double dz;
	double dt;
};

struct extrn_data // data structure for external problem
{
	char *cmdline; // command line for external execution (TODO we should allow for more than one commands)
	int ntpl; // number of template files
	char **fn_tpl; // template filename
	char **fn_out; // model input filename associated with the template file
	int nins; // number of instruction files
	char **fn_ins; // instruction filename
	char **fn_obs; // model output filename associated with the instruction file
};

struct gsens_data // global sensitivity analysis data structure
{
	double **var_a_lhs;	// sample a for global sensitivity analysis
	double **var_b_lhs;	// sample b for global sensitivity analysis
	double f_hat_0;		// total output mean
	double *f_a;		// sample a phis
	double *f_b;		// sample b phis
	double **fmat_a;	// matrix of phis - rows \theta_i^a; column \theta_{~i}^b
	double **fmat_b;	// matrix of phis - rows \theta_i^b; column \theta_{~i}^a
	double D_hat_t;		// total output variance
	double *D_hat; 		// component output variance (\hat{D}_i)
	double *D_hat_n; 	// not component output variance (\hat{D}_{~i})
	double ep;          // absolute first moment
};

struct anal_data
{
	int num_param;
	int debug;
	double xe; // x coordinate; needed only for the functions during integration (can be a subclass)
	double ye; // y coordinate; needed only for the functions during integration (can be a subclass)
	double ze; // z coordinate; needed only for the functions during integration (can be a subclass)
	double te; // t coordinate; needed only for the functions during integration (can be a subclass)
	double var[NUM_ANAL_PARAMS]; // optimized model parameters; needed only for the functions during integration (can be a subclass)
};

// mads.c
int optimize_lm( struct opt_data *op ); // LM (Levenberg-Marquardt) optimization
int optimize_pso( struct opt_data *op ); // PSO optimization
int eigen( struct opt_data *op, double *f_x, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar ); // Eigen analysis
void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op, int debug ); // Random sampling
void print_results( struct opt_data *op, int verbosity ); // Print final results
void save_final_results( char *filename, struct opt_data *op, struct grid_data *gd ); // Save final results
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var );
void ave_sorted( double data[], int n, double *ave, double *ep );
char *timestamp(); // create time stamp
char *datestamp(); // create date stamp
// mads_func.c
int func_extrn( double *x, void *data, double *f );
int func_extrn_write( int ieval, double *x, void *data );
int func_extrn_read( int ieval, void *data, double *f );
int func_extrn_check_read( int ieval, void *data );
int func_extrn_r( double *x, void *data, double *f );
int func_intrn( double *x, void *data, double *f );
void func_levmar( double *x, double *f, int m, int n, void *data );
void func_dx_levmar( double *x, double *f, double *jacobian, int m, int n, void *data );
int func_dx( double *x, double *f_x, void *data, double *jacobian );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );
// mads_mem.c
char **char_matrix( int maxCols, int maxRows );
float **float_matrix( int maxCols, int maxRows );
double **double_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
void zero_double_matrix( double **matrix, int maxCols, int maxRows );
void *malloc_check( const char *what, size_t n );
char *white_trim( char *x );
void white_skip( char **s );
// mads_io.c
int parse_cmd( char *buf, struct calc_data *cd );
int load_problem( char *filename, int argn, char *argv[], struct opt_data *op );
int save_problem( char *filename, struct opt_data *op );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc2( char *filename, char *filename2, struct opt_data *op );
void compute_btc( char *filename, struct opt_data *op );
// mads_io_external.c
int load_pst( char *filename, struct opt_data *op );
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );

#endif /* MADS_H */
