// MADS: Model Analyses & Decision Support (v1.1) 2011
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

#ifndef MADS_H
#define MADS_H

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

enum IO_TYPE {IO_TEXT = 0, IO_YAML, IO_XML };
enum PROBLEM_TYPE {UNKNOWN = -3, CHECK, CREATE, FORWARD, CALIBRATE, LOCALSENS, EIGEN, MONTECARLO, GLOBALSENS, ABAGUS, INFOGAP, POSTPUA, GLUE, BAYES };
enum CALIBRATION_TYPE {SIMPLE, PPSD, IGPD, IGRND};
enum GSA_TYPE {SOBOL, SALTELLI, MOAT};
enum OBJFUNC_TYPE {SSR = 0, SSDR, SSD0, SSDX, SSDA, SCR };
enum SOLUTION_TYPE {TEST = -2, EXTERNAL = -1, POINT = 0, PLANE = 1, PLANE3D = 2, BOX = 3, GAUSSIAN2D = 4, GAUSSIAN3D = 5, POINT_TRIANGLE_TIME = 6 };
enum LEVY_TYPE {NO_LEVY = 0, FULL_LEVY = 1, SYM_LEVY = 2};
#define NUM_ANAL_PARAMS_SOURCE 9
enum SOURCE_PARAM_TAGS {SOURCE_X = 0, SOURCE_Y, SOURCE_Z, SOURCE_DX, SOURCE_DY, SOURCE_DZ, FLUX, TIME_INIT, TIME_END };
#define NUM_ANAL_PARAMS_AQUIFER 17
enum AQUIFER_PARAM_TAGS { POROSITY = NUM_ANAL_PARAMS_SOURCE, RF, LAMBDA, FLOW_ANGLE, VX, VY, VZ, AX, AY, AZ, TSCALE_DISP, TSCALE_ADV, TSCALE_REACT, ALPHA, BETA, NLC0, NLC1 };

#define SHELL "/usr/bin/env tcsh -f -c"

int (*func_global)( double *x, void *data, double *f ); // global pointer to the model evaluation func (external or internal)
void tprintf( char const *fmt, ... );
char version_id[80];
extern const char *gitversion;
#ifdef MAIN
FILE *mads_output = NULL;
int quiet = 0;
#else
FILE *mads_output;
int quiet;
#endif

#define iswhite(c) ((c)== ' ' || (c)=='\t' || (c)=='\n' || (c)=='\r' )

#define COMPARE_EPSILON pow( FLT_EPSILON, (double) 1/2 ) // EPSILON FOR BOUND COMPARISON

struct opt_data // TODO class MADS (in C++)
{
	double phi; // current objective function (can be applied as termination criteria)
	int success; // success (can be applied as termination criteria)
	int global_success; // global success for TRIBES and SQUADS (can be applied as termination criteria)
	int counter; // current run counter (Monte Carlo / paranoid / retry / igrd / igpd ... )
	char *root; // problem name (filename root)
	char *filename; // problem filename
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
	struct anal_data *ad; // parameter subclass for the internal analytical simulations
	struct source_data *sd; // source parameters
	struct aquifer_data *qd; // aquifer parameters
	// TODO model of the MADS functions should be part of this class
};

struct calc_data // calculation parameters; TODO some of the flags can be boolean type
{
	int quit;
	int ioml; // YAML/XML input / output format
	int problem_type; // problem type: forward, calibration, ...
	int analysis_type; // calibration type: simple, igpd, ...
	int gsa_type; // global sensitivity analysis type: sobol, saltelli, moat, ...
	int paranoid; // paranoid calibration
	int num_sources; // number of contaminant sources (internal solutions)
	int num_source_params; // number of parameters for contaminant sources
	int num_aquifer_params; // number of aquifer parameters
	int *solution_type; // external / internal (box, ... )
	int levy;  // Levy dispersion YES/NO
	int nlmo; // number of LM calls
	int nreal; // number of realizations 
	int ireal; // execution of specific realization (case)
	int nretries; // number of paranoid retries
	int retry_ind; // retry evaluation counter
	int init_particles; // number of tribes (squads, tribes)
	int lm_niter; // number of iterations
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
	int parallel_type; // type of parallelization
	bool posix; // Write/Execute in a POSIX thread
	bool omp; // OpenMP parallelization
	int omp_threads; // number of threads for parallel execution
	int num_proc; // number of processors for parallel execution
	int proc_per_task; // number of processors per external job task for srun
	int lm_num_parallel_lambda; // number of parallel lambda runs in the case of Levenberg-Marquardt optimization (<= num_proc)
	int restart; // flag restart for parallel jobs
	bool bin_restart;
	int njob; // number of parallel jobs
	int energy; // starting energy for pssa particles
	int disp_tied;
	int disp_scaled;
	int time_step; // 1 => End time == Time step
	int save;
	int pargen;
	int obs_int;
	double lm_factor;
	int lm_acc; // to accelerate or npt LM
	int lm_indir;
	double lm_mu;
	int lm_nu;
	int lm_num_lambda_searches;
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
	char *restart_container; // directory name or zip filename with the restart files
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
	int tdebug;
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
	int nAnalParam; // global number of parameter for internal analytical solutions
	int nOptParam; // number of optimized parameters
	int nFlgParam; // number of special (flagged) parameters
	int nFixParam; // number of fixed parameters
	int nExpParam; // number of parameters with computational expressions (i.e. tied parameters)
	int nIgnParam; // number of ignored parameters
	int *var_index; // parameter index array
	char **var_name; // parameter identifier (name)
	char **var_id; // short parameter identifier (name) for the case of internal models only
	double *var; // parameter value (initial/final)
	int *var_opt; // parameter flag
	int *var_log; // flag for log transformation
	double *var_dx; // parameter increase/decrease step (discretization/derivatives)
	double *var_min; // parameter min value
	double *var_max; // parameter max value
	double *var_init_min; // parameter min value for initialization
	double *var_init_max; // parameter max value for initialization
	double *var_range; // parameter range: range = max - min
	double *var_current; // parameter value (current)
	double *var_best; // parameter value (current best)
	double *var_truth; // true parameter value
	int *param_expressions_index; // index array for parameters with computational expressions (i.e. tied parameters)
	void **param_expression; // math expressions for "tied" parameters
	gsl_vector *var_current_gsl; // current model parameters as GSL vector
};

struct regul_data // data structure for regularization terms
{
	int nRegul; // number of regularization terms
	void **regul_expression; // math expressions for the regularization terms
	char **regul_id; // regularization identifier (name)
	double *regul_target; // regularization value (target)
	double *regul_weight; // regularization weight
	int *regul_log; // flag for log transformation
	double *regul_min; // regularization min
	double *regul_max; // regularization max
	int regul_nMap; // number of mapped parameters and observations
	char **regul_map_id; // mapping of parameter and observation id's
	double *regul_map_val; // mapping of parameter and observation values
};

struct obs_data // data structure for observation data (EXTERNAL PROBLEM)
{
	int nTObs; // total number of observations and regularization terms; nGObs = nObs + nRegul
	int nObs; // number of observations (internal problem: nObs = nCObs; external problem: nObs > nCObs (all observations) )
	int nCObs; // total number of calibration targets observations with weight greater than zero
	int nPred; // number of performance criterion prediction
	int include_predictions;
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
	double *obs_scale; //gives the scale for the distribution
	double *obs_location; //gives the location for the distribution
	double *obs_alpha; //gives the stability index for the distribution
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
	int    *nWellObs; // number of observations at the wells
	double **obs_time; // observation time
	double **obs_target; // observation value (target)
	double **obs_weight; // observation weight
	int **obs_log; // log transformation
	double **obs_min; // observation min
	double **obs_max; // observation max
	double **obs_scale; //gives the scale of the distribution
	double **obs_location; //gives the location of the distribution
	double **obs_alpha; //gives the stability index for the distribution
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

struct anal_data
{
	int num_param;
	int debug;
	int scaling_dispersion;
	int time_step;
	double xe; // x coordinate; needed only for the functions during integration (can be a subclass)
	double ye; // y coordinate; needed only for the functions during integration (can be a subclass)
	double ze; // z coordinate; needed only for the functions during integration (can be a subclass)
	double te; // t coordinate; needed only for the functions during integration (can be a subclass)
	double *var; // optimized model parameters; needed only for the functions during integration (can be a subclass)
};

struct source_data
{
	int num_param;
	char **param_id;
	char **param_name;
};

struct aquifer_data
{
	int num_param;
	char **param_id;
	char **param_name;
};

struct class_data
{
	int proc_problem;
	int proc_solution;
	int proc_parameters;
	int proc_executable;
	int proc_temp_ins;
	int proc_grid;
	int proc_time;
	char **class_id;
};

// mads.c
int optimize_lm( struct opt_data *op ); // LM (Levenberg-Marquardt) optimization
int optimize_pso( struct opt_data *op ); // PSO optimization
int eigen( struct opt_data *op, double *f_x, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar ); // Eigen analysis
void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op, int debug ); // Random sampling
void print_results( struct opt_data *op, int verbosity ); // Print final results
void save_results( int final, char *filename, struct opt_data *op, struct grid_data *gd ); // Save final results
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var );
void ave_sorted( double data[], int n, double *ave, double *ep );
char *timestamp(); // create time stamp
char *datestamp(); // create date stamp
void mads_quits( char *root );
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
int load_problem_text( char *filename, int argn, char *argv[], struct opt_data *op );
int save_problem_text( char *filename, struct opt_data *op );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc2( char *filename, char *filename2, struct opt_data *op );
void compute_btc( char *filename, struct opt_data *op );
// mads_io_external.c
int load_pst( char *filename, struct opt_data *op );
int check_ins_obs( int nobs, char **obs_id, int *check, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, int *check, char *fn_in_t, char *fn_in_d, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
//io.c
int Ftest( char *filename );
#else
extern int quiet;
extern int mads_output_file;
#endif /* MADS_H */
