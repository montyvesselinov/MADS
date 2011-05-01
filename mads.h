#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

enum PROBLEM_TYPE {UNKNOWN = -2, CREATE, FORWARD, CALIBRATE, LOCALSENS, EIGEN, MONTECARLO, GLOBALSENS, ABAGUS, INFOGAP, POSTPUA };
enum CALIBRATION_TYPE {SIMPLE, PPSD, IGPD, IGRND};
enum OBJFUNC_TYPE {SSR = 0, SSDR, SSD0, SSDA, SCR };
enum SOLUTION_TYPE {TEST = -2, EXTERNAL = -1, POINT = 0, PLANE = 1, PLANE3D = 2, BOX = 3 };
enum PARAM_TAGS {SOURCE_X = 0, SOURCE_Y, SOURCE_Z, SOURCE_DX, SOURCE_DY, SOURCE_DZ, C0, TIME_INIT, TIME_END, POROSITY, KD, LAMBDA, FLOW_ANGLE, VX, VY, VZ, AX, AY, AZ };

int (*func)( double *x, void *data, double *f ); // global pointer to the model evaluation func (external or internal)

struct opt_data // TODO class MADS (in C++)
{
	double phi; // current objective function (can be applied as termination criteria)
	int success; // success (can be applied as termination criteria)
	int counter; // current run counter (Monte Carlo / paranoid / retry / igrd / igpd ... )
	char *root; // problem name (filename root)
	char *label; // problem label (for output file generation)
	char *datetime_stamp; // date & time of the simulation
	FILE *f_ofe; // runtime output file with current best objective function
	struct param_data *pd; // parameters subclass
	struct obs_data *od; // observations subclass
	struct obs_data *preds; // predictions subclass (special observations)
	struct well_data *wd; // well-data subclass
	struct calc_data *cd; // calculation parameters subclass
	struct grid_data *gd; // grid subclass to compute model predictions
	struct extrn_data *ed; // parameter subclass for external simulations
	// model of the MADS functions should be part of this class
};

struct calc_data // calculation parameters; TODO some of the flags can be boolean type
{
	int problem_type; // problem type: forward, calibration, ...
	int calib_type; // calibration type: simple, igpd, ...
	int paranoid; // paranoid calibration
	int solution_type; // external / internal (box, ... )
	int nlmo; // number of LM calls
	int nreal; // number of realizations 
	int ireal; // execution of specific realization (case)
	int nretries; // number of paranoid retries
	int retry_ind; // retry evaluation counter
	int ntribe; // number of tribes (squads, tribes)
	int niter; // number of iterations
	int neval; // current number of evaluations (can be applied as termination criteria)
	int maxeval; // maximum number of evaluations (termination criteria)
	int standalone; // flag standalone LM run (yes/no)
	int seed; //random seed
	int test_func; // test function id
	int test_func_dim; // test function dimensionality
	int num_proc; // number of processors for parallel execution
	int restart; // flag restart for parallel jobs
	int njob; // number of parallel jobs
	int energy; // starting energy for pssa particles
	double lmfactor;
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
	double phi_cutoff; // objective function cutoff value (termination criteria)
	double sindx; // discretization step / increments for model parameters in the sin transformed space
	int debug; // various debug / verbosity levels
	int fdebug;
	int ldebug;
	int pdebug;
	int mdebug;
	int odebug;
	int tpldebug;
	int insdebug;
	int pardebug;
	int leigen; // flag for eigen analysis
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
	int *var_index; // parameter index array
	char **var_id; // parameter identifier (name)
	double *var; // parameter value (initial/final)
	int *var_opt; // parameter flag
	int *var_log; // flag for log transformation
	double *var_dx; // parameter increase/decrease step (discretization/derivatives)
	double *var_min; // parameter min value
	double *var_max; // parameter max value
	double *var_range; // parameter range: range = max - min
	double *var_current; // parameter value (current)
	double *var_best; // parameter value (current best)
	gsl_vector *var_current_gsl; // current model parameters as GSL vector
};

struct obs_data // data structure for observation data (EXTERNAL PROBLEM)
{
	int nObs; // number of observation for calibration: nObs = nTObs - nPreds
	int nTObs; // total number of observations: nTObs = nObs + nPreds
	int nPreds; // number of performance criterion prediction: nPreds = nTObs - nObs
	char **obs_id; // observation identifier (name) 
	char **preds_id; // performance criterion identifier (name)
	double *obs_target; // observation value (target)
	double *pred_crit; // value of prediction criterion
	double *obs_current; // current model predicted observation; NOTE: redundant
	int *obs_log; // flag for log transformation
	int *well_index; // well index for observations
	int *time_index; // time index for observations
	int *pwell_index; // well index for prediction
	int *ptime_index; // time index for prediction
	double *obs_weight; // observation weight
	double *obs_min; // observation min
	double *obs_max; // observation max
	double *res; // current residual: res = obs_target - obs_current
	gsl_vector *obs_current_gsl; // current model predicted observation as GSL vector; NOTE: redundant
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

// mads.c
int optimize_lm( struct opt_data *op ); // LM (Levenberg-Marquardt) optimization
int optimize_pso( struct opt_data *op ); // PSO optimization
int eigen( struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar ); // Eigen analysis
void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op ); // Random sampling
void print_results( struct opt_data *op ); // Print final results
void save_results( char *filename, struct opt_data *op, struct grid_data *gd ); // Save final results
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
void compute_btc2( char *filename, struct opt_data *op );
void compute_btc( char *filename, struct opt_data *op );
// pesting/pesting.c
int load_pst( char *filename, struct opt_data *op );
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );

