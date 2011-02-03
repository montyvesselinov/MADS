#include <gsl/gsl_vector.h>
#include <stdio.h>

enum PROBLEM_TYPE {UNKNOWN = -2, CREATE, FORWARD, CALIBRATE, LOCALSENS, EIGEN, MONTECARLO, GLOBALSENS, ABAGUS, POSTPUA };
enum CALIBRATION_TYPE {SIMPLE, PPSD, IGPD, IGRND};
enum OBJFUNC_TYPE {SSR = 0, SSDR, SSD0, SSDA, SCR };
enum SOLUTION_TYPE {TEST = -2, EXTERNAL = -1, POINT = 0, PLANE = 1, PLANE3D = 2, BOX = 3 };
enum PARAM_TAGS {SOURCE_X = 0, SOURCE_Y, SOURCE_Z, SOURCE_DX, SOURCE_DY, SOURCE_DZ, C0, TIME_INIT, TIME_END, POROSITY, KD, LAMBDA, FLOW_ANGLE, VX, VY, VZ, AX, AY, AZ };

struct problem_data
{
	int narg;
	char **arg;
};

int (*func)( double *x, void *data, double *f );

struct calc_data
{
	int problem_type;
	int calib_type;
	int paranoid;
	int solution_type;
	int nlmo;
	int nreal;
	int nretries;
	int ntribe;
	int niter;
	int maxeval;
	int eval;
	int standalone;
	int seed;
	int test;
	int dim;
	int num_proc;
	int njob;
	int ncase;
	char *solution_id;
	char *opt_method;
	char *smp_method;
	char *paran_method;
	char **paral_dirs;
	char **paral_hosts;
	char *infile; //! old results file from pssa to be read in to initialize kdtree
	double phi_cutoff;
	double sindx;
	int debug;
	int fdebug;
	int ldebug;
	int pdebug;
	int mdebug;
	int odebug;
	int leigen;
	int sintrans;
	int plogtrans;
	int ologtrans;
	int oweight;
	int pderiv;
	int oderiv;
	int tpldebug;
	int insdebug;
	int objfunc;
	int check_success;
	int compute_phi;
	int energy; //! starting energy for pssa particles
	double xe; /** needed only for the functions during integration */
	double ye;
	double ze;
	double te;
	double *var;
};

struct param_data
{
	int nParam;
	int nOptParam;
	int nFlgParam;
	int *var_index;
	char **var_id;
	double *var;
	int *var_opt;
	int *var_log;
	double *var_dx;
	double *var_max;
	double *var_min;
	double *var_range;
	double *var_current;
	double *var_best;
	gsl_vector *var_current_gsl;
};

struct obs_data
{
	int nObs;
	int nTObs;
	char **obs_id;
	double *obs_target;
	double *obs_current;
	int *obs_log;
	int *well_index;
	int *time_index;
	double *obs_weight;
	double *obs_max;
	double *obs_min;
	double *res;
	gsl_vector *obs_current_gsl;
};

struct well_data
{
	int    nW;
	char   **id;
	double *x;
	double *y;
	double *z1;
	double *z2;
	double *xa;
	double *ya;
	double *za1;
	double *za2;
	int    *nWellObs;
	double **obs_time;
	double **obs_target;
	int **obs_log;
	double **obs_weight;
	double **obs_min;
	double **obs_max;
};

struct grid_data
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

struct opt_data
{
	double phi;
	int success;
	int s;
	char *root;
	FILE *f_ofe;
	struct param_data *pd;
	struct obs_data *od;
	struct well_data *wd;
	struct calc_data *cd;
	struct extrn_data *ed;
};

struct extrn_data
{
	char *cmdline;
	int ntpl;
	char **fn_tpl;
	char **fn_out;
	int nins;
	char **fn_ins;
	char **fn_obs;
};

//! global sensitivity analysis data structure
struct gsens_data
{
	double **var_a_lhs;	//! sample a for global sensitivity analysis
	double **var_b_lhs;	//! sample b for global sensitivity analysis
	double f_hat_0;		//! total output mean
	double *f_a;		//! sample a phis
	double *f_b;		//! sample b phis
	double **fmat_a;	//! matrix of phis - rows \theta_i^a; column \theta_{~i}^b
	double **fmat_b;	//! matrix of phis - rows \theta_i^b; column \theta_{~i}^a
	double D_hat_t;		//! total output variance
	double *D_hat; 		//! component output variance (\hat{D}_i)
	double *D_hat_n; 	//! not component output variance (\hat{D}_{~i})
	double ep;		//! absolute first moment
};
