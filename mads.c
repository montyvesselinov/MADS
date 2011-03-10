#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include "levmar-2.5/levmar.h"
#include "mads.h"

#define FIT(i) gsl_vector_get(solver->x, i)
#define MAX(X,Y) ( ((X) > (Y)) ? (X) : (Y) )

/* Functions here */
int optimize_lm( struct opt_data *op );
int eigen( struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar );
int optimize_pso( struct opt_data *op );
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var );
void ave_sorted( double data[], int n, double *ave, double *ep );
void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op );
void print_results( struct opt_data *op );
void save_results( char *filename, struct opt_data *op, struct grid_data *gd );
char *timestamp();
char *datestamp();

/* Functions elsewhere */
int pssa( struct opt_data *op );
int postpua( struct opt_data *op );
int pso_tribes( struct opt_data *op );
int pso_std( struct opt_data *op );
int parse_cmd( char *buf, struct calc_data *cd );
int load_problem( char *filename, int argn, char *argv[], struct opt_data *op, struct calc_data *cd, struct param_data *pd, struct obs_data *od, struct well_data *wd, struct grid_data *gd, struct extrn_data *ed );
int load_pst( char *filename, struct calc_data *cd, struct param_data *pd, struct obs_data *od, struct extrn_data *ed );
int save_problem( char *filename, struct calc_data *cd, struct param_data *pd, struct obs_data *od, struct well_data *wd, struct grid_data *gd, struct extrn_data *ed );
int check_ins_obs( int nobs, char **obs_id, double *obs, char *fn_in_t, int debug );
int check_par_tpl( int npar, char **par_id, double *par, char *fn_in_t, int debug );
int ins_obs( int nobs, char **obs_id, double *obs, double *check, char *fn_in_t, char *fn_in_d, int debug );
int par_tpl( int npar, char **par_id, double *par, char *fn_in_t, char *fn_out, int debug );
void compute_grid( char *filename, struct calc_data *cd, struct grid_data *gd );
void compute_btc( char *filename, struct opt_data *od, struct grid_data *gd );
void compute_btc2( char *filename, struct opt_data *od, struct grid_data *gd );
int Ftest( char *filename );
FILE *Fread( char *filename );
FILE *Fwrite( char *filename );
FILE *Fappend( char *filename );
char *Fdatetime( char *filename, int debug );
time_t Fdatetime_t( char *filename, int debug );
int lm_opt( int func( double x[], void *data, double f[] ), int func_dx( double *x, double *f_x, void *data, double *jacobian ), void *data,
			int nObs, int nParam, int nsig, double eps, double delta, int max_eval, int max_iter,
			int iopt, double parm[], double x[], double *phi, double f[],
			double jacobian[], int nJacobian, double jacTjac[], int *infer );
int zxssqch( int func( double x[], void *, double f[] ), void *func_data,
			 int m, int n, int nsig, double eps, double delta, int maxfn,
			 int iopt, double parm[], double x[], double *phi, double f[],
			 double xjac[], int ixjac, double xjtj[], int *infer );
int lm_gsl( gsl_vector *x, struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *covar );
double epsilon();
void lhs_imp_dist( int nvar, int npoint, int d, int *seed, double x[] );
void lhs_center( int nvar, int npoint, int *seed, double x[] );
void lhs_edge( int nvar, int npoint, int *seed, double x[] );
void lhs_random( int nvar, int npoint, int *seed, double x[] );
void smp_random( int nvar, int npoint, int *seed, double x[] );
int get_seed( );
double **double_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
void zero_double_matrix( double **matrix, int maxCols, int maxRows );
int func_gsl_dx( const gsl_vector *x, void *data, gsl_matrix *J );
int func_gsl_xdx( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J );
double func_gsl_deriv( double x, void *data );
int func_gsl_deriv_dx( const gsl_vector *x, void *data, gsl_matrix *J );
char **char_matrix( int maxCols, int maxRows );
int create_mprun_dirs( int nDir, char **dirs );
int delete_mprun_dirs( int nDir, char **dirs );
int func_extrn( double *x, void *data, double *f );
int func_intrn( double *x, void *data, double *f );
void func_levmar( double *x, double *f, int m, int n, void *data );
void func_dx_levmar( double *x, double *f, double *jacobian, int m, int n, void *data );
int func_dx( double *x, double *f_x, void *data, double *jacobian );
double func_solver( double x, double y, double z1, double z2, double t, void *data );
double func_solver1( double x, double y, double z, double t, void *data );
int func_extrn_write( int ieval, double *x, void *data );
int func_extrn_exec_serial( int ieval, void *data );
int func_extrn_read( int ieval, void *data, double *f );
void Transform( double *v, void *data, double *vt );
void DeTransform( double *v, void *data, double *vt );
int delete_mprun_outputs( int nDir, void *data );
int mprun( int nJob, void *data );
char *dir_hosts( void *data, char *timedate_stamp );
char *white_trim( char *x );

int main( int argn, char *argv[] )
{
	int i, j, k, ier, npar, status, status_global, status_all, count, debug_level, predict = 0, compare, eval_total, bad_data;
	double c, err, phi, phi_min, *orig_params, *opt_params, *var_lhs, *var_a_lhs, *var_b_lhs;
	struct calc_data cd;
	struct param_data pd;
	struct obs_data od;
	struct well_data wd;
	struct extrn_data ed;
	struct grid_data gd;
	struct opt_data op;
	struct gsens_data gs;
	char filename[80], root[80], extension[80], buf[255], *dot, *cwd;
	int ( *f )( struct opt_data * op );
	char *host, *nodelist, *hostlist, *proclist, *lsblist, *beowlist;
	FILE *in, *out, *out2;
	time_t t1, t2, dt;
	pid_t pid;
	struct tm *ptr_ts;
	t1 = time( NULL );
	op.datetime_stamp = datestamp();
	op.pd = &pd;
	op.od = &od;
	op.wd = &wd;
	op.cd = &cd;
	op.ed = &ed;
	cd.eval = 0;
	cd.nlmo = 0;
	cd.standalone = 1; // LM variable
	cd.pderiv = cd.oderiv = -1; // Do not compute derivatives
	op.phi = HUGE;
	op.success = 0;
	op.f_ofe = NULL;
	printf( "MADS: Model Analyses & Decision Support (v1.1) 2011\n" );
	printf( "---------------------------------------------------\n" );
	if( argn < 2 )
	{
		printf( "USAGE: %s problem_name [ keywords | options ]       OR\n", argv[0] );
		printf( "       %s MADS_input_file [ keywords | options ]    OR\n", argv[0] );
		printf( "       %s PEST_input_file [ keywords | options ] (MADS is compatible with PEST control, template and instruction files)\n\n", argv[0] );
		printf( "problem_name:         name of the solved problem; MADS_input_file named problem_name.mads is expected\n" );
		printf( "MADS_input_file:      problem input file in MADS format (typically *.mads)\n" );
		printf( "PEST_input_file:      problem input file in PEST format (PEST control file; *.pst)\n\n" );
		printf( "keywords & options (can be provided in any order):\n\n" );
		printf( "problem type keywords:\n" );
		printf( "   create             - create calibration problem based on provided model parameters\n" );
		printf( "   forward            - forward model run\n" );
		printf( "   calibrate          - calibration run [default]\n" );
		printf( "   montecarlo         - Monte Carlo analysis\n" );
		printf( "   gsens              - global sensitivity analysis\n" );
		printf( "   lsens              - local sensitivity analysis (standalone or at the end of the calibration)\n" );
		printf( "   eigen              - local eigensystem analysis (standalone or at the end of the calibration)\n" );
		printf( "   abagus             - Agent-Based Global Uncertainty & Sensitivity Analysis (ABAGUS)\n" );
		printf( "   postpua            - predictive uncertainty analysis of sampling results (currently works with ABAGUS output files only)\n" );
		printf( "\ncalibration method keywords (select one):\n" );
		printf( "   single             - single calibration using initial guesses provided in the input file [default]\n" );
		printf( "   igrnd              - sequential calibrations using a set of random initial values (number of realizations defined by real=X)\n" );
		printf( "   igpd               - sequential calibrations using a set of discretized initial values (discretization defined in the input file)\n" );
		printf( "   ppsd               - sequential calibrations using partial parameter space discretization (PPSD) method\n" );
		printf( "\ncalibration termination criteria:\n" );
		printf( "   eval=[integer]     - terminate calibration if functional evaluations exceed the predefined value [default eval=5000]\n" );
		printf( "   cutoff=[real]      - terminate calibration if the objective function is below the cutoff value [default cutoff=0]\n" );
		printf( "   success            - terminate calibration if model predictions are within predefined calibration ranges\n" );
		printf( "\ncalibration options:\n" );
		printf( "   retry=[integer]    - number of optimization retries [default retry=0]\n" );
		printf( "   iter=[integer]     - number of Levenberg-Marquardt iterations [default iter=50]\n" );
		printf( "   tribe=[integer]    - number of APSO tribes [default tribe=10+sqrt(Number_of_parameters)]\n" );
		printf( "   leigen|eigen       - eigen analysis of the final optimized solution\n" );
		printf( "\noptimization method (opt=[string]; various combinations are possible, e.g. pso_std_lm_gsl):\n" );
		printf( "   opt=lm             - Local Levenberg-Marquardt optimization [default]\n" );
		printf( "   opt=lm_levmar      - Local Levenberg-Marquardt optimization using LEVMAR library\n" );
		printf( "   opt=lm_gsl         - Local Levenberg-Marquardt optimization using GSL library\n" );
		printf( "   opt=pso            - Global Particle Swarm optimization (default Standard2006)\n" );
		printf( "   opt=apso           - Global Adaptive Particle Swarm optimization (default TRIBES)\n" );
		printf( "   opt=swarm          - Global Particle Swarm optimization Standard2006 (also opt=pso_std)\n" );
		printf( "   opt=tribes         - Global Particle Swarm optimization TRIBES (also opt=pso_tribes)\n" );
		printf( "   opt=squads         - SQUADS: Adaptive hybrid optimization using coupled local and global optimization techniques\n" );
		printf( "\nsampling method (smp=[string] OR paran=[string] for paranoid LM analysis using multiple retries):\n" );
		printf( "   smp=olhs           - Optimal Latin Hyper Cube sampling [default] (if real = 1 RANDOM; if real > IDLHS; if real > 500 LHS)\n" );
		printf( "   smp=lhs            - Latin Hyper Cube sampling (LHS)\n" );
		printf( "   smp=idlhs          - Improved Distance Latin Hyper Cube sampling (IDLHS)\n" );
		printf( "   smp=random         - Random sampling\n" );
		printf( "\nsampling options:\n" );
		printf( "   real=[integer]     - number of random realizations / samples [default real=100]\n" );
		printf( "   case=[integer]     - execute a single case from all the realizations / samples (applied in PPSD, IGDP, IGRND, MONTECARLO)\n" );
		printf( "   seed=[integer]     - random seed value [randomly generated by default]\n" );
		printf( "\nobjective function functional form options (select one):\n" );
		printf( "   ssr                - sum of the squared residuals [default]\n" );
		printf( "   ssd0               - sum of the squared discrepancies\n" );
		printf( "   ssda               - sum of the squared discrepancies and absolute residuals\n" );
		printf( "   ssdr               - sum of the squared discrepancies and squared residuals\n" );
		printf( "\ntransformation of parameter space and observations:\n" );
		printf( "   nosin              - Sin transformation of optimized parameters is not applied [parameters are sin transformed by default]\n" );
		printf( "   plog=[-1,0,1]      - Log transformation of all optimized parameters is enforced (1) or disabled (0)\n" );
		printf( "                        [default plog=-1; log transformation is explicitly defined for each parameter in the input file]\n" );
		printf( "   olog=[-1,0,1]      - Log transformation of all the observations (simulated and measured) is enforced (1) or disabled (0)\n" );
		printf( "                        [default olog=-1; log transformation is explicitly defined for each observation in the input file]\n" );
		printf( "   oweight=[-1,0,1,2] - Weights for all the observation residuals are defined:\n" );
		printf( "                        0 = zero weight, 1 = unit weight, 2 = weight reversely proportional to observation\n" );
		printf( "                        [default oweight=-1; weights for each observation are explicitly defined in the input file]\n" );
		printf( "\nparallelization (parallelization environment and available resources are internally detected by default):\n" );
		printf( "   np=[integer]       - Number of requested parallel jobs [optional]\n" );
		printf( "   rstfile=[string]   - name of existing ZIP restart file to be used (created by previous Parallel MADS run) [optional]\n" );
		printf( "   restart=[integer]  - restart=1 (default; automatic restart if possible); restart=0 (force no restart); restart=2 (force restart)\n" );
		printf( "                        by default the analyses will be restarted automatically (restart=1)\n" );
		printf( "\nABAGUS (Agent-Based Global Uncertainty & Sensitivity Analysis) options:\n" );
		printf( "   infile=[string]    - name of previous results file to be used to initialize Kd-tree [default=NULL]\n" );
		printf( "   energy=[integer]   - initial energy for particles [default energy=10000]\n" );
		printf( "\nbuild-in analytical solutions:\n" );
		printf( "   point              - point contaminant source in 3D flow domain\n" );
		printf( "   plane              - areal contaminant source in 3D flow domain\n" );
		printf( "   box                - brick contaminant source in 3D flow domain\n" );
		printf( "\nbuild-in test problems for global optimization / uncertainty-quantification techniques (local techniques will not work):\n" );
		printf( "   test=[integer]     - test problem ID [default=0]:\n" );
		printf( "                           0: Parabola (Sphere)\n" );
		printf( "                           1: De Jong's Function #4\n" );
		printf( "                           2: Griewank\n" );
		printf( "                           3: Rosenbrock\n" );
		printf( "                           4: Step\n" );
		printf( "                           6: Foxholes 2D\n" );
		printf( "                           7: Polynomial fitting\n" );
		printf( "                           8: Alpine function (Clerc's Function #1)\n" );
		printf( "                           9: Rastrigin\n" );
		printf( "                          10: Ackley\n" );
		printf( "                          13: 2D Tripod function\n" );
		printf( "                          17: Krishna Kumar\n" );
		printf( "                          18: Eason 2D\n" );
		printf( "   dim=[integer]      - dimensionality of parameter space for the test problem (fixed for some of the problems) [default=2]\n" );
		printf( "\ndebugging / verbose levels:\n" );
		printf( "   debug=[0-5]        - general debugging [default debug=0]\n" );
		printf( "   fdebug=[0-5]       - model evaluation debugging [default fdebug=0]\n" );
		printf( "   ldebug=[0-3]       - Levenberg-Marquardt optimization debugging [default ldebug=0]\n" );
		printf( "   pdebug=[0-3]       - Particle Swarm optimization debugging [default pdebug=0]\n" );
		printf( "   mdebug=[0-3]       - Random sampling debugging [default mdebug=0]\n" );
		printf( "   odebug=[0-1]       - Record objective function progress in a file with extension \'ofe\' [default odebug=0]\n" );
		printf( "   tpldebug=[0-3]     - Debug the writing of external files [default tpldebug=0]\n" );
		printf( "   insdebug=[0-3]     - Debug the reading of external files [default insdebug=0]\n" );
		printf( "   pardebug=[0-3]     - Debug the parallel execution [default pardebug=0]\n" );
		printf( "\nExamples: (code WELLS can be obtained at www.ees.lanl.gov/staff/monty/codes/wells)\n" );
		printf( "   mads a01 test=3 opt=pso igrnd real=1 (no input files are needed for execution)\n" );
		printf( "   mads example/contamination/s01 ldebug (file s01.mads is located in example/contamination)\n" );
		printf( "   mads example/contamination/s02 ldebug (file s02.mads is located in example/contamination)\n" );
		printf( "   mads example/contamination/s01 ldebug igrnd real=1\n" );
		printf( "   mads example/contamination/s01 seed=1549170842 success igrnd real=1\n" );
		printf( "   mads example/contamination/s01 opt=squads seed=1549170842 eigen success pdebug igrnd real=1\n" );
		printf( "   mads example/contamination/s01 opt=pso seed=1549170842 eigen success igrnd real=1\n" );
		printf( "   mads w01 np=2 ldebug pardebug=2 (files associated with problem w01 are located in example/wells)\n" );
		printf( "\nComparisons (code PEST can be obtained at http://www.sspa.com/pest/):\n" );
		printf( "   mads s02 ldebug (file s02.mads is located in example/contamination)\n" );
		printf( "   pest s02pest (file s02pest.pst is located in example/contamination)\n" );
		printf( "   mads w01 ldebug (files associated with problem w01 are located in example/wells)\n" );
		printf( "   pest w01pest (files associated with problem w01pest are located in example/wells)\n" );
		printf( "   pest w02pest (files associated with problem w02pest are located in example/wells)\n" );
		printf( "\nFor additional information:\n" );
		printf( "   web:   www.ees.lanl.gov/staff/monty/codes/mads\n" );
		printf( "   email: Velimir Vesselinov (monty) <vvv@lanl.gov>\n" );
		exit( 1 );
	}
	else
		printf( "Velimir Vesselinov (monty) <vvv@lanl.gov>\nhttp://www.ees.lanl.gov/staff/monty/mads.html\n\n" );
	if( cd.debug ) printf( "Argument[1]: %s\n", argv[1] );
	strcpy( root, argv[1] ); // Defined problem name (root)
	dot = strrchr( root, '.' );
	if( dot != NULL && dot[1] != '/' )
	{
		strcpy( filename, argv[1] );
		strcpy( extension, &dot[1] );
		dot[0] = 0;
	}
	else
	{
		sprintf( filename, "%s.mads", argv[1] );
		extension[0] = 0;
	}
	printf( "Input file name: %s\n", filename );
	cd.time_infile = Fdatetime_t( filename, 0 );
	cd.datetime_infile = Fdatetime( filename, 0 );
	if( cd.debug ) printf( "Root: %s", root );
	op.root = root;
	op.s = 0;
	if( cd.debug && extension[0] != 0 ) printf( " Extension: %s\n", extension );
	printf( "\n" );
	sprintf( buf, "%s.running", op.root ); // File named root.running is used to prevent simultaneous execution of multiple problems
	if( Ftest( buf ) == 0 ) // If file already exists quit ...
	{
		printf( "Potentially another MADS run is currently performed for problem \'%s\' since file %s exists!\n", op.root, buf );
		printf( "If there is no other MADS run, delete %s to execute (sorry for the inconvenience)!\n", buf );
		exit( 1 );
	}
	sprintf( buf, "touch %s.running", op.root ); system( buf ); // Create a file named root.running to prevent simultaneous execution of multiple problems
	/*
	 *  Read input data
	 */
	if( strcasecmp( extension, "pst" ) == 0 ) // PEST Problem
	{
		printf( "PEST problem:\n" );
		load_pst( filename, &cd, &pd, &od, &ed );
		if( cd.opt_method[0] == 0 ) { strcpy( cd.opt_method, "lm" ); cd.calib_type = SIMPLE; cd.problem_type = CALIBRATE; }
		cd.solution_type = EXTERNAL;
		func = func_extrn;
		buf[0] = 0;
		for( i = 2; i < argn; i++ ) { strcat( buf, " " ); strcat( buf, argv[i] ); }
		if( parse_cmd( buf, &cd ) == -1 ) { sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
	}
	else // MADS Problem
	{
		if(( ier = load_problem( filename, argn, argv, &op, &cd, &pd, &od, &wd, &gd, &ed ) ) <= 0 )
		{
			printf( "MADS quits! Data input problem!\nExecute \'mads\' without any arguments to check the acceptable command-line keywords and options.\n" );
			if( ier == 0 )
			{
				sprintf( filename, "%s-error.mads", op.root );
				save_problem( filename, &cd, &pd, &od, &wd, &gd, &ed );
				printf( "\nMADS problem file named %s-error.mads is created to debug.\n", op.root );
			}
			sprintf( buf, "rm -f %s.running", op.root ); system( buf ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			exit( 0 );
		}
		if( cd.solution_type == EXTERNAL ) func = func_extrn;
		else func = func_intrn;
	}
	/*
	 *  Check for parallel environment
	 */
	cd.paral_hosts = NULL;
	hostlist = NULL;
	if(( nodelist = getenv( "NODELIST" ) ) != NULL )
	{
		if( cd.debug ) printf( "\nParallel environment is detected (environmental variable NODELIST is defined)\n" );
		if( cd.debug ) printf( "Node list %s\n", nodelist );
		hostlist = nodelist;
	}
	if(( beowlist = getenv( "BEOWULF_JOB_MAP" ) ) != NULL )
	{
		if( cd.debug ) printf( "\nParallel environment is detected (environmental variable BEOWULF_JOB_MAP is defined)\n" );
		if( cd.debug ) printf( "Node list %s\n", beowlist );
		hostlist = beowlist;
	}
	if(( lsblist = getenv( "LSB_HOSTS" ) ) != NULL )
	{
		if( cd.debug ) printf( "\nParallel environment is detected (environmental variable LSB_HOSTS is defined)\n" );
		if( cd.debug ) printf( "Node list %s\n", lsblist );
		hostlist = lsblist;
		if(( proclist = getenv( "LSB_MCPU_HOSTS" ) ) != NULL && cd.debug ) printf( "LSB_MCPU_HOSTS Processors list %s\n", proclist );
	}
	if( hostlist != NULL )
	{
		if( cd.debug == 0 ) printf( "\nParallel environment is detected.\n" );
		if(( host = getenv( "HOSTNAME" ) ) == NULL ) host = getenv( "HOST" );
		printf( "Host: %s\n", host );
		k = strlen( hostlist );
		i = count = 0;
		printf( "Nodes:" );
		while( i <= k )
		{
			sscanf( &hostlist[i], "%s", buf );
			printf( " \'%s\'", buf );
			i += strlen( buf ) + 1;
			count++;
		}
		printf( "\n" );
		printf( "Number of available nodes for parallel execution: %i\n", count );
		if( count <= 0 )
			printf( "ERROR: There is problem with the description of execution nodes!\n" );
		else
		{
			cd.num_proc = count;
			cd.paral_hosts = char_matrix( cd.num_proc, 95 );
			k = strlen( hostlist );
			i = 0;
			j = 0;
			while( i <= k )
			{
				sscanf( &hostlist[i], "%s", cd.paral_hosts[j] );
				i += strlen( cd.paral_hosts[j++] ) + 1;
			}
		}
	}
	else if( cd.num_proc > 1 )
	{
		printf( "\nLocal parallel execution is requested using %d processors (np=%d)\n", cd.num_proc, cd.num_proc );
		cwd = getenv( "OSTYPE" ); printf( "OS type: %s\n", cwd );
		if( strncasecmp( cwd, "darwin", 6 ) == 0 )
			system( "\\rm -f num_proc >& /dev/null; ( sysctl hw.logicalcpu | cut -d : -f 2 ) > num_proc" ); // MAC OS
		else
			system( "\\rm -f num_proc >& /dev/null; ( cat /proc/cpuinfo | grep processor | wc -l ) > num_proc" ); // LINUX
		in = Fread( "num_proc" );
		fscanf( in, "%d", &k );
		fclose( in );
		system( "\\rm -f num_proc >& /dev/null" );
		printf( "Number of local processors available for parallel execution: %i\n", k );
		if( k < cd.num_proc ) printf( "WARNING: Number of requested processors exceeds the available resources!\n" );
	}
	if( cd.num_proc > 1 ) // Parallel job
	{
		pid = getpid();
		if( cd.debug ) printf( "Parent ID [%d]\n", pid );
		cwd = getenv( "PWD" ); dot = strrchr( cwd, '/' );
		cd.mydir = &dot[1];
		if( cd.debug ) printf( "Working directory: %s (%s)\n", cwd, cd.mydir );
		cd.mydir_hosts = dir_hosts( &op, op.datetime_stamp ); // Directories for parallel execution have unique name based on the execution time
	}
	op.label = ( char * ) malloc( 10 * sizeof( char ) ); op.label[0] = 0;
	if(( orig_params = ( double * ) malloc( pd.nParam * sizeof( double ) ) ) == NULL ) { printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf );  exit( 1 ); }
	/*
	 *  Problem based on external model
	 */
	if( cd.solution_type == EXTERNAL ) // Check the files for external execution
	{
		if( cd.debug || cd.tpldebug || cd.insdebug ) printf( "Checking the instruction and template files for errors ...\n" );
		bad_data = 0;
		for( i = 0; i < pd.nParam; i++ ) orig_params[i] = ( double ) - 1;
		for( i = 0; i < ed.ntpl; i++ ) // Check template files ...
			if( check_par_tpl( pd.nParam, pd.var_id, orig_params, ed.fn_tpl[i], cd.tpldebug ) == -1 )
				bad_data = 1;
		for( i = 0; i < pd.nParam; i++ )
		{
			if( orig_params[i] < 0 )
			{
				printf( "ERROR: Model parameter \'%s\' is not represented in the template file(s)!\n", pd.var_id[i] );
				bad_data = 1;
			}
			else if( orig_params[i] > 1.5 )
			{
				printf( "WARNING: Model parameter \'%s\' is represented more than once (%d) in the template file(s)!\n", pd.var_id[i], ( int ) orig_params[i] );
				bad_data = 1;
			}
		}
		for( i = 0; i < od.nObs; i++ ) od.obs_current[i] = ( double ) - 1;
		for( i = 0; i < ed.nins; i++ )
			if( check_ins_obs( od.nObs, od.obs_id, od.obs_current, ed.fn_ins[i], cd.insdebug ) == -1 ) // Check instruction files.
				bad_data = 1;
		for( i = 0; i < od.nObs; i++ )
		{
			if( od.obs_current[i] < 0 )
			{
				printf( "ERROR: Observation \'%s\' is not defined in the instruction file(s)!\n", od.obs_id[i] );
				bad_data = 1;
			}
			else if( od.obs_current[i] > 1.5 )
				printf( "WARNING: Observation \'%s\' is defined more than once (%d) in the instruction file(s)! Arithmetic average will be computed!\n", od.obs_id[i], ( int ) od.obs_current[i] );
		}
		if( bad_data )
		{
			sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			system( buf );
			exit( -1 );
		}
	}
	/*
	 *  Check for restart conditions
	 */
	printf( "\nExecution date & time stamp: %s\n", op.datetime_stamp ); // Stamp will be applied to name / rename various output files
	if( cd.solution_type == EXTERNAL && cd.num_proc > 1 )
	{
		if( cd.restart == 1 ) // Restart by default
		{
			strcpy( buf, filename ); // Temporarily preserve the input file name
			sprintf( filename, "%s.restart_%s.zip", op.root, cd.datetime_infile );
			strcpy( cd.restart_zip_file, filename );
			if( Ftest( cd.restart_zip_file ) != 0 ) { if( cd.pardebug ) printf( "ZIP file (%s) with restart information is not available.\n", cd.restart_zip_file ); cd.restart = 0; }
			else
			{
				dt = cd.time_infile - Fdatetime_t( cd.restart_zip_file, 0 ); // time_infile - time_zipfile ...
				if( dt >= 0 ) { if( cd.pardebug ) printf( "No restart: the zip file (%s) with restart information is older than the MADS input file (%s)\n(restart can be enforced using \'restart=-1\' or \'rstfile=%s\')\n", cd.restart_zip_file, buf, cd.restart_zip_file ); cd.restart = 0; } // No restart
				else cd.restart = 1; // Attempt restart
			}
			if( cd.restart )
				printf( "DEFAULT Restart: zip file %s is consistent with date/time stamp of the MADS input file\n(to ignore the zip file, either delete the zip file, or use keyword \'restart=0\' ... \n", filename );
		}
		else if( cd.restart == -1 ) // Forced restart
		{
			sprintf( filename, "%s.restart_%s.zip", op.root, cd.datetime_infile );
			if( cd.restart_zip_file[0] == 0 ) strcpy( cd.restart_zip_file, filename );
			if( Ftest( cd.restart_zip_file ) != 0 ) { printf( "Restart is requested but a zip file (%s) with restart information is not available.\n", cd.restart_zip_file ); cd.restart = 0; }
			else printf( "FORCED Restart: using zip file %s ...\n", cd.restart_zip_file );
		}
		if( cd.restart )
		{
			printf( "MADS  input  file \'%40s\' last modified on %s\n", buf, Fdatetime( buf, 0 ) );
			sprintf( filename, "%s.results", op.root ); if( Ftest( filename ) != 0 ) printf( "MADS results file \'%40s\' last modified on %s\n", filename, Fdatetime( filename, 0 ) );
			printf( "MADS restart file \'%40s\' last modified on %s\n", cd.restart_zip_file, Fdatetime( cd.restart_zip_file, 0 ) );
			printf( "ZIP file (%s) with restart information is unzipped ... \n", cd.restart_zip_file );
			sprintf( buf, "rm -fR ../%s* %s.restart_info; unzip -u -: %s >& /dev/null", cd.mydir_hosts, op.root, cd.restart_zip_file ); // the input file name was temporarily in buf; not any more ...
			system( buf );
			sprintf( filename, "%s.restart_info", op.root );
			in = Fread( filename );
			fgets( buf, 255, in );
			white_trim( buf );
			cd.mydir_hosts = dir_hosts( &op, buf ); // Directories for parallel execution have unique name based on the old execution time (when restart files were created)
			fclose( in );
			printf( "Date & time stamp of the previous run: %s\n", buf );
		}
		// Preserve the existing restart zip file
		if( Ftest( cd.restart_zip_file ) == 0 )
		{
			if( cd.pardebug ) printf( "Previous restart file (%s) exists!\n", cd.restart_zip_file );
			if( cd.restart ) sprintf( buf, "cp %s %s.restart_%s_%s.zip >& /dev/null", cd.restart_zip_file, op.root, cd.datetime_infile, Fdatetime( cd.restart_zip_file, 0 ) );  // Copy if restart
			else sprintf( buf, "mv %s %s.restart_%s_%s.zip >& /dev/null", cd.restart_zip_file, op.root, cd.datetime_infile, Fdatetime( cd.restart_zip_file, 0 ) );  // Move if no restart
			system( buf );
		}
		if( cd.restart == 0 )
		{
			sprintf( filename, "%s.restart_info", op.root );
			out = Fwrite( filename );
			fprintf( out, "%s\n", op.datetime_stamp );
			for( i = 0; i < argn; i++ )
				fprintf( out, "%s ", argv[i] );
			fprintf( out, "\n" );
			fclose( out );
			sprintf( buf, "zip %s %s.restart_info >& /dev/null", cd.restart_zip_file, op.root );
			system( buf );
		}
	}
//
// ------------------------ IGRND
//
	if( cd.problem_type == CALIBRATE && cd.calib_type == IGRND ) /* Calibration analysis using random initial guessed */
	{
		strcpy( op.label, "igrnd" );
		if(( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		printf( "\nSEQUENTIAL CALIBRATION using random initial guesses for model parameters (realizations = %d):\n", cd.nreal );
		if( pd.nFlgParam != 0 ) { printf( "Only flagged parameters are randomized\n" ); npar = pd.nFlgParam; }
		else if( pd.nOptParam != 0 ) { printf( "No flagged parameters; all optimizable parameters are randomized\n" ); npar = pd.nOptParam; }
		else { printf( "No flagged or optimizable parameters; all parameters are randomized\n" ); npar = pd.nParam; }
		if(( var_lhs = ( double * ) malloc( npar * cd.nreal * sizeof( double ) ) ) == NULL )
		{
			printf( "Not enough memory!\n" );
			sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
			system( buf );
			exit( 1 );
		}
		if( cd.seed < 0 ) { cd.seed *= -1; printf( "Imported seed: %d\n", cd.seed ); }
		else if( cd.seed == 0 ) { printf( "New " ); cd.seed = get_seed(); }
		else printf( "Current seed: %d\n", cd.seed );
		printf( "Random sampling (variables %d; realizations %d) using ", npar, cd.nreal );
		sampling( npar, cd.nreal, &cd.seed, var_lhs, &op );
		printf( "done.\n" );
		if( cd.mdebug )
		{
			sprintf( filename, "%s.igrnd_set", op.root );
			out = Fwrite( filename );
			for( count = 0; count < cd.nreal; count ++ )
			{
				for( k = 0; k < npar; k++ )
					fprintf( out, "%.15g ", var_lhs[k+count*npar] );
				fprintf( out, "\n" );
			}
			fclose( out );
			printf( "Random sampling set saved in %s.igrnd_set\n", op.root );
		}
		for( i = 0; i < pd.nParam; i++ )
			orig_params[i] = pd.var[i]; // Save original initial values for all parameters
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) f = optimize_lm; // Define optimization method: LM
		else f = optimize_pso; // Define optimization method: PSO
		// File management
		sprintf( filename, "%s.igrnd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.igrnd.zip %s.igrnd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.igrnd.zip %s.igrnd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.igrnd.zip %s.igrnd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.igrnd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.igrnd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		sprintf( filename, "%s.igrnd-opt=%s_eval=%d_real=%d", op.root, cd.opt_method, cd.maxeval, cd.nreal );
		out2 = Fwrite( filename );
		if( pd.nOptParam == 0 )
			printf( "WARNING: No parameters to optimize! Forward runs performed instead (ie Monte Carlo analysis)\n" );
		phi_min = HUGE;
		status_global = eval_total = 0;
		if( cd.ncase != 0 ) k = cd.ncase - 1; // applied if execution of a specific realization is requested (ncase)
		else k = 0;
		printf( "\n" );
		for( count = k; count < cd.nreal; count++ )
		{
			cd.eval = 0;
			fprintf( out, "%d : init var", count + 1 );
			printf( "\nRandom set #%d: ", count + 1 );
			if( cd.mdebug || cd.nreal == 1 ) printf( "\n" );
			op.s = count + 1;
			for( k = i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] == 2 || ( pd.var_opt[i] == 1 && pd.nFlgParam == 0 ) )
				{
					pd.var[i] = var_lhs[k+count*npar] * pd.var_range[i] + pd.var_min[i];
					if( pd.var_log[i] )
					{
						if( cd.mdebug || cd.nreal == 1 ) printf( "%s %.15g\n", pd.var_id[i], pow( 10, pd.var[i] ) );
						fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
					}
					else
					{
						if( cd.mdebug || cd.nreal == 1 ) printf( "%s %.15g\n", pd.var_id[i], pd.var[i] );
						fprintf( out, " %.15g", pd.var[i] );
					}
					k++;
				}
				else pd.var[i] = orig_params[i];
			if( pd.nOptParam > 0 ) status = f( &op ); // Optimize
			else
			{
				for( i = 0; i < pd.nParam; i++ )
					pd.var[i] = var_lhs[i+count*npar] * pd.var_range[i] + pd.var_min[i];
				if( cd.mdebug ) { printf( "Forward run ... \n" ); debug_level = cd.fdebug; cd.fdebug = 3; }
				cd.compute_phi = 1;
				func( pd.var, &op, od.res ); // pd.var is a dummy variable because pd.nOptParam == 0
				cd.compute_phi = 0;
				if( cd.mdebug ) cd.fdebug = debug_level;
			}
			if( cd.mdebug || cd.nreal == 1 )
			{
				printf( "\n" );
				print_results( &op );
			}
			status_global += op.success;
			if( op.phi < phi_min ) { phi_min = op.phi; for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]]; }
			if( cd.pdebug || cd.ldebug ) printf( "\n" ); // extra new line if the optimization process is debugged
			fprintf( out2, "%g %d %d\n", op.phi, op.success, cd.eval );
			fflush( out2 );
			fprintf( out, " : OF %g status %d : final var", op.phi, op.success );
			for( i = 0; i < pd.nOptParam; i++ ) // Print only optimized parameters (including flagged); ignore fixed parameters
			{
				k = pd.var_index[i];
				if( pd.var_log[k] ) fprintf( out, " %.15g", pow( 10, pd.var[k] ) );
				else fprintf( out, " %.15g", pd.var[k] );
			}
			fprintf( out, "\n" );
			fflush( out );
			if( op.success && cd.nreal > 1 ) save_results( "igrnd", &op, &gd );
			else                             save_results( "igrnd", &op, &gd );
			eval_total += cd.eval;
			if( cd.ncase != 0 ) break;
		}
		op.s = 0;
		free( var_lhs );
		cd.eval = eval_total; // provide the correct number of total evaluations
		op.phi = phi_min;
		for( i = 0; i < pd.nOptParam; i++ ) opt_params[i] = pd.var[pd.var_index[i]] = pd.var_current[i] = pd.var_best[i]; // get the best estimate
		printf( "Minimum objective function: %g\n", phi_min );
		printf( "Total number of evaluations = %d\n", eval_total );
		fprintf( out, "Minimum objective function: %g\n", phi_min );
		if( cd.nreal > 1 )
		{
			if( status_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions within calibration ranges!\n", cd.nreal );
			else printf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", status_global, cd.nreal, ( double ) status_global / cd.nreal );
		}
		fprintf( out2, "min %g\n", phi_min );
		fprintf( out2, "status %d\n", status_global );
		fprintf( out, "Number of the sequential calibration runs producing predictions within calibration ranges = %d\n", status_global );
		fprintf( out2, "eval %d\n", eval_total );
		fprintf( out, "Number of evaluations = %d\n", eval_total );
		fclose( out ); fclose( out2 );
		printf( "Results are saved in %s.igrnd.results and %s.igrnd-opt=%s_eval=%d_real=%d\n", op.root, op.root, cd.opt_method, cd.maxeval, cd.nreal );
		printf( "Repeat the run producing the best results and save them ...\n" );
		if( cd.debug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
		Transform( opt_params, &op, opt_params );
		cd.compute_phi = 1;
		func( opt_params, &op, od.res );
		cd.compute_phi = 0;
		if( cd.debug ) cd.fdebug = debug_level;
		free( opt_params );
		save_results( "", &op, &gd );
	}
//
// ------------------------ IGPD
//
	if( cd.problem_type == CALIBRATE && cd.calib_type == IGPD ) /* Calibration analysis using discretized initial guesses */
	{
		strcpy( op.label, "igpd" );
		if(( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		printf( "\nSEQUENTIAL CALIBRATION using discretized initial guesses for model parameters:\n" );
		if( pd.nFlgParam == 0 )
			printf( "WARNING: No flagged parameters! Discretization of the initial guesses cannot be performed! Forward run will be performed instead.\n" );
		if( pd.nOptParam == 0 )
			printf( "WARNING: No parameters to optimize! Forward run will be performed instead.\n" );
		for( i = 0; i < pd.nParam; i++ ) orig_params[i] = pd.var[i]; // Save original initial values for all parameters
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) f = optimize_lm; // Define optimization method: LM
		else f = optimize_pso; // Define optimization method: PSO
		// File management
		sprintf( filename, "%s.igpd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.igpd.zip %s.igpd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.igpd.zip %s.igpd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.igpd.zip %s.igpd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.igpd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.igpd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		sprintf( filename, "%s.igpd-opt=%s_eval=%d_real=%d", op.root, cd.opt_method, cd.maxeval, cd.nreal );
		out2 = Fwrite( filename );
		for( i = 0; i < pd.nParam; i++ )
			if( pd.var_opt[i] == 2 )
				orig_params[i] = pd.var_min[i];
		phi_min = HUGE;
		count = eval_total = 0;
		do
		{
			cd.eval = 0;
			count++;
			if( cd.ncase == 0 || cd.ncase == count )
			{
				fprintf( out, "%d : init var", count ); // counter
				printf( "SEQUENTIAL CALIBRATION #%d: ", count );
				op.s = count;
				if( cd.debug == 0 ) printf( "\n" );
				for( i = 0; i < pd.nParam; i++ )
				{
					pd.var[i] = orig_params[i];
					if( pd.var_opt[i] == 2 ) // Print flagged parameters
					{
						printf( "%s %g\n", pd.var_id[i], orig_params[i] );
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, orig_params[i] ) );
						else fprintf( out, " %.15g", orig_params[i] );
					}
				}
				if( pd.nOptParam > 0 ) status = f( &op ); /* Optimize */
				else
				{
					if( cd.debug ) { printf( "Forward run ... \n" ); debug_level = cd.fdebug; cd.fdebug = 3; }
					cd.compute_phi = 1;
					func( pd.var, &op, od.res ); // pd.var is dummy because pd.nOptParam == 0
					cd.compute_phi = 0;
					if( cd.debug ) cd.fdebug = debug_level;
				}
				eval_total += cd.eval;
				if( cd.debug )
				{
					printf( "\n" );
					fprintf( out, " : OF %g status %d : final var ", op.phi, op.success );
				}
				else printf( "Objective function: %g Success: %d\n", op.phi, op.success );
				if( op.phi < phi_min ) { phi_min = op.phi; for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]]; }
				fprintf( out, " %d", op.success );
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 ) // Print only optimized parameters (including flagged); ignore fixed parameters
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				fprintf( out, "\n" );
				fflush( out );
				if( op.success ) save_results( "igpd", &op, &gd );
				if( cd.ncase != 0 ) break;
			}
			if( pd.nFlgParam == 0 || pd.nOptParam == 0 ) break;
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] == 2 )
				{
					if( orig_params[i] < pd.var_max[i] ) { orig_params[i] += pd.var_dx[i]; break; }
					else orig_params[i] = pd.var_min[i];
				}
			if( i == pd.nParam ) break;
		}
		while( 1 );
		op.s = 0;
		cd.eval = eval_total; // provide the correct number of total evaluations
		printf( "Total number of evaluations = %d\n", eval_total );
		op.phi = phi_min;
		for( i = 0; i < pd.nOptParam; i++ ) opt_params[i] = pd.var[pd.var_index[i]] = pd.var_current[i] = pd.var_best[i];
		printf( "Minimum objective function: %g\n", phi_min );
		printf( "Total number of evaluations = %d\n", eval_total );
		fprintf( out, "Minimum objective function: %g\n", phi_min );
		if( cd.nreal > 1 )
		{
			if( status_global == 0 ) printf( "None of the %d sequential calibration runs produced predictions within calibration ranges!\n", cd.nreal );
			else printf( "Number of the sequential calibration runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", status_global, cd.nreal, ( double ) status_global / cd.nreal );
		}
		fprintf( out2, "min %g\n", phi_min );
		fprintf( out2, "status %d\n", status_global );
		fprintf( out, "Number of the sequential calibration runs producing predictions within calibration ranges = %d\n", status_global );
		fprintf( out2, "eval %d\n", eval_total );
		fprintf( out, "Number of evaluations = %d\n", eval_total );
		fclose( out ); fclose( out2 );
		printf( "Results are saved in %s.igpd.results\n", op.root );
		printf( "Repeat the run producing the best results ...\n" );
		if( cd.debug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
		Transform( opt_params, &op, opt_params );
		cd.compute_phi = 1;
		func( opt_params, &op, od.res );
		cd.compute_phi = 0;
		if( cd.debug ) cd.fdebug = debug_level;
		free( opt_params );
		save_results( "", &op, &gd );
	}
//
// ------------------------ PPSD
//
	if( cd.problem_type == CALIBRATE && cd.calib_type == PPSD ) /* Calibration analysis using discretized parameters */
	{
		strcpy( op.label, "ppsd" );
		printf( "\nSEQUENTIAL CALIBRATION using partial parameter-space discretization (PPSD):\n" );
		if( pd.nFlgParam == 0 )
			printf( "WARNING: No flagged parameters! Discretization of the initial guesses cannot be performed!\n" );
		if( pd.nOptParam == 0 )
			printf( "WARNING: No parameters to optimize! Forward runs performed instead\n" );
		for( i = 0; i < pd.nParam; i++ ) orig_params[i] = pd.var[i]; // Save original initial values for all parameters
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) f = optimize_lm; // Define optimization method: LM
		else f = optimize_pso; // Define optimization method: PSO
		sprintf( filename, "%s.igpd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.ppsd.zip %s.ppsd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.ppsd.zip %s.ppsd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.ppsd.zip %s.ppsd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.ppsd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.ppsd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		for( i = 0; i < pd.nParam; i++ )
			if( pd.var_opt[i] == 2 ) cd.var[i] = pd.var_min[i];
		phi_min = HUGE;
		count = eval_total = 0;
		do
		{
			cd.eval = 0;
			count++;
			if( cd.ncase == 0 || cd.ncase == count )
			{
				fprintf( out, "%d : ", count );
				printf( "\nDISCRETIZED CALIBRATION #%d:", count );
				op.s = count;
				if( cd.debug ) printf( "\nDiscretized parameters:\n" );
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] == 2 ) // Print flagged parameters
					{
						if( cd.debug ) printf( "%s %g\n", pd.var_id[i], cd.var[i] );
						fprintf( out, "%g ", cd.var[i] );
					}
					else pd.var[i] = orig_params[i]; // these are the true original parameters
				if( pd.nOptParam > 0 ) status = f( &op ); /* Optimize */
				else
				{
					if( cd.debug ) { printf( "Forward run ... \n" ); debug_level = cd.fdebug; cd.fdebug = 3; }
					cd.compute_phi = 1;
					func( pd.var, &op, od.res ); // pd.var is dummy because pd.nOptParam == 0
					cd.compute_phi = 0;
					if( cd.debug ) cd.fdebug = debug_level;
				}
				eval_total += cd.eval;
				if( cd.debug )
				{
					printf( "\n" );
					print_results( &op );
				}
				else printf( "Objective function: %g Status: %d\n", op.phi, op.success );
				fprintf( out, " : OF %g status %d : final var", op.phi, op.success );
				if( op.phi < phi_min ) { phi_min = op.phi; for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]]; }
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] == 1 ) // Print only optimized parameters; ignore fixed and flagged parameters
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				fprintf( out, "\n" );
				fflush( out );
				if( op.success ) save_results( "ppsd", &op, &gd );
				if( cd.ncase != 0 ) break;
			}
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] == 2 )
				{
					if( cd.var[i] < pd.var_max[i] ) { cd.var[i] += pd.var_dx[i]; break; }
					else cd.var[i] = pd.var_min[i];
				}
			if( i == pd.nParam ) break;
		}
		while( 1 );
		cd.eval = eval_total; // provide the correct number of total evaluations
		printf( "Total number of evaluations = %d\n", eval_total );
		op.s = 0;
		fclose( out );
		printf( "Results are saved in %s.ppsd.results\n", op.root );
	}
	//
	// ------------------------ MONTECARLO
	//
	if( cd.problem_type == MONTECARLO ) /* Monte Carlo analysis */
	{
		strcpy( op.label, "mcrnd" );
		if(( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		npar = pd.nOptParam;
		if(( var_lhs = ( double * ) malloc( npar * cd.nreal * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		printf( "\nMonte Carlo analysis using latin-hyper cube sampling:\n" );
		if( cd.seed < 0 ) { cd.seed *= -1; printf( "Imported seed: %d\n", cd.seed ); }
		else if( cd.seed == 0 ) { printf( "New " ); cd.seed = get_seed(); }
		else printf( "Current seed: %d\n", cd.seed );
		printf( "Random sampling (variables %d; realizations %d) using ", npar, cd.nreal );
		sampling( npar, cd.nreal, &cd.seed, var_lhs, &op );
		printf( "done.\n" );
		if( cd.mdebug )
		{
			sprintf( filename, "%s.mcrnd_set", op.root );
			out = Fwrite( filename );
			for( count = 0; count < cd.nreal; count ++ )
			{
				for( k = 0; k < npar; k++ )
					fprintf( out, "%.15g ", var_lhs[k+count*npar] );
				fprintf( out, "\n" );
			}
			fclose( out );
			printf( "Random sampling set saved in %s.mcrnd_set\n", op.root );
		}
		sprintf( filename, "%s.mcrnd.zip", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.mcrnd.zip %s.mcrnd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		sprintf( buf, "zip -m %s.mcrnd.zip %s.mcrnd-[0-9]*.* >& /dev/null", op.root, op.root ); system( buf );
		sprintf( buf, "mv %s.mcrnd.zip %s.mcrnd_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
		sprintf( filename, "%s.mcrnd.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.mcrnd_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		cd.compute_phi = 1;
		status_global = 0;
		phi_min = HUGE;
		if( cd.ncase != 0 ) k = cd.ncase - 1;
		else k = 0;
		if( cd.solution_type == EXTERNAL && cd.num_proc > 1 && k == 0 ) // Parallel job:
		{
			if( cd.debug ) printf( "Parallel execution of external jobs ...\n" );
			for( count = 0; count < cd.nreal; count ++ ) // Write all the files
			{
				fprintf( out, "%d : ", count + 1 ); // counter
				if( cd.mdebug ) printf( "\n" );
				printf( "Random set #%d: ", count + 1 );
				fflush( stdout );
				for( i = 0; i < pd.nOptParam; i++ )
				{
					k = pd.var_index[i];
					opt_params[i] = pd.var[k] = var_lhs[i+count*npar] * pd.var_range[k] + pd.var_min[k];
				}
				if( cd.mdebug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
				Transform( opt_params, &op, opt_params );
				func_extrn_write( count + 1, opt_params, &op );
				printf( "external model input file(s) generated ...\n" );
				if( cd.mdebug ) cd.fdebug = debug_level;
				if( cd.mdebug )
				{
					printf( "\nRandom parameter values:\n" );
					for( i = 0; i < pd.nOptParam; i++ )
						if( pd.var_log[pd.var_index[i]] == 0 ) printf( "%s %g\n", pd.var_id[pd.var_index[i]], pd.var[pd.var_index[i]] );
						else printf( "%s %g\n", pd.var_id[pd.var_index[i]], pow( 10, pd.var[pd.var_index[i]] ) );
				}
			}
			if( cd.pardebug > 4 )
			{
				for( count = 0; count < cd.nreal; count ++ ) // Perform all the runs in serial model (for testing)
				{
					printf( "Execute model #%d ... ", count + 1 );
					fflush( stdout );
					func_extrn_exec_serial( count + 1, &op );
					printf( "done!\n" );
				}
			}
			else if( mprun( cd.nreal, &op ) < 0 ) // Perform all the runs in parallel
			{
				printf( "ERROR: there is a problem with the parallel execution!\n" );
				sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
				system( buf );
			}
			for( count = 0; count < cd.nreal; count ++ ) // Read all the files
			{
				printf( "Model results #%d: ", count + 1 );
				bad_data = 0;
				bad_data = func_extrn_read( count + 1, &op, od.res );
				if( bad_data )
				{
					sprintf( buf, "rm -f %s.running", op.root ); // Delete a file named root.running to prevent simultaneous execution of multiple problems
					system( buf );
					exit( -1 );
				}
				if( cd.mdebug > 1 ) { printf( "\n" ); print_results( &op ); }
				else if( cd.mdebug )
				{
					printf( "Objective function: %g Success: %d\n", op.phi, status_all );
					if( status_all ) printf( "All the predictions are within calibration ranges!\n" );
					else printf( "At least one of the predictions is outside calibration ranges!\n" );
				}
				else
					printf( "Objective function: %g Success = %d\n", op.phi, status_all );
				if( op.phi < phi_min ) { phi_min = op.phi; for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]]; }
				if( status_all ) status_global++;
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 )
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				fprintf( out, " OF %g status %d\n", op.phi, status_all );
				fflush( out );
				if( status_all ) save_results( "mcrnd", &op, &gd );
			}
		}
		else // Serial job
		{
			for( count = k; count < cd.nreal; count ++ )
			{
				fprintf( out, "%d : ", count + 1 ); // counter
				if( cd.mdebug ) printf( "\n" );
				printf( "Random set #%d: ", count + 1 );
				fflush( stdout );
				for( i = 0; i < pd.nOptParam; i++ )
				{
					k = pd.var_index[i];
					opt_params[i] = pd.var[k] = var_lhs[i+count*npar] * pd.var_range[k] + pd.var_min[k];
				}
				if( cd.mdebug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
				Transform( opt_params, &op, opt_params );
				func( opt_params, &op, od.res );
				if( cd.mdebug ) cd.fdebug = debug_level;
				if( cd.mdebug )
				{
					printf( "\nRandom parameter values:\n" );
					for( i = 0; i < pd.nOptParam; i++ )
						if( pd.var_log[pd.var_index[i]] == 0 ) printf( "%s %g\n", pd.var_id[pd.var_index[i]], pd.var[pd.var_index[i]] );
						else printf( "%s %g\n", pd.var_id[pd.var_index[i]], pow( 10, pd.var[pd.var_index[i]] ) );
				}
				if( cd.mdebug > 1 ) { printf( "\nPredicted calibration targets:\n" ); print_results( &op ); }
				else if( cd.mdebug )
				{
					printf( "Objective function: %g Success: %d\n", op.phi, status_all );
					if( status_all ) printf( "All the predictions are within calibration ranges!\n" );
					else printf( "At least one of the predictions is outside calibration ranges!\n" );
				}
				else
					printf( "Objective function: %g Success = %d\n", op.phi, status_all );
				if( op.phi < phi_min ) { phi_min = op.phi; for( i = 0; i < pd.nOptParam; i++ ) pd.var_best[i] = pd.var[pd.var_index[i]]; }
				if( status_all ) status_global++;
				for( i = 0; i < pd.nParam; i++ )
					if( pd.var_opt[i] >= 1 )
					{
						if( pd.var_log[i] ) fprintf( out, " %.15g", pow( 10, pd.var[i] ) );
						else fprintf( out, " %.15g", pd.var[i] );
					}
				fprintf( out, " OF %g status %d\n", op.phi, status_all );
				fflush( out );
				if( status_all ) save_results( "mcrnd", &op, &gd );
				if( cd.ncase != 0 ) break;
			}
		}
		op.s = 0;
		free( var_lhs );
		fclose( out );
		op.phi = phi_min;
		for( i = 0; i < pd.nOptParam; i++ ) opt_params[i] = pd.var[pd.var_index[i]] = pd.var_current[i] = pd.var_best[i];
		printf( "Results are saved in %s.mcrnd.results\n", op.root );
		printf( "Minimum objective function: %g\n", phi_min );
		if( status_global == 0 ) printf( "None of the Monte-Carlo runs produced predictions within calibration ranges!\n" );
		else printf( "Number of Monte-Carlo runs producing predictions within calibration ranges = %d (out of %d; success ratio %g)\n", status_global, cd.nreal, ( double ) status_global / cd.nreal );
		printf( "Repeat the Monte-Carlo run producing the best results ...\n" );
		if( cd.debug ) { debug_level = cd.fdebug; cd.fdebug = 3; }
		Transform( opt_params, &op, opt_params );
		func( opt_params, &op, od.res );
		if( cd.debug ) cd.fdebug = debug_level;
		free( opt_params );
		save_results( "", &op, &gd );
		cd.compute_phi = 0;
	}
	//
	// ------------------------ GLOBALSENS
	//
	if( cd.problem_type == GLOBALSENS ) // Global sensitivity analysis run
	{
		strcpy( op.label, "gsens" );
		double fhat, fhat2, *phis_full, *phis_half;
		int n_sub; //! number of samples for subsets a and b
		//		gsl_qrng *q = gsl_qrng_alloc( gsl_qrng_sobol, pd.nOptParam );
		n_sub = cd.nreal / 2;	// set to half of user specified reals
		if(( opt_params = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Temporary variable to store cd.nreal phis
		if(( phis_full = ( double * ) malloc( cd.nreal * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Temporary variable to store m_sub phis
		if(( phis_half = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Temporary variable to store random sample a
		if(( var_a_lhs = ( double * ) malloc( pd.nOptParam * n_sub * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Sample a phis
		if(( gs.f_a = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Sample b phis
		if(( gs.f_b = ( double * ) malloc( n_sub * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Temporary variable to store random sample b
		if(( var_b_lhs = ( double * ) malloc( pd.nOptParam * n_sub * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// matrices to store lhs samples
		gs.var_a_lhs = double_matrix( n_sub, pd.nOptParam );
		gs.var_b_lhs = double_matrix( n_sub, pd.nOptParam );
		// Matrices to store phis with different combinations of parameters from samples a and b
		if(( gs.fmat_a = double_matrix( pd.nOptParam, n_sub ) ) == NULL )
			{ printf( "Error creating 3D matrix\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		if(( gs.fmat_b = double_matrix( pd.nOptParam, n_sub ) ) == NULL )
			{ printf( "Error creating 3D matrix\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 0 ); }
		// Vector of variances for individual component contribution
		if(( gs.D_hat = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		// Vector of variances for total component contribution
		if(( gs.D_hat_n = ( double * ) malloc( pd.nOptParam * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		printf( "\nGlobal sensitivity analysis using random sampling:\n" );
		// Create samples
		if( cd.seed < 0 ) { cd.seed *= -1; printf( "Imported seed: %d\n", cd.seed ); }
		else if( cd.seed == 0 ) { printf( "New " ); cd.seed = get_seed(); }
		else printf( "Current seed: %d\n", cd.seed );
		printf( "Random sampling set 1 (variables %d; realizations %d) using ", pd.nOptParam, cd.nreal );
		sampling( pd.nOptParam, n_sub, &cd.seed, var_a_lhs, &op );
		printf( "done.\n" );
		printf( "Random sampling set 2 (variables %d; realizations %d) using ", pd.nOptParam, cd.nreal );
		sampling( pd.nOptParam, n_sub, &cd.seed, var_b_lhs, &op );
		printf( "done.\n" );
		// Create samples using Sobol's quasi-random sequence
		/*		for( count = 0; count < n_sub; count++ )
				{
					double v[ pd.nOptParam ];
					gsl_qrng_get( q, v);
					for( i = 0; i < pd.nOptParam; i++ )
					{
						k = pd.var_index[i];
						gs.var_a_lhs[count][i] = v[i] * pd.var_range[k] + pd.var_min[k];
					}
				}

				for( count = 0; count < n_sub; count++ )
				{
					double v[ pd.nOptParam ];
					gsl_qrng_get( q, v);
					for( i = 0; i < pd.nOptParam; i++ )
					{
						k = pd.var_index[i];
						gs.var_b_lhs[count][i] = v[i] * pd.var_range[k] + pd.var_min[k];
					}
				}*/
		// Copy temp lhs vectors to matrices
		for( count = 0; count < n_sub; count++ )
			for( i = 0; i < pd.nOptParam; i++ )
			{
				k = pd.var_index[i];
				gs.var_a_lhs[count][i] = var_a_lhs[i+count*pd.nOptParam] * pd.var_range[k] + pd.var_min[k];
				gs.var_b_lhs[count][i] = var_b_lhs[i+count*pd.nOptParam] * pd.var_range[k] + pd.var_min[k];
			}
		free( var_a_lhs );
		free( var_b_lhs );
		// Output samples to files
		if( cd.mdebug )
		{
			sprintf( filename, "%s.gsens.zip", op.root );
			if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s.gsens.zip %s.gsens_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
			sprintf( buf, "zip -m %s.gsens.zip %s.gsens_set_* >& /dev/null", op.root, op.root ); system( buf );
			sprintf( buf, "mv %s.gsens.zip %s.gsens_%s.zip >& /dev/null", op.root, op.root, Fdatetime( filename, 0 ) ); system( buf );
			sprintf( filename, "%s.gsens_set_a", op.root ); out = Fwrite( filename );
			sprintf( filename, "%s.gsens_set_b", op.root ); out2 = Fwrite( filename );
			for( count = 0; count < n_sub; count ++ )
			{
				for( k = 0; k < pd.nOptParam; k++ )
				{
					fprintf( out, "%.15g ", gs.var_a_lhs[count][k] );
					fprintf( out2, "%.15g ", gs.var_b_lhs[count][k] );
				}
				fprintf( out, "\n" );
				fprintf( out2, "\n" );
			}
			fclose( out );
			fclose( out2 );
			printf( "Random sampling sets a and b saved in %s.mcrnd_set_a and %s.mcrnd_set_b\n", op.root, op.root );
		}
		sprintf( filename, "%s.gsens.results", op.root );
		if( Ftest( filename ) == 0 ) { sprintf( buf, "mv %s %s.gsens_%s.results >& /dev/null", filename, op.root, Fdatetime( filename, 0 ) ); system( buf ); }
		out = Fwrite( filename );
		// Accumulate phis into fhat and fhat2 for total output mean and variance
		fhat = fhat2 = 0;
		cd.compute_phi = 1;
		printf( "Computing phis to calculate total output mean and variance...\n" );
		// Compute sample a phis
		for( count = 0; count < n_sub; count ++ )
		{
			for( i = 0; i < pd.nOptParam; i++ )
			{
				k = pd.var_index[i];
				opt_params[i] = pd.var[k] = gs.var_a_lhs[count][i];
			}
			Transform( opt_params, &op, opt_params );
			func( opt_params, &op, od.res );
			// Sum phi and phi^2
			fhat += op.phi;
			fhat2 += pow( op.phi, 2 );
			// Save sample a phis
			gs.f_a[count] = op.phi;
			phis_full[count] = op.phi;
			// save to results file
			fprintf( out, "%d : ", count + 1 ); // counter
			fprintf( out, "%g :", op.phi );
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] >= 1 )
					fprintf( out, " %.15g", pd.var[i] );
			fprintf( out, "\n" );
			fflush( out );
		}
		// Compute sample b phis
		for( count = 0; count < n_sub; count ++ )
		{
			for( i = 0; i < pd.nOptParam; i++ )
			{
				k = pd.var_index[i];
				opt_params[i] = pd.var[k] = gs.var_b_lhs[count][i];
			}
			Transform( opt_params, &op, opt_params );
			func( opt_params, &op, od.res );
			// Sum phi and phi^2
			fhat += op.phi;
			fhat2 += pow( op.phi, 2 );
			// Save sample b phis
			gs.f_b[count] = op.phi;
			phis_full[ n_sub + count ] = op.phi;
			// save to results file
			fprintf( out, "%d : ", n_sub + count ); // counter
			fprintf( out, "%g :", op.phi );
			for( i = 0; i < pd.nParam; i++ )
				if( pd.var_opt[i] >= 1 )
					fprintf( out, " %.15g", pd.var[i] );
			fprintf( out, "\n" );
			fflush( out );
		}
		fclose( out );
		printf( "Global Sensitivity MC results are saved in %s.gsens.results\n", op.root );
		// Calculate total output mean and variance based on sample a
		gs.f_hat_0 = fhat / ( 2 * n_sub );
		gs.D_hat_t = fhat2 / ( 2 * n_sub ) - gs.f_hat_0;
		printf( "Total output mean = %g\n", gs.f_hat_0 );
		printf( "Total output variance = %g\n", gs.D_hat_t );
		gs.f_hat_0 = gsl_stats_mean( phis_full, 1, cd.nreal );
		gs.D_hat_t = gsl_stats_variance( phis_full, 1, cd.nreal );
		printf( "Total output mean = %g\n", gs.f_hat_0 );
		printf( "Total output variance = %g\n", gs.D_hat_t );
		gs.f_hat_0 = gs.D_hat_t = 0.0;
		ave_sorted( phis_full, cd.nreal, &gs.f_hat_0, &gs.ep );
		printf( "Total output mean = %g abs 1st moment = %g\n", gs.f_hat_0, gs.ep );
		var_sorted( phis_full, phis_full, cd.nreal, gs.f_hat_0, gs.ep, &gs.D_hat_t );
		printf( "Total output variance = %g\n", gs.D_hat_t );
		/*		// Subtract f_hat_0 from phis and recalculate total output variance
				fhat2 = 0;
				for( count = 0; count < n_sub; count++ )
				{
					gs.f_a[count] -= gs.f_hat_0;
					gs.f_b[count] -= gs.f_hat_0;
					phis_full[ count ] = gs.f_a[count];
					phis_full[ n_sub + count ] = gs.f_b[count];
					fhat2 += pow( gs.f_a[count], 2 );
					fhat2 += pow( gs.f_b[count], 2 );
				}
				gs.D_hat_t = fhat2 / (2 * n_sub);
				printf( "Total output variance = %g\n", gs.D_hat_t );
				gs.D_hat_t = gsl_stats_variance( phis_full, 1, cd.nreal );
				printf( "Total output variance = %g\n", gs.D_hat_t );
		*/		free( phis_full );
		// Collect matrix of phis for fmat_a
		printf( "Computing phis for calculation of individual output variances:\n" );
		fflush( stdout );
		for( i = 0; i < pd.nOptParam; i++ )
		{
			printf( "Parameter %d...\n", i + 1 );
			for( count = 0; count < n_sub; count ++ )
			{
				for( j = 0; j < pd.nOptParam; j++ )
				{
					k = pd.var_index[j];
					if( i == j ) // then select from sample a
						opt_params[j] = pd.var[k] = gs.var_a_lhs[count][j];
					else // else select from sample b
						opt_params[j] = pd.var[k] = gs.var_b_lhs[count][j];
				}
				Transform( opt_params, &op, opt_params );
				func( opt_params, &op, od.res );
				// Save phi to fmat_a
				gs.fmat_a[i][count] = op.phi;
			}
		}
		// Collect matrix of phis for fmat_b
		printf( "Computing phis for calculation of individual plus interaction output variances:\n" );
		for( i = 0; i < pd.nOptParam; i++ )
		{
			printf( "Parameter %d...\n", i + 1 );
			for( count = 0; count < n_sub; count ++ )
			{
				for( j = 0; j < pd.nOptParam; j++ )
				{
					k = pd.var_index[j];
					if( i == j ) // then select from sample b
						opt_params[j] = pd.var[k] = gs.var_b_lhs[count][j];
					else // else select from sample a
						opt_params[j] = pd.var[k] = gs.var_a_lhs[count][j];
				}
				Transform( opt_params, &op, opt_params );
				func( opt_params, &op, od.res );
				// Save phi to fmat_b
				gs.fmat_b[i][count] = op.phi;
			}
		}
		printf( "done.\n" );
		// Calculate individual and interaction output variances
		for( i = 0; i < pd.nOptParam; i++ )
		{
			fhat2 = 0;
			for( j = 0; j < n_sub; j++ )
			{
				fhat2 += ( gs.f_a[j] * gs.fmat_a[i][j] );
				phis_half[ j ] = ( gs.f_a[j] * gs.fmat_a[i][j] );
			}
			gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			gs.D_hat[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			gs.D_hat[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_a[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			var_sorted( gs.f_a, gs.fmat_a[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat[i] );
			printf( "hat{D}_%d = %g\n", i, gs.D_hat[i] );
			//gs.D_hat[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
			fhat2 = 0;
			for( j = 0; j < n_sub; j++ )
			{
				fhat2 += ( gs.f_a[j] * gs.fmat_b[i][j] );
				phis_half[ j ] = ( gs.f_a[j] * gs.fmat_b[i][j] );
			}
			gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			gs.D_hat_n[i] = gsl_stats_mean( phis_half, 1, n_sub ) - pow( gs.f_hat_0, 2 );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			gs.D_hat_n[i] = gsl_stats_covariance_m( gs.f_a, 1, gs.fmat_b[i], 1, n_sub, gs.f_hat_0, gs.f_hat_0 );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			var_sorted( gs.f_a, gs.fmat_b[i], n_sub, gs.f_hat_0, gs.ep, &gs.D_hat_n[i] );
			printf( "hat{D}_n%d = %g\n", i, gs.D_hat_n[i] );
			//gs.D_hat_n[i] = ( fhat2 / n_sub ) - pow( gs.f_hat_0, 2 );
		}
		// Print sensitivity indices
		printf( "\nParameter sensitivity indices:\n" );
		printf( "parameter individual interaction\n" );
		for( i = 0; i < pd.nOptParam; i++ ) printf( "%d %g %g\n", i + 1, gs.D_hat[i] / gs.D_hat_t, 1 - ( gs.D_hat_n[i] / gs.D_hat_t ) );
		printf( "\n" );
		free( opt_params ); free( phis_half ); free( gs.f_a ); free( gs.f_b ); free( gs.D_hat ); free( gs.D_hat_n );
		free_matrix(( void ** ) gs.var_a_lhs, n_sub );
		free_matrix(( void ** ) gs.var_b_lhs, n_sub );
		free_matrix(( void ** ) gs.fmat_a, pd.nOptParam );
		free_matrix(( void ** ) gs.fmat_b, pd.nOptParam );
	}
//
// ------------------------ SIMPLE CALIBRATION
//
	if( cd.problem_type == CALIBRATE && cd.calib_type == SIMPLE ) /* Inverse analysis */
	{
		if( cd.nretries > 1 ) printf( "\nPARANOID CALIBRATION using a series of random initial guesses:\n" );
		else printf( "\nSINGLE CALIBRATION: single optimization based on initial guesses provided in the input file:\n" );
		if( strncasecmp( cd.opt_method, "lm", 2 ) == 0 ) f = optimize_lm; // Define optimization method: LM
		else f = optimize_pso; // Define optimization method: PSO
		for( i = 0; i < pd.nParam; i++ ) cd.var[i] = pd.var[i]; // Set all the initial values
		status = f( &op ); // Optimize
		if( status == 0 )
			{ printf( "ERROR: Optimization did not start! Optimization method mismatch!\n" ); sprintf( buf, "rm -f %s.running", op.root ); system( buf ); exit( 1 ); }
		sprintf( filename, "%s-rerun.mads", op.root );
		save_problem( filename, &cd, &pd, &od, &wd, &gd, &ed );
		printf( "\n" );
		print_results( &op );
		save_results( "", &op, &gd );
		if( od.nObs < od.nTObs )
		{
			predict = 1; // Produce outputs for predictions that are not calibration targets (see below PREDICT)
			sprintf( filename, "%s.results", op.root );
			out = Fappend( filename ); // Reopen results file
			printf( "\nModel predictions that are not calibration targets:\n" );
			fprintf( out, "\nModel predictions that are not calibration targets:\n" );
		}
		else predict = 0; // There are no observations that are not calibration targets
	}
	strcpy( op.label, "" ); // No labels needed below
//
// ------------------------ EIGEN || LOCALSENS
//
	if( cd.problem_type == EIGEN || cd.problem_type == LOCALSENS )
		eigen( &op, NULL, NULL ); // Eigen or sensitivity analysis run
//
//------------------------- ABAGUS
//
	if( cd.problem_type == ABAGUS ) // Particle swarm sensitivity analysis run
		status = pssa( &op ); // Optimize
//
//------------------------ POSTPUA
//
	if( cd.problem_type == POSTPUA ) // Predictive uncertainty analysis of sampling results
		status = postpua( &op );
//
// ------------------------ FORWARD
//
	if( cd.problem_type == FORWARD ) // Forward run
	{
		sprintf( filename, "%s.results", op.root );
		out = Fwrite( filename );
		fprintf( out, "Model parameter values:\n" );
		for( i = 0; i < pd.nParam; i++ )
			fprintf( out, "%s %g\n", pd.var_id[i], pd.var[i] );
		printf( "\nModel predictions (forward run; no calibration):\n" );
		fprintf( out, "\nModel predictions (forward run; no calibration):\n" );
		sprintf( filename, "%s.forward", op.root );
		out2 = Fwrite( filename );
	}
//
// ------------------------ CREATE
//
	if( cd.problem_type == CREATE ) // Create a MADS file based on a forward run
		printf( "\nModel predictions (forward run; no calibration):\n" );
//
// ------------------------ PREDICT
//
	if( cd.problem_type == FORWARD || cd.problem_type == CREATE || predict )
	{
		status_all = 1;
		compare = 0;
		phi = 0;
		if( cd.solution_type == EXTERNAL )
			for( i = 0; i < od.nObs; i++ )
			{
				if( cd.problem_type == CALIBRATE && od.obs_weight[i] != 0 ) continue;
				compare = 1;
				c = od.obs_current[i];
				err = od.obs_target[i] - c;
				phi += ( err * err ) * od.obs_weight[i];
				if( c < od.obs_min[i] || c > od.obs_max[i] ) { status_all = 0; status = 0; }
				else status = 1;
				if( od.nObs < 50 || ( i < 20 || i > od.nObs - 20 ) ) printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", od.obs_id[i], od.obs_target[i], c, err, err * od.obs_weight[i], status, od.obs_min[i], od.obs_max[i] );
				if( od.nObs > 50 && i == 21 ) printf( "...\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", od.obs_id[i], od.obs_target[i], c, err, err * od.obs_weight[i], status, od.obs_min[i], od.obs_max[i] );
				else od.obs_target[i] = c; // Save computed values as calibration targets
			}
		else if( cd.solution_type != TEST )
			for( i = 0; i < wd.nW; i++ )
				for( j = 0; j < wd.nWellObs[i]; j++ )
				{
					if( cd.problem_type == CALIBRATE && wd.obs_weight[i][j] > DBL_EPSILON ) continue;
					compare = 1;
					c = func_solver( wd.x[i], wd.y[i], wd.z1[i], wd.z2[i], wd.obs_time[i][j], &cd );
					err = wd.obs_target[i][j] - c;
					if( cd.problem_type != CALIBRATE ) phi += ( err * err ) * wd.obs_weight[i][j];
					else phi += ( err * err );
					if( c < wd.obs_min[i][j] || c > wd.obs_max[i][j] ) { status_all = 0; status = 0; }
					else status = 1;
					if( cd.problem_type != CALIBRATE )
						printf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err * wd.obs_weight[i][j], status, wd.obs_min[i][j], wd.obs_max[i][j] );
					else
						printf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err, status, wd.obs_min[i][j], wd.obs_max[i][j] );
					if( cd.problem_type != CREATE ) fprintf( out, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d\n", wd.id[i], wd.obs_time[i][j], wd.obs_target[i][j], c, err, err * wd.obs_weight[i][j], status );
					else wd.obs_target[i][j] = c; // Save computed values as calibration targets
					if( cd.problem_type == FORWARD ) fprintf( out2, "%s(%g) %g\n", wd.id[i], wd.obs_time[i][j], c ); // Forward run
				}
		cd.eval++;
		if( compare )
		{
			op.phi = phi;
			printf( "Objective function: %g Success: %d\n", op.phi, status_all );
			if( cd.problem_type != CREATE ) fprintf( out, "Objective function = %g Success: %d\n", op.phi, status_all );
			if( status_all )
			{
				printf( "All the predictions are within acceptable ranges!\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "All the predictions are within acceptable ranges!\n" );
			}
			else
			{
				printf( "At least one of the predictions is outside acceptable ranges!\n" );
				if( cd.problem_type != CREATE ) fprintf( out, "At least one of the predictions is outside acceptable ranges!\n" );
			}
		}
		else    printf( "No calibration targets!\n" );
		fclose( out );
	}
	if( cd.problem_type == FORWARD ) fclose( out2 );
	if( cd.problem_type == FORWARD || cd.problem_type == CALIBRATE )
	{
		if( gd.min_t > 0 && cd.solution_type != TEST )
		{
			printf( "\nCompute breakthrough curves at all the wells ...\n" );
			sprintf( filename, "%s.btc2", op.root );
			compute_btc2( filename, &op, &gd );
//			sprintf( filename, "%s.btc", root );
//			compute_btc( filename, &op, &gd );
		}
		if( gd.time > 0 && cd.solution_type != TEST )
		{
			printf( "\nCompute spatial distribution of predictions at t = %g ...\n", gd.time );
			sprintf( filename, "%s.vtk", op.root );
			compute_grid( filename, &cd, &gd );
		}
		sprintf( filename, "%s.phi", op.root );
		out2 = Fwrite( filename );
		fprintf( out2, "%g\n", op.phi ); // Write phi in a separate file
		fclose( out2 );
	}
	if( cd.problem_type == CREATE ) /* Create a file with calibration targets equal to the model predictions */
	{
		cd.problem_type = CALIBRATE;
		sprintf( filename, "%s-truth.mads", op.root );
		save_problem( filename, &cd, &pd, &od, &wd, &gd, &ed );
		printf( "\nMADS problem file named %s-truth.mads is created; modify the file if needed\n\n", op.root );
	}
	free( orig_params );
	// Finalize the run
	t2 = time( NULL );
	dt = t2 - t1;
	if( dt > 86400 ) printf( "Simulation time = %g days\n", (( double ) dt / 86400 ) );
	else if( dt > 3600 ) printf( "Simulation time = %g hours\n", (( double ) dt / 3600 ) );
	else if( dt > 60 ) printf( "Simulation time = %g minutes\n", (( double ) dt / 60 ) );
	else printf( "Simulation time = %ld seconds\n", dt );
	printf( "Functional evaluations = %d\n", cd.eval );
	if( cd.problem_type == CALIBRATE ) printf( "Levenberg-Marquardt optimizations = %d\n", cd.nlmo );
	if( dt > 0 )
	{
		c = cd.eval / dt;
		if( c < (( double ) 1 / 86400 ) ) printf( "Functional evaluations per day = %g\n", c * 86400 );
		else if( c < (( double ) 1 / 3600 ) ) printf( "Functional evaluations per hour = %g\n", c * 3600 );
		else if( c < (( double ) 1 / 60 ) ) printf( "Functional evaluations per minute = %g\n", c * 60 );
		else printf( "Functional evaluations per second = %g\n", c );
	}
	ptr_ts = gmtime( &t1 );
	printf( "Execution  started  on %s", asctime( ptr_ts ) );
	ptr_ts = gmtime( &t2 );
	printf( "Execution completed on %s", asctime( ptr_ts ) );
	printf( "Execution date & time stamp: %s\n", op.datetime_stamp );
	sprintf( buf, "rm -f %s.running", op.root ); system( buf );
	exit( 0 );
}

int optimize_pso( struct opt_data *op )
{
	if( op->cd->debug ) printf( "\nParticle-Swarm Optimization:" );
	if( strncasecmp( op->cd->opt_method, "pso", 3 ) == 0 || strcasestr( op->cd->opt_method, "std" ) != NULL || strncasecmp( op->cd->opt_method, "swarm", 5 ) == 0 )
	{
		if( op->cd->debug ) printf( " Standard (2006)\n" );
		pso_std( op );
	}
	else
	{
		if( op->cd->debug ) printf( " TRIBES\n" );
		pso_tribes( op );
	}
	if( op->cd->debug )
	{
		printf( "\n------------------------- Optimization Results:\n" );
		print_results( op );
	}
	if( op->cd->leigen ) eigen( op, NULL, NULL ); // Execute eigen analysis of the final results
	return( 1 );
}

int optimize_lm( struct opt_data *op )
{
	double phi, phi_min;
	double *opt_params, *res, *x_c;
	int   nsig, maxfn, maxiter, iopt, infer, ier, debug, standalone;
	int   i, j, k, debug_level, count, count_set, npar;
	double opt_parm[4], *jacobian, *jacTjac, *covar, *work, eps, delta, *var_lhs;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	char buf[80];
	if( op->od->nObs == 0 ) { printf( "ERROR: Number of observations is equal to zero! Levenberg-Marquardt Optimization cannot be performed!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if( op->pd->nOptParam == 0 ) { printf( "ERROR: Number of optimized model parameters is equal to zero! Levenberg-Marquardt Optimization cannot be performed!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	gsl_matrix *gsl_jacobian = gsl_matrix_alloc( op->od->nObs, op->pd->nOptParam );
	gsl_matrix *gsl_covar = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam );
	gsl_vector *gsl_opt_params = gsl_vector_alloc( op->pd->nOptParam );
	if(( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if(( x_c = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if(( res = ( double * ) malloc( op->od->nObs * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	debug = op->cd->debug;
	standalone = op->cd->standalone;
	if( standalone && !debug ) { printf( "Levenberg-Marquardt Optimization ... " ); fflush( stdout ); }
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]];
	if( op->cd->paranoid )
	{
		npar = op->pd->nOptParam;
		if(( var_lhs = ( double * ) malloc( npar * op->cd->nretries * sizeof( double ) ) ) == NULL )
			{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
		if( op->cd->seed < 0 ) { op->cd->seed *= -1; printf( "Imported seed: %d\n", op->cd->seed ); }
		else if( op->cd->seed == 0 ) { printf( "New " ); op->cd->seed = get_seed(); }
		else printf( "Current seed: %d\n", op->cd->seed );
		printf( "Random sampling for paranoid optimization (variables %d; realizations %d) using ", npar, op->cd->nretries );
		if( op->cd->paran_method[0] != 0 ) { strcpy( buf, op->cd->smp_method ); strcpy( op->cd->smp_method, op->cd->paran_method ); }
		sampling( npar, op->cd->nretries, &op->cd->seed, var_lhs, op );
		if( op->cd->paran_method[0] != 0 ) strcpy( op->cd->smp_method, buf );
		printf( "done.\n" );
		count = count_set = 0;
	}
	phi_min = HUGE;
	do // BEGIN Paranoid loop
	{
		if( op->cd->maxeval <= op->cd->eval ) { if( debug || op->cd->paranoid ) printf( "Maximum number of evaluations is exceeded (%d <= %d)!\n", op->cd->maxeval, op->cd->eval ); break; }
		op->cd->nlmo++;
		if( op->cd->paranoid )
		{
			count++;
			if( op->cd->calib_type == IGRND && count == 1 )
				printf( "CALIBRATION %d: initial guesses from the provided IGRND random set: ", count );
			else
			{
				printf( "CALIBRATION %d: initial guesses from internal paranoid random set #%d: ", count, count_set + 1 );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					k = op->pd->var_index[i];
					opt_params[i] = var_lhs[i+count_set*npar] * op->pd->var_range[k] + op->pd->var_min[k];
					if( debug )
					{
						if( op->pd->var_log[k] ) printf( "%s %.15g\n", op->pd->var_id[k], pow( 10, opt_params[i] ) );
						else printf( "%s %.15g\n", op->pd->var_id[k], opt_params[i] );
					}
				}
				count_set++;
			}
			if( debug ) printf( "\n" );
			fflush( stdout );
		}
		if( standalone ) Transform( opt_params, op, opt_params ); // Transform if standalone; do not tranform is part of PSO run
		if( debug && standalone )
		{
			printf( "\n-------------------- Initial state:\n" );
			op->cd->pderiv = op->cd->oderiv = -1;
			debug_level = op->cd->fdebug; op->cd->fdebug = 3;
			op->cd->compute_phi = 1;
			func( opt_params, op, res );
			op->cd->compute_phi = 0;
			op->cd->fdebug = debug_level;
		}
		// LM optimization ...
		if( strcasestr( op->cd->opt_method, "mon" ) != NULL || strcasestr( op->cd->opt_method, "chav" ) != NULL ) // Monty/Chavo versions
		{
			if( debug && standalone ) printf( "\nLevenberg-Marquardt Optimization:\n" );
			if(( jacobian = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->od->nObs ) ) == NULL )
				{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
			if(( jacTjac = ( double * ) malloc( sizeof( double ) * (( op->pd->nOptParam + 1 ) * op->pd->nOptParam / 2 ) ) ) == NULL )
				{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
			iopt = 2; /*    iopt=0 Brown's algorithm without strict descent
			       		iopt=1 strict descent and default values for input vector parm
					iopt=2 strict descent with user parameter choices in input vector parm */
			opt_parm[0] = 10; /* initial value of the Marquardt parameter */
			opt_parm[1] = 2.0; /* scaling factor used to modify the Marquardt parameter */
			opt_parm[2] = 1e30; /* upper bound for increasing the Marquardt parameter */
			opt_parm[3] = 500; /* value for indicating when central differencing is to be used for calculating the jacobian */
			/* First convergence criterion */
			nsig = 8; /* parameter estimates agree to nsig digits */
			/* Second convergence criterion */
			eps = epsilon(); /* two successive iterations the residual sum of squares estimates have relative difference less than or equal to eps */
			/* Third convergence criterion */
			delta = 0; /* norm of the approximate gradient is less than or equal to delta */
			maxfn = op->cd->maxeval - op->cd->eval; /* maximum number of function evaluations; remove the number of evaluation already performed */
			maxiter = op->cd->niter; /* maximum number of iterations */
			if( strcasestr( op->cd->opt_method, "chav" ) != NULL ) ier = zxssqch( func, op, op->od->nObs, op->pd->nOptParam, nsig, eps, delta, maxfn, iopt, opt_parm, opt_params, &phi, res, jacobian, op->od->nObs, jacTjac, &infer ); // Chavo's version
			else ier = lm_opt( func, func_dx, op, op->od->nObs, op->pd->nOptParam, nsig, eps, delta, maxfn, maxiter, iopt, opt_parm, opt_params, &phi, res, jacobian, op->od->nObs, jacTjac, &infer ); // Monty's version
			for( k = i = 0; i < op->pd->nOptParam; i++ )
				for( j = 0; j < op->od->nObs; j++, k++ )
					gsl_matrix_set( gsl_jacobian, j, i, jacobian[k] );
			gsl_multifit_covar( gsl_jacobian, 0.0, gsl_covar );
			free( jacTjac ); free( jacobian );
			op->phi = phi;
		}
		else if( strcasestr( op->cd->opt_method, "gsl" ) != NULL ) // GSL version of LM
		{
			if( debug && standalone ) printf( "\nLevenberg-Marquardt Optimization using GSL library:\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
				gsl_vector_set( gsl_opt_params, i, opt_params[i] );
			lm_gsl( gsl_opt_params, op, gsl_jacobian, gsl_covar );
			for( i = 0; i < op->pd->nOptParam; i++ )
				opt_params[i] = gsl_vector_get( gsl_opt_params, i );
			phi = op->phi;
		}
		else // DEFAULT LevMar version of LM
		{
			if(( covar = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->pd->nOptParam ) ) == NULL )
				{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
			// LM_DIF_WORKSZ(m,n) = 4*n+4*m + n*m + m*m
			if(( work = ( double * ) malloc( sizeof( double ) * LM_DIF_WORKSZ( op->pd->nOptParam, op->od->nObs ) ) ) == NULL )
				{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
			for( i = 0; i < op->od->nObs; i++ ) res[i] = 0;
			jacobian = work + op->pd->nOptParam + 2 * op->od->nObs;
			if( debug && standalone ) printf( "\nLevenberg-Marquardt Optimization using LevMar library:\n" );
			else if( op->cd->ldebug ) printf( "\n" );
			opts[0] = 1e-3; opts[1] = 1E-5; opts[2] = 1E-5;
			opts[3] = op->cd->phi_cutoff;
			if( op->cd->sintrans == 0 ) opts[4] = 0.1; // Forward difference; Central difference if negative; DO NOT USE CENTRAL DIFFERENCE
			else opts[4] = op->cd->sindx;
			while( op->cd->maxeval > op->cd->eval )
			{
				// Levmar has no termination creteria based on the number of functional evaluations or number of jacobian evaluations
				if( opts[4] > 0 ) maxiter = ( double )( op->cd->maxeval - op->cd->eval ) / ( op->pd->nOptParam + 10 ) + 1; // Forward derivatives
				else              maxiter = ( double )( op->cd->maxeval - op->cd->eval ) / ( 2 * op->pd->nOptParam + 10 ) + 1; // Central derivatives
				if( maxiter > op->cd->niter ) maxiter = op->cd->niter;
				maxiter *= 10; // Assuming about 10 lambda searches per iteration
				if( strcasestr( op->cd->opt_method, "dif" ) != NULL ) ier = dlevmar_dif( func_levmar, opt_params, res, op->pd->nOptParam, op->od->nObs, maxiter, opts, info, work, covar, op );
				else ier = dlevmar_der( func_levmar, func_dx_levmar, opt_params, res, op->pd->nOptParam, op->od->nObs, maxiter, opts, info, work, covar, op );
				if( info[6] == 4 || info[6] == 5 ) { opts[0] *= 10; if( op->cd->ldebug ) printf( "Rerun with larger initial lambda (%g)\n", opts[0] ); }
				else break;
			}
			if( op->cd->ldebug > 1 )
			{
				printf( "Levenberg-Marquardt Optimization completed after %g iteration (reason %g) (returned value %d)\n", info[5], info[6], ier );
				printf( "initial phi %g final phi %g ||J^T e||_inf %g ||Dp||_2 %g mu/max[J^T J]_ii %g\n", info[0], info[1], info[2], info[3], info[4] );
				printf( "function evaluation %g jacobian evaluations %g linear systems solved %g\n", info[7], info[8], info[9] );
			}
			for( k = j = 0; j < op->od->nObs; j++ )
				for( i = 0; i < op->pd->nOptParam; i++ )
					gsl_matrix_set( gsl_jacobian, j, i, jacobian[k++] );
			for( i = 0; i < op->pd->nOptParam; i++ )
				for( j = 0; j < op->pd->nOptParam; j++ )
					gsl_matrix_set( gsl_covar, i, j, covar[i * op->pd->nOptParam + j] );
			op->phi = phi = info[1];
			free( work ); free( covar );
		}
		if( standalone ) // if LM is stand alone (not part of PSO run)
		{
			if( debug )
			{
				printf( "\n------------------------- Final state:\n" );
				op->cd->pderiv = op->cd->oderiv = -1;
				debug_level = op->cd->fdebug; op->cd->fdebug = 3;
				op->cd->compute_phi = 1;
				func( opt_params, op, op->od->res ); // opt_params are already transformed
				op->cd->compute_phi = 0;
				op->cd->fdebug = debug_level;
			}
			else
			{
				// Make a Forward run with the best results
				op->cd->compute_phi = 1;
				func( opt_params, op, op->od->res ); // opt_params are already transformed
				op->cd->compute_phi = 0;
			}
			DeTransform( opt_params, op, x_c );
			for( i = 0; i < op->pd->nOptParam; i++ )
				op->pd->var[op->pd->var_index[i]] = x_c[i]; // Save the obtained results
			if( debug )
			{
				printf( "\n------------------------- LM Optimization Results:\n" );
				print_results( op );
			}
		}
		else // if LM is part of PSO run
			for( i = 0; i < op->pd->nOptParam; i++ )
				op->pd->var[op->pd->var_index[i]] = opt_params[i];
		if( op->cd->paranoid )
		{
			if( op->phi < phi_min ) { phi_min = op->phi; for( i = 0; i < op->pd->nOptParam; i++ ) op->pd->var_best[i] = op->pd->var[op->pd->var_index[i]]; }
			if( debug == 0 ) printf( "Objective function: %g Success %d\n", op->phi, op->success );
			if( phi_min < op->cd->phi_cutoff )
				{ printf( "Calibration objective function is below the cutoff value after %d random initial guess runs\n", count ); break; }
			if( op->cd->check_success && op->success )
				{ printf( "Calibration within calibration ranges after %d random initial guess runs\n", count ); break; }
			else if( count == op->cd->nretries )
				{ printf( "Calibration attempts terminated after %d random initial guess runs\n", count ); break; }
		}
		else break; // Quit if not Paranoid run
	}
	while( 1 ); // END Paranoid loop
	if( op->cd->paranoid ) // Recompute for the best results
	{
		op->phi = phi_min;
		for( i = 0; i < op->pd->nOptParam; i++ )
			op->pd->var[op->pd->var_index[i]] = op->pd->var_best[i];
		Transform( op->pd->var_best, op, opt_params );
		op->cd->compute_phi = 1;
		func( opt_params, op, op->od->res );
		op->cd->compute_phi = 0;
	}
	if(( op->cd->leigen || debug ) && standalone ) eigen( op, gsl_jacobian, gsl_covar );  // Eigen analysis
	if( op->cd->paranoid ) free( var_lhs );
	free( opt_params ); free( x_c ); free( res );
	gsl_matrix_free( gsl_jacobian ); gsl_matrix_free( gsl_covar ); gsl_vector_free( gsl_opt_params );
	if( standalone && !debug && op->cd->problem_type != CALIBRATE ) printf( "\n" );
	return( 1 );
}

int eigen( struct opt_data *op, gsl_matrix *gsl_jacobian, gsl_matrix *gsl_covar )
{
	FILE *out;
	double phi, dof, stddev_scale, gf;
	double *opt_params, *x_u, *x_d, *stddev, *jacobian;
	int   debug, compute_covar, compute_jacobian;
	int   i, j, k, ier, debug_level, status;
	double eps;
	char filename[200], buf[20];
	gsl_vector *gsl_opt_params = gsl_vector_alloc( op->pd->nOptParam );
	gsl_matrix *eigenvec = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam );
	gsl_vector *eigenval = gsl_vector_alloc( op->pd->nOptParam );
	gsl_eigen_symmv_workspace *eigenwork = gsl_eigen_symmv_alloc( op->pd->nOptParam );
	if(( jacobian = ( double * ) malloc( sizeof( double ) * op->pd->nOptParam * op->od->nObs ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	compute_jacobian = compute_covar = 0;
	if( gsl_jacobian == NULL ) { gsl_jacobian = gsl_matrix_alloc( op->od->nObs, op->pd->nOptParam ); compute_jacobian = 1; }
	if(( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if(( x_u = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if(( x_d = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	if(( stddev = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 1 ); }
	debug = op->cd->debug;
	for( i = 0; i < op->pd->nOptParam; i++ )
		opt_params[i] = op->pd->var[op->pd->var_index[i]];
	Transform( opt_params, op, opt_params );
	for( i = 0; i < op->pd->nOptParam; i++ )
		gsl_vector_set( gsl_opt_params, i, opt_params[i] );
	printf( "\nEigen analysis ...\n" );
	op->cd->pderiv = op->cd->oderiv = -1;
	if( compute_jacobian )
	{
		/*		op->pd->var_current_gsl = gsl_vector_alloc( op->pd->nOptParam );
				op->od->obs_current_gsl = gsl_vector_alloc(op->od->nObs);
				func_gsl_deriv_dx( gsl_opt_params, op, gsl_jacobian ); // Using GSL function
				gsl_vector_free( op->od->obs_current_gsl );
				gsl_vector_free( op->pd->var_current_gsl ); */
		func_gsl_dx( gsl_opt_params, op, gsl_jacobian ); // Compute Jacobian using forward difference
		func_dx( opt_params, NULL, op, jacobian ); // Compute Jacobian using forward difference
	}
	if( debug )
	{
		printf( "Analyzed state:\n" );
		debug_level = op->cd->fdebug; op->cd->fdebug = 3;
	}
	op->cd->compute_phi = 1;
	func( opt_params, op, op->od->res );
	op->cd->compute_phi = 0;
	phi = op->phi;
	if( debug )
	{
		printf( "\nJacobian matrix\n" ); // Print Jacobian
		printf( "%-25s :", "Observations" );
		for( k = 0; k < op->od->nObs; k++ )
		{
			if( op->od->nObs < 30 || ( k < 10 || k > op->od->nObs - 10 ) ) printf( " %s", op->od->obs_id[k] );
			if( op->od->nObs > 30 && k == 11 ) printf( " ..." );
		}
		printf( "\n" );
		for( k = i = 0; i < op->pd->nOptParam; i++ )
		{
			printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->od->nObs; j++ )
			{
				eps = gsl_matrix_get( gsl_jacobian, j, i );
				if( fabs( eps ) > 1e3 ) sprintf( buf, " %6.0e", eps );
				else sprintf( buf, " %6.2f", eps );
				if( op->od->nObs < 30 || ( j < 10 || j > op->od->nObs - 10 ) ) printf( " %s", buf );
				if( op->od->nObs > 30 && j == 11 ) printf( " ..." );
			}
			printf( "\n" );
			/*
			printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->od->nObs; j++, k++ )
			{
				eps = jacobian[k];
				if( fabs( eps ) > 1e3 ) sprintf( buf, " %6.0e", eps );
				else sprintf( buf, " %6.2f", eps );
				if( op->od->nObs < 30 || ( j < 10 || j > op->od->nObs - 10 ) ) printf( " %s", buf );
				if( op->od->nObs > 30 && j == 11 ) printf( " ..." );
			}
			printf( "\n" );
			*/
		}
	}
	sprintf( filename, "%s", op->root );
	if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
	if( op->s > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->s );
	strcat( filename, ".jacobian" );
	out = Fwrite( filename );
	fprintf( out, "%-25s :", "Parameters" );
	for( i = 0; i < op->pd->nOptParam; i++ )
		fprintf( out, " \"%s\"", op->pd->var_id[op->pd->var_index[i]] );
	fprintf( out, "\n" );
	for( j = 0; j < op->od->nObs; j++ )
	{
		fprintf( out, "%-25s :", op->od->obs_id[j] );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fprintf( out, " %g", gsl_matrix_get( gsl_jacobian, j, i ) );
		fprintf( out, "\n" );
	}
	fclose( out );
	printf( "Jacobian matrix stored (%s)\n", filename );
	if( op->cd->problem_type == EIGEN || op->cd->leigen )
	{
		if( gsl_covar == NULL ) { gsl_covar = gsl_matrix_alloc( op->pd->nOptParam, op->pd->nOptParam ); compute_covar = 1; }
		if( compute_covar ) // Standalone eigen analysis
		{
			ier = gsl_multifit_covar( gsl_jacobian, 0.0, gsl_covar );
			if( ier != GSL_SUCCESS ) { printf( "Problem computing covariance matrix!\n" ); ier = 1; }
			else ier = 0;
		}
		if( debug )
		{
			printf( "\nCovariance matrix\n" );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					printf( " %7.0e", gsl_matrix_get( gsl_covar, i, j ) );
				printf( "\n" );
			}
		}
		sprintf( filename, "%s", op->root );
		if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
		if( op->s > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->s );
		strcat( filename, ".covariance" );
		out = Fwrite( filename );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
			for( j = 0; j < op->pd->nOptParam; j++ )
				fprintf( out, " %g", gsl_matrix_get( gsl_covar, i, j ) );
			fprintf( out, "\n" );
		}
		fclose( out );
		printf( "Covariance matrix stored (%s)\n", filename );
		for( ier = 0, i = 0; i < op->pd->nOptParam; i++ )
		{
			stddev[i] = sqrt( gsl_matrix_get( gsl_covar, i, i ) ); // compute standard deviations before the covariance matrix is destroyed by eigen functions
			if( stddev[i] < DBL_EPSILON ) ier = 1;
		}
		if( ier == 0 )
		{
			if( debug )
			{
				printf( "\nCorrelation matrix\n" );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					printf( "%-25s : ", op->pd->var_id[op->pd->var_index[i]] );
					for( j = 0; j < op->pd->nOptParam; j++ )
						printf( " %6.3f", gsl_matrix_get( gsl_covar, i, j ) / ( stddev[i] * stddev[j] ) );
					printf( "\n" );
				}
			}
			sprintf( filename, "%s", op->root );
			if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
			if( op->s > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->s );
			strcat( filename, ".correlation" );
			out = Fwrite( filename );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					fprintf( out, " %g", gsl_matrix_get( gsl_covar, i, j ) / ( stddev[i] * stddev[j] ) );
				fprintf( out, "\n" );
			}
			fclose( out );
			printf( "Correlation matrix stored (%s)\n", filename );
			// GSL_COVAR is destroyed during eigen computation
			gsl_eigen_symmv( gsl_covar, eigenval, eigenvec, eigenwork );
			if( debug )
			{
				printf( "\nEigenvectors (sorted by absolute values of eigenvalues)\n" );
				gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_ABS_ASC );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
					for( j = 0; j < op->pd->nOptParam; j++ )
						printf( " %6.3f", gsl_matrix_get( eigenvec, i, j ) );
					printf( "\n" );
				}
				printf( "%-25s :", "Eigenvalues" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					printf( " %6.0e", gsl_vector_get( eigenval, i ) );
				printf( "\n" );
				printf( "\nEigenvectors (sorted by eigenvalues)\n" );
				gsl_eigen_symmv_sort( eigenval, eigenvec, GSL_EIGEN_SORT_VAL_ASC );
				for( i = 0; i < op->pd->nOptParam; i++ )
				{
					printf( "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
					for( j = 0; j < op->pd->nOptParam; j++ )
						printf( " %6.3f", gsl_matrix_get( eigenvec, i, j ) );
					printf( "\n" );
				}
				printf( "%-25s :", "Eigenvalues" );
				for( i = 0; i < op->pd->nOptParam; i++ )
					printf( " %6.0e", gsl_vector_get( eigenval, i ) );
				printf( "\n" );
			}
			sprintf( filename, "%s", op->root );
			if( op->label[0] != 0 ) sprintf( filename, "%s.%s", filename, op->label );
			if( op->s > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->s );
			strcat( filename, ".eigen" );
			out = Fwrite( filename );
			for( i = 0; i < op->pd->nOptParam; i++ )
			{
				fprintf( out, "%-25s :", op->pd->var_id[op->pd->var_index[i]] );
				for( j = 0; j < op->pd->nOptParam; j++ )
					fprintf( out, " %g", gsl_matrix_get( eigenvec, i, j ) );
				fprintf( out, "\n" );
			}
			fprintf( out, "%-25s :", "Eigenvalues" );
			for( i = 0; i < op->pd->nOptParam; i++ )
				fprintf( out, " %g", gsl_vector_get( eigenval, i ) );
			fprintf( out, "\n" );
			fclose( out );
			printf( "Eigen vactors and eigen values stored (%s)\n", filename );
		}
		else
			printf( "Correlation matrix and eigen vectors cannot be computed!\n" );
		dof = op->od->nObs - op->pd->nOptParam;
		stddev_scale = sqrt( phi / dof );
		gf = phi / dof;
		printf( "\nObtained fit is " );
		if( gf > 200 ) printf( "not very good (chi^2/dof = %g > 200)\n", gf );
		else printf( "relatively good (chi^2/dof = %g < 200)\n", gf );
		printf( "\nOptimized parameters:\n" );
		if( debug ) printf( "Transformed space (applied during optimization):\n" );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			if( op->cd->sintrans == 1 ) opt_params[i] = asin( sin( opt_params[i] ) );
			stddev[i] *= stddev_scale;
			x_u[i] = opt_params[i] + ( double ) 3 * stddev[i];
			x_d[i] = opt_params[i] - ( double ) 3 * stddev[i];
			status = 0;
			if( op->cd->sintrans == 1 )
			{
				if( x_d[i] < -M_PI / 2 ) { status = 1; x_d[i] = -M_PI / 2; }
				if( x_u[i] > M_PI / 2 ) { status = 1; x_u[i] = M_PI / 2; }
			}
			else
			{
				if( x_d[i] < op->pd->var_min[k] ) { status = 1; x_d[i] = op->pd->var_min[k]; }
				if( x_u[i] > op->pd->var_max[k] ) { status = 1; x_u[i] = op->pd->var_max[k]; }
			}
			if( debug )
			{
				printf( "%-40s : %12g stddev %12g (%12g - %12g)", op->pd->var_id[k], opt_params[i], stddev[i], x_d[i], x_u[i] );
				if( status ) printf( " Estimated ranges are constrained by prior uncertainty bounds\n" );
				else printf( "\n" );
			}
		}
		DeTransform( x_u, op, x_u );
		DeTransform( x_d, op, x_d );
		if( debug ) printf( "Untransformed space:\n" );
		for( i = 0; i < op->pd->nOptParam; i++ )
		{
			k = op->pd->var_index[i];
			printf( "%-40s : ", op->pd->var_id[k] );
			if( op->pd->var_log[k] == 0 ) printf( "%12g stddev %12g (%12g - %12g)", op->pd->var[k], stddev[i], x_d[i], x_u[i] );
			else printf( "%12g stddev %12g (%12g - %12g)", pow( 10, op->pd->var[k] ), stddev[k], pow( 10, x_d[i] ), pow( 10, x_u[i] ) );
			status = 0;
			if( x_d[i] <= op->pd->var_min[k] ) { status = 1; x_d[i] = op->pd->var_min[k]; }
			if( x_u[i] >= op->pd->var_max[k] ) { status = 1; x_u[i] = op->pd->var_max[k]; }
			if( status ) printf( " Estimated ranges are constrained by prior uncertainty bounds\n" );
			else printf( "\n" );
		}
	}
	free( opt_params ); free( stddev ); free( x_u ); free( x_d ); free( jacobian );
	gsl_vector_free( gsl_opt_params );
	gsl_matrix_free( eigenvec ); gsl_vector_free( eigenval ); gsl_eigen_symmv_free( eigenwork );
	if( compute_jacobian ) gsl_matrix_free( gsl_jacobian );
	if( compute_covar ) gsl_matrix_free( gsl_covar );
	return( 1 );
}

int postpua( struct opt_data *op )
{
	FILE *fl, *outfl;
	double *opt_params, of;
	char buf[80], filename[80];
	int i, n;
	if( op->cd->infile[0] == 0 ) { printf( "\nInfile must be specified for postpua run\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 0 );}
	fl = fopen( op->cd->infile, "r" );
	if( fl == NULL ) { printf( "\nError opening %s\n", op->cd->infile ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 0 ); }
	sprintf( filename, "%s.pua", op->root );
	outfl = fopen( filename , "w" );
	if( outfl == NULL ) { printf( "\nError opening %s\n", filename ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 0 ); }
	printf( "\nComputing predictions for %s...", op->cd->infile );
	fflush( stdout );
	if(( opt_params = ( double * ) malloc( op->pd->nOptParam * sizeof( double ) ) ) == NULL )
		{ printf( "Not enough memory!\n" ); sprintf( buf, "rm -f %s.running", op->root ); system( buf ); exit( 0 ); }
	fgets( buf, sizeof buf, fl ); // Skip header
	fprintf( outfl, "Number       OF           " );
	for( i = 0; i < op->od->nObs; i++ )
		fprintf( outfl, " %-12s", op->od->obs_id[i] );
	fprintf( outfl, "\n" );
	while( fscanf( fl, "%d %lf", &n, &of ) > 0 )
	{
		fprintf( outfl, "%-12d %-12lf ", n, of );
		for( i = 0; i < op->pd->nOptParam; i++ )
			fscanf( fl, "%lf", &opt_params[i] );
		fscanf( fl, " \n" );
		func( opt_params, op, op->od->res );
		for( i = 0; i < op->od->nObs; i++ )
			fprintf( outfl, " %-12g", op->od->obs_current[i] );
		fprintf( outfl, "\n" );
	}
	fclose( fl );
	fclose( outfl );
	printf( "Done\n" );
	printf( "Results written to %s\n\n", filename );
	return 1;
}

// Modified from Numerical Recipes in C: The Art of Scientific Computing (ISBN 0-521-43108-5)
// corrected three-pass algorithm to minimize roundoff error in variance
void var_sorted( double data[], double datb[], int n, double ave, double ep, double *var )
{
	int j;
	double dev2[n];
	// First pass to calculate mean
//	s = 0.0;
//	for( j = 0; j < n; j++ ) s += data[j];
//	*ave = s/n;
	// Second pass to calculate absolute deviations
	for( j = 0; j < n; j++ )
		dev2[j] = ( data[j] - ave ) * ( datb[j] - ave );
	// Sort devs
	gsl_sort( dev2, 1, n );
	// Third pass to calculate first (absolute) and second moments
	*var = 0.0;
	for( j = 0; j < n; j++ )
		*var += dev2[j];
	*var = ( *var - ep * ep / n ) / ( n - 1 );
}

void ave_sorted( double data[], int n, double *ave, double *ep )
{
	int j;
	double s, dev[n];
	// First pass to calculate mean
	s = 0.0;
	for( j = 0; j < n; j++ ) s += data[j];
	*ave = s / n;
	// Second pass to calculate absolute deviations
	for( j = 0; j < n; j++ )
		dev[j] = data[j] - *ave;
	// Sort devs
	gsl_sort( dev, 1, n );
	// Third pass to calculate first (absolute) moment
	*ep = 0.0;
	for( j = 0; j < n; j++ )
		*ep += dev[j];
}

void sampling( int npar, int nreal, int *seed, double var_lhs[], struct opt_data *op )
{
	printf( "%s\n", op->cd->smp_method );
	if( nreal == 1 || strncasecmp( op->cd->smp_method, "random", 6 ) == 0 )
	{
		printf( "Pure random sampling method ... " );
		fflush( stdout );
		smp_random( npar, nreal, seed, var_lhs );
	}
	else if(( nreal <= 500 && op->cd->smp_method[0] == 0 ) || strncasecmp( op->cd->smp_method, "idlhs", 5 ) == 0 )
	{
		printf( "Improved Distances LHS method " );
		if( strncasecmp( op->cd->smp_method, "idlhs", 5 ) != 0 ) printf( "( real < 500 ) " );
		printf( "... " );
		fflush( stdout );
		lhs_imp_dist( npar, nreal, 5, seed, var_lhs );
	}
	else if(( nreal > 500 && op->cd->smp_method[0] == 0 ) || strncasecmp( op->cd->smp_method, "lhs", 3 ) == 0 )
	{
		printf( "Standard LHS method " );
		if( strncasecmp( op->cd->smp_method, "lhs", 3 ) != 0 ) printf( "( real > 500 ) " );
		printf( "... " );
		fflush( stdout );
		lhs_random( npar, nreal, seed, var_lhs );
	}
}

void print_results( struct opt_data *op )
{
	int i, j, k, status, status_all;
	double c, err;
	status_all = 1;
	printf( "Model parameters:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->pd->var_log[k] == 0 ) printf( "%s %g\n", op->pd->var_id[k], op->pd->var[k] );
		else printf( "%s %g\n", op->pd->var_id[k], pow( 10, op->pd->var[k] ) );
	}
	if( op->od->nObs > 0 ) printf( "\nCalibration targets:\n" );
	if( op->cd->solution_type == EXTERNAL )
		for( i = 0; i < op->od->nObs; i++ )
		{
			if( op->od->obs_weight[i] == 0 ) continue;
			c = op->od->obs_current[i];
			err = op->od->obs_target[i] - c;
			if( c < op->od->obs_min[i] || c > op->od->obs_max[i] ) { status_all = 0; status = 0; }
			else status = 1;
			if( op->od->nObs < 50 || ( i < 20 || i > op->od->nObs - 20 ) )
				printf( "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], status, op->od->obs_min[i], op->od->obs_max[i] );
			if( op->od->nObs > 50 && i == 21 ) printf( "...\n" );
		}
	else
	{
		for( k = 0, i = 0; i < op->wd->nW; i++ )
			for( j = 0; j < op->wd->nWellObs[i]; j++ )
			{
				if( op->wd->obs_weight[i][j] == 0 ) continue;
				c = op->od->obs_current[k++];
				err = op->wd->obs_target[i][j] - c;
				if( c < op->wd->obs_min[i][j] || c > op->wd->obs_max[i][j] ) { status_all = 0; status = 0; }
				else status = 1;
				printf( "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], status, op->wd->obs_min[i][j], op->wd->obs_max[i][j] );
			}
	}
	op->success = status_all;
	printf( "Objective function: %g Success: %d \n", op->phi, op->success );
	if( status_all ) printf( "All the predictions are within calibration ranges!\n" );
	else printf( "At least one of the predictions is outside calibration ranges!\n" );
	printf( "Number of function evaluations = %d\n", op->cd->eval );
}

void save_results( char *label, struct opt_data *op, struct grid_data *gd )
{
	FILE *out, *out2;
	int i, j, k, status, status_all;
	double c, err;
	char filename[80], f[80];
	status_all = 1;
	sprintf( filename, "%s", op->root );
	if( label[0] != 0 ) sprintf( filename, "%s.%s", filename, label );
	if( op->s > 0 && op->cd->nreal > 1 ) sprintf( filename, "%s-%08d", filename, op->s );
	strcpy( f, filename );
	strcat( filename, ".results" );
	out = Fwrite( filename );
	strcpy( filename, f );
	strcat( filename, ".residuals" );
	out2 = Fwrite( filename );
	fprintf( out, "Optimized parameter values:\n" );
	for( i = 0; i < op->pd->nOptParam; i++ )
	{
		k = op->pd->var_index[i];
		if( op->pd->var_log[k] == 0 ) fprintf( out, "%s %g\n", op->pd->var_id[k], op->pd->var[k] );
		else fprintf( out, "%s %g\n", op->pd->var_id[k], pow( 10, op->pd->var[k] ) );
	}
	if( op->od->nObs > 0 ) fprintf( out, "\nOptimized calibration targets:\n" );
	if( op->cd->solution_type == EXTERNAL )
		for( i = 0; i < op->od->nObs; i++ )
		{
			if( op->od->obs_weight[i] != 0 )
			{
				c = op->od->obs_current[i];
				err = op->od->obs_target[i] - c;
				if( c < op->od->obs_min[i] || c > op->od->obs_max[i] ) { status_all = 0; status = 0; }
				else status = 1;
			}
			else status = 0;
			fprintf( out, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], status, op->od->obs_min[i], op->od->obs_max[i] );
			fprintf( out2, "%-20s:%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->od->obs_id[i], op->od->obs_target[i], c, err, err * op->od->obs_weight[i], status, op->od->obs_min[i], op->od->obs_max[i] );
		}
	else
	{
		for( k = 0, i = 0; i < op->wd->nW; i++ )
			for( j = 0; j < op->wd->nWellObs[i]; j++ )
			{
				if( op->wd->obs_weight[i][j] != 0 )
				{
					c = op->od->obs_current[k++];
					err = op->wd->obs_target[i][j] - c;
					if( c < op->wd->obs_min[i][j] || c > op->wd->obs_max[i][j] ) { status_all = 0; status = 0; }
					else status = 1;
				}
				else status = 0;
				fprintf( out, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], status, op->wd->obs_min[i][j], op->wd->obs_max[i][j] );
				fprintf( out2, "%-10s(%5g):%12g - %12g = %12g (%12g) success %d range %12g - %12g\n", op->wd->id[i], op->wd->obs_time[i][j], op->wd->obs_target[i][j], c, err, err * op->wd->obs_weight[i][j], status, op->wd->obs_min[i][j], op->wd->obs_max[i][j] );
			}
	}
	op->success = status_all;
	fprintf( out, "Objective function: %g Success: %d \n", op->phi, op->success );
	if( status_all ) fprintf( out, "All the predictions are within calibration ranges!\n" );
	else fprintf( out, "At least one of the predictions is outside calibration ranges!\n" );
	fprintf( out, "Number of function evaluations = %d\n", op->cd->eval );
	if( op->cd->seed > 0 ) fprintf( out, "Seed = %d\n", op->cd->seed );
	fclose( out ); fclose( out2 );
	if( gd->min_t > 0 )
	{
		printf( "\nCompute breakthrough curves at all the wells ..." );
		fflush( stdout );
		sprintf( filename, "%s.btc", f );
		compute_btc2( filename, op, gd );
		printf( "done.\n" );
	}
	if( gd->time > 0 )
	{
		printf( "\nCompute spatial distribution of predictions at t = %g ...", gd->time );
		fflush( stdout );
		sprintf( filename, "%s.vtk", f );
		compute_grid( filename, op->cd, gd );
		printf( "done.\n" );
	}
}

char *timestamp()
{
	time_t raw_time;
	struct tm *ptr_ts;
	char *datetime;
	datetime = malloc( 10 * sizeof( char ) );
	time( &raw_time );
	ptr_ts = localtime( &raw_time );
//	printf( "%s\n", asctime( ptr_ts ) );
	sprintf( datetime, "%02d:%02d:%02d", ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
	return( datetime );
}

char *datestamp()
{
	time_t raw_time;
	struct tm *ptr_ts;
	char *datetime;
	datetime = malloc( 16 * sizeof( char ) );
	time( &raw_time );
	ptr_ts = localtime( &raw_time );
//	printf( "%s\n", asctime( ptr_ts ) );
	sprintf( datetime, "%4d%02d%02d-%02d%02d%02d", ptr_ts->tm_year + 1900, ptr_ts->tm_mon + 1, ptr_ts->tm_mday, ptr_ts->tm_hour, ptr_ts->tm_min, ptr_ts->tm_sec );
	return( datetime );
}
