MADS is an open-source code designed as an integrated high-performance computational framework performing a wide range of model-based analyses:
* Sensitivity Analysis
* Parameter Estimation
* Model Inversion and Calibration
* Uncertainty Quantification
* Model Selection and Averaging
* Decision Support

MADS utilizes adaptive rules and techniques which allows the analyses to be performed with minimum user input.
The code provides a series of alternative algorithms to perform each type of model analyses.
The code allows for coupled model parameters and regularization terms that are internally computed based on user-defined mathematical expressions.

For more information check out he MADS website at http://mads.lanl.gov

Licensing: GPLv3: http://www.gnu.org/licenses/gpl-3.0.html

Download:  
* http://gitlab.com/monty/mads  
* http://github.com/montyvesselinov/MADS  
* http://bitbucket.org/monty_vesselinov/mads  

Clone GIT repositories:  
* git clone git@gitlab.com:monty/mads.git  
* git clone git@github.com:montyvesselinov/MADS.git  
* git clone git@bitbucket.org:monty_vesselinov/mads.git  

Required third party libraries (TPL's):  
* GSL: http://www.gnu.org/s/gsl  
* LAPACK: http://www.netlib.org/lapack  
* BLAS: http://www.netlib.org/blas  
* MATHEVAL: http://www.gnu.org/software/libmatheval  
* YAML: http://pyyaml.org/wiki/LibYAML; http://www.yaml.org  
* GLIB: http://developer.gnome.org/glib  
  
The required TPL's can be downloaded at http://gitlab.com/mads/mads-tpls  
  
The code has been tested on Apple MAC OS X, Linux (RHEL / CentOS / Fedora / Ubuntu / Debian) and Cygwin/Microsoft Windows.  
  
Checkout the 'readme' file for installation instructions.  
  
Compilation:  
* `make` OR
* `cmake -f CMakeLists.txt` (cmake version 3.1 is required)
  
Verification:  
* `make verify` (all test problems listed below)  
* `make verify-internal` (internal test functions)  
* `make verify-forward` (forward contaminant transport simulations)  
* `make verify-contaminant` (inverse contaminant transport analyses)  
* `make verify-multistart1` verify-multistart2 (multi-start inverse problems using random initial parameter guesses)  
* `make verify-external` verify-external-short (inverse problems using the code WELLS; http://wells.lanl.gov)  
* `make verify-parallel` (parallel inverse analysis)  
* `make verify-sa` (global sensitivity analysis)  
  
Examples:  
* `make examples`  
* additional examples can be found in directory 'example' (check the 'readme' files in the directory 'example')  
* see also http://mads.lanl.gov/#examples and http://mads.lanl.gov/#demos  
  
Comparisons:
* test probelems are available
* comparisons with the code PEST are available as well (in subdirectories of 'example', check the 'readme' files)  
* see also http://mads.lanl.gov/#comparisons  
  
Manual: http://mads.lanl.gov/#manual  
  
USAGE:
*      `mads problem_name    [ keywords | options ]`  OR
*      `mads MADS_input_file [ keywords | options ]`  OR
*      `mads PEST_input_file [ keywords | options ]`  (MADS is compatible with PEST control, template and instruction files)

where:
* `problem_name`:         name of the solved problem; `MADS_input_file` named `problem_name.mads` is expected
* `MADS_input_file`:      problem input file in MADS format (typically `*.mads`)
* `PEST_input_file`:      problem input file in PEST format (PEST control file; `*.pst`)

keywords & options (can be provided in any order):

problem type keywords:
*  check              - check model setup and input/output files (no model execution)
*  create             - create calibration problem based on provided model parameters (single model execution to compute calibration targets)
*  forward            - forward model run (no calibration); recommended for model setup checking
*  calibrate          - calibration run [default if no problem type keyword is provided]
*  montecarlo         - Monte-Carlo analysis
*  gsens              - global sensitivity analysis
*  glue               - Generalized Likelihood Uncertainty Estimation
*  lsens              - local sensitivity analysis (standalone or at the end of the calibration)
*  eigen              - local eigensystem analysis (standalone or at the end of the calibration)
*  abagus             - Agent-Based Global Uncertainty & Sensitivity Analysis
*  infogap            - Info-gap decision analysis
*  postpua            - predictive uncertainty analysis of sampling results (currently for abagus PSSA files only)
*  bayes              - Bayesian parameter sampling using DREAM (DiffeRential Evolution Adaptive Metropolis; Vrugt et al. 2009)

calibration method keywords (select one):
*  single             - single calibration using initial guesses provided in the input file [default]
*  igrnd              - sequential calibrations using a set of random initial values (number of realizations defined by real=X)
*  igpd               - sequential calibrations using a set of discretized initial values (discretization defined in the input file)
*  ppsd               - sequential calibrations using partial parameter space discretization (PPSD) method (discretization defined in the input file)

calibration termination criteria:
*  eval=[integer]     - functional evaluations exceed the predefined value [default eval=5000]
*  cutoff=[real]      - objective function is below the cutoff value [default cutoff=0]
*  obsrange           - model predictions are within predefined calibration ranges
*  obserror=[real]    - model predictions are within a predefined absolute error [default obserror=0.1]
*  parerror=[real]    - model parameters are within a predefined absolute error from their known 'true' values [default parerror=0.1]

user-enforced termination:
*  problem_name.quit  - if file with this name exists in the running directory, MADS terminates as soon as possible.
*  problem_name.stop  - if file with this name exists in the running directory, MADS terminates after saving intermediate results.

optimization method (opt=[string]; various combinations are possible, e.g. pso_std_lm_gsl):
*  opt=lm             - Local Levenberg-Marquardt optimization [default]
*  opt=lm_levmar      - Local Levenberg-Marquardt optimization using LEVMAR library
*  opt=lm_gsl         - Local Levenberg-Marquardt optimization using GSL library
*  opt=lm_ms          - Local Multi-Start (Multi-Try) Levenberg-Marquardt (MSLM) optimization using multiple random initial guesses
*  opt=pso            - Global Particle Swarm optimization (default Standard2006; http://clerc.maurice.free.fr/pso)
*  opt=apso           - Global Adaptive Particle Swarm optimization (default TRIBES)
*  opt=swarm          - Global Particle Swarm optimization Standard2006 (also opt=pso_std)
*  opt=tribes         - Global Particle Swarm optimization TRIBES-D (Clerc 2004; http://clerc.maurice.free.fr/pso)
*  opt=squads         - SQUADS: Adaptive hybrid optimization using coupled local and global optimization techniques

general calibration/optimization options:
*  retry=[integer]    - number of optimization retries [default retry=0]
*  particles=[integer]- number of initial particles or tribes [default `particles=10+2*sqrt(Number_of_parameters)`]
*  lmeigen|eigen      - eigen analysis of the intermediate / final optimized solution

Levenberg-Marquardt optimization options:
*  lmfactor=[double]  - multiplier applied to compute when to initiate LM searches within SQUADS algorithm [default lmfactor=1.0]
*  lmindirect         - Indirect computation of LM alpha coefficient [default DIRECT/Delayed gratification computation]
*  lmmu=[double]      - LM alpha multiplier for direct computation of LM alpha coefficient when OF decreases [default lmmu=0.1]
*  lmnu=[integer]     - LM alpha multiplier for direct computation of LM alpha coefficient when OF increases [default lmnu=10]
*  lmaccel            - LM geodesic acceleration as proposed by Transtrum et al (2011) [default NO acceleration]
*  lmratio=[double]   - LM acceleration velocity ratio for recomputing the Jacobian [default lmratio=(3/4)^2]
*  lmh=[double]       - LM acceleration multiplier [default lmh=0.1]
*  lmiter=[integer]   - number of LM iterations [default computed internally based on number of evaluations]
*  lmnlam=[integer]   - Maximum number of linear solves (lambda searches) after each Jacobian estimate [default lmnlam=10]
*  lmnlamof=[integer] - Number of acceptable linear solves (lambda searches) with similar OF's during LM optimization [default lmnlamof=3]
*  lmnjacof=[integer] - Number of acceptable jacobian iterations with similar OF's during LM optimization [default lmnjacof=5]
*  lmerror=[double]   - LM convergence error [default lmerror=1e-5]

sampling method (smp=[string] OR mslm=[string] for Multi-Start Levenberg-Marquardt (MSLM) analysis using multiple retries):
*  smp=olhs           - Optimal Latin Hyper Cube sampling [default] (if real=1 RANDOM; if real>1 IDLHS; if real>500 LHS)
*  smp=lhs            - Latin Hyper Cube sampling (LHS)
*  smp=idlhs          - Improved Distributed Latin Hyper Cube sampling (IDLHS)
*  smp=random         - Random sampling

sampling options:
*  real=[integer]     - number of random realizations / samples [default real=100]
*  case=[integer]     - execute a single case from all the realizations / samples (applied in PPSD, IGDP, IGRND, MONTECARLO)
*  seed=[integer]     - random seed value [randomly generated by default]

objective function functional form options (select one; regularizaiton terms can be added separately):
*  ssr                - sum of the squared residuals [default]
*  ssd0               - sum of the squared discrepancies
*  ssdx               - sum of the squared discrepancies increased to get in the bounds
*  ssda               - sum of the squared discrepancies and absolute residuals
*  ssdr               - sum of the squared discrepancies and squared residuals

transformation of parameter and observation properties:
*  nosin              - Sin transformation of optimized parameters is not applied [parameters are sin transformed by default]
*  sindx              - Parameter space step for numerical derivatives of sin transformed parameters
*                       [default sindx=1e-7 for internal problems; sindx=0.1 for external problems]
*  lindx              - Parameter space step for numerical derivatives of not transformed parameters [default lindx=0.001]
*  pardx              - Parameter space step for parameter space discretization [default pardx=0.1]
*  plog=[-1,0,1]      - Log transformation of all optimized parameters is enforced (1) or disabled (0)
*                       [default plog=-1; log transformation is explicitly defined for each parameter in the input file]
*  olog=[-1,0,1]      - Log transformation of all the observations (simulated and measured) is enforced (1) or disabled (0)
*                       [default olog=-1; log transformation is explicitly defined for each observation in the input file]
*  oweight=[-1,0,1,2] - Weights for all the observation residuals are defined:
*                       0 = zero weight, 1 = unit weight, 2 = weight reversely proportional to observation
*                       [default oweight=-1; weights for each observation are explicitly defined in the input file]
*  obsdomain=[float]  - observation space domain size [default provided in the MADS input file]
*  obsstep=[float]    - observation space domain step to explore info-gap observation uncertainty [default ignored]

parallelization (available parallelization resources are internally detected by default; supported: Slurm, Moab, BProc, LSB, PBS, Beowulf, OpenMPI, OpenMP, POSIX, ...):
*  np=[integer]       - Number of requested parallel jobs [optional]
*  mpi                - Use OpenMPI parallelization [optional]
*  posix[=integer]    - POSIX parallel threading [optional]; number of threads can be defined as well
*  omp[=integer]      - OpenMP parallel threading [optional]; number of threads can be defined as well
*  ppt=[integer]      - Number of processors per an external model task [optional]
*  nplambda=[integer] - Number of requested parallel lambda runs in the case of Levenberg-Marquardt optimization [optional; nplambda <= np]
*  rstfile=[string]   - name of existing ZIP restart file to be used (created by previous Parallel MADS run) [optional]
*  rstdir=[string]    - name of existing restart directory to be used (created by previous Parallel MADS run) [optional]
*  restart=[integer]  - restart=1 (default; automatic restart if possible); restart=0 (force no restart); restart=-1 (force restart) by default the analyses will be restarted automatically (restart=1)
*  bin_restart        - restart information is stored in binary files (by default, model outputs are added in a zip file for restart)

ABAGUS (Agent-Based Global Uncertainty & Sensitivity Analysis) options:
*  infile=[string]    - name of previous results file to be used to initialize Kd-tree [default=NULL]
*  energy=[integer]   - initial energy for particles [default energy=10000]

options for the build-in analytical solutions:
*  point              - point contaminant source in 3D flow domain
*  plane              - plane contaminant source in 3D flow domain
*  box                - brick contaminant source in 3D flow domain
*  obs_int=[1,2]      - concentration integration along observation well screens (1 - mid point; 2 - two end points [default=1]
*  disp_tied          - lateral and vertical transverse dispersivities are fractions of the longitudinal dispersivity
*  disp_scaled        - longitudinal dispersivity is scaled with the travel distance
*                       if disp_scaled and disp_tied are applied, longitudinal dispersivity is scaled and
*                       lateral and vertical transverse dispersivities are fractions of the scaled longitudinal dispersivity
*  disp_scaled=2      - longitudinal, lateral and vertical transverse dispersivities are scaled with the travel distance
*  time_step          - parameter "end time" is representing the period (dt) within which the source is active

build-in test problems for optimization / uncertainty-quantification techniques (local and global methods):
*  test=[integer]     - test problem ID [default=1]:
	*                          1: Parabola (Sphere)
	*                          2: Griewank
	*                          3: Rosenbrock
	*                          4: De Jong's Function #4
	*                          5: Step
	*                          6: Alpine function (Clerc's Function #1)
	*                          7: Rastrigin
	*                          8: Krishna Kumar
	*                          9: Tripod function 2D
	*                         10: Shekel's Foxholes 2D
	*                         11: Shekel's Foxholes 5D
	*                         12: Shekel's Foxholes 10D
	*                         20: Shekel's Foxholes 2D (alternative; global methods only)
	*                         21: Polynomial fitting (global methods only)
	*                         22: Ackley (global methods only)
	*                         23: Eason 2D (global methods only)
	*                         31: Rosenbrock (2D simplified alternative)
	*                         32: Griewank (alternative)
	*                         33: Rosenbrock (alternative with d*(d-1) observations
	*                         34: Powell's Quadratic
	*                         35: Booth
	*                         36: Beale
	*                         37: Parsopoulos
	*                         Curve-fitting test functions:
	*                         40: Sin/Cos test function (2 parameters)
	*                         41: Sin/Cos test function (4 parameters)
	*                         42: Sin/Cos test function (2 parameters; simplified)
	*                         43: Exponential Data Fitting I (5 parameters)
	*                         44: Exponential Data Fitting II (11 parameters)
*  dim=[integer]      - dimensionality of parameter space for the test problem (fixed for some of the problems) [default=2]
*  npar=[integer]     - number of model parameters for the test problem (fixed for some of the problems) [default=2]
*  nobs=[integer]     - number of observations for the test problem (fixed for some of the problems) [default=2]
*  pardomain=[float]  - parameter space domain size [default pardomain=100]

debugging / verbose levels:
*  force | f          - enforce running even if a file *.running exists (the file *.running is to prevent overlapping executions)
*  quiet | q          - no screen output (all the screen output is saved in a file with extension mads_output
*  debug=[0-5]        - general debugging [default debug=0]
*  fdebug=[0-5]       - model evaluation debugging [default fdebug=0]
*  ldebug=[0-3]       - Levenberg-Marquardt optimization debugging [default ldebug=0]
*  pdebug=[0-3]       - Particle Swarm optimization debugging [default pdebug=0]
*  mdebug=[0-3]       - Random sampling debugging [default mdebug=0]
*  odebug=[0-1]       - Record objective function progress in a file with extension 'ofe' [default odebug=0]
*  tdebug=[0-1]       - Output process times for various tasks [default tdebug=0]
*  tpldebug=[0-3]     - Debug the writing of external files [default tpldebug=0]
*  insdebug=[0-3]     - Debug the reading of external files [default insdebug=0]
*  pardebug=[0-3]     - Debug the parallel execution [default pardebug=0]

pre-/post-processing:
*  yaml               - Reads/writes/converts MADS files in YAML format
*  xml                - Reads/writes/converts MADS files in XML format
*  resultsfile=[file] - Post process results saved in a previously generated resultsfile (e.g. using 'igrnd' analysis)
*  resultscase=[int]  - Post process specific case saved in resultsfile (if resultscase<0, first abs(resultscase) cases)
*  cutoff=[real]      - Post process all cases saved in resultsfile with objective function below the cutoff value
*  obsrange           - Post process all cases saved in resultsfile with model predictions within predefined calibration ranges
*  pargen             - Generate MADS input files for independent multi-processor runs [default pargen=0]
*  save               - Save MADS input/output files for successful parameter sets [default save=0]

Examples:
*  mads a01 test=2 opt=lm eigen igrnd real=1 (no input files are needed for execution)
*  mads a01 test=2 opt=squads igrnd real=1
*  mads a01 test=2 abagus cutoff=0.1 eval=100000 (collect solutions of Griewank function below phi cutoff)
*  mads a01 test=3 abagus cutoff=20 eval=100000  (collect solutions of Rosenbrock function below phi cutoff)
*  mads a01 test=111 gsens dim=8 real=10000 smp=lhs pardomain=0.5 seed=15176048200 (Sobol's global sensitivity analysis)
*  mads examples/contamination/s01 ldebug lmeigen (file s01.mads is located in examples/contamination)
*  mads examples/contamination/s01 ldebug lmeigen igrnd real=1
*  mads examples/contamination/s01 seed=1549170842 obsrange igrnd real=1
*  mads examples/contamination/s01 opt=squads seed=1549170842 eigen obsrange pdebug igrnd real=1
*  mads examples/contamination/s01 opt=pso seed=1549170842 eigen obsrange igrnd real=1
*  mads examples/contamination/s01-flagged ppsd (Partial Parameter Space Discretization)
*  mads examples/contamination/s01-flagged igpd (Initial Guesses based on Discreetly distributed model parameters)
*  mads examples/contamination/s01 igrnd real=10 (Random Initial Guesses; all parameters)
*  mads examples/contamination/s01-flagged igrnd real=10 (Random Initial Guesses; only flagged parameters)
*  mads examples/contamination/s01 monte real=10 (Monte Carlo Analysis)
*  cd examples/wells; mads w01.mads igrnd real=1 seed=501228648 eigen

Parallel Levenberg-Marquardt optimization:
*  cd examples/wells-short
*  mads w01parallel.mads restart=0 np=2 ldebug pardebug=2 (Parallel optimization using 2 processors)
*  mads w01parallel.mads restart=0 np=11 nplambda=11
*  mads w01parallel.mads restart=0 np=3 nplambda=3 lmnlam=21 lmnlamof=12 (if a small number of processors is used, lmnlam & lmnlamof should be increased)

Comparisons between local and global methods:
*  mads a01 test=3 opt=lm     igrnd real=1000 cutoff=1e-3 (Levenberg-Marquardt optimization)
*  mads a01 test=3 opt=lm_ms  igrnd real=1000 cutoff=1e-3 (Multi-Start Levenberg-Marquardt optimization)
*  mads a01 test=3 opt=swarm  igrnd real=1000 cutoff=1e-3 (Particle Swarm optimization Standard2006)
*  mads a01 test=3 opt=tribes igrnd real=1000 cutoff=1e-3 (Particle Swarm optimization TRIBES-D)
*  mads a01 test=3 opt=squads igrnd real=1000 cutoff=1e-3 (Adaptive hybrid optimization Squads)

Comparisons with PEST (http://www.sspa.com/pest/):
*  cd examples/contamination
*  mads s02 lmeigen                  (file s02.mads is located in examples/contamination)
*  pest s02pest                   (file s02pest.pst is located in examples/contamination)
*  mads w01 lmeigen     (files associated with problem w01 are located in examples/wells)
*  pest w01pest     (files associated with problem w01pest are located in examples/wells)
*  pest w02pest     (files associated with problem w02pest are located in examples/wells)

For additional information:
*  web:   http://mads.lanl.gov -:- http://www.ees.lanl.gov/staff/monty/codes/mads
*  repo:  http://gitlab.com/monty/mads
*  git:   git clone git@gitlab.com:monty/mads.git
*  email: Velimir V Vesselinov (monty) vvv@lanl.gov -:- velimir.vesselinov@gmail.com
