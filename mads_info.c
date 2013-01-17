/*
 * mads_info.c
 *
 *  Created on: Mar 4, 2012
 *      Author: monty
 */

#include <stdio.h>

void mads_info()
{
	printf( "MADS is an open-source code designed as an integrated high-performance computational framework performing\n" );
	printf( "a wide range of model-based analyses: Sensitivity Analysis (SA), Parameter Estimation (PE), Model Inversion and Calibration,\n" );
	printf( "Uncertainty Quantification (UQ), Model Selection and Averaging, and Decision Support (DS).\n\n" );
	printf( "MADS utilizes adaptive rules and techniques which allows the analyses to be performed with minimum user input. The code\n" );
	printf( "provides a series of alternative algorithms to perform each type of model analyses. The code allows for coupled model\n" );
	printf( "parameters and regularization terms that are internally computed based on user-defined mathematical expressions.\n\n" );
	printf( "USAGE: mads problem_name [ keywords | options ]       OR\n" );
	printf( "       mads MADS_input_file [ keywords | options ]    OR\n" );
	printf( "       mads PEST_input_file [ keywords | options ] (MADS is compatible with PEST control, template and instruction files)\n\n" );
	printf( "problem_name:         name of the solved problem; MADS_input_file named problem_name.mads is expected\n" );
	printf( "MADS_input_file:      problem input file in MADS format (typically *.mads)\n" );
	printf( "PEST_input_file:      problem input file in PEST format (PEST control file; *.pst)\n\n" );
	printf( "keywords & options (can be provided in any order):\n\n" );
	printf( "problem type keywords:\n" );
	printf( "   check              - check model setup and input/output files (no model execution)\n" );
	printf( "   create             - create calibration problem based on provided model parameters (single model execution to compute calibration targets)\n" );
	printf( "   forward            - forward model run (no calibration); recommended for model setup checking\n" );
	printf( "   calibrate          - calibration run [default if no problem type keyword is provided]\n" );
	printf( "   montecarlo         - Monte-Carlo analysis\n" );
	printf( "   gsens              - global sensitivity analysis\n" );
	printf( "   glue               - Generalized Likelihood Uncertainty Estimation\n" );
	printf( "   lsens              - local sensitivity analysis (standalone or at the end of the calibration)\n" );
	printf( "   eigen              - local eigensystem analysis (standalone or at the end of the calibration)\n" );
	printf( "   abagus             - Agent-Based Global Uncertainty & Sensitivity Analysis\n" );
	printf( "   infogap            - Info-gap decision analysis\n" );
	printf( "   postpua            - predictive uncertainty analysis of sampling results (currently for abagus PSSA files only)\n" );
	printf( "\ncalibration method keywords (select one):\n" );
	printf( "   single             - single calibration using initial guesses provided in the input file [default]\n" );
	printf( "   igrnd              - sequential calibrations using a set of random initial values (number of realizations defined by real=X)\n" );
	printf( "   igpd               - sequential calibrations using a set of discretized initial values (discretization defined in the input file)\n" );
	printf( "   ppsd               - sequential calibrations using partial parameter space discretization (PPSD) method\n" );
	printf( "\ncalibration termination criteria:\n" );
	printf( "   eval=[integer]     - functional evaluations exceed the predefined value [default eval=5000]\n" );
	printf( "   cutoff=[real]      - objective function is below the cutoff value [default cutoff=0]\n" );
	printf( "   obsrange           - model predictions are within predefined calibration ranges\n" );
	printf( "   obserror=[real]    - model predictions are within a predefined absolute error [default obserror=0.1]\n" );
	printf( "   parerror=[real]    - model parameters are within a predefined absolute error from their known 'true' values [default parerror=0.1]\n" );
	printf( "\noptimization method (opt=[string]; various combinations are possible, e.g. pso_std_lm_gsl):\n" );
	printf( "   opt=lm             - Local Levenberg-Marquardt optimization [default]\n" );
	printf( "   opt=lm_levmar      - Local Levenberg-Marquardt optimization using LEVMAR library\n" );
	printf( "   opt=lm_gsl         - Local Levenberg-Marquardt optimization using GSL library\n" );
	printf( "   opt=lm_ms          - Local Multi-Start (Multi-Try) Levenberg-Marquardt (MSLM) optimization using multiple random initial guesses\n" );
	printf( "   opt=pso            - Global Particle Swarm optimization (default Standard2006; http://clerc.maurice.free.fr/pso)\n" );
	printf( "   opt=apso           - Global Adaptive Particle Swarm optimization (default TRIBES)\n" );
	printf( "   opt=swarm          - Global Particle Swarm optimization Standard2006 (also opt=pso_std)\n" );
	printf( "   opt=tribes         - Global Particle Swarm optimization TRIBES-D (Clerc 2004; http://clerc.maurice.free.fr/pso)\n" );
	printf( "   opt=squads         - SQUADS: Adaptive hybrid optimization using coupled local and global optimization techniques\n" );
	printf( "\ngeneral calibration/optimization options:\n" );
	printf( "   retry=[integer]    - number of optimization retries [default retry=0]\n" );
	printf( "   particles=[integer]- number of initial particles or tribes [default particles=10+2*sqrt(Number_of_parameters)]\n" );
	printf( "   lmeigen|eigen      - eigen analysis of the intermediate / final optimized solution\n" );
	printf( "\nLevenberg-Marquardt optimization options:\n" );
	printf( "   lmfactor=[double]  - multiplier applied to compute when to initiate LM searches within SQUADS algorithm [default lmfactor=1]\n" );
	printf( "   lmindirect         - Indirect computation of LM alpha coefficient [default DIRECT/Delayed gratification computation]\n" );
	printf( "   lmmu=[double]      - LM alpha multiplier for direct computation of LM alpha coefficient when OF decreases [default lmmu=0.1]\n" );
	printf( "   lmnu=[integer]     - LM alpha multiplier for direct computation of LM alpha coefficient when OF increases [default lmnu=10]\n" );
	printf( "   lmaccel            - LM geodesic acceleration as proposed by Transtrum et al (2011) [default NO acceleration]\n" );
	printf( "   lmratio=[double]   - LM acceleration velocity ratio for recomputing the Jacobian [default lmratio=(3/4)^2]\n" );
	printf( "   lmh=[double]       - LM acceleration multiplier [default lmh=0.1]\n" );
	printf( "   lmiter=[integer]   - number of LM iterations [default computed internally based on number of evaluations]\n" );
	printf( "   lmerror=[double]   - LM convergence error [default lmerror=1e-5]\n" );
	printf( "   lmnlamof=[integer] - Number of acceptable linear solves (lambda searches) with similar OF's during LM optimization [default lmnlamof=3]\n" );
	printf( "   lmnjacof=[integer] - Number of acceptable jacobian iterations with similar OF's during LM optimization [default lmnjacof=4]\n" );
	printf( "\nsampling method (smp=[string] OR mslm=[string] for Multi-Start Levenberg-Marquardt (MSLM) analysis using multiple retries):\n" );
	printf( "   smp=olhs           - Optimal Latin Hyper Cube sampling [default] (if real = 1 RANDOM; if real > IDLHS; if real > 500 LHS)\n" );
	printf( "   smp=lhs            - Latin Hyper Cube sampling (LHS)\n" );
	printf( "   smp=idlhs          - Improved Distributed Latin Hyper Cube sampling (IDLHS)\n" );
	printf( "   smp=random         - Random sampling\n" );
	printf( "\nsampling options:\n" );
	printf( "   real=[integer]     - number of random realizations / samples [default real=100]\n" );
	printf( "   case=[integer]     - execute a single case from all the realizations / samples (applied in PPSD, IGDP, IGRND, MONTECARLO)\n" );
	printf( "   seed=[integer]     - random seed value [randomly generated by default]\n" );
	printf( "\nobjective function functional form options (select one; regularizaiton terms can be added separately):\n" );
	printf( "   ssr                - sum of the squared residuals [default]\n" );
	printf( "   ssd0               - sum of the squared discrepancies\n" );
	printf( "   ssdx               - sum of the squared discrepancies increased to get in the bounds\n" );
	printf( "   ssda               - sum of the squared discrepancies and absolute residuals\n" );
	printf( "   ssdr               - sum of the squared discrepancies and squared residuals\n" );
	printf( "\ntransformation of parameter and observation properties:\n" );
	printf( "   nosin              - Sin transformation of optimized parameters is not applied [parameters are sin transformed by default]\n" );
	printf( "   sindx              - Parameter space step for numerical derivatives of sin transformed parameters\n" );
	printf( "                        [default sindx=1e-7 for internal problems; sindx=0.1 for external problems]\n" );
	printf( "   lindx              - Parameter space step for numerical derivatives of not transformed parameters [default lindx=0.001]\n" );
	printf( "   pardx              - Parameter space step for parameter space discretization [default pardx=0.1]\n" );
	printf( "   plog=[-1,0,1]      - Log transformation of all optimized parameters is enforced (1) or disabled (0)\n" );
	printf( "                        [default plog=-1; log transformation is explicitly defined for each parameter in the input file]\n" );
	printf( "   olog=[-1,0,1]      - Log transformation of all the observations (simulated and measured) is enforced (1) or disabled (0)\n" );
	printf( "                        [default olog=-1; log transformation is explicitly defined for each observation in the input file]\n" );
	printf( "   oweight=[-1,0,1,2] - Weights for all the observation residuals are defined:\n" );
	printf( "                        0 = zero weight, 1 = unit weight, 2 = weight reversely proportional to observation\n" );
	printf( "                        [default oweight=-1; weights for each observation are explicitly defined in the input file]\n" );
	printf( "   obsdomain=[float]  - observation space domain size [default provided in the MADS input file]\n" );
	printf( "   obsstep=[float]    - observation space domain step to explore info-gap observation uncertainty [default ignored]\n" );
	printf( "\nparallelization (parallelization environment and available resources are internally detected by default):\n" );
	printf( "   np=[integer]       - Number of requested parallel jobs [optional]\n" );
	printf( "   rstfile=[string]   - name of existing ZIP restart file to be used (created by previous Parallel MADS run) [optional]\n" );
	printf( "   restart=[integer]  - restart=1 (default; automatic restart if possible); restart=0 (force no restart); restart=2 (force restart)\n" );
	printf( "                        by default the analyses will be restarted automatically (restart=1)\n" );
	printf( "\nABAGUS (Agent-Based Global Uncertainty & Sensitivity Analysis) options:\n" );
	printf( "   infile=[string]    - name of previous results file to be used to initialize Kd-tree [default=NULL]\n" );
	printf( "   energy=[integer]   - initial energy for particles [default energy=10000]\n" );
	printf( "\noptions for the build-in analytical solutions:\n" );
	printf( "   point              - point contaminant source in 3D flow domain\n" );
	printf( "   plane              - areal contaminant source in 3D flow domain\n" );
	printf( "   box                - brick contaminant source in 3D flow domain\n" );
	printf( "   disp_tied          - lateral and vertical transverse dispersivities are fractions of the longitudinal dispersivity\n" );
	printf( "   disp_scaled        - longitudinal dispersivity is scaled with the travel distance\n" );
	printf( "                        if disp_scaled and disp_tied are applied, longitudinal dispersivity is scaled and\n" );
	printf( "                        lateral and vertical transverse dispersivities are fractions of the scaled longitudinal dispersivity\n" );
	printf( "   disp_scaled=2      - longitudinal, lateral and vertical transverse dispersivities are scaled with the travel distance\n" );
	printf( "\nbuild-in test problems for optimization / uncertainty-quantification techniques (local and global methods):\n" );
	printf( "   test=[integer]     - test problem ID [default=1]:\n" );
	printf( "                           1: Parabola (Sphere)\n" );
	printf( "                           2: Griewank\n" );
	printf( "                           3: Rosenbrock\n" );
	printf( "                           4: De Jong's Function #4\n" );
	printf( "                           5: Step\n" );
	printf( "                           6: Alpine function (Clerc's Function #1)\n" );
	printf( "                           7: Rastrigin\n" );
	printf( "                           8: Krishna Kumar\n" );
	printf( "                           9: Tripod function 2D\n" );
	printf( "                          10: Shekel's Foxholes 2D\n" );
	printf( "                          11: Shekel's Foxholes 5D\n" );
	printf( "                          12: Shekel's Foxholes 10D\n" );
	printf( "                          20: Shekel's Foxholes 2D (alternative; global methods only)\n" );
	printf( "                          21: Polynomial fitting (global methods only)\n" );
	printf( "                          22: Ackley (global methods only)\n" );
	printf( "                          23: Eason 2D (global methods only)\n" );
	printf( "                          31: Rosenbrock (2D simplified alternative)\n" );
	printf( "                          32: Griewank (alternative)\n" );
	printf( "                          33: Rosenbrock (alternative with d*(d-1) observations\n" );
	printf( "                          34: Powell's Quadratic\n" );
	printf( "                          35: Booth\n" );
	printf( "                          36: Beale\n" );
	printf( "                          37: Parsopoulos\n" );
	printf( "                          Curve-fitting test functions:\n" );
	printf( "                          40: Sin/Cos test function (2 parameters)\n" );
	printf( "                          41: Sin/Cos test function (4 parameters)\n" );
	printf( "                          42: Sin/Cos test function (2 parameters; simplified)\n" );
	printf( "                          43: Exponential Data Fitting I (5 parameters)\n" );
	printf( "                          44: Exponential Data Fitting II (11 parameters)\n" );
	printf( "   dim=[integer]      - dimensionality of parameter space for the test problem (fixed for some of the problems) [default=2]\n" );
	printf( "   npar=[integer]     - number of model parameters for the test problem (fixed for some of the problems) [default=2]\n" );
	printf( "   nobs=[integer]     - number of observations for the test problem (fixed for some of the problems) [default=2]\n" );
	printf( "   pardomain=[float]  - parameter space domain size [default pardomain=100]\n" );
	printf( "\ndebugging / verbose levels:\n" );
	printf( "   quiet              - no screen output (all the output is saved in a file with extension mads_output\n" );
	printf( "   debug=[0-5]        - general debugging [default debug=0]\n" );
	printf( "   fdebug=[0-5]       - model evaluation debugging [default fdebug=0]\n" );
	printf( "   ldebug=[0-3]       - Levenberg-Marquardt optimization debugging [default ldebug=0]\n" );
	printf( "   pdebug=[0-3]       - Particle Swarm optimization debugging [default pdebug=0]\n" );
	printf( "   mdebug=[0-3]       - Random sampling debugging [default mdebug=0]\n" );
	printf( "   odebug=[0-1]       - Record objective function progress in a file with extension \'ofe\' [default odebug=0]\n" );
	printf( "   tpldebug=[0-3]     - Debug the writing of external files [default tpldebug=0]\n" );
	printf( "   insdebug=[0-3]     - Debug the reading of external files [default insdebug=0]\n" );
	printf( "   pardebug=[0-3]     - Debug the parallel execution [default pardebug=0]\n" );
	printf( "\npre-/post-processing:\n" );
	printf( "   resultsfile        - Post process results saved in resultsfile\n" );
	printf( "   resultscase        - Post process specific case saved in resultsfile (if resultscase<0, first abs(resultscase) cases)\n" );
	printf( "   cutoff=[real]      - Post process all cases saved in resultsfile with objective function below the cutoff value\n" );
	printf( "   obsrange           - Post process all cases saved in resultsfile with model predictions within predefined calibration ranges\n" );
	printf( "   pargen             - Generate MADS input files for independent multi-processor runs [default pargen=0]\n" );
	printf( "   save               - Save MADS input/output files for successful parameter sets [default save=0]\n" );
	printf( "\nExamples:\n" );
	printf( "   mads a01 test=2 opt=lm eigen igrnd real=1 (no input files are needed for execution)\n" );
	printf( "   mads a01 test=2 opt=squads igrnd real=1\n" );
	printf( "   mads a01 test=2 abagus cutoff=0.1 eval=100000 (collect solutions of Griewank function below phi cutoff)\n" );
	printf( "   mads a01 test=3 abagus cutoff=20 eval=100000  (collect solutions of Rosenbrock function below phi cutoff)\n" );
	printf( "   mads example/contamination/s01 ldebug lmeigen (file s01.mads is located in example/contamination)\n" );
	printf( "   mads example/contamination/s01 ldebug lmeigen igrnd real=1\n" );
	printf( "   mads example/contamination/s01 seed=1549170842 obsrange igrnd real=1\n" );
	printf( "   mads example/contamination/s01 opt=squads seed=1549170842 eigen obsrange pdebug igrnd real=1\n" );
	printf( "   mads example/contamination/s01 opt=pso seed=1549170842 eigen obsrange igrnd real=1\n" );
	printf( "   mads example/contamination/s01-flagged ppsd (Partial Parameter Space Discretization)\n" );
	printf( "   mads example/contamination/s01-flagged igpd (Initial Guesses based on Discreetly distributed model parameters)\n" );
	printf( "   mads example/contamination/s01 igrnd real=10 (Random Initial Guesses; all parameters)\n" );
	printf( "   mads example/contamination/s01-flagged igrnd real=10 (Random Initial Guesses; only flagged parameters)\n" );
	printf( "   mads example/contamination/s01 monte real=10 (Monte Carlo Analysis)\n" );
	printf( "   mads w01 np=2 ldebug pardebug=2 (Parallel optimization using 2 processors; files associated with problem w01 are located in example/wells)\n" );
	printf( "             (for w01 example, code WELLS can be obtained at http://wells.lanl.gov)\n" );
	printf( "\nComparisons between local and global methods:\n" );
	printf( "   mads a01 test=3 opt=lm     igrnd real=1000 cutoff=1e-3 (Levenberg-Marquardt optimization)\n" );
	printf( "   mads a01 test=3 opt=lm_ms  igrnd real=1000 cutoff=1e-3 (Multi-Start Levenberg-Marquardt optimization)\n" );
	printf( "   mads a01 test=3 opt=swarm  igrnd real=1000 cutoff=1e-3 (Particle Swarm optimization Standard2006)\n" );
	printf( "   mads a01 test=3 opt=tribes igrnd real=1000 cutoff=1e-3 (Particle Swarm optimization TRIBES-D)\n" );
	printf( "   mads a01 test=3 opt=squads igrnd real=1000 cutoff=1e-3 (Adaptive hybrid optimization Squads)\n" );
	printf( "\nComparisons with PEST (http://www.sspa.com/pest/):\n" );
	printf( "   mads s02 lmeigen                  (file s02.mads is located in example/contamination)\n" );
	printf( "   pest s02pest                   (file s02pest.pst is located in example/contamination)\n" );
	printf( "   mads w01 lmeigen     (files associated with problem w01 are located in example/wells)\n" );
	printf( "   pest w01pest     (files associated with problem w01pest are located in example/wells)\n" );
	printf( "   pest w02pest     (files associated with problem w02pest are located in example/wells)\n" );
	printf( "               (for w01 comparison, code WELLS can be obtained at http://wells.lanl.gov)\n" );
	printf( "\nFor additional information:\n" );
	printf( "   web:   http://mads.lanl.gov -:- http://www.ees.lanl.gov/staff/monty/codes/mads\n" );
	printf( "   email: Velimir Vesselinov (monty) vvv@lanl.gov -:- velimir.vesselinov@gmail.com\n" );
}
