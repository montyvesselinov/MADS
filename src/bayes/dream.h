// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
// Brianeisha Eure
// Leif Zinn-Bjorkman
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

#ifndef DREAM_H_
#define DREAM_H_

enum PPCR_TYPE { PPCR_UPDATE };
enum BOUND_HANDLING_TYPE { BH_REFLECT, BH_BOUND, BH_FOLD };
enum INIT_POPULATION_TYPE { IP_LHS_BASED };

#ifdef __cplusplus
extern "C"
#endif
struct MCMC *get_posterior_parameter_samples( struct opt_data *od );
#ifdef __cplusplus
extern "C"
#endif
void free_mcmc( struct MCMC *mcmc );
#ifdef __cplusplus
extern "C"
#endif
void tprintf( char const *fmt, ... );
#ifdef __cplusplus
extern "C"
#endif
void symmetric_astable_pdf_interp( double x, double alpha, double gamma, double lambda, double *val );

struct MCMC  //introduce a data structure called MCMC
{
	int n;                         // Dimension of the problem (number of parameters to be estimated)
	int ndraw;                     // Maximum number of function evaluations
	double parallelUpdate;         // Fraction of parallel direction updates
	int seq;                       // Number of Markov Chains / sequences
	int seq_length;                // Sequence length
	int xx;                        // Z array length
	int np3;                       // Dimension of the problem (number of parameters to be estimated) + 3
	int DEpairs;                   // Number of chain pairs to generate candidate points
	double gamma;                  // Kurtosis parameter Bayesian Inference Scheme
	int nCR;                       // Number of crossover values used
	int m0;                        // Initial size of Z
	int m;                         // How many elements are actually stored in z
	int nzee;                      // Length of the array z
	int nzoff;
	int k;                         // Thinning parameter for appending X to Z
	double eps;                    // Perturbation for ergodicity
	int steps;                     // Number of steps before calculating convergence diagnostics
	double **CR;
	int ppCR;
	int reduced_sample_collection;
	int BoundHandling;           // boundary handling
	int save_in_memory;          // save in memory or not
	int save_in_file;            // save in file or not
	int InitPopulation;
	int verbose;
	int reduced_seq_interval;
	int reduced_seq_length;
	int qcov;
	//z is double[nzee][np2] containing the sampled parameters and in the last two columns their posterior likelihood and log-likelihood
	//the number of rows that actually contain samples is m, so use only up to z[m]
	double **z;
	//	double qcorr;
};

struct Range
{
	double *minn;                 // minimum parameter ranges
	double *maxn;                 // maximum parameter ranges
};

struct Measure
{
	double *MeasData;           // Define the Measured Streamflow
	double Sigma;
	int N;
};

struct out
{
	int *nEval;
	double *Acceptance_Rate;
	double **CR;
	int *R_stat_iter;
	double **R_stat;
};
#endif /* DREAM_H_ */
