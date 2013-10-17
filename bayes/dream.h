/*
 * dream.h
 *
 *  Created on: Mar 18, 2012
 *      Author: monty
 */

#ifndef DREAM_H_
#define DREAM_H_

enum PPCR_TYPE { PPCR_UPDATE };
enum BOUND_HANDLING_TYPE { BH_REFLECT, BH_BOUND, BH_FOLD };
enum INIT_POPULATION_TYPE { IP_LHS_BASED };

struct MCMC  //introduce a data structure called MCMC
{
	int n;                         // Dimension of the problem (number of parameters to be estimated)
	int ndraw;                     // Maximum number of function evaluations
	double parallelUpdate;         // Fraction of parallel direction updates
	int seq;                       // Number of Markov Chains / sequences
	int seq_length;                // Sequence length
	int xx;                        // Z array length
	int np2;                       // Dimension of the problem (number of parameters to be estimated) + 2
	int DEpairs;                   // Number of chain pairs to generate candidate points
	double gamma;                  // Kurtosis parameter Bayesian Inference Scheme
	int nCR;                       // Number of crossover values used
	int m0;                        // Initial size of Z
	int m;
	int nzee;
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
#ifdef __cplusplus
extern "C"
#endif
struct MCMC *get_posterior_parameter_samples( struct opt_data *od );
#ifdef __cplusplus
extern "C"
#endif
void free_mcmc( struct MCMC *mcmc );
#endif /* DREAM_H_ */
