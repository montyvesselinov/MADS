/*
 * dream.h
 *
 *  Created on: Mar 18, 2012
 *      Author: monty
 */

#ifndef DREAM_H_
#define DREAM_H_

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
	string ppCR;
	string reduced_sample_collection;
	string BoundHandling;           // boundary handling
	string save_in_memory;          // save in memory or not
	string save_in_file;            // save in file or not
	string InitPopulation;
	int verbose;
	int reduced_seq_interval;
	int reduced_seq_length;
	int qcov;
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
