/* -----------------DREAM with sampling from past and snooker updates: DREAM_ZS---------------
 */

// Different test examples from SIAM paper
// example 1: n-dimensional Gaussian distribution
// example 2: multivariate student t distribution
// example 3: n-dimensional banana shaped Gaussian distribution
// example 4: n-dimensional multimodal mixture distribution
// example 5: real-world example using hymod model

#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include <limits>

using namespace std;

#include "dream.h"
#include "../mads.h"

// function headers
struct MCMC *get_posterior_parameter_samples( struct opt_data *od );
void dream_zs( struct MCMC *MCMCPar, Range ParRange, Measure Measurement, struct opt_data *od, int option,
			   double ***Reduced_Seq, double **Zee, out &output, double **Zinit,
			   double *ModPred, double *p, double *log_p,
			   double ***Sequences, double **pCR,
			   double **xold, double *p_xold, double *log_p_xold,
			   int *RandArray, double **xnew, double **CRS,
			   int **DEversion, double *alpha_s, double **Table_JumpRate,
			   int *accept, double *p_xnew, double *log_p_xnew,
			   double **lCR, double *delta_tot, double *R_stat, double **lCRnew );
void LHSU( Range ParRange, struct MCMC *MCMCPar, double **Zinit );
void comp_likelihood( double **x, struct MCMC *MCMCPar, Measure Measurement, struct opt_data *od, double *ModPred, double *p, double *log_p );
void InitVariables( struct MCMC *MCMCPar, int &nEval, int &output_teller, int &reduced_seq_teller,
					double **Table_JumpRate, double ***Reduced_Seq, out &output,
					double ***Sequences, double **pCR, double **CRS, double **lCR );
void CalcCbWb( struct MCMC *MCMCPar, double &Cb, double &Wb );
void GenCR( struct MCMC *MCMCPar, double **lCR, double **pCR, double **CRS );
void multrnd( int nn, int m, int k, double **lCR, double **Y, double **pCR );
void Gelman( struct MCMC *MCMCPar, double *R_stat, double ***Sequences, int start_loc, int end_loc );
void GelmanCum( struct MCMC *MCMCPar, double *R_stat, double **X2 );
void GetLocation( double **X2, struct MCMC *MCMCPar, double **xold, double *p_xold, double *log_p_xold );
void randsample( int n, int k, int *ro, int verbose );
void offde( double **xold, double **Zoff, struct MCMC *MCMCPar, string Update, double **Table_JumpRate, Range ParRange, double **xnew,
			double *alpha_s, double **CRS, int **DEversion );
void DEStrategy( struct MCMC *MCMCPar, int **DEversion );
void ReflectBounds( double **xnew, int nIndivs, int nDim, Range ParRange );
void metrop( double **xnew, double *p_xnew, double *log_p_x, double **xold, double *p_xold, double *log_p_xold,
			 double *alpha_s, Measure Measurement, struct MCMC *MCMCPar, double option, double **newgen, int *accept );
void CalcDelta( struct MCMC *MCMCPar, double *delta_normX, double **CRS, double *delta_tot );
void AdaptpCR( struct MCMC *MCMCPar, double *delta_tot, double **lCR, double **pCR );

double **double_matrix( int maxCols, int maxRows );
int **int_matrix( int maxRows, int maxCols );
void free_matrix( void **matrix, int maxCols );
int cum_seq_size = 0;
double **meanseqsum;
double **meanseqsum2;
double **varseqsum;

//main rundream code
extern "C" struct MCMC *get_posterior_parameter_samples( struct opt_data *od )
{
	struct MCMC *MCMCPar;
	struct Measure Measurement;
	struct Range ParRange;
	struct out output;
	int option = 2;
	int i, k;
	int verbose = 0;
	srand( time( 0 ) ); // initialize seed "randomly"
	// srand( 1 );
	// Problem specific parameter settings
	MCMCPar = (struct MCMC *) malloc( sizeof( struct MCMC ) );
	MCMCPar->n = od->pd->nOptParam;
	MCMCPar->ndraw = od->cd->nreal;                  // Maximum number of function evaluations
	MCMCPar->parallelUpdate = 0.9;           // Fraction of parallel direction updates
	MCMCPar->seq_length = MCMCPar->ndraw;
	MCMCPar->np2 = MCMCPar->n + 2;
	// Recommended parameter settings
	MCMCPar->z = NULL;
	MCMCPar->seq = 5;                        // Number of Markov Chains / sequences
	MCMCPar->DEpairs = 1;                    // Number of chain pairs to generate candidate points
	MCMCPar->CR = double_matrix( MCMCPar->DEpairs, MCMCPar->seq );
	MCMCPar->gamma = 0;                      // Kurtosis parameter Bayesian Inference Scheme
	MCMCPar->nCR = MCMCPar->seq;              // Number of crossover values used
	MCMCPar->m0 = 10 * MCMCPar->n;            // Initial size of Z
	MCMCPar->m = MCMCPar->m0;
	MCMCPar->steps = 10;                     // Number of steps before calculating convergence diagnostics
	MCMCPar->k = 10;                         // Thinning parameter for appending X to Z
	MCMCPar->eps = 0.05;                     // Perturbation for ergodicity
	MCMCPar->nzee = MCMCPar->m0 + ( MCMCPar->seq * ( MCMCPar->ndraw - MCMCPar->m0 ) ) / ( MCMCPar->seq * MCMCPar->k ) + 2 * MCMCPar->seq;
	MCMCPar->nzoff = 2 * MCMCPar->DEpairs * MCMCPar->seq;
	// -----------------------------------------------------------------------------------------------------------------------
	MCMCPar->ppCR = PPCR_UPDATE;                   // Adaptive tuning of crossover values
	// -----------------------------------------------------------------------------------------------------------------------
	// --------------------------------------- Added for reduced sample storage ----------------------------------------------
	MCMCPar->reduced_sample_collection = 0; // Thinned sample collection?
	MCMCPar->reduced_seq_interval = 100;
	MCMCPar->reduced_seq_length = MCMCPar->seq_length / MCMCPar->reduced_seq_interval; // Every Tth sample is collected
	// -----------------------------------------------------------------------------------------------------------------------
	// What type of initial sampling
	MCMCPar->InitPopulation = IP_LHS_BASED;
	// Define the boundary handling
	MCMCPar->BoundHandling = BH_REFLECT;
	// Save in memory or not
	MCMCPar->save_in_memory = 0;
	MCMCPar->save_in_file = 0;
	MCMCPar->verbose = verbose;
	// Give the parameter ranges (minimum and maximum values)
	ParRange.minn = ( double * ) malloc( sizeof( double ) * MCMCPar->n );
	ParRange.maxn = ( double * ) malloc( sizeof( double ) * MCMCPar->n );
	for( i = 0; i < MCMCPar->n; i++ )
	{
		k = od->pd->var_index[i];
		ParRange.minn[i] = od->pd->var_min[k];
		ParRange.maxn[i] = od->pd->var_max[k];
	}
	if( verbose > 4 )
	{
		for( int i = 0; i < MCMCPar->n; i++ )
			printf( "range %g %g\n", ParRange.minn[i], ParRange.maxn[i] );
	}
	//Put the obs_data in the Measure data structure
	Measurement.MeasData = new double[od->od->nCObs];
	Measurement.N = od->od->nCObs;
	for( i = 0, k = 0; i < od->od->nTObs; i++ )
	{
		if( od->od->obs_weight[i] >= 0 )
		{
			Measurement.MeasData[k] = od->od->obs_target[i];
			k++;
		}
	}
	// make ModPred pointer to pass back and forth from functions
	double *ModPred = new double[Measurement.N];
	// make Zinit pointer to pass back and forth from functions
	double **Zinit = double_matrix( MCMCPar->m0 + MCMCPar->seq, MCMCPar->n );;
	// make 3dArray Sequences pointer to pass back and forth from functions
	double ***Sequences = NULL;
	if( MCMCPar->save_in_memory )
	{
		Sequences = (double ***) malloc( sizeof( double ** ) * MCMCPar->seq );
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			Sequences[row] = double_matrix( MCMCPar->seq_length, MCMCPar->np2 );
		}
	}
	// make 3dArray Sequences pointer to pass back and forth from functions
	double ***Reduced_Seq = NULL;
	if( MCMCPar->reduced_sample_collection )
	{
		Reduced_Seq = new double **[MCMCPar->seq];
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			Reduced_Seq[row] = double_matrix( MCMCPar->reduced_seq_length, MCMCPar->np2 );
		}
	}
	// make 2dArray lCR pointer to pass back and forth from functions
	double **lCR = double_matrix( 1, MCMCPar->nCR );;
	// make 2dArray lCRnew pointer to pass back and forth from functions
	double **lCRnew = double_matrix( 1, MCMCPar->nCR );;
	// make 2dArray log_p pointer to pass back and forth from functions
	double **pCR = double_matrix( 1, MCMCPar->nCR );
	// make 1dArray delta_tot pointer to pass back and forth from functions
	double *delta_tot = new double [MCMCPar->nCR];;
	// make 2dArray xnew pointer to pass back and forth from functions
	double **xnew = double_matrix( MCMCPar->seq, MCMCPar->n );
	// make 2dArray xold pointer to pass back and forth from functions
	double **xold = double_matrix( MCMCPar->seq, MCMCPar->n );
	// make 2dArray p pointer to pass back and forth from functions
	double *p = new double[MCMCPar->seq];;
	// make 2dArray log_p pointer to pass back and forth from functions
	double *log_p = new double[MCMCPar->seq];;
	// make 2dArray p_xold pointer to pass back and forth from functions
	double *p_xold = new double[MCMCPar->seq];;
	// make 2dArray log_p_xold pointer to pass back and forth from functions
	double *log_p_xold = new double[MCMCPar->seq];
	// make 1dArray RandArray pointer to pass back and forth from functions
	int *RandArray = new int [MCMCPar->nzoff]; // MCMCPar->nzoff * MCMCPar->DEpairs * MCMCPar->seq
	// make 2dArray CRS pointer to pass back and forth from functions
	double **CRS = double_matrix( MCMCPar->nCR, MCMCPar->steps );
	// make 2dArray DEversion pointer to pass back and forth from functions
	int **DEversion = int_matrix( MCMCPar->seq, MCMCPar->DEpairs );
	// make 1dArray alpha_s pointer to pass back and forth from functions
	double *alpha_s = new double [MCMCPar->seq];;
	// make 2dArray xnew pointer to pass back and forth from functions
	double *p_xnew = new double[MCMCPar->seq];;
	// make 2dArray log_p_xnew pointer to pass back and forth from functions
	double *log_p_xnew = new double[MCMCPar->seq];;
	// make 1dArray accept pointer to pass back and forth from functions
	int *accept = new int [MCMCPar->seq];
	// make 1dArray R_stat pointer to pass back and forth from functions
	double *R_stat = new double [MCMCPar->n];
	// make 2dArray Table_JumpRate pointer to pass back and forth from functions
	double **Table_JumpRate = double_matrix( MCMCPar->n, MCMCPar->DEpairs );
	double **Zee= double_matrix( MCMCPar->nzee, MCMCPar->np2 );
	// output arrays
	int outlines = MCMCPar->ndraw / MCMCPar->seq + MCMCPar->seq;
	output.nEval = new int[ outlines ];
	output.Acceptance_Rate = new double[ outlines ];
	output.CR = double_matrix( outlines, MCMCPar->nCR );
	output.R_stat = double_matrix( outlines, MCMCPar->n );
	meanseqsum = double_matrix( MCMCPar->seq, MCMCPar->n );
	meanseqsum2 = double_matrix( MCMCPar->seq, MCMCPar->n );
	varseqsum = double_matrix( MCMCPar->seq, MCMCPar->n );
	for( int row = 0; row < MCMCPar->seq; row++ )
		for( int col = 0; col < MCMCPar->n; col++ )
			meanseqsum[row][col] = meanseqsum2[row][col] = varseqsum[row][col] = 0;
	// Run the distributed DREAM algorithm with sampling from past
	dream_zs( MCMCPar, ParRange, Measurement, od, option, Reduced_Seq, Zee, output, Zinit,
			  ModPred, p, log_p, Sequences, pCR, xold, p_xold, log_p_xold, RandArray, xnew, CRS, DEversion, alpha_s, Table_JumpRate,
			  accept, p_xnew, log_p_xnew, lCR, delta_tot, R_stat, lCRnew );
	// There are a lot of frees and deletes that still need to be implemented
	delete Measurement.MeasData;
	delete ModPred;
	free_matrix( (void **) Zinit, MCMCPar->m0 + MCMCPar->seq );
	if( MCMCPar->save_in_memory )
	{
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			free_matrix( (void **) Sequences[row], MCMCPar->seq_length );
		}
		free( Sequences );
	}
	if( MCMCPar->reduced_sample_collection )
	{
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			free_matrix( (void **) Reduced_Seq[row], MCMCPar->reduced_seq_length );
		}
		delete Reduced_Seq;
	}
	free_matrix( (void **) lCR, 1 );
	free_matrix( (void **) lCRnew, 1);
	free_matrix( (void **) pCR, 1 );
	delete delta_tot;
	free_matrix( (void **) xnew, MCMCPar->seq );
	free_matrix( (void **) xold, MCMCPar->seq );
	delete p;
	delete log_p;
	delete p_xold;
	delete log_p_xold;
	delete RandArray;
	free_matrix( (void **) CRS, MCMCPar->nCR );
	free_matrix( (void **) DEversion, MCMCPar->seq );
	delete alpha_s;
	delete p_xnew;
	delete log_p_xnew;
	delete accept;
	delete R_stat;
	free_matrix( (void **) Table_JumpRate, MCMCPar->n );
	//don't free Zee, because we're sending it back to the calling function
	MCMCPar->z = Zee;//this is how we send the samples back to the caller
	free( ParRange.minn );
	free( ParRange.maxn );
	delete output.nEval;
	delete output.Acceptance_Rate;
	free_matrix( (void **) output.CR, MCMCPar->nCR );
	free_matrix( (void **) output.R_stat, MCMCPar->n );
	free_matrix( (void **) meanseqsum, MCMCPar->n );
	free_matrix( (void **) meanseqsum2, MCMCPar->n );
	free_matrix( (void **) varseqsum, MCMCPar->n );
	cout << "done." << endl;
	return MCMCPar;
}

//function DREAM_zs
void dream_zs( struct MCMC *MCMCPar, Range ParRange, Measure Measurement, struct opt_data *od, int option,
			   double ***Reduced_Seq, double **Zee, out &output, double **Zinit,
			   double *ModPred, double *p, double *log_p,
			   double ***Sequences, double **pCR,
			   double **xold, double *p_xold, double *log_p_xold,
			   int *RandArray, double **xnew, double **CRS,
			   int **DEversion, double *alpha_s, double **Table_JumpRate,
			   int *accept, double *p_xnew, double *log_p_xnew,
			   double **lCR, double *delta_tot, double *R_stat, double **lCRnew )
{
	ofstream outfile;
	outfile.open( "dream-rstat.dat" );
	int nEvaluations = 0;
	int output_teller, reduced_seq_teller;
	int iloc = 0;
	// Step 0: Initialize variables
	InitVariables( MCMCPar, nEvaluations, output_teller, reduced_seq_teller, Table_JumpRate, Reduced_Seq,
				   output, Sequences, pCR, CRS, lCR );
	// Step 1: Sample MCMCPar->m0 points in the parameter space and store in Z
	if( MCMCPar->InitPopulation == IP_LHS_BASED )
		LHSU( ParRange, MCMCPar, Zinit ); // Latin hypercube sampling
	if( MCMCPar->verbose > 5 )
	{
		printf( "Zinit\n" );
		for( int row = 0; row < MCMCPar->m0; row++ )
		{
			for( int col = 0; col < MCMCPar->n; col++ )
				printf( " %g", Zinit[row][col] );
			printf( "\n" );
		}
	}
	double **Zoff;
	Zoff = double_matrix( MCMCPar->nzoff, MCMCPar->np2 );
	for( int row = MCMCPar->m0; row < MCMCPar->seq; row++ )
		for( int col = 0; col < MCMCPar->np2; col++ )
			Zee[row][col] = 0;
	// Define initial MCMCPar->m0 rows of Z to be initial sample -- posterior density is not needed and thus not evaluated!!
	for( int row = 0; row < MCMCPar->m0; row++ )
		for( int col = 0; col < MCMCPar->n; col++ )
			Zee[row][col] = Zinit[row][col];
	// Define initial population from last MCMCPar->seq samples of Zinit
	double **X2;
	X2 = double_matrix( MCMCPar->seq, MCMCPar->np2 ); // MCMCPar->np2 = MCMCPar->n or 5+2
	if( 1 ) // real random case
	{
		for( int row = 0; row < MCMCPar->seq; row++ )
			for( int col = 0; col < MCMCPar->n; col++ )
				X2[row][col] = Zinit[row + MCMCPar->m0][col];
	}
	else // test case
	{
		double **ex;
		ex = double_matrix( MCMCPar->seq, MCMCPar->n );
		ex[0][0] = 232.1355; // use "ex" for testing onlt
		ex[0][1] = 0.5516;
		ex[0][2] =  0.1130;
		ex[0][3] = 0.0420;
		ex[0][4] =  0.5226;
		ex[1][0] = 84.9467;
		ex[1][1] = 1.4588;
		ex[1][2] =  0.9871;
		ex[1][3] = 0.0828;
		ex[1][4] = 0.6742;
		ex[2][0] = 318.9574;
		ex[2][1] = 0.3821;
		ex[2][2] =  0.1466;
		ex[2][3] = 0.0591;
		ex[2][4] = 0.8211;
		for( int row = 0; row < MCMCPar->seq; row++ )
			for( int col = 0; col < MCMCPar->n; col++ )
				X2[row][col] = ex[row][col];
		free_matrix( ( void ** ) ex, MCMCPar->seq );
	}
	// Calculate posterior density associated with each value of X
	comp_likelihood( X2, MCMCPar, Measurement, od, ModPred, p, log_p );
	// Append X with information about posterior density (or transformation thereof) make it into X2
	for( int row = 0; row < MCMCPar->seq; row++ )
	{
		X2[row][MCMCPar->n] = p[row];
		X2[row][MCMCPar->n + 1] = log_p[row];
	}
	if( MCMCPar->verbose > 5 )
	{
		cout << "X2 init" << endl;
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			for( int col = 0; col < MCMCPar->np2; col++ )
				cout << " " << X2[row][col];
			cout << endl;
		}
	}
	// Initialize the sequences
	if( MCMCPar->save_in_memory )
	{
		for( int k = 0; k < MCMCPar->seq; k++ )
		{
			for( int j = 0; j < MCMCPar->np2; j++ )
			{
				//if( Sequences[k][0] == NULL ) printf( "seq null\n" );
				//printf( "%d, %d\n", k, j );
				Sequences[k][0][j] = X2[k][j];
			}
		}
	}
	GelmanCum( MCMCPar, R_stat, X2 );
	int iloc_2;
	if( MCMCPar->reduced_sample_collection )
		iloc_2 = 0; // Reduced sample collection
	for( int i = 0; i < MCMCPar->nCR; i++ )
		delta_tot[i] = 0;
	// Save N_CR in memory and initialize sum_p2
	output.nEval[0] = nEvaluations;
	for( int i = 0; i < MCMCPar->nCR; i++ )
		output.CR[0][i] = pCR[0][i];
	// Sequences [width = MCMCPar->seq] [rows = iloc] [columns = MCMCPar->n]
	// Compute the R-statistic of Gelman and Rubin
	Gelman( MCMCPar, R_stat, Sequences, 0, iloc );
	for( int col = 0; col < MCMCPar->n; col++ )
		output.R_stat[0][col] = R_stat[col];
	outfile << output.nEval[0];
	for( int col = 0; col < MCMCPar->n; col++ )
		outfile  << " " << output.R_stat[0][col];
	cout << endl;
	outfile << endl;
	if( MCMCPar->save_in_file )
	{
		ofstream outFile;
		char filename[100];
		for( int k = 0; k < MCMCPar->seq; k++ )
		{
			sprintf( filename, "dream_sequence_%04d.dat", k + 1 );
			remove( filename );
		}
	}
	double U;
	string Update;
	// Move prior population to posterior population ...
	int total_accept = 0;
	while( nEvaluations < MCMCPar->ndraw ) // Main Loop
	{
		if( MCMCPar->verbose > 5 ) cout << "\nEvaluations " << nEvaluations << endl;
		// Initialize total accepted realizations;
		int iter_accept = 0;
		// Loop a number of times before calculating convergence diagnostic, etc.
		for( int gen_number = 0; gen_number < MCMCPar->steps; gen_number++ ) // MCMC loop for MCMCPar->steps
		{
			// Initialize teller
			if( MCMCPar->reduced_sample_collection ) reduced_seq_teller++;
			// Define the current locations and associated posterior densities
			if( MCMCPar->verbose ) cout << "GetLocation ..." << endl;
			GetLocation( X2, MCMCPar, xold, p_xold, log_p_xold );
			if( MCMCPar->verbose > 5 )
			{
				printf( "xold\n" );
				for( int row = 0; row < MCMCPar->seq; row++ )
				{
					for( int i = 0; i < MCMCPar->n; i++ )
						printf( " %g", xold[row][i] );
					printf( "\n" );
				}
				printf( "p_xold\n" );
				for( int row = 0; row < MCMCPar->seq; row++ )
					printf( " %g (%g)", p_xold[row], log_p_xold[row] );
				printf( "\n" );
			}
			if( MCMCPar->m < MCMCPar->nzoff )
			{
				// The number of elements of Z is not sufficient
				cout << "size of Z not sufficient to generate offspring with selected MCMCPar->m0, MCMCPar->seq, and MCMCPar->DEpairs" << endl;
			}
			else
			{
				// Without replacement draw rows from Z for proposal creation
				if( MCMCPar->m > MCMCPar->nzee )
				{
					printf( "ERROR: Memory problem. Increase MCMCPar->nzee (%d). MCMCPar->m = %d\n", MCMCPar->nzee, MCMCPar->m );
					exit( 1 );
				}
				if( MCMCPar->verbose ) cout << "Randsample ..." << endl;
				randsample( MCMCPar->m, MCMCPar->nzoff, RandArray, MCMCPar->verbose );
				for( int row = 0; row < MCMCPar->nzoff; row++ )
					for( int col = 0; col < MCMCPar->n; col++ )
						Zoff[row][col] = Zee[RandArray[row]][col];
			}
			if( MCMCPar->verbose > 5 )
			{
				printf( "Zoff\n" );
				for( int row = 0; row < MCMCPar->nzoff; row++ )
				{
					for( int col = 0; col < MCMCPar->n; col++ )
						printf( " %g", Zoff[row][col] );
					printf( "\n" );
				}
			}
			// First generate a random number between 0 and 1
			U = ( double ) rand() / RAND_MAX;
			// Determine to do parallel direction or snooker update depending on random number U
			if( U <= MCMCPar->parallelUpdate ) Update = "Parallel_Direction_Update";
			else                              Update = "Snooker_Update";
			// Generate candidate points (proposal) in each chain using either snooker or parallel direction update
			if( MCMCPar->verbose ) cout << "offde using " << Update << " ..." << endl;
			offde( xold, Zoff, MCMCPar, Update, Table_JumpRate, ParRange, xnew, alpha_s, CRS, DEversion );
			// Compute the likelihood of each proposal in each chain
			if( MCMCPar->verbose ) cout << "comp_likelihood ..." << endl;
			comp_likelihood( xnew, MCMCPar, Measurement, od, ModPred, p_xnew, log_p_xnew );
			// Update number of evaluations
			nEvaluations += MCMCPar->seq;
			if( MCMCPar->verbose > 5 )
			{
				printf( "xnew\n" );
				for( int row = 0; row < MCMCPar->seq; row++ )
				{
					for( int i = 0; i < MCMCPar->n; i++ )
						printf( " %g", xnew[row][i] );
					printf( "\n" );
				}
				printf( "p_xnew\n" );
				for( int row = 0; row < MCMCPar->seq; row++ )
					printf( " %g (%g)", p_xnew[row], log_p_xnew[row] );
				printf( "\n" );
				printf( "xold\n" );
				for( int row = 0; row < MCMCPar->seq; row++ )
				{
					for( int i = 0; i < MCMCPar->n; i++ )
						printf( " %g", xold[row][i] );
					printf( "\n" );
				}
				printf( "p_xold\n" );
				for( int row = 0; row < MCMCPar->seq; row++ )
					printf( " %g (%g)", p_xold[row], log_p_xold[row] );
				printf( "\n" );
			}
			// p_xnew = p_xnew(:,1); TODO check do we need this ..
			// Apply the acceptance/rejectance rule
			if( MCMCPar->verbose ) cout << "metrop ..." << endl;
			metrop( xnew, p_xnew, log_p_xnew, xold, p_xold, log_p_xold, alpha_s, Measurement, MCMCPar, option, X2, accept );
			// How many candidate points have been accepted -- for Acceptance Rate
			int cur_iter_accept = 0;
			for( int col = 0; col < MCMCPar->seq; col++ )
				cur_iter_accept += accept[col];
			iter_accept += cur_iter_accept;
			if( MCMCPar->verbose > 5 )
			{
				cout << "X2 updated" << endl;
				for( int row = 0; row < MCMCPar->seq; row++ )
				{
					for( int col = 0; col < MCMCPar->np2; col++ )
						cout << " " << X2[row][col];
					cout << endl;
				}
			}
			// if( strcmp( MCMCPar->save_in_file.c_str(), "No" ) == 0 && strcmp( MCMCPar->reduced_sample_collection.c_str(), "No" ) == 0 )
			if( 1 ) // it should be always done because it is cheap ...
			{
				GelmanCum( MCMCPar, R_stat, X2 );
				/*printf( "Evals %6d Acceptance Rate %7.4g Accepted %4d R_stat", nEvaluations, ( double ) cur_iter_accept / MCMCPar->seq, cur_iter_accept );
				for( int col = 0; col < MCMCPar->n; col++ )
					printf( " %7.3g", R_stat[col] );
				printf( "\n" );
				*/
			}
			// Check whether to add to sequence or to only store current point
			if( MCMCPar->save_in_memory )
			{
				iloc++; // Define Sequences idx based on iloc
				if( iloc < MCMCPar->seq_length )
				{
					if( MCMCPar->verbose ) cout << "Augment Sequences ..." << endl;
					for( int i = 0; i < MCMCPar->np2; i++ )
						for( int j = 0; j < MCMCPar->seq; j++ )
							Sequences[j][iloc][i] = X2[j][i]; // Update the location of the chains
					if( MCMCPar->verbose > 15 )
					{
						Gelman( MCMCPar, R_stat, Sequences, 0, iloc + 1 ); // this is double check the accuracy of GelmanCum()
						printf( "FEvals %6d R_stat", nEvaluations );
						for( int col = 0; col < MCMCPar->n; col++ )
							printf( " %7.3g", R_stat[col] );
						printf( "\n" );
						for( int k = 0; k < MCMCPar->seq; k++ )
						{
							cout << "Sequence " << k + 1	 << endl;
							for( int i = 0; i < MCMCPar->np2; i++ )
							{
								cout << "Param " << i + 1;
								for( int j = 0; j <= iloc; j++ )
									cout << " " << Sequences[k][j][i];
								cout << endl;
							}
						}
					}
				}
				else
				{
					printf( "No memory to save Sequences! Increase MCMCPar->seq_length ... \n" );
					exit( 1 );
				}
			}
			// Check whether to store in external files
			if( MCMCPar->save_in_file )
			{
				ofstream outFile;
				char filename[100];
				for( int k = 0; k < MCMCPar->seq; k++ )
				{
					sprintf( filename, "dream_sequence_%04d.dat", k + 1 );
					outFile.open( filename, ios::app );
					for( int i = 0; i < MCMCPar->n; i++ )
						outFile << " " << X2[k][i];
					outFile << endl;
					outFile.close();
				}
			}
			// Check whether to store a reduced sample
			if( MCMCPar->reduced_sample_collection )
			{
				if( reduced_seq_teller == MCMCPar->reduced_seq_interval )
				{
					if( MCMCPar->verbose ) cout << "Augment Reduced Sequences ..." << endl;
					reduced_seq_teller = 0; // reset new_teller
					if( iloc_2 < MCMCPar->reduced_seq_length )
						for( int i = 0; i < MCMCPar->np2; i++ )
							for( int j = 0; j < MCMCPar->seq; j++ )
								Reduced_Seq[j][iloc_2][i] = X2[j][i]; // Reduced sample collection
					else
					{
						printf( "No memory to save Reduced_Seq! Increase MCMCPar->reduced_seq_length ... \n" );
						exit( 1 );
					}
					iloc_2++; // Update iloc_2
				}
			}
			double delta_normX[MCMCPar->seq], r[MCMCPar->n];
			double sum, average, tmp;
			int row, col;
			// Compute squared jumping distance for each CR value
			if( MCMCPar->ppCR == PPCR_UPDATE )
			{
				if( MCMCPar->verbose ) cout << "Update CR ..." << endl;
				for( col = 0; col < MCMCPar->n; col++ )
				{
					// Calculate the mean of X
					for( row = 0, sum = 0; row < MCMCPar->seq; row++ )
						sum += X2[row][col];
					average = sum / MCMCPar->seq;
					// Calculate the standard deviation of each dimension of X
					for( row = 0, sum = 0; row < MCMCPar->seq; row++ )
					{
						tmp = X2[row][col] - average;
						sum += tmp * tmp;
					}
					r[col] = sqrt( sum / ( MCMCPar->seq - 1 ) );
				}
				if( MCMCPar->verbose > 5 )
				{
					cout << "R ";
					for( col = 0; col < MCMCPar->n; col++ )
						cout << " " << r[col];
					cout << endl;
				}
				// Compute the Euclidean distance between new X and old X
				for( row = 0; row < MCMCPar->seq; row++ )
				{
					for( col = 0, sum = 0; col < MCMCPar->n; col++ )
					{
						tmp = ( xold[row][col] - X2[row][col] ) / r[col];
						sum += tmp * tmp;
					}
					delta_normX[row] = sum;
				}
				// Use this information to update sum_p2 to update N_CR
				if( MCMCPar->verbose > 3 ) cout << "CalcDelta ..." << endl;
				CalcDelta( MCMCPar, delta_normX, CRS, delta_tot );
			}
			// Check whether to append X to Z
			if( ( ( gen_number + 1 ) % MCMCPar->k ) == 0 )
			{
				// Append X to Z
				if( MCMCPar->verbose ) cout << "Augment Z ..." << endl;
				if( MCMCPar->m + MCMCPar->seq <= MCMCPar->nzee )
					for( row = 0; row < MCMCPar->seq; row++ )
						for( col = 0; col < MCMCPar->np2; col++ )
							Zee[MCMCPar->m + row][col] = X2[row][col];
				else
				{
					printf( "No memory to save Zee Sequences (%d>%d) Increase MCMCPar->nzee ... \n", MCMCPar->m + MCMCPar->seq, MCMCPar->nzee );
					exit( 1 );
				}
				MCMCPar->m += MCMCPar->seq; // Update MCMCPar->m
				if( MCMCPar->verbose > 5 ) // MCMCPar->verbose > 5
				{
					printf( "Zee appended\n" );
					for( row = 0; row < MCMCPar->m; row++ )
					{
						for( col = 0; col < MCMCPar->np2; col++ )
							printf( "%g ", Zee[row][col] );
						printf( "\n" );
					}
				}
			}
			// Update Iteration
			/* This code appears to do nothing
			double RMSE[MCMCPar->seq];
			for( int row = 0; row < MCMCPar->seq; row++ )
				RMSE[row] = sqrt( X2[row][MCMCPar->n - 1] / Measurement.N );
			*/
		}
		if( MCMCPar->verbose ) cout << "\nDone with MCMC steps ..." << endl;
		if( MCMCPar->verbose > 5 )
		{
			cout << "Sequences after MCMC steps" << endl;
			for( int k = 0; k < MCMCPar->seq; k++ )
			{
				cout << "Sequence " << k + 1	 << endl;
				for( int i = 0; i < MCMCPar->n; i++ )
				{
					cout << "Param " << i + 1;
					for( int j = 0; j <= iloc; j++ )
						cout << " " << Sequences[k][j][i];
					cout << endl;
				}
			}
		}
		// Store Important Diagnostic information -- Acceptance Rate
		output.nEval[output_teller] = nEvaluations;
		output.Acceptance_Rate[output_teller] = ( double ) iter_accept / ( MCMCPar->steps * MCMCPar->seq );
		// Store Important Diagnostic information -- Probability of individual crossover values
		for( int col = 0; col < MCMCPar->nCR; col++ )
			output.CR[output_teller][col] = pCR[0][col];
		// Check whether to update individual pCR values
		if( nEvaluations <= ( 0.1 * MCMCPar->ndraw ) )
			if( MCMCPar->ppCR == PPCR_UPDATE )
				AdaptpCR( MCMCPar, delta_tot, lCR, pCR ); // Update pCR values
		// Generate CR values based on current pCR values
		if( MCMCPar->verbose ) cout << "GenCR ..." << endl;
		GenCR( MCMCPar, lCRnew, pCR, CRS );
		for( int col = 0; col < MCMCPar->seq; col++ )
			lCR[0][col] += lCRnew[0][col];
		int start_loc;
		int end_loc;
		// Calculate Gelman and Rubin Convergence Diagnostic
		if( MCMCPar->save_in_memory )
		{
			// Use the full sequence
			if( iloc / 2 > 1 ) start_loc = iloc / 2 - 1;
			else start_loc = 0;
			end_loc = iloc + 1; // to use < end_loc instead of <=iloc; this is the number of locations
			// Compute the R-statistic using 50% burn-in from Sequences
			if( MCMCPar->verbose ) cout << "Gelman Convergence Diagnostic using Full Sequences ... (" << start_loc  << "," << end_loc << ")" << endl;
			Gelman( MCMCPar, R_stat, Sequences, 0, end_loc );
			for( int col = 0; col < MCMCPar->n; col++ )
				output.R_stat[output_teller][col] = R_stat[col];
		}
		else if( MCMCPar->reduced_sample_collection )
		{
			// Use the reduced sequence
			if( iloc_2 / 2 > 1 ) start_loc = iloc_2 / 2;
			else start_loc = 0;
			end_loc = iloc_2; // to use < end_loc instead of <=iloc; this is the number of locations
			// Compute the R-statistic using 50% burn-in from Reduced_Seq
			if( MCMCPar->verbose ) cout << "Gelman Convergence Diagnostic using Reduced Sequences ... (" << start_loc  << "," << end_loc << ")" << endl;
			Gelman( MCMCPar, R_stat, Reduced_Seq, start_loc, end_loc );
			for( int col = 0; col < MCMCPar->n; col++ )
				output.R_stat[output_teller][col] = R_stat[col];
		}
		else // R_stat are already computed using compute using GelmanCum()
			for( int col = 0; col < MCMCPar->n; col++ )
				output.R_stat[output_teller][col] = R_stat[col];
		if( MCMCPar->save_in_memory || MCMCPar->reduced_sample_collection ) // If GelmanCum() only no output needed
		{
			// produce screen output
			printf( "Evals %6d Acceptance Rate %7.4g Accepted %4d R_stat", output.nEval[output_teller], output.Acceptance_Rate[output_teller], iter_accept );
			for( int col = 0; col < MCMCPar->n; col++ )
				printf( " %7.3g", output.R_stat[output_teller][col] );
			cout << " (full)" << endl;
		}
		// produce file output
		outfile << output.nEval[output_teller];
		for( int col = 0; col < MCMCPar->n; col++ )
			outfile << " " << output.R_stat[output_teller][col];
		outfile << endl;
		outfile.flush();
		// Update the teller
		output_teller++;
		total_accept += iter_accept;
	}
	printf( "Total Accepted %5d\n", total_accept );
	outfile.close();
}

// Latin Hypercube sampling
void LHSU( Range ParRange, struct MCMC *MCMCPar, double **Zinit )
{
	// Declare and initialize variables
	double xmin[MCMCPar->n];
	double xmax[MCMCPar->n];
	// min, max, nsample
	for( int i = 0; i < MCMCPar->n; i++ )
	{
		xmin[i] = ParRange.minn[i];
		xmax[i] = ParRange.maxn[i];
	}
	int nsample = MCMCPar->m0 + MCMCPar->seq;
	int nvar = MCMCPar->n;
	double ran[nsample][nvar];
	// Initialize array ran with random numbers
	for( int row = 0; row < nsample; row++ )
		for( int col = 0; col < nvar; col++ )
			ran[row][col] = ( double ) rand() / RAND_MAX; //random number between 0 & 1
	double idx[nsample];
	double P[nsample];
	double sum[nsample][nvar];
	// Initialize array s  to 0
	for( int i = 0; i < nsample; i++ )
		for( int j = 0; j < nvar; j++ )
			Zinit[i][j] = 0;
	//random permutation for idx
	for( int i = 0; i < nsample; ++i )
	{
		int j = rand() % ( i + 1 );
		idx[i] = idx[j];
		idx[j] = i + 1;
	}
	for( int j = 0; j < nvar; j++ )
	{
		for( int i = 0; i < nsample; i++ )
			P[i] = ( idx[i] - ran[i][j] ) / nsample;
		// Now fill s
		for( int i = 0; i < nsample; i++ )
		{
			sum[i][j] = P[i] * ( xmax[j] - xmin[j] ); //broken into 2 parts
			Zinit[i][j] = xmin[j] + sum[i][j];        // fill the jth column
		}
	}
}

// This function computes the likelihood for each value of x
void comp_likelihood( double **x, struct MCMC *MCMCPar, Measure Measurement, struct opt_data *od, double *ModPred, double *p, double *log_p )
{
	int ii, i, k;
	// Loop over the individual parameter combinations of x
	for( ii = 0; ii < MCMCPar->seq; ii++ )
	{
		log_p[ii] = 0;
		// Call model to generate simulated data
		func_global( x[ii], od, ModPred );
		//Note that after calling func_intrn, ModPred does not contain the model predictions.
		//ModPred contains the errors times the weights
		//The model predictions are contained in od->od->obs_current[k]
		for( i = 0, k = 0; i < od->od->nTObs; i++ )
		{
			if( od->od->obs_weight[i] >= 0 )
			{
				//TODO: Implement full covariance stuff
				//In these calculations, the weight is basically 1 / (2 * the variance in the measurement error)
				log_p[ii] -= od->od->obs_weight[i] * ( Measurement.MeasData[k] - od->od->obs_current[i] ) * ( Measurement.MeasData[k] - od->od->obs_current[i] );
				k++;
			}
			p[ii] = exp( log_p[ii] );
		}
	}
}


// function initializes important variables for use in the algorithm
void InitVariables( struct MCMC *MCMCPar, int &nEval, int &output_teller, int &reduced_seq_teller, double **Table_JumpRate, double ***Reduced_Seq,
					out &output, double ***Sequences, double **pCR, double **CRS, double **lCR )
{
	// Calculate the parameters in the exponential power density function of Box and Tiao (1973)
	double Cb, Wb;
	CalcCbWb( MCMCPar, Cb, Wb );
	double ones[3];
	// Define the crossover values as geometrical series
	for( int col = 0; col < MCMCPar->nCR; col++ )
		ones[col] = ( double ) 1 * ( 1 / MCMCPar->nCR ); //.3333
	MCMCPar->CR[0][0] = ones[0];
	for( int i = 1; i < MCMCPar->nCR; i++ )
		MCMCPar->CR[0][i] = MCMCPar->CR[0][i - 1] + ones[i]; // cumulative sum (cumsum)
	// Derive the number of elements in the output file
	double Nelem;
	Nelem = floor( ( double ) MCMCPar->ndraw / MCMCPar->seq ) + 1;
	double q;
	q = floor( ( double ) Nelem / MCMCPar->steps );
	int cr;
	cr = static_cast<int>( q ); //change from double to integer
	// Initialize output information -- N_CR
	for( int row = 0; row < cr; row++ )
		for( int col = 0; col < MCMCPar->nCR; col++ )
			output.CR[row][col] = 0;
	// Initialize output information -- Acceptance_Rate
	for( int row = 0; row < cr; row++ )
		output.Acceptance_Rate[row] = 0;
	output.nEval[0] = MCMCPar->seq;
	output.Acceptance_Rate[0] = -1;
	// Initialize output information -- R statistic
	for( int i = 0; i < cr; i++ )
		for( int j = 0; j < MCMCPar->n; j++ )
			output.R_stat[i][j] = 0;
	// Calculate multinomial probabilities of each of the nCR CR values
	for( int i = 0; i < 1; i++ )
		for( int j = 0; j < MCMCPar->nCR; j++ )
			pCR[i][j] = ( double ) 1 * ( 1 / MCMCPar->nCR );
	// Calculate the actual CR values based on p
	GenCR( MCMCPar, lCR, pCR, CRS );
	// Check what to save in memory
	if( MCMCPar->save_in_memory )
	{
		// Initialize Sequences with zeros
		for( int i = 0; i < MCMCPar->seq; i++ )
			for( int j = 0; j < MCMCPar->seq_length; j++ )
				for( int k = 0; k < MCMCPar->np2; k++ )
					Sequences[i][j][k] = 0;
	}
	if( MCMCPar->reduced_sample_collection )
	{
		// Initialize Reduced Sequences with zeros
		reduced_seq_teller = 0; // reduced sequence count
		for( int i = 0; i < MCMCPar->seq; i++ )
			for( int j = 0; j < MCMCPar->n; j++ )
				for( int k = 0; k < MCMCPar->np2; k++ )
					Reduced_Seq[i][j][k] = 0;
	}
	// Generate the Table with JumpRates (dependent on number of dimensions and number of pairs
	for( int zz = 0; zz < MCMCPar->DEpairs; zz++ )
		for( int i = 0; i < MCMCPar->n; i++ )
			Table_JumpRate[i][zz] = ( double ) 2.38 / sqrt( 2 * ( zz + 1 ) * ( i + 1 ) );
	nEval = MCMCPar->seq; // Initialize nEval
	output_teller = 0; // output count
}

// This function calculates the parameters for the exponential power density
void CalcCbWb( struct MCMC *MCMCPar, double &Cb, double &Wb )
{
	// Equation [20] paper by Thiemann et al. WRR 2001, Vol 37, No 10, 2521-2535
	double A1, A2;
	A1 = tgamma( ( double ) 3 * ( MCMCPar->gamma + 1 )  / 2 );
	A2 = tgamma( ( double )( MCMCPar->gamma + 1 ) / 2 );
	Cb = pow( A1 / A2, ( double ) 1 / ( MCMCPar->gamma + 1 ) );
	Wb = sqrt( A1 / ( pow( A2, 1.5 ) ) * MCMCPar->gamma + 1 );
}

// This function generates CR values based on current probabilities
void GenCR( struct MCMC *MCMCPar, double **lCR, double **pCR, double **CRS )
{
	int m = 1;
	double L2[MCMCPar->nCR];
	double **Y;
	Y = double_matrix( m, MCMCPar->nCR );
	int nn = MCMCPar->seq * MCMCPar->steps;
	// How many candidate points for each crossover value?
	multrnd( nn, m, MCMCPar->nCR, lCR, Y, pCR ); // TODO why 1 ?!
	if( MCMCPar->verbose > 5 )
	{
		printf( "lCR" );
		for( int col = 0; col < MCMCPar->nCR; col++ )
			printf( " %g", lCR[0][col] );
		printf( "\n" );
	}
	for( int col = 0; col < MCMCPar->nCR; col++ )
		L2[col] = lCR[0][col];
	for( int col = 1; col < MCMCPar->nCR; col++ )
		L2[col] += L2[col - 1]; // cumsum TODO L2 not used !?
	if( MCMCPar->verbose > 5 )
	{
		printf( "L2" );
		for( int col = 0; col < MCMCPar->nCR; col++ )
			printf( " %g", L2[col] );
		printf( "\n" );
	}
	int RR[MCMCPar->seq * MCMCPar->steps];
	int idx[MCMCPar->seq * MCMCPar->steps];
	// Then select which candidate points are selected with what CR random permutation for r
	for( int i = 0; i < ( MCMCPar->seq * MCMCPar->steps ); ++i ) // random permutation
	{
		int j = rand() % ( i + 1 );
		RR[i] = RR[j];
		RR[j] = i;
	}
	if( MCMCPar->verbose > 5 )
	{
		printf( "RR" );
		for( int col = 0; col < ( MCMCPar->seq * MCMCPar->steps ); col++ )
			printf( " %i", RR[col] );
		printf( "\n" );
	}
	// Define start and end
	int i_start = 0;
	int i_end = MCMCPar->np2;
	double CRR[MCMCPar->seq * MCMCPar->steps];
	// Then generate CR values for each chain
	for( int zz = 0; zz < MCMCPar->nCR; zz++ )
	{
		if( i_start == ( MCMCPar->seq * MCMCPar->steps ) / 2 - 1 )
		{
			i_end = MCMCPar->seq * MCMCPar->steps;
			for( int start = i_start; start < i_end; start++ )
				idx[start] = RR[start];
		}
		// Assign these indices MCMCPar->CR(zz)
		for( int start = i_start; start < i_end; start++ )
			idx[start] = RR[start];
		for( int i = i_start; i < i_end; i++ )
			CRR[idx[i]] = ( double )( zz + 1 ) / MCMCPar->nCR;
		// new start and end
		i_start = i_end;
		i_end += MCMCPar->np2;
	}
	// Now reshape CR
	int counts = 0;
	while( counts < ( MCMCPar->seq * MCMCPar->steps ) )
	{
		for( int i = 0; i < MCMCPar->steps; i++ )
			for( int j = 0; j < MCMCPar->nCR; j++ )
				CRS[j][i] = CRR[counts++];
	}
	if( MCMCPar->verbose > 5 )
	{
		cout << "CR\n";
		for( int i = 0; i < MCMCPar->steps; i++ )
		{
			for( int j = 0; j < MCMCPar->nCR; j++ )
				cout << " " << CRS[j][i];
			cout << endl;
		}
	}
}

void multrnd( int nn, int m, int k, double **lCR, double **Y, double **pCR )
{
	double P = 0;
	for( int i = 0; i < m; i++ )
		for( int j = 0; j < k; j++ )
			P += pCR[i][j];
	int o[m][nn];
	double s[m][nn];
	double r[m][nn];
	int countP = 0;
	for( int i = 0; i < m; i++ )
	{
		s[i][0] = pCR[i][0] / P;
		for( int col = 0; col < nn; col++ )
		{
			o[i][col] = 1;
			r[i][col] = ( double ) rand() / RAND_MAX;
		}
		for( int col = 1; col < k; col++ )
			s[i][col] = pCR[i][col] / P + s[i][col - 1]; //cumsum
		for( int col = 0; col < k; col++ )
			countP++;
	}
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < countP; j++ )
			for( int col = 0; col < nn; col++ )
				if( r[i][col] > s[i][j] )
					o[i][col]++;
		for( int c = 0; c < countP; c++ )
			lCR[i][c] = 0;
		for( int j = 1; j <= countP; j++ )
			for( int col = 0; col < nn; col++ )
				if( o[i][col] == j )
					lCR[i][j - 1]++;
	}
	for( int i = 0; i < m; i++ )
		for( int col = 0; col < countP; col++ )
			Y[i][col] = lCR[i][col] / nn;
}


// Calculates the R-statistic convergence diagnostic
/* ----------------------------------------------------
   For more information please refer to: Gelman, A. and D.R. Rubin, 1992.
   Inference from Iterative Simulation Using Multiple Sequences,
   Statistical Science, Volume 7, Issue 4, 457-472.
   ----------------------------------------------------*/
void Gelman( struct MCMC *MCMCPar, double *R_stat, double ***Sequences, int start_loc, int end_loc )
{
	int verbose;
	verbose = MCMCPar->verbose > 5;
	int seq_size = end_loc - start_loc;
	if( verbose ) printf( "Gelman: size %d range %d %d\n", seq_size, start_loc, end_loc ); //MCMCPar->verbose
	if( seq_size < 2 ) // Set the R-statistic to a large value
	{
		for( int i = 0; i < MCMCPar->n; i++ )
			R_stat[i] = HUGE_VAL;
		return;
	}
	// Compute R-statistics
	// Write the current Sequences
	int i, k, j;
	if( verbose > 2 )
	{
		for( k = 0; k < MCMCPar->seq; k++ )
		{
			for( i = 0; i < MCMCPar->n; i++ )
			{
				for( j = start_loc; j < end_loc; j++ )
					printf( " %g", Sequences[k][j][i] );
				printf( "\n" );
			}
			printf( "\n" );
		}
		printf( "\n" );
	}
	double sum;
	double **meanSeq;
	meanSeq = double_matrix( MCMCPar->seq, MCMCPar->n );
	// Step 1: Determine the sequence means
	if( verbose ) cout << "mean sequence" << endl;
	for( k = 0; k < MCMCPar->seq; k++ )
	{
		for( i = 0; i < MCMCPar->n; i++ )
		{
			for( j = start_loc, sum = 0; j < end_loc; j++ )
				sum += Sequences[k][j][i];
			meanSeq[k][i] = sum / seq_size;
			if( verbose ) cout << " " << meanSeq[k][i]; // MCMCPar->verbose
		}
		if( verbose ) cout << endl;
	}
	if( verbose ) cout << endl;
	double mean[MCMCPar->n];
	// Calculate the mean
	if( verbose ) cout << "mean parameter" << endl;
	for( i = 0; i < MCMCPar->n; i++ )
	{
		for( k = 0, sum = 0; k < MCMCPar->seq; k++ )
			sum += meanSeq[k][i];
		mean[i] = sum / MCMCPar->seq;
		if( verbose ) cout << " " << mean[i];
	}
	if( verbose ) cout << endl;
	double stdDev_loc[MCMCPar->n];
	// Step 1: Determine the variance between the sequence means
	double tmp;
	if( verbose ) cout << "stddev parameter" << endl;
	for( i = 0; i < MCMCPar->n; i++ )
	{
		for( k = 0, sum = 0; k < MCMCPar->seq; k++ )
		{
			tmp = meanSeq[k][i] - mean[i];
			sum += tmp * tmp;
		}
		stdDev_loc[i] = sum * seq_size / ( MCMCPar->seq - 1 );
		if( verbose ) cout << " " << stdDev_loc[i];
	}
	if( verbose ) cout << endl;
	// Step 2: Compute the variance of the various sequences
	double W[MCMCPar->n];
	double varSeq[MCMCPar->seq][MCMCPar->n];
	for( k = 0; k < MCMCPar->seq; k++ )
	{
		for( i = 0; i < MCMCPar->n; i++ )
		{
			for( j = start_loc, sum = 0; j < end_loc; j++ )
			{
				tmp = Sequences[k][j][i] - meanSeq[k][i];
				sum += tmp * tmp;
				//				printf( "s%d %g\n ", j, sum );
			}
			varSeq[k][i] = sum / ( seq_size - 1 );
		}
	}
	if( verbose )
	{
		cout << "var sequence\n";
		for( k = 0; k < MCMCPar->seq; k++ )
		{
			for( i = 0; i < MCMCPar->n; i++ )
				cout << varSeq[k][i] << " ";
			cout << endl;
		}
	}
	// Step 3: Calculate the average of the within sequence variances
	for( i = 0; i < MCMCPar->n; i++ )
	{
		for( k = 0, sum = 0; k < MCMCPar->seq; k++ )
			sum += varSeq[k][i];
		W[i] = sum / MCMCPar->seq;
	}
	/* Step 3: Estimate the target mean
	mu = mean(meanSeq); */
	double sigma;
	for( i = 0; i < MCMCPar->n; i++ )
	{
		sigma = W[i] * ( ( double )( seq_size - 1 ) / seq_size ) + stdDev_loc[i] * ( ( double ) 1.0 / seq_size );  // Step 4: Estimate the target variance
		R_stat[i] = ( double )( MCMCPar->seq + 1.0 ) / MCMCPar->seq * sigma / W[i] - ( double )( seq_size - 1.0 ) / MCMCPar->seq / seq_size;   // Step 5: Compute the R-statistic
	}
	if( verbose )
	{
		cout << "R_stat " << endl;
		for( i = 0; i < MCMCPar->n; i++ )
			cout << " " << R_stat[i];
		cout << endl;
	}
}

// Calculates the R-statistic convergence diagnostic
/* ----------------------------------------------------
   For more information please refer to: Gelman, A. and D.R. Rubin, 1992.
   Inference from Iterative Simulation Using Multiple Sequences,
   Statistical Science, Volume 7, Issue 4, 457-472.
   ----------------------------------------------------*/
void GelmanCum( struct MCMC *MCMCPar, double *R_stat, double **X2 )
{
	int verbose;
	int i, k;
	verbose = MCMCPar->verbose > 5;
	cum_seq_size++;
	if( verbose ) printf( "Gelman Cummulative: size %d\n", cum_seq_size ); //MCMCPar->verbose
	if( cum_seq_size < 2 )
	{
		for( i = 0; i < MCMCPar->n; i++ )
			R_stat[i] = HUGE_VAL; // Set the R-statistic to a large value
		for( k = 0; k < MCMCPar->seq; k++ )
			for( i = 0; i < MCMCPar->n; i++ )
				meanseqsum[k][i] = meanseqsum2[k][i] += X2[k][i];
		return;
	}
	// Compute R-statistics
	// Write the current Sequences
	double **meanSeq;
	meanSeq = double_matrix( MCMCPar->seq, MCMCPar->n );
	// Step 1: Determine the sequence means
	if( verbose ) cout << "mean sequence" << endl;
	for( k = 0; k < MCMCPar->seq; k++ )
	{
		for( i = 0; i < MCMCPar->n; i++ )
		{
			meanseqsum[k][i] += X2[k][i];
			meanSeq[k][i] = meanseqsum[k][i] / cum_seq_size;
			if( verbose ) cout << " " << meanSeq[k][i]; // MCMCPar->verbose
		}
		if( verbose ) cout << endl;
	}
	if( verbose ) cout << endl;
	double sum;
	double mean[MCMCPar->n];
	// Calculate the mean
	if( verbose ) cout << "mean parameter" << endl;
	for( i = 0; i < MCMCPar->n; i++ )
	{
		for( k = 0, sum = 0; k < MCMCPar->seq; k++ )
			sum += meanSeq[k][i];
		mean[i] = sum / MCMCPar->seq;
		if( verbose ) cout << " " << mean[i];
	}
	if( verbose ) cout << endl;
	double stdDev_loc[MCMCPar->n];
	// Step 1: Determine the variance between the sequence means
	double tmp;
	if( verbose ) cout << "stddev parameter" << endl;
	for( i = 0; i < MCMCPar->n; i++ )
	{
		for( k = 0, sum = 0; k < MCMCPar->seq; k++ )
		{
			tmp = meanSeq[k][i] - mean[i];
			sum += tmp * tmp;
		}
		stdDev_loc[i] = sum * cum_seq_size / ( MCMCPar->seq - 1 );
		if( verbose ) cout << " " << stdDev_loc[i];
	}
	if( verbose ) cout << endl;
	// Step 2: Compute the variance of the various sequences
	double varSeq[MCMCPar->seq][MCMCPar->n];
	for( k = 0; k < MCMCPar->seq; k++ )
	{
		for( i = 0; i < MCMCPar->n; i++ )
		{
			tmp = X2[k][i] - meanseqsum2[k][i];
			meanseqsum2[k][i] += tmp / cum_seq_size;
			varseqsum[k][i] += tmp * ( X2[k][i] - meanseqsum2[k][i] );
			varSeq[k][i] = varseqsum[k][i] / ( cum_seq_size - 1 );
		}
	}
	if( verbose )
	{
		cout << "var sequence\n";
		for( k = 0; k < MCMCPar->seq; k++ )
		{
			for( i = 0; i < MCMCPar->n; i++ )
				cout << varSeq[k][i] << " ";
			cout << endl;
		}
	}
	// Step 3: Calculate the average of the within sequence variances
	double W[MCMCPar->n];
	for( i = 0; i < MCMCPar->n; i++ )
	{
		for( k = 0, sum = 0; k < MCMCPar->seq; k++ )
			sum += varSeq[k][i];
		W[i] = sum / MCMCPar->seq;
	}
	/* Step 3: Estimate the target mean
	mu = mean(meanSeq); */
	double sigma;
	for( i = 0; i < MCMCPar->n; i++ )
	{
		sigma = W[i] * ( ( double )( cum_seq_size - 1 ) / cum_seq_size ) + stdDev_loc[i] * ( ( double ) 1.0 / cum_seq_size );  // Step 4: Estimate the target variance
		R_stat[i] = ( double )( MCMCPar->seq + 1.0 ) / MCMCPar->seq * sigma / W[i] - ( double )( cum_seq_size - 1.0 ) / MCMCPar->seq / cum_seq_size;   // Step 5: Compute the R-statistic
	}
	if( verbose )
	{
		cout << "R_stat " << endl;
		for( i = 0; i < MCMCPar->n; i++ )
			cout << " " << R_stat[i];
		cout << endl;
	}
}

// Extracts the current location and density of the chain
void GetLocation( double **X2, struct MCMC *MCMCPar, double **xold, double *p_xold, double *log_p_xold )
{
	// First get the current location
	for( int row = 0; row < MCMCPar->seq; row++ )
		for( int col = 0; col < MCMCPar->n; col++ )
			xold[row][col] = X2[row][col];
	// Then get the current density
	for( int row = 0; row < MCMCPar->seq; row++ )
		p_xold[row] = X2[row][MCMCPar->n];
	// Then get the current log density
	for( int row = 0; row < MCMCPar->seq; row++ )
		log_p_xold[row] = X2[row][MCMCPar->n + 1];
}

// RANDSAMPLE Random sampling, without replacement
void randsample( int n, int k, int *ro, int verbose )
{
	double q;
	int i, a, rp[n];
	double sumx;
	if( ( 4 * k ) > n )
	{
		for( i = 0; i < n; ++i ) //random permutation
			rp[i] = i;
		for( i = 0; i < n; ++i )
		{
			int j = rand() % ( i + 1 );
			rp[i] = rp[j];
			rp[j] = i;
		}
		for( i = 0; i < k; i++ )
			ro[i] = rp[i];
	}
	/*If the sample is a small fraction of the population, a full
	    sort is wasteful.  Repeatedly sample with replacement until
	    there are k unique values.*/
	else
	{
		for( i = 0; i < n; i++ )
			rp[i] = 0; // flags
		sumx = 0;
		while( sumx < k )
		{
			for( i = 0; i < k; i++ )
			{
				// sample w/replacement
				q = ( ( double ) rand() / RAND_MAX ) * n;
				q = ceil( q ) - 1;
				a = static_cast<int>( q );
				ro[i] = a;
				rp[a] = 1;
				// count how many unique elements so far
				sumx += rp[a];
			}
		};
		if( verbose > 15 )
		{
			printf( "rp\n" );
			for( i = 0; i < n; i++ )
				printf( " %d", rp[i] );
			printf( "\n" );
			printf( "ro\n" );
			for( i = 0; i < k; i++ )
				printf( " %d", ro[i] );
			printf( "\n" );
		}
	}
	if( verbose > 5 )
	{
		printf( "random sample: " );
		for( i = 0; i < k; i++ )
			printf( " %d", ro[i] );
		printf( "\n" );
	}
	/* a scalar loop version
	     x = 1:n;
	     n = n:(-1):(n-k+1);
	     y = zeros(1,k);
	     j = ceil(n .* rand(1,k));
	     for i = 1:k
	         y(i) = x(j(i));
	         x(j(i)) = x(n(i));
	 */
}

// Generates offspring using METROPOLIS HASTINGS monte-carlo markov chain
void offde( double **xold, double **Zoff, struct MCMC *MCMCPar, string Update, double **Table_JumpRate, Range ParRange,
			double **xnew, double *alpha_s, double **CRS, int **DEversion )
{
	double xjump[MCMCPar->seq][MCMCPar->n]; // TODO replace MCMCPar->seq with MCMCPar->nCR ?!
	double noise_x[MCMCPar->seq][MCMCPar->n];
	double D[MCMCPar->seq][MCMCPar->n];
	int rr[MCMCPar->seq][4];
	int it[MCMCPar->seq][MCMCPar->n];
	int count = 0;
	double jump[MCMCPar->n];
	double Gamma = 0;
	// Determine how many pairs to use for each jump in each chain
	DEStrategy( MCMCPar, DEversion );
	if( MCMCPar->verbose > 5 )
	{
		printf( "DEversion offde\n" );
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			for( int i = 0; i < MCMCPar->DEpairs; i++ )
				printf( " %d", DEversion[row][i] );
			printf( "\n" );
		}
	}
	// Generate uniform random numbers for each chain to determine which dimension to update
	for( int i = 0; i < MCMCPar->seq; i++ )
		for( int j = 0; j < MCMCPar->n; j++ )
			D[i][j] = ( double ) rand() / RAND_MAX;    //random numbers between 0 and 1
	// Generate noise to ensure ergodicity for each individual chain
	for( int i = 0; i < MCMCPar->seq; i++ )
		for( int j = 0; j < MCMCPar->n; j++ )
			noise_x[i][j] = MCMCPar->eps * ( ( double ) rand() * 2 / RAND_MAX - 1 );
	// Initialize the delta update to zero
	for( int i = 0; i < MCMCPar->seq; i++ )
		for( int j = 0; j < MCMCPar->n; j++ )
			xjump[i][j] = 0;
	if( strcmp( Update.c_str(), "Parallel_Direction_Update" ) == 0 ) // PARALLEL DIRECTION UPDATER
	{
		if( MCMCPar->verbose > 3 ) cout << "offde Parallel_Direction_Update" << endl;
		// Define which points of Zoff to use to generate jumps
		rr[0][0] = 0;
		rr[0][1] = rr[0][0] + DEversion[0][0] - 1;
		rr[0][2] = rr[0][1] + 1;
		rr[0][3] = rr[0][2] + DEversion[0][0] - 1;
		// Do this for each chain
		for( int qq = 1; qq < MCMCPar->seq; qq++ )
		{
			// Define rr to be used for population evolution
			rr[qq][0] = rr[qq - 1][3] + 1;
			rr[qq][1] = rr[qq][0] + DEversion[qq][0] - 1;
			rr[qq][2] = rr[qq][1] + 1;
			rr[qq][3] = rr[qq][2] + DEversion[qq][0] - 1;
		}
		// Each chain evolves using information from other chains to create offspring
		for( int qq = 0; qq < MCMCPar->seq; qq++ )
		{
			// ------------ WHICH DIMENSIONS TO UPDATE? USE CROSSOVER ----------
			count = 0;
			for( int i = 0; i < MCMCPar->n; i++ )
			{
				if( MCMCPar->verbose > 5 ) printf( "c %g %g\n", D[qq][i], 1 - CRS[qq][0] );
				if( D[qq][i] > ( double ) 1 - CRS[qq][0] ) it[qq][i] = i; // TODO CRS[qq][0] does not change after the optimization is initiated; not good!
				else                                   { it[qq][i] = MCMCPar->n; count++; } // Ignore this direction
			}
			if( MCMCPar->verbose > 5 )
			{
				printf( "it[qq] " );
				for( int k = 0; k < MCMCPar->n; k++ )
					printf( " %i", it[qq][k] );
				printf( "\n" );
			}
			// if( it[qq][0] == MCMCPar->n && it[qq][1] == MCMCPar->n && it[qq][2] == MCMCPar->n && it[qq][3] == MCMCPar->n && it[qq][4] == MCMCPar->n )
			if( count == MCMCPar->n ) // If all the directions are ignored ...
			{
				for( int i = 0; i < MCMCPar->n; ++i ) // random permutation
				{
					int j = rand() % ( i + 1 );
					it[qq][i] = it[qq][j];
					it[qq][j] = i;
				}
				for( int k = 0; k < MCMCPar->n; k++ )
					it[qq][k] = it[qq][0];
				if( MCMCPar->verbose > 5 )
				{
					printf( "it[qq] new " );
					for( int k = 0; k < MCMCPar->n; k++ )
						printf( " %i", it[qq][k] );
					printf( "\n" );
				}
			}
			// -----------------------------------------------------------------
			// Generate a random number between 0 and 1
			double U = ( double ) rand() / RAND_MAX;
			// Select the appropriate JumpRate and create a jump
			if( U < 0.8 )
			{
				if( MCMCPar->verbose > 5 ) cout << "offde not done ul " << qq << " / " << MCMCPar->seq << endl;
				// Select the JumpRate (dependent of NrDim and number of pairs)
				int NrDim = count - 1;
				if( NrDim < 0 ) NrDim = 0;
				Gamma = Table_JumpRate[NrDim][0];
				// Produce the difference of the pairs used for population evolution
				if( MCMCPar->verbose > 5 ) cout << "offde not done ul " << rr[qq][0] << " - " << rr[qq][3] << endl;
				for( int j = 0; j < MCMCPar->n; j++ )
					jump[j] = Zoff[rr[qq][0]][j] - Zoff[rr[qq][3]][j];
				if( MCMCPar->verbose > 6 )
				{
					printf( "jump " );
					for( int j = 0; j < MCMCPar->n; j++ )
						printf( " %g", jump[j] );
					printf( "\n" );
					printf( "Zoff[rr[qq][0]][j] " );
					for( int j = 0; j < MCMCPar->n; j++ )
						printf( " %g", Zoff[rr[qq][0]][j] );
					printf( "\n" );
					printf( "Zoff[rr[qq][3]][j] " );
					for( int j = 0; j < MCMCPar->n; j++ )
						printf( " %g", Zoff[rr[qq][3]][j] );
					printf( "\n" );
				}
				// Then fill update the dimension
				for( int i = 0; i < MCMCPar->n; i++ )
				{
					if( it[qq][i] != MCMCPar->n ) xjump[qq][i] = Gamma * ( noise_x[qq][it[qq][i]] + 1 ) * jump[it[qq][i]];
					else                         xjump[qq][i] = 0;
				}
				if( MCMCPar->verbose > 5 ) cout << "offde done ul " << qq << " / " << MCMCPar->seq << endl;
			}
			else
			{
				if( MCMCPar->verbose > 5 ) cout << "offde not done ug " << qq << " / " << MCMCPar->seq << endl;
				// Full space jump
				Gamma = 1;
				CRS[qq][0] = 1;
				// Compute delta from one pair
				if( MCMCPar->verbose > 5 ) cout << "offde not done ug " << rr[qq][0] << " - " << rr[qq][3] << endl;
				for( int j = 0; j < MCMCPar->n; j++ )
					jump[j] = Zoff[rr[qq][0]][j] - Zoff[rr[qq][3]][j];
				// Now jumprate to facilitate jumping from one mode to the other in all dimensions
				for( int j = 0; j < MCMCPar->n; j++ )
					xjump[qq][j] = Gamma * jump[j];
				if( MCMCPar->verbose > 5 ) cout << "offde done ug " << qq << " / " << MCMCPar->seq << endl;
			}
			count = 0;
		}
	}
	//declare variables
	double zR1[MCMCPar->n];
	double zR2[MCMCPar->n];
	int r[MCMCPar->nzoff];
	int t[MCMCPar->nzoff - 2];
	double z[MCMCPar->seq][MCMCPar->n];
	double F[MCMCPar->seq][MCMCPar->n];
	double Fnew[MCMCPar->seq][MCMCPar->n];
	double zP[MCMCPar->seq][MCMCPar->n];	string Update2;
	if( strcmp( Update.c_str(), "Snooker_Update" ) == 0 ) // SNOOKER UPDATER
	{
		if( MCMCPar->verbose > 3 ) cout << "offde Snooker_Update" << endl;
		// Determine the number of rows of Zoff
		// Define rr
		int k = 0;
		for( int i = 0; i < MCMCPar->seq; i++ )
			for( int j = 0; j < 2; j++ )
				rr[i][j] = k++;
		// Define JumpRate -- uniform rand number between 1.2 and 2.2
		Gamma = ( double ) rand() / RAND_MAX + 1.2;
		// Loop over the individual chains
		for( int qq = 0; qq < MCMCPar->seq; qq++ )
		{
			// Define which points of Zoff z_r1, z_r2
			for( int j = 0; j < MCMCPar->n; j++ )
			{
				zR1[j] = Zoff[rr[qq][0]][j];
				zR2[j] = Zoff[rr[qq][1]][j];
			}
			// Now select z from Zoff; z cannot be zR1 and zR2
			for( int i = 0; i < MCMCPar->nzoff; i++ )
				r[i] = i;
			r[rr[qq][0]] = 0;
			r[rr[qq][1]] = 0;
			for( int i = 0; i < MCMCPar->nzoff; i++ )
				if( r[i] > 0 )
					r[i] = i;
			for( int i = 0; i < ( MCMCPar->nzoff - 2 ); ++i ) // random permutation
			{
				int j = rand() % ( i + 1 );
				t[i] = t[j];
				t[j] = i;
			}
			// Define z
			for( int col = 0; col < MCMCPar->n; col++ )
			{
				int sum = t[0] + 2; // TODO why 2?
				sum = r[sum];
				z[qq][col] = Zoff[sum][col];
			}
			// Define projection vector x(qq) - z
			for( int col = 0; col < MCMCPar->n; col++ )
				F[qq][col] = xold[qq][col] - z[qq][col];
			for( int i = 0; i < MCMCPar->n; i++ )
				Fnew[qq][i] = F[qq][i] * F[qq][i];
			double D = 0;
			for( int i = 0; i < MCMCPar->n; i++ )
				D += Fnew[qq][i];
			// Orthogonally project of zR1 and zR2 onto F
			for( int col = 0; col < MCMCPar->n; col++ )
				zP[qq][col] = F[qq][col] * ( zR1[col] - zR2[col] * F[qq][col] ) / D;
			// And define the jump
			for( int i = 0; i < MCMCPar->n; i++ )
				xjump[qq][i] = Gamma * zP[qq][i];
			// Update CR because we only consider full dimensional updates
			CRS[qq][1] = 1;
		}
		Update2 = "Snooker_Update";
	}
	// Now propose new x
	for( int qq = 0; qq < MCMCPar->seq; qq++ )
		for( int j = 0; j < MCMCPar->n; j++ )
			xnew[qq][j] = xold[qq][j] + xjump[qq][j];
	if( MCMCPar->verbose > 5 )
	{
		printf( "xjump offde\n" );
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			for( int i = 0; i < MCMCPar->n; i++ )
				printf( " %g", xjump[row][i] );
			printf( "\n" );
		}
	}
	// Define alpha_s
	if( strcmp( Update.c_str(), "Snooker_Update" ) == 0 )
	{
		double sum, sum1, tmp;
		double power = ( double )( MCMCPar->n - 1 ) / 2;
		// Determine Euclidean distance
		for( int qq = 0; qq < MCMCPar->seq; qq++ )
		{
			sum = sum1 = 0;
			for( int i = 0; i < MCMCPar->n; i++ )
			{
				tmp = xold[qq][i] - z[qq][i];
				sum1 += tmp * tmp;
				tmp = xnew[qq][i] - z[qq][i];
				sum += tmp * tmp;
			}
			alpha_s[qq] = pow( sum / sum1, power );
		}
	}
	else
		for( int i = 0; i < MCMCPar->seq; i++ )
			alpha_s[i] = 1;
	if( MCMCPar->verbose > 8 )
	{
		printf( "xnew pre bound\n" );
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			for( int i = 0; i < MCMCPar->n; i++ )
				printf( " %g", xnew[row][i] );
			printf( "\n" );
		}
	}
	// Do boundary handling -- what to do when points fall outside bound
	if( MCMCPar->BoundHandling == BH_REFLECT )
	{
		if( MCMCPar->verbose ) cout << "ReflectBounds" << endl;
		ReflectBounds( xnew, MCMCPar->seq, MCMCPar->n, ParRange );
	}
	else if( MCMCPar->BoundHandling == BH_BOUND )
	{
		if( MCMCPar->verbose ) cout << "Bound not done" << endl;
		//           SetToBounds(xnew, MCMCPar->seq, MCMCPar->n, ParRange);
	}
	else if( MCMCPar->BoundHandling == BH_FOLD )
	{
		if( MCMCPar->verbose ) cout << "Fold not done" << endl;
		//           FoldBounds(xnew, MCMCPar->seq, MCMCPar->n, ParRange);
	}
	if( MCMCPar->verbose > 8 )
	{
		printf( "xnew post bound\n" );
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			for( int i = 0; i < MCMCPar->n; i++ )
				printf( " %g", xnew[row][i] );
			printf( "\n" );
		}
		cout << "offde done." << endl;
	}
}

// Determine which sequences to evolve with what DE strategy
void DEStrategy( struct MCMC *MCMCPar, int **DEversion )
{
	double p_pair[MCMCPar->DEpairs + 1];
	double Z[MCMCPar->seq];
	int z = 0;
	// Determine probability of selecting a given number of pairs;
	for( int j = 0; j < MCMCPar->DEpairs; j++ )
		p_pair[j] = ( double )( 1.0 / MCMCPar->DEpairs ) * 1;
	for( int j = 0; j < MCMCPar->DEpairs; j++ )
		p_pair[j + 1] = p_pair[j];  //cumsum of p_pair                  ???? need ????
	p_pair[0] = 0;
	// Generate a random number between 0 and 1
	for( int i = 0; i < MCMCPar->seq; i++ )
		Z[i] = ( double ) rand() / RAND_MAX;    //random numbers between 0 and 1
	// Select number of pairs
	for( int qq = 0; qq < MCMCPar->seq; qq++ )
	{
		if( Z[qq] > p_pair[1] ) z = qq;
		else                    z = p_pair[1];
		DEversion[qq][0] = z;
	}
}
// Checks the bounds of the parameters
void ReflectBounds( double **xnew, int nIndivs, int nDim, Range ParRange )
{
	int i, qq;
	// Now check whether points are within bond
	for( i = 0; i < nDim; i++ )
		for( qq = 0; qq < nIndivs; qq++ )
			if( xnew[qq][i] < ParRange.minn[i] ) // TODO this might not be the best apporach
				xnew[qq][i] = 2 * ParRange.minn[i] - xnew[qq][i]; // reflect in min
	// Now check whether points are within bond
	for( i = 0; i < nDim; i++ )
		for( qq = 0; qq < nIndivs; qq++ )
			if( xnew[qq][i] > ParRange.maxn[i] )
				xnew[qq][i] = 2 * ParRange.maxn[i] - xnew[qq][i]; // reflect in max
	// Now double check if all elements are within bounds
	for( i = 0; i < nDim; i++ )
		for( qq = 0; qq < nIndivs; qq++ )
			if( xnew[qq][i] < ParRange.minn[i] )
				xnew[qq][i] = ParRange.minn[i] + ( double ) rand() / RAND_MAX * ( ParRange.maxn[i] - ParRange.minn[i] );
	// Now double check if all elements are within bounds
	for( i = 0; i < nDim; i++ )
		for( qq = 0; qq < nIndivs; qq++ )
			if( xnew[qq][i] > ParRange.maxn[i] )
				xnew[qq][i] = ParRange.minn[i] + ( double ) rand() / RAND_MAX * ( ParRange.maxn[i] - ParRange.minn[i] );
}

// Metropolis rule for acceptance or rejection
void metrop( double **xnew, double *p_xnew, double *log_p_x, double **xold, double *p_xold, double *log_p_xold,
			 double *alpha_s, Measure Measurement, struct MCMC *MCMCPar, double option, double **newgen, int *accept )
{
	// Calculate the number of Chains
	double alpha[MCMCPar->seq];
	double alp = 0;
	double asp = 0;
	// First set newgen to the old positions in X
	for( int row = 0; row < MCMCPar->seq; row++ )
		for( int col = 0; col < MCMCPar->n; col++ )
			newgen[row][col] = xold[row][col];
	for( int row = 0; row < MCMCPar->seq; row++ )
	{
		newgen[row][MCMCPar->n] = p_xold[row];
		newgen[row][MCMCPar->n + 1] = log_p_xold[row];
	}
	// And initialize accept with zeros
	for( int i = 0; i < MCMCPar->seq; i++ )
		accept[i] = 0;
	// -------------------- Now check the various options ----------------------
	if( option == 1 )
		for( int row = 0; row < MCMCPar->seq; row++ )
			alpha[row] = ( p_xnew[row] / p_xold[row] );
	if( option == 2 || option == 4 ) // Lnp probability evaluation
		for( int row = 0; row < MCMCPar->seq; row++ )
			alpha[row] = exp( log_p_x[row] - log_p_xold[row] );
	if( option == 3 ) // SSE probability evaluation
	{
		alp = ( double ) - Measurement.N * ( MCMCPar->gamma + 1 ) / 2;
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			alpha[row] = pow( p_xnew[row] / p_xold[row], alp );
			if( MCMCPar->verbose > 5 ) printf( "alp %g %g %g\n", alpha[row], p_xnew[row], p_xold[row] );
		}
	}
	if( option == 5 ) // Similar as 3 but now weighted with Measurement.Sigma
	{
		for( int row = 0; row < MCMCPar->seq; row++ )
		{
			asp = Measurement.Sigma * Measurement.Sigma;
			alp = ( double ) - 0.5 * ( -p_xnew[row] + p_xold[row] ) / asp;
			alpha[row] = exp( alp );
			// signs are different because we write -SSR
		}
	}
	if( option == 6 )
	{
		alp = ( double ) - Measurement.N * ( MCMCPar->gamma + 1 ) / 2;
		for( int row = 0; row < MCMCPar->seq; row++ )
			alpha[row] = pow( p_xnew[row] / p_xold[row], alp );
	}
	// -------------------------------------------------------------------------
	for( int row = 0; row < MCMCPar->seq; row++ )
		alpha[row] *= alpha_s[row];
	double Z[MCMCPar->seq];
	// Generate random numbers
	for( int row = 0; row < MCMCPar->seq; row++ )
		Z[row] = ( double ) rand() / RAND_MAX;
	// Find which alpha's are greater than Z
	for( int row = 0; row < MCMCPar->seq; row++ )
	{
		if( MCMCPar->verbose > 5 ) printf( "compare %g %g\n", alpha[row], Z[row] );
		if( alpha[row] > Z[row] )
		{
			if( MCMCPar->verbose > 5 ) printf( "accept sequence %d\n", row );
			// And update these chains
			for( int i = 0; i < MCMCPar->n; i++ )
				newgen[row][i] = xnew[row][i];
			newgen[row][MCMCPar->n] = p_xnew[row];
			newgen[row][MCMCPar->n + 1] = log_p_x[row];
			// And indicate that these chains have been accepted
			accept[row] = 1;
		}
		else accept[row] = 0;
	}
}

// Calculate total normalized Euclidean distance for each crossover value
void CalcDelta( struct MCMC *MCMCPar, double *delta_normX, double **CRS, double *delta_tot )
{
	int idx = 0;
	delta_tot[0] = 0;
	double sum = 0;
	// Derive sum_p2 for each different CR value
	for( int zz = 1; zz < MCMCPar->nCR; zz++ )
		for( int col = 0; col < MCMCPar->steps; col++ )
		{
			for( int row = 0; row < MCMCPar->seq; row++ )
			{
				// Find which chains are updated with zz/MCMCPar->nCR
				if( CRS[row][col] - ( double ) zz / MCMCPar->nCR < numeric_limits<double>::epsilon( ) )
					idx = row; //should be row or col not row and col
				sum += delta_normX[idx]; //TODO strange
			}
			// Add the normalized squared distance tot the current delta_tot;
			delta_tot[zz] += sum;
		}
}

// Updates the probabilities of the various crossover values
void AdaptpCR( struct MCMC *MCMCPar, double *delta_tot, double **lCR, double **pCR )
{
	double sum = 0;
	for( int col = 0; col < MCMCPar->nCR; col++ )
		sum += delta_tot[col];
	// Adapt pCR using information from averaged normalized jumping distance
	for( int col = 0; col < MCMCPar->nCR; col++ )
		pCR[0][col] = MCMCPar->seq * ( delta_tot[col] / lCR[0][col] / sum );
	sum = 0;
	// Normalize pCR
	for( int col = 0; col < MCMCPar->nCR; col++ )
		sum += pCR[0][col];
	for( int col = 0; col < MCMCPar->nCR; col++ )
		pCR[0][col] = pCR[0][col] / sum;
}

double **double_matrix( int maxRows, int maxCols )
{
	double **matrix;
	//matrix = new double* [maxRows];
	matrix = (double **) malloc( sizeof( double * ) * maxRows );
	for( int row = 0; row < maxRows; row++ )
	{
		//matrix[row] = new double[maxCols];
		matrix[row] = (double *) malloc( sizeof( double ) * maxCols );
	}
	return( matrix );
}

int **int_matrix( int maxRows, int maxCols )
{
	int **matrix;
	//matrix = new int* [maxRows];
	matrix = (int **) malloc( sizeof( int * ) * maxRows );
	for( int row = 0; row < maxRows; row++ )
	{
		//matrix[row] = new int[maxCols];
		matrix[row] = (int *) malloc( sizeof( int ) * maxCols );
	}
	return( matrix );
}

void free_matrix( void **matrix, int maxCols )
{
	int i;
	for( i = 0; i < maxCols; i++ )
		free( matrix[i] );
	free( matrix );
}

extern "C" void free_mcmc( struct MCMC *mcmc )
{
	free_matrix( (void **) mcmc->z, mcmc->nzee );
	free_matrix( (void **) mcmc->CR, mcmc->DEpairs );
	free( mcmc );
}
