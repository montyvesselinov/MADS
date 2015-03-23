// MADS: Model Analyses & Decision Support (v.1.1.14) 2013
//
// Velimir V Vesselinov (monty), vvv@lanl.gov, velimir.vesselinov@gmail.com
// Dan O'Malley, omalled@lanl.gov
// Dylan Harp, dharp@lanl.gov
//
// http://mads.lanl.gov
// http://www.ees.lanl.gov/staff/monty/codes/mads
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

#include <sys/types.h>  /* include this before any other sys headers */
#include <sys/wait.h>   /* header for waitpid() and various macros */
#include <signal.h>     /* header for signal functions */
#include <stdio.h>      /* header for fprintf() */
#include <unistd.h>     /* header for vfork() */
#include <strings.h>    /* header for strcpy() */
#include <string.h>     /* header for strcpy() */
#include <stdlib.h>     /* header for exit() and malloc() */
#include <time.h>       /* header for time() and ctime() */

#include "../mads.h"

static void handler( int sig );

#define MAXATTEMPTS 5

int nProc, nKids, nHosts, debug;
int *kidstatus;
pid_t *kidids;

/* Functions elsewhere */
char **char_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
char *timestamp();
int func_extrn_check_read( int ieval, void *data );

int mprunwrite( int nJob, void *data, double *var_mat[], double *phi, double *f[] )
{
	struct opt_data *p = ( struct opt_data * )data;
	struct sigaction act;
	int w, i, j, ieval, type, cJob, nFailed, child, child1, wait, job_wait, done, next, refork = 0, refresh, destroy, rerun, rJob, *kidattempt, *skip_job;
	pid_t pid, return_fork;
	char *exec_name, **kidhost, **kiddir, **rerundir, dir[1025], buf[1025], *atime;
	double *opt_params;
	if( ( opt_params = ( double * ) malloc( p->pd->nOptParam * sizeof( double ) ) ) == NULL ) { tprintf( "Not enough memory!\n" ); return( 0 ); }
	if( p->cd->num_proc <= 1 ) { tprintf( "\nERROR: Number of available processors is 1; cannot parallelize!\n" ); return( -1 ); }
	if( p->cd->pardebug > 3 )
	{
		if( p->cd->paral_hosts == NULL ) tprintf( "\nWARNING: Local runs using %d processors! No parallel hosts!\n", p->cd->num_proc );
		else tprintf( "Parallel runs using %d hosts!\n", p->cd->num_proc );
	}
	type = p->cd->parallel_type;
	nProc = nHosts = p->cd->num_proc; // Number of processors/hosts available initially
	exec_name = p->ed->cmdline; // Executable / Execution command line
	ieval = p->cd->neval; // Current number of model evaluations
	if( p->cd->pardebug )
	{
		if( nJob > 1 ) tprintf( "Parallel writing of %d jobs using %d processors ... ", nJob, nProc );
		else           tprintf( "Parallel writing of 1 job ... " );
		tprintf( "\n" );
	}
	if( ( skip_job = ( int * ) malloc( nJob * sizeof( int ) ) ) == NULL ) tprintf( "ERROR!!!\n" );
	for( i = 0; i < nJob; i++ )
		skip_job[i] = 0;
	// for( w = 0; w < nJob; w++ )
	// tprintf( "%d: init skip_job %d\n", w, skip_job[w] );
	if( p->cd->restart ) // Check for already computed jobs (smart restart)
	{
		done = 0;
		for( i = 0; i < nJob; i++ )
		{
			done += skip_job[i] = func_extrn_check_read( ieval + i + 1, p );
			if( p->cd->pardebug > 1 )
			{
				if( skip_job[i] == 1 ) tprintf( "Job %d is already completed; it will be skipped!\n", ieval + i + 1 );
				else tprintf( "Job %d will be executed!\n", ieval + i + 1 );
			}
		}
		if( done == nJob ) // All the jobs will be skipped
		{
			free( skip_job );
			if( p->cd->pardebug ) tprintf( "Restart: All %d jobs are already completed!\n", nJob );
			return( 1 );
		}
		else if( done > 0 )
			tprintf( "WARNING Restart: %d jobs out of %d will be skipped because it appears to be already completed!\n", done, nJob );
	}
	debug = ( p->cd->pardebug > 3 ) ? 1 : 0; // Debug level
	kidids = ( pid_t * ) malloc( nProc * sizeof( pid_t ) ); // ID's of external jobs
	kidstatus = ( int * ) malloc( nProc * sizeof( int ) );  // Status of external jobs
	kidattempt = ( int * ) malloc( nProc * sizeof( int ) ); // Number of attempts to execute each external job
	kiddir = char_matrix( nProc, 1025 ); // Directories for external jobs
	rerundir = char_matrix( nProc, 1025 ); // Rerun directories for external jobs
	for( w = 0; w < nProc; w++ )
	{
		kidids[w] = 0;
		kidstatus[w] = kidattempt[w] = -1;
	}
	if( type == 0 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ ) strcpy( kidhost[i], "local" );
	}
	else if( type == 1 )
		kidhost = p->cd->paral_hosts; // List of processors/hosts
	else if( type == 2 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ )
			strcpy( kidhost[i], "slurm" );
	}
	/*
	for( w = 0; w < nProc; w++ )
		tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
	for( w = 0; w < nJob; w++ )
		tprintf( "skip_job %d\n", skip_job[w] );
	*/
	act.sa_handler = handler; // POSIX process handler
	sigemptyset( &act.sa_mask );
	act.sa_flags = 0;
	if( sigaction( SIGCHLD, &act, NULL ) < 0 )
	{
		tprintf( "sigaction failed!!!\n" );
		free( skip_job );
		return( -1 );
	}
	nFailed = 0; nKids = 0; cJob = 0; rJob = 0; wait = 0; done = 0;
	while( 1 ) // Main loop
	{
//		if( rJob >= nJob || nProc <= 0 )
		if( rJob > nJob || nProc <= 0 )
		{
			tprintf( "ERROR: None of the processors is responding properly! Parallel execution fails!\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
			act.sa_handler = SIG_DFL; // Default handler
			sigaction( SIGCHLD, &act, NULL );
			free( opt_params );
			free( ( void * ) kidids );
			free( ( void * ) kidstatus );
			free( ( void * ) kidattempt );
			free( skip_job );
			free_matrix( ( void ** ) kiddir, nProc );
			free_matrix( ( void ** ) rerundir, nProc );
			if( type != 1 ) free_matrix( ( void ** ) kidhost, nProc );
			return( -1 );
		}
		job_wait = 1;
		/*
		for( w = 0; w < nProc; w++ )
			tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
		for( w = 0; w < nJob; w++ )
			tprintf( "skip_job %d\n", skip_job[w] );
		*/
		if( rJob > 0 && nKids < ( nHosts - nFailed ) )
		{
			if( p->cd->pardebug > 1 )
			{
				for( i = 0; i < nHosts; i++ )
					tprintf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				tprintf( "Free processor to rerun: " );
			}
			for( j = 0; j < nHosts; j++ )
				if( kidattempt[j] <= 0 )
				{
					j++;
					if( p->cd->pardebug > 1 ) tprintf( "Processor %d\n", j );
					job_wait = 0;
				}
			if( p->cd->pardebug > 1 ) tprintf( "NONE! %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
		}
		if( job_wait )
		{
			// double check
			j = nKids;
			w = 0;
			for( i = 0; i < nHosts; i++ )
			{
				if( kidstatus[i] == 1 && kidids[i] > 0 )
				{
					if( p->cd->pardebug > 1 ) tprintf( "Processor %i [%d] : status = %d %d", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
					if( kidids[i] > 0 && kill( kidids[i], 0 ) == 0 )
					{
						if( p->cd->pardebug > 1 ) tprintf( " is running!\n" );
						w++;
					}
					else if( errno == ESRCH )
					{
						if( p->cd->pardebug > 1 ) tprintf( " does not exist!\n" );
						tprintf( "ERROR: Processor %i [%d] : status = %d %d does not exist\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
						kidstatus[i] = -1;
						kidids[i] = 0;
						j--;
					}
					else
					{
						if( p->cd->pardebug > 1 ) tprintf( " there is a problem to determine the status!\n" );
					}
				}
			}
			if( j < nKids ) { tprintf( "ERROR: %d missing jobs!\n", nKids - j ); }
			if( w < nKids ) { tprintf( "ERROR: job count mismatch (%d < %d); number of running jobs is %d!\n", w, nKids, w ); nKids = w; }
			if( !done )
			{
				if( nKids >= nProc )
				{
					wait = 1;
					if( p->cd->pardebug > 2 ) tprintf( "Waiting ... %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
					else if( p->cd->pardebug > 1 ) tprintf( "Waiting ...\n" );
					sigsuspend( &act.sa_mask );
				}
			}
			else
			{
				if( wait != 2 && p->cd->pardebug > 1 ) tprintf( "All the jobs are started!\n" );
				int nKidsnewcount = 0; for( j = 0; j < nHosts; j++ ) if( kidstatus[j] == 1 && kidattempt[j] > 0 ) nKidsnewcount++;
				if( nKidsnewcount != nKids ) { tprintf( "Kids count mismatch (%d != %d)\n", nKidsnewcount, nKids ); nKids = nKidsnewcount; }
				if( nKids > 0 )
				{
					wait = 2;
					if( p->cd->pardebug > 2 ) tprintf( "Waiting ... %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
					else if( p->cd->pardebug > 1 ) tprintf( "Waiting ... (%d)\n", nKids );
					sigsuspend( &act.sa_mask );
				}
				else
				{
					if( p->cd->pardebug > 2 ) tprintf( "All the jobs are done!\n" );
					break;
				}
			}
			for( j = 0; j < nHosts; ) if( kidstatus[j++] != 1 ) break; // find a kid with status != 1
			if( j > nHosts )
			{
				if( p->cd->pardebug > 1 )
				{
					tprintf( "All the processors are busy! kids = %d\n", nKids );
					for( i = 0; i < nHosts; i++ )
						tprintf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				}
				continue;
			}
		}
		child1 = j; // available kid with status != 1
		child  = j - 1;
		atime = timestamp();
		if( p->cd->pardebug > 1 ) tprintf( "Processor %i [%s:%d] %s : ", child1, kidhost[child], kidids[child], atime );
		destroy = refresh = 0;
		if( rJob > 0 ) { rerun = 1; next = 0; }
		else           { rerun = 0; next = 1; }
		if( kidstatus[child] == 0 ) // kid suspended
		{
			if( done && !rerun ) { if( p->cd->pardebug > 1 ) tprintf( "suspended; will be killed! " ); kiddir[child][0] = 0; destroy = 1; }
			else { if( p->cd->pardebug > 1 ) tprintf( "suspended; new job " ); refork = 0; refresh = 1; }
		}
		else if( kidstatus[child] == -1 ) // kid finished
		{
			if( done && !rerun )
			{
				kiddir[child][0] = 0;
				kidhost[child][0] = 0;
				kidstatus[child] = 1;
				kidattempt[child] = 0;
				kidids[child] = 0;
				if( p->cd->pardebug > 1 ) tprintf( "finished!\n" );
				continue;
			}
			else
			{
				nKids++;
				refork = 1;  // refork; do not refresh
				refresh = 0;
				kidids[child] = 0;
				if( kidattempt[child] == -1 ) { if( p->cd->pardebug > 1 ) tprintf( "Initializing " ); }
				else                          { if( p->cd->pardebug > 1 ) tprintf( "finished; Starting " ); }
				kidattempt[child] = 1;
			}
		}
		else if( kidstatus[child] == -2 ) // kid killed
		{
			if( kiddir[child][0] != 0 )
			{
				if( p->cd->pardebug > 1 ) tprintf( "killed; " );
				if( kidattempt[child] >= MAXATTEMPTS ) // ignore the processor
				{
					nFailed++;
					nProc--;
					if( p->cd->pardebug > 1 ) tprintf( "(child %d on %s ignored)\n", child1, kidhost[child] );
					strcpy( rerundir[rJob], kiddir[child] );
					rJob++;
					kiddir[child][0] = 0;
					kidhost[child][0] = 0;
					kidstatus[child] = 1;
					kidids[child] = 0;
					tprintf( "WARNING: The number of currently used processors is decreased: %d of %d!\n", nProc, nHosts );
					continue;
				}
				else // try to restart
				{
					nKids++;
					kidattempt[child]++;
					if( p->cd->pardebug > 1 ) tprintf( "restart (attempt %d of %d)", kidattempt[child], MAXATTEMPTS );
					strcpy( dir, kiddir[child] );
					rerun = 0;
					next = 0;
					refork = 1;
					refresh = 0;
					kidids[child] = 0;
					sleep( kidattempt[child] );
				}
			}
			else { if( p->cd->pardebug > 1 ) tprintf( "killed internally!\n" ); continue; }
		}
		if( destroy && kidids[child] > 0 )
		{
			if( kidids[child] > 0 ) kill( kidids[child] * -1, SIGCONT );
			else tprintf( "ERROR: something is wrong destroying process %d\n!", child1 );
			if( p->cd->pardebug > 1 ) tprintf( "\n" );
			continue;
		}
		if( rerun )
		{
			next = 0;
			rJob--;
			strcpy( dir, rerundir[rJob] );
		}
		else if( next ) // Next job
		{
			cJob++;
			if( cJob >= nJob ) done = 1;
			if( cJob > nJob ) continue;
			if( p->cd->pardebug ) tprintf( "Job %d", cJob );
			if( skip_job[cJob - 1] == 1 )
			{
				nKids--;
				if( p->cd->pardebug ) tprintf( " skipped because it is already completed!\n" );
				continue;
			}
			else
			{
				sprintf( dir, "../%s_%08d", p->cd->mydir_hosts, ieval + cJob ); // Name of directory for parallel runs
				strcpy( kiddir[child], dir );
			}
		}
		if( p->cd->pardebug > 1 ) tprintf( " : writing in \'%s\'\n", dir );
		else if( p->cd->pardebug ) { if( nJob > 1 ) tprintf( " ...\n" ); else tprintf( " ... " ); }
		if( refork )
		{
			if( kidhost[child][0] == 0 ) continue;
			if( cJob > nJob ) continue;
			if( ( return_fork = fork() ) == 0 )
			{
				/*
				for( w = 0; w < nProc; w++ )
					tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
				for( w = 0; w < nJob; w++ )
					tprintf( "skip_job %d\n", skip_job[w] );
				*/
				// TODO sleep below can be important; try to run without; if it fails, comment it out
				struct timespec tspec, tspec2;
				tspec.tv_sec = 0;
				tspec.tv_nsec = 500L;
				nanosleep( &tspec, &tspec2 );
				pid = getpid();
				setpgid( pid, pid );
				if( p->cd->pardebug > 3 ) tprintf( "Forked Process %i [%s:%d] : writing in \'%s\'\n", child1, kidhost[child], pid, dir );
				if( p->cd->debug || p->cd->mdebug || p->cd->pardebug > 3 ) tprintf( "Parallel writing of the model output files for case %d ...\n", ieval + cJob );
				int ii;
				for( ii = 0; ii < p->pd->nOptParam; ii++ )
				{
					int kk = p->pd->var_index[ii];
					opt_params[ii] = p->pd->var[kk] = var_mat[cJob - 1][ii];
					tprintf( "%s %.12g\n", p->pd->var_name[kk], p->cd->var[kk] );
				}
				func_extrn_write( ieval + cJob, opt_params, p );
				// execlp( "/usr/bin/env", "/usr/bin/env", "tcsh", "-f", "-c", "sleep 4", ( char * ) 0 );
				// sleep( 1 );
				_exit( 7 );
			}
			if( return_fork > 0 )
			{
				kidids[child] = return_fork;
				kidstatus[child] = 1;
				strcpy( kiddir[child], dir );
			}
			else
			{
				nKids--;
				tprintf( "WARNING: fork failed!!!\n" );
			}
		}
		else if( refresh )
		{
			strcpy( kiddir[child], dir );
			tprintf( "Refreshed Process %i [%s:%d] : writing in \'%s\'\n", child1, kidhost[child], kidids[child], kiddir[child] );
			sleep( 1 );
			if( kidids[child] > 0 ) kill( kidids[child] * -1, SIGCONT );
			else tprintf( "ERROR: something is wrong destroying process\n!" );
		}
	}
	/*
	for( w = 0; w < nProc; w++ )
		tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
	for( w = 0; w < nJob; w++ )
		tprintf( "skip_job %d\n", skip_job[w] );
	*/
	tprintf( "Done:\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
	act.sa_handler = SIG_DFL; // Default handler
	sigaction( SIGCHLD, &act, NULL );
	if( opt_params != NULL ) free( opt_params );
	if( kidids != NULL ) free( ( void * ) kidids );
	if( kidstatus != NULL ) free( ( void * ) kidstatus );
	if( kidattempt != NULL ) free( ( void * ) kidattempt );
	if( skip_job != NULL ) free( skip_job );
	if( kiddir != NULL ) free_matrix( ( void ** ) kiddir, nProc );
	if( rerundir != NULL ) free_matrix( ( void ** ) rerundir, nProc );
	if( kidhost != NULL && type != 1 ) free_matrix( ( void ** ) kidhost, nProc );
	return( 1 );
}

int mprunread( int nJob, void *data, double *var_mat[], double *phi, double *f[] )
{
	struct opt_data *p = ( struct opt_data * )data;
	struct sigaction act;
	int w, i, j, ieval, type, cJob, nFailed, child, child1, wait, job_wait, done, next, refork = 0, refresh, destroy, rerun, rJob, *kidattempt, *skip_job;
	pid_t pid, return_fork;
	char *exec_name, **kidhost, **kiddir, **rerundir, dir[1025], buf[1025], *atime;
	if( p->cd->num_proc <= 1 ) { tprintf( "\nERROR: Number of available processors is 1; cannot parallelize!\n" ); return( -1 ); }
	if( p->cd->pardebug > 3 )
	{
		if( p->cd->paral_hosts == NULL ) tprintf( "\nWARNING: Local runs using %d processors! No parallel hosts!\n", p->cd->num_proc );
		else tprintf( "Parallel runs using %d hosts!\n", p->cd->num_proc );
	}
	type = p->cd->parallel_type;
	nProc = nHosts = p->cd->num_proc; // Number of processors/hosts available initially
	exec_name = p->ed->cmdline; // Executable / Execution command line
	ieval = p->cd->neval - nJob; // Current number of model evaluations
	if( p->cd->pardebug )
	{
		if( nJob > 1 ) tprintf( "Parallel reading of %d jobs using %d processors ... ", nJob, nProc );
		else           tprintf( "Parallel reading of 1 job ... " );
		tprintf( "\n" );
	}
	if( ( skip_job = ( int * ) malloc( nJob * sizeof( int ) ) ) == NULL ) tprintf( "ERROR!!!\n" );
	for( i = 0; i < nJob; i++ )
		skip_job[i] = 0;
	// for( w = 0; w < nJob; w++ )
	// tprintf( "%d: init skip_job %d\n", w, skip_job[w] );
	if( p->cd->restart ) // Check for already computed jobs (smart restart)
	{
		done = 0;
		for( i = 0; i < nJob; i++ )
		{
			done += skip_job[i] = func_extrn_check_read( ieval + i + 1, p );
			if( p->cd->pardebug > 1 )
			{
				if( skip_job[i] == 1 ) tprintf( "Job %d is already completed; it will be skipped!\n", ieval + i + 1 );
				else tprintf( "Job %d will be executed!\n", ieval + i + 1 );
			}
		}
		if( done == nJob ) // All the jobs will be skipped
		{
			free( skip_job );
			if( p->cd->pardebug ) tprintf( "Restart: All %d jobs are already completed!\n", nJob );
			return( 1 );
		}
		else if( done > 0 )
			tprintf( "WARNING Restart: %d jobs out of %d will be skipped because it appears to be already completed!\n", done, nJob );
	}
	debug = ( p->cd->pardebug > 3 ) ? 1 : 0; // Debug level
	kidids = ( pid_t * ) malloc( nProc * sizeof( pid_t ) ); // ID's of external jobs
	kidstatus = ( int * ) malloc( nProc * sizeof( int ) );  // Status of external jobs
	kidattempt = ( int * ) malloc( nProc * sizeof( int ) ); // Number of attempts to execute each external job
	kiddir = char_matrix( nProc, 1025 ); // Directories for external jobs
	rerundir = char_matrix( nProc, 1025 ); // Rerun directories for external jobs
	for( w = 0; w < nProc; w++ )
	{
		kidids[w] = 0;
		kidstatus[w] = kidattempt[w] = -1;
	}
	if( type == 0 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ ) strcpy( kidhost[i], "local" );
	}
	else if( type == 1 )
		kidhost = p->cd->paral_hosts; // List of processors/hosts
	else if( type == 2 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ )
			strcpy( kidhost[i], "slurm" );
	}
	/*
	for( w = 0; w < nProc; w++ )
		tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
	for( w = 0; w < nJob; w++ )
		tprintf( "skip_job %d\n", skip_job[w] );
	*/
	act.sa_handler = handler; // POSIX process handler
	sigemptyset( &act.sa_mask );
	act.sa_flags = 0;
	if( sigaction( SIGCHLD, &act, NULL ) < 0 )
	{
		tprintf( "sigaction failed!!!\n" );
		free( skip_job );
		return( -1 );
	}
	nFailed = 0; nKids = 0; cJob = 0; rJob = 0; wait = 0; done = 0;
	while( 1 ) // Main loop
	{
//		if( rJob >= nJob || nProc <= 0 )
		if( rJob > nJob || nProc <= 0 )
		{
			tprintf( "ERROR: None of the processors is responding properly! Parallel execution fails!\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
			act.sa_handler = SIG_DFL; // Default handler
			sigaction( SIGCHLD, &act, NULL );
			free( ( void * ) kidids ); free( ( void * ) kidstatus ); free( ( void * ) kidattempt );
			free( skip_job );
			free_matrix( ( void ** ) kiddir, nProc );
			free_matrix( ( void ** ) rerundir, nProc );
			if( type != 1 ) free_matrix( ( void ** ) kidhost, nProc );
			return( -1 );
		}
		job_wait = 1;
		/*
		for( w = 0; w < nProc; w++ )
			tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
		for( w = 0; w < nJob; w++ )
			tprintf( "skip_job %d\n", skip_job[w] );
		*/
		if( rJob > 0 && nKids < ( nHosts - nFailed ) )
		{
			if( p->cd->pardebug > 1 )
			{
				for( i = 0; i < nHosts; i++ )
					tprintf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				tprintf( "Free processor to rerun: " );
			}
			for( j = 0; j < nHosts; j++ )
				if( kidattempt[j] <= 0 )
				{
					j++;
					if( p->cd->pardebug > 1 ) tprintf( "Processor %d\n", j );
					job_wait = 0;
				}
			if( p->cd->pardebug > 1 ) tprintf( "NONE! %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
		}
		if( job_wait )
		{
			// double check
			j = nKids;
			w = 0;
			for( i = 0; i < nHosts; i++ )
			{
				if( kidstatus[i] == 1 && kidids[i] > 0 )
				{
					if( p->cd->pardebug > 1 ) tprintf( "Processor %i [%d] : status = %d %d", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
					if( kidids[i] > 0 && kill( kidids[i], 0 ) == 0 )
					{
						if( p->cd->pardebug > 1 ) tprintf( " is running!\n" );
						w++;
					}
					else if( errno == ESRCH )
					{
						if( p->cd->pardebug > 1 ) tprintf( " does not exist!\n" );
						tprintf( "ERROR: Processor %i [%d] : status = %d %d does not exist\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
						kidstatus[i] = -1;
						kidids[i] = 0;
						j--;
					}
					else
					{
						if( p->cd->pardebug > 1 ) tprintf( " there is a problem to determine the status!\n" );
					}
				}
			}
			if( j < nKids ) { tprintf( "ERROR: %d missing jobs!\n", nKids - j ); }
			if( w < nKids ) { tprintf( "ERROR: job count mismatch (%d < %d); number of running jobs is %d!\n", w, nKids, w ); nKids = w; }
			if( !done )
			{
				if( nKids >= nProc )
				{
					wait = 1;
					if( p->cd->pardebug > 2 ) tprintf( "Waiting ... %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
					else if( p->cd->pardebug > 1 ) tprintf( "Waiting ...\n" );
					sigsuspend( &act.sa_mask );
				}
			}
			else
			{
				if( wait != 2 && p->cd->pardebug > 1 ) tprintf( "All the jobs are started!\n" );
				int nKidsnewcount = 0; for( j = 0; j < nHosts; j++ ) if( kidstatus[j] == 1 && kidattempt[j] > 0 ) nKidsnewcount++;
				if( nKidsnewcount != nKids ) { tprintf( "Kids count mismatch (%d != %d)\n", nKidsnewcount, nKids ); nKids = nKidsnewcount; }
				if( nKids > 0 )
				{
					wait = 2;
					if( p->cd->pardebug > 2 ) tprintf( "Waiting ... %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
					else if( p->cd->pardebug > 1 ) tprintf( "Waiting ... (%d)\n", nKids );
					sigsuspend( &act.sa_mask );
				}
				else
				{
					if( p->cd->pardebug > 2 ) tprintf( "All the jobs are done!\n" );
					break;
				}
			}
			for( j = 0; j < nHosts; ) if( kidstatus[j++] != 1 ) break; // find a kid with status != 1
			if( j > nHosts )
			{
				if( p->cd->pardebug > 1 )
				{
					tprintf( "All the processors are busy! kids = %d\n", nKids );
					for( i = 0; i < nHosts; i++ )
						tprintf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				}
				continue;
			}
		}
		child1 = j; // available kid with status != 1
		child  = j - 1;
		atime = timestamp();
		if( p->cd->pardebug > 1 ) tprintf( "Processor %i [%s:%d] %s : ", child1, kidhost[child], kidids[child], atime );
		destroy = refresh = 0;
		if( rJob > 0 ) { rerun = 1; next = 0; }
		else           { rerun = 0; next = 1; }
		if( kidstatus[child] == 0 ) // kid suspended
		{
			if( done && !rerun ) { if( p->cd->pardebug > 1 ) tprintf( "suspended; will be killed! " ); kiddir[child][0] = 0; destroy = 1; }
			else { if( p->cd->pardebug > 1 ) tprintf( "suspended; new job " ); refork = 0; refresh = 1; }
		}
		else if( kidstatus[child] == -1 ) // kid finished
		{
			if( done && !rerun )
			{
				kiddir[child][0] = 0;
				kidhost[child][0] = 0;
				kidstatus[child] = 1;
				kidattempt[child] = 0;
				kidids[child] = 0;
				if( p->cd->pardebug > 1 ) tprintf( "finished!\n" );
				continue;
			}
			else
			{
				nKids++;
				refork = 1;  // refork; do not refresh
				refresh = 0;
				kidids[child] = 0;
				if( kidattempt[child] == -1 ) { if( p->cd->pardebug > 1 ) tprintf( "Initializing " ); }
				else                          { if( p->cd->pardebug > 1 ) tprintf( "finished; Starting " ); }
				kidattempt[child] = 1;
			}
		}
		else if( kidstatus[child] == -2 ) // kid killed
		{
			if( kiddir[child][0] != 0 )
			{
				if( p->cd->pardebug > 1 ) tprintf( "killed; " );
				if( kidattempt[child] >= MAXATTEMPTS ) // ignore the processor
				{
					nFailed++;
					nProc--;
					if( p->cd->pardebug > 1 ) tprintf( "(child %d on %s ignored)\n", child1, kidhost[child] );
					strcpy( rerundir[rJob], kiddir[child] );
					rJob++;
					kiddir[child][0] = 0;
					kidhost[child][0] = 0;
					kidstatus[child] = 1;
					kidids[child] = 0;
					tprintf( "WARNING: The number of currently used processors is decreased: %d of %d!\n", nProc, nHosts );
					continue;
				}
				else // try to restart
				{
					nKids++;
					kidattempt[child]++;
					if( p->cd->pardebug > 1 ) tprintf( "restart (attempt %d of %d)", kidattempt[child], MAXATTEMPTS );
					strcpy( dir, kiddir[child] );
					rerun = 0;
					next = 0;
					refork = 1;
					refresh = 0;
					kidids[child] = 0;
					sleep( kidattempt[child] );
				}
			}
			else { if( p->cd->pardebug > 1 ) tprintf( "killed internally!\n" ); continue; }
		}
		if( destroy && kidids[child] > 0 )
		{
			if( kidids[child] > 0 ) kill( kidids[child] * -1, SIGCONT );
			else tprintf( "ERROR: something is wrong destroying process %d\n!", child1 );
			if( p->cd->pardebug > 1 ) tprintf( "\n" );
			continue;
		}
		if( rerun )
		{
			next = 0;
			rJob--;
			strcpy( dir, rerundir[rJob] );
		}
		else if( next ) // Next job
		{
			cJob++;
			if( cJob >= nJob ) done = 1;
			if( cJob > nJob ) continue;
			if( p->cd->pardebug ) tprintf( "Job %d", cJob );
			if( skip_job[cJob - 1] == 1 )
			{
				nKids--;
				if( p->cd->pardebug ) tprintf( " skipped because it is already completed!\n" );
				continue;
			}
			else
			{
				sprintf( dir, "../%s_%08d", p->cd->mydir_hosts, ieval + cJob ); // Name of directory for parallel runs
				strcpy( kiddir[child], dir );
			}
		}
		if( p->cd->pardebug > 1 ) tprintf( " : reading in \'%s\'\n", dir );
		else if( p->cd->pardebug ) { if( nJob > 1 ) tprintf( " ...\n" ); else tprintf( " ... " ); }
		if( refork )
		{
			if( kidhost[child][0] == 0 ) continue;
			if( cJob > nJob ) continue;
			ieval++;
			if( ( return_fork = fork() ) == 0 )
			{
				/*
				for( w = 0; w < nProc; w++ )
					tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
				for( w = 0; w < nJob; w++ )
					tprintf( "skip_job %d\n", skip_job[w] );
				*/
				// TODO sleep below can be important; try to run without; if it fails, comment it out
				struct timespec tspec, tspec2;
				tspec.tv_sec = 0;
				tspec.tv_nsec = 500L;
				nanosleep( &tspec, &tspec2 );
				pid = getpid();
				setpgid( pid, pid );
				if( p->cd->pardebug > 3 ) tprintf( "Forked Process %i [%s:%d] : reading in \'%s\'\n", child1, kidhost[child], pid, dir );
				if( p->cd->debug || p->cd->mdebug ) tprintf( "Parallel reading the model output files for case %d ... \n", ieval );
				if( p->cd->pardebug > 3 )
				{
					tprintf( "iteration %d : ", ieval ); // counter
					int kk, ii;
					for( ii = 0; ii < p->pd->nOptParam; ii++ ) // re
					{
						kk = p->pd->var_index[ii];
						tprintf( "%s\n", p->pd->var_name[kk] );
						tprintf( "%s %.12g\n", p->pd->var_name[kk], var_mat[cJob - 1][ii] );
					}
				}
				if( func_extrn_read( ieval, p, p->od->res ) ) return( 0 );
				if( p->cd->debug || p->cd->mdebug ) tprintf( "Parallel reading OF = %g for case %d ... \n", p->phi, ieval );
				if( phi != NULL ) phi[cJob - 1] = p->phi;
				if( f != NULL )
					for( j = 0; j < p->od->nTObs; j++ )
						f[cJob - 1][j] = p->od->res[j];
				// sleep( 1 );
				// execlp( "/usr/bin/env", "/usr/bin/env", "tcsh", "-f", "-c", "sleep 4", ( char * ) 0 );
				_exit( 7 );
			}
			if( return_fork > 0 )
			{
				kidids[child] = return_fork;
				kidstatus[child] = 1;
				strcpy( kiddir[child], dir );
			}
			else
			{
				nKids--;
				tprintf( "WARNING: fork failed!!!\n" );
			}
		}
		else if( refresh )
		{
			strcpy( kiddir[child], dir );
			tprintf( "Refreshed Process %i [%s:%d] : reading in \'%s\'\n", child1, kidhost[child], kidids[child], kiddir[child] );
			sleep( 1 );
			if( kidids[child] > 0 ) kill( kidids[child] * -1, SIGCONT );
			else tprintf( "ERROR: something is wrong destroying process\n!" );
		}
	}
	/*
	for( w = 0; w < nProc; w++ )
		tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
	for( w = 0; w < nJob; w++ )
		tprintf( "skip_job %d\n", skip_job[w] );
	*/
	tprintf( "Done:\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
	act.sa_handler = SIG_DFL; // Default handler
	sigaction( SIGCHLD, &act, NULL );
	free( ( void * ) kidids );
	free( ( void * ) kidstatus );
	free( ( void * ) kidattempt );
	free( skip_job );
	free_matrix( ( void ** ) kiddir, nProc );
	free_matrix( ( void ** ) rerundir, nProc );
	if( type != 1 ) free_matrix( ( void ** ) kidhost, nProc );
	return( 1 );
}


int mprun( int nJob, void *data )
{
	struct opt_data *p = ( struct opt_data * )data;
	struct sigaction act;
	int w, i, j, ieval, type, cJob, nFailed, child, child1, wait, job_wait, done, next, refork = 0, refresh, destroy, rerun, rJob, *kidattempt, *skip_job;
	pid_t pid, return_fork;
	char *exec_name, **kidhost, **kiddir, **rerundir, dir[1025], buf[1025], *atime;
	if( p->cd->num_proc <= 1 ) { tprintf( "\nERROR: Number of available processors is 1; cannot parallelize!\n" ); return( -1 ); }
	if( p->cd->pardebug > 3 )
	{
		if( p->cd->paral_hosts == NULL ) tprintf( "\nWARNING: Local runs using %d processors! No parallel hosts!\n", p->cd->num_proc );
		else tprintf( "Parallel runs using %d hosts!\n", p->cd->num_proc );
	}
	type = p->cd->parallel_type;
	nProc = nHosts = p->cd->num_proc; // Number of processors/hosts available initially
	exec_name = p->ed->cmdline; // Executable / Execution command line
	ieval = p->cd->neval; // Current number of model evaluations
	if( p->cd->pardebug )
	{
		if( nJob > 1 ) tprintf( "Parallel execution of %d jobs using %d processors ... ", nJob, nProc );
		else           tprintf( "Parallel execution of 1 job ... " );
		tprintf( "\n" );
	}
	if( ( skip_job = ( int * ) malloc( nJob * sizeof( int ) ) ) == NULL ) tprintf( "ERROR!!!\n" );
	for( i = 0; i < nJob; i++ )
		skip_job[i] = 0;
	// for( w = 0; w < nJob; w++ )
	// tprintf( "%d: init skip_job %d\n", w, skip_job[w] );
	if( p->cd->restart ) // Check for already computed jobs (smart restart)
	{
		done = 0;
		for( i = 0; i < nJob; i++ )
		{
			done += skip_job[i] = func_extrn_check_read( ieval + i + 1, p );
			if( p->cd->pardebug > 1 )
			{
				if( skip_job[i] == 1 ) tprintf( "Job %d is already completed; it will be skipped!\n", ieval + i + 1 );
				else tprintf( "Job %d will be executed!\n", ieval + i + 1 );
			}
		}
		if( done == nJob ) // All the jobs will be skipped
		{
			p->cd->neval += nJob;
			free( skip_job );
			if( p->cd->pardebug ) tprintf( "Restart: All %d jobs are already completed!\n", nJob );
			return( 1 );
		}
		else if( done > 0 )
			tprintf( "WARNING Restart: %d jobs out of %d will be skipped because it appears to be already completed!\n", done, nJob );
	}
	debug = ( p->cd->pardebug > 3 ) ? 1 : 0; // Debug level
	kidids = ( pid_t * ) malloc( nProc * sizeof( pid_t ) ); // ID's of external jobs
	kidstatus = ( int * ) malloc( nProc * sizeof( int ) );  // Status of external jobs
	kidattempt = ( int * ) malloc( nProc * sizeof( int ) ); // Number of attempts to execute each external job
	kiddir = char_matrix( nProc, 1025 ); // Directories for external jobs
	rerundir = char_matrix( nProc, 1025 ); // Rerun directories for external jobs
	for( w = 0; w < nProc; w++ )
	{
		kidids[w] = 0;
		kidstatus[w] = kidattempt[w] = -1;
	}
	if( type == 0 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ ) strcpy( kidhost[i], "local" );
	}
	else if( type == 1 )
		kidhost = p->cd->paral_hosts; // List of processors/hosts
	else if( type == 2 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ )
			strcpy( kidhost[i], "slurm" );
	}
	/*
	for( w = 0; w < nProc; w++ )
		tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
	for( w = 0; w < nJob; w++ )
		tprintf( "skip_job %d\n", skip_job[w] );
	*/
	act.sa_handler = handler; // POSIX process handler
	sigemptyset( &act.sa_mask );
	act.sa_flags = 0;
	if( sigaction( SIGCHLD, &act, NULL ) < 0 )
	{
		tprintf( "sigaction failed!!!\n" );
		free( skip_job );
		return( -1 );
	}
	nFailed = 0; nKids = 0; cJob = 0; rJob = 0; wait = 0; done = 0;
	while( 1 ) // Main loop
	{
//		if( rJob >= nJob || nProc <= 0 )
		if( rJob > nJob || nProc <= 0 )
		{
			tprintf( "ERROR: None of the processors is responding properly! Parallel execution fails!\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
			act.sa_handler = SIG_DFL; // Default handler
			sigaction( SIGCHLD, &act, NULL );
			free( ( void * ) kidids ); free( ( void * ) kidstatus ); free( ( void * ) kidattempt );
			free( skip_job );
			free_matrix( ( void ** ) kiddir, nProc );
			free_matrix( ( void ** ) rerundir, nProc );
			if( type != 1 ) free_matrix( ( void ** ) kidhost, nProc );
			return( -1 );
		}
		job_wait = 1;
		/*
		for( w = 0; w < nProc; w++ )
			tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
		for( w = 0; w < nJob; w++ )
			tprintf( "skip_job %d\n", skip_job[w] );
		*/
		if( rJob > 0 && nKids < ( nHosts - nFailed ) )
		{
			if( p->cd->pardebug > 1 )
			{
				for( i = 0; i < nHosts; i++ )
					tprintf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				tprintf( "Free processor to rerun: " );
			}
			for( j = 0; j < nHosts; j++ )
				if( kidattempt[j] <= 0 )
				{
					j++;
					if( p->cd->pardebug > 1 ) tprintf( "Processor %d\n", j );
					job_wait = 0;
				}
			if( p->cd->pardebug > 1 ) tprintf( "NONE! %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
		}
		if( job_wait )
		{
			// double check
			j = nKids;
			w = 0;
			for( i = 0; i < nHosts; i++ )
			{
				if( kidstatus[i] == 1 && kidids[i] > 0 )
				{
					if( p->cd->pardebug > 1 ) tprintf( "Processor %i [%d] : status = %d %d", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
					if( kidids[i] > 0 && kill( kidids[i], 0 ) == 0 )
					{
						if( p->cd->pardebug > 1 ) tprintf( " is running!\n" );
						w++;
					}
					else if( errno == ESRCH )
					{
						if( p->cd->pardebug > 1 ) tprintf( " does not exist!\n" );
						tprintf( "ERROR: Processor %i [%d] : status = %d %d does not exist\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
						kidstatus[i] = -1;
						kidids[i] = 0;
						j--;
					}
					else
					{
						if( p->cd->pardebug > 1 ) tprintf( " there is a problem to determine the status!\n" );
					}
				}
			}
			if( j < nKids ) { tprintf( "ERROR: %d missing jobs!\n", nKids - j ); }
			if( w < nKids ) { tprintf( "ERROR: job count mismatch (%d < %d); number of running jobs is %d!\n", w, nKids, w ); nKids = w; }
			if( !done )
			{
				if( nKids >= nProc )
				{
					wait = 1;
					if( p->cd->pardebug > 2 ) tprintf( "Waiting ... %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
					else if( p->cd->pardebug > 1 ) tprintf( "Waiting ...\n" );
					sigsuspend( &act.sa_mask );
				}
			}
			else
			{
				if( wait != 2 && p->cd->pardebug > 1 ) tprintf( "All the jobs are started!\n" );
				int nKidsnewcount = 0; for( j = 0; j < nHosts; j++ ) if( kidstatus[j] == 1 && kidattempt[j] > 0 ) nKidsnewcount++;
				if( nKidsnewcount != nKids ) { tprintf( "Kids count mismatch (%d != %d)\n", nKidsnewcount, nKids ); nKids = nKidsnewcount; }
				if( nKids > 0 )
				{
					wait = 2;
					if( p->cd->pardebug > 2 ) tprintf( "Waiting ... %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
					else if( p->cd->pardebug > 1 ) tprintf( "Waiting ... (%d)\n", nKids );
					sigsuspend( &act.sa_mask );
				}
				else
				{
					if( p->cd->pardebug > 2 ) tprintf( "All the jobs are done!\n" );
					break;
				}
			}
			for( j = 0; j < nHosts; ) if( kidstatus[j++] != 1 ) break; // find a kid with status != 1
			if( j > nHosts )
			{
				if( p->cd->pardebug > 1 )
				{
					tprintf( "All the processors are busy! kids = %d\n", nKids );
					for( i = 0; i < nHosts; i++ )
						tprintf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				}
				continue;
			}
		}
		child1 = j; // available kid with status != 1
		child  = j - 1;
		atime = timestamp();
		if( p->cd->pardebug > 1 ) tprintf( "Processor %i [%s:%d] %s : ", child1, kidhost[child], kidids[child], atime );
		destroy = refresh = 0;
		if( rJob > 0 ) { rerun = 1; next = 0; }
		else           { rerun = 0; next = 1; }
		if( kidstatus[child] == 0 ) // kid suspended
		{
			if( done && !rerun ) { if( p->cd->pardebug > 1 ) tprintf( "suspended; will be killed! " ); kiddir[child][0] = 0; destroy = 1; }
			else { if( p->cd->pardebug > 1 ) tprintf( "suspended; new job " ); refork = 0; refresh = 1; }
		}
		else if( kidstatus[child] == -1 ) // kid finished
		{
			if( done && !rerun )
			{
				kiddir[child][0] = 0;
				kidhost[child][0] = 0;
				kidstatus[child] = 1;
				kidattempt[child] = 0;
				kidids[child] = 0;
				if( p->cd->pardebug > 1 ) tprintf( "finished!\n" );
				continue;
			}
			else
			{
				nKids++;
				refork = 1;  // refork; do not refresh
				refresh = 0;
				kidids[child] = 0;
				if( kidattempt[child] == -1 ) { if( p->cd->pardebug > 1 ) tprintf( "qizing " ); }
				else                          { if( p->cd->pardebug > 1 ) tprintf( "finished; Starting " ); }
				kidattempt[child] = 1;
			}
		}
		else if( kidstatus[child] == -2 ) // kid killed
		{
			if( kiddir[child][0] != 0 )
			{
				if( p->cd->pardebug > 1 ) tprintf( "killed; " );
				if( kidattempt[child] >= MAXATTEMPTS ) // ignore the processor
				{
					nFailed++;
					nProc--;
					if( p->cd->pardebug > 1 ) tprintf( "(child %d on %s ignored)\n", child1, kidhost[child] );
					strcpy( rerundir[rJob], kiddir[child] );
					rJob++;
					kiddir[child][0] = 0;
					kidhost[child][0] = 0;
					kidstatus[child] = 1;
					kidids[child] = 0;
					tprintf( "WARNING: The number of currently used processors is decreased: %d of %d!\n", nProc, nHosts );
					continue;
				}
				else // try to restart
				{
					nKids++;
					kidattempt[child]++;
					if( p->cd->pardebug > 1 ) tprintf( "restart (attempt %d of %d)", kidattempt[child], MAXATTEMPTS );
					strcpy( dir, kiddir[child] );
					rerun = 0;
					next = 0;
					refork = 1;
					refresh = 0;
					kidids[child] = 0;
					sleep( kidattempt[child] );
				}
			}
			else { if( p->cd->pardebug > 1 ) tprintf( "killed internally!\n" ); continue; }
		}
		if( destroy && kidids[child] > 0 )
		{
			if( kidids[child] > 0 ) kill( kidids[child] * -1, SIGCONT );
			else tprintf( "ERROR: something is wrong destroying process %d\n!", child1 );
			if( p->cd->pardebug > 1 ) tprintf( "\n" );
			continue;
		}
		if( rerun )
		{
			next = 0;
			rJob--;
			strcpy( dir, rerundir[rJob] );
		}
		else if( next ) // Next job
		{
			cJob++;
			if( cJob >= nJob ) done = 1;
			if( cJob > nJob ) continue;
			if( p->cd->pardebug ) tprintf( "Job %d", cJob );
			if( skip_job[cJob - 1] == 1 )
			{
				nKids--;
				if( p->cd->pardebug ) tprintf( " skipped because it is already completed!\n" );
				continue;
			}
			else
			{
				sprintf( dir, "../%s_%08d", p->cd->mydir_hosts, ieval + cJob ); // Name of directory for parallel runs
				strcpy( kiddir[child], dir );
			}
		}
		if( p->cd->pardebug > 1 ) tprintf( " : \'%s\' in \'%s\'\n", exec_name, dir );
		else if( p->cd->pardebug ) { if( nJob > 1 ) tprintf( " ...\n" ); else tprintf( " ... " ); }
		if( refork )
		{
			if( kidhost[child][0] == 0 ) continue;
			if( cJob > nJob ) continue;
			if( ( return_fork = vfork() ) == 0 )
			{
				/*
				for( w = 0; w < nProc; w++ )
					tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
				for( w = 0; w < nJob; w++ )
					tprintf( "skip_job %d\n", skip_job[w] );
				*/
				// TODO sleep below can be important; try to run without; if it fails, comment it out
				struct timespec tspec, tspec2;
				tspec.tv_sec = 0;
				tspec.tv_nsec = 500L;
				nanosleep( &tspec, &tspec2 );
				pid = getpid();
				setpgid( pid, pid );
				sprintf( buf, "cd %s; %s", dir, exec_name );
				if( p->cd->pardebug > 3 ) tprintf( "Forked Process %i [%s:%d] : \'%s\' in \'%s\'\n", child1, kidhost[child], pid, exec_name, dir );
				if( type == 1 )      execlp( "bpsh", "bpsh", kidhost[child], "/usr/bin/env", "tcsh", "-f", "-c", buf, ( char * ) 0 );
				else if( type == 2 ) execlp( "srun", "srun", "--exclusive", "-N1", "-n1", "/usr/bin/env", "tcsh", "-f", "-c", buf, ( char * ) 0 );
				else                 execlp( "/usr/bin/env", "/usr/bin/env", "tcsh", "-f", "-c", buf, ( char * ) 0 );
				_exit( 7 );
			}
			if( return_fork > 0 )
			{
				kidids[child] = return_fork;
				kidstatus[child] = 1;
				strcpy( kiddir[child], dir );
			}
			else
			{
				nKids--;
				tprintf( "WARNING: fork failed!!!\n" );
			}
		}
		else if( refresh )
		{
			strcpy( kiddir[child], dir );
			tprintf( "Refreshed Process %i [%s:%d] : \'%s\' in \'%s\'\n", child1, kidhost[child], kidids[child], exec_name, kiddir[child] );
			sleep( 1 );
			if( kidids[child] > 0 ) kill( kidids[child] * -1, SIGCONT );
			else tprintf( "ERROR: something is wrong destroying process\n!" );
		}
	}
	/*
	for( w = 0; w < nProc; w++ )
		tprintf( "%i: kidids %10i kidstatus %d kidattempt %d kiddir %s rerundir %s kidhost %s\n", w, kidids[w], kidstatus[w], kidattempt[w], kiddir[w], rerundir[w], kidhost[w] );
	for( w = 0; w < nJob; w++ )
		tprintf( "skip_job %d\n", skip_job[w] );
	*/
	tprintf( "Done:\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
	act.sa_handler = SIG_DFL; // Default handler
	sigaction( SIGCHLD, &act, NULL );
	free( ( void * ) kidids );
	free( ( void * ) kidstatus );
	free( ( void * ) kidattempt );
	free( skip_job );
	free_matrix( ( void ** ) kiddir, nProc );
	free_matrix( ( void ** ) rerundir, nProc );
	if( type != 1 ) free_matrix( ( void ** ) kidhost, nProc );
	p->cd->neval += nJob;
	return( 1 );
}

static void handler( int sig )
{
	pid_t pid;
	int status, i, child, child1 = 0;
	while( ( pid = waitpid( ( pid_t ) - 1, &status, WNOHANG ) ) > 0 )
	{
		for( i = 0; i < nHosts; )
			if( kidids[i++] == pid )
			{
				child1 = i;
				break;
			}
		if( debug ) tprintf( "PROCESS %d [%d] : ", child1, pid );
		child = child1 - 1;
		if( WIFEXITED( status ) )
		{
			nKids--;
			if( debug ) tprintf( "finished (exit %d)\n", WEXITSTATUS( status ) );
			if( WEXITSTATUS( status ) == 1 ) kidstatus[child] = -2;
			else kidstatus[child] = -1;
			kidids[child] = 0;
		}
		else if( WIFSIGNALED( status ) )
		{
			nKids--;
			if( debug ) tprintf( "killed (signal %d)\n", WTERMSIG( status ) );
			kidstatus[child] = -2;
			kidids[child] = 0;
		}
		else { if( debug ) tprintf( "EXIT STATUS IS NOT UNDERSTOOD?!?\n" ); }
	}
}
