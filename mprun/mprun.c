#include <sys/types.h>  /* include this before any other sys headers */
#include <sys/wait.h>   /* header for waitpid() and various macros */
#include <signal.h>     /* header for signal functions */
#include <stdio.h>      /* header for fprintf() */
#include <unistd.h>     /* header for fork() */
#include <strings.h>    /* header for strcpy() */
#include <string.h>     /* header for strcpy() */
#include <stdlib.h>     /* header for exit() and malloc() */
#include <time.h>       /* header for time() and ctime() */

#include "../mads.h"

static void handler( int sig );

#define MAXATTEMPTS 5

volatile int nProc, nKids, nHosts, *kidstatus, debug;
volatile pid_t *kidids;

/* Functions elsewhere */
char **char_matrix( int maxCols, int maxRows );
void free_matrix( void **matrix, int maxCols );
char *timestamp();
int func_extrn_check_read( int ieval, void *data );

int mprun( int nJob, void *data )
{
	struct opt_data *p = ( struct opt_data * )data;
	struct sigaction act;
	int i, j, ieval, type, cJob, nFailed, child, child1, wait, job_wait, done, next, refork, refresh, destroy, rerun, rJob, *kidattempt, *skip_job;
	pid_t pid, return_fork;
	char *exec_name, **kidhost, **kiddir, **rerundir, dir[255], buf[255], *atime;
	if( p->cd->num_proc <= 1 ) { printf( "\nERROR: Number of available processors is 1; cannot parallelize!\n" ); return( -1 ); }
	else if( p->cd->paral_hosts == NULL ) { if( p->cd->pardebug > 3 ) printf( "\nWARNING: Local runs using %d processors! No parallel hosts!\n", p->cd->num_proc ); type = 0; }
	else { if( p->cd->pardebug ) printf( "Parallel runs using %d hosts!\n", p->cd->num_proc ); type = 1; }
	nProc = nHosts = p->cd->num_proc; // Number of processors/hosts available initially
	exec_name = p->ed->cmdline; // Executable / Execution command line
	ieval = p->cd->neval; // Current number of model evaluations
	kidhost = p->cd->paral_hosts; // List of processors/hosts
	if( nJob > 1 )
	{
		printf( "Parallel execution of %d jobs using %d processors ... ", nJob, nProc );
		if( p->cd->pardebug ) printf( "\n" );
	}
	else if( p->cd->pardebug ) printf( "Parallel execution of 1 job ...\n" );
	skip_job = ( int * ) malloc( nJob * sizeof( int ) );
	if( p->cd->restart ) // Check for already computed jobs (smart restart)
	{
		done = 0;
		for( i = 0; i < nJob; i++ )
		{
			done += skip_job[i] = func_extrn_check_read( ieval + i + 1, p );
			if( p->cd->pardebug > 1 )
			{
				if( skip_job[i] == 1 ) printf( "Job %d will be skipped!\n", ieval + i + 1 );
				else printf( "Job %d will be executed!\n", ieval + i + 1 );
			}
		}
		if( done == nJob ) // All the jobs will be skipped
		{
			p->cd->neval += nJob;
			free( skip_job );
			return( 1 );
		}
	}
	debug = ( p->cd->pardebug > 3 ) ? 1 : 0; // Debug level
	act.sa_handler = handler; // POSIX process handler
	sigemptyset( &act.sa_mask );
	act.sa_flags = 0;
	if( sigaction( SIGCHLD, &act, NULL ) < 0 )
	{
		printf( "sigaction failed!!!\n" );
		free( skip_job );
		return( -1 );
	}
	kidids = ( pid_t * ) malloc( nProc * sizeof( pid_t ) ); memset(( pid_t * ) kidids, ( pid_t ) 0, nProc * sizeof( pid_t ) ); // ID's of external jobs
	kidstatus = ( int * ) malloc( nProc * sizeof( int ) ); memset(( int * ) kidstatus, ( int ) - 1, nProc * sizeof( int ) ); // Status of external jobs
	kidattempt = ( int * ) malloc( nProc * sizeof( int ) ); memset(( int * ) kidattempt, ( int ) - 1, nProc * sizeof( int ) ); // Number of attempts to execute each external job
	kiddir = char_matrix( nProc, 95 ); // Directories for external jobs
	rerundir = char_matrix( nProc, 95 ); // Rerun directories for external jobs
	if( type == 0 )
	{
		kidhost = char_matrix( nProc, 95 );
		for( i = 0; i < nProc; i++ ) strcpy( kidhost[i], "local" );
	}
	nFailed = 0; nKids = 0; cJob = 0; rJob = 0; wait = 0; done = 0;
	while( 1 ) // Main loop
	{
//		if( rJob >= nJob || nProc <= 0 )
		if( rJob > nJob || nProc <= 0 )
		{
			printf( "None of the processors is responding properly! Parallel execution fails!\nrJob = %d nJob = %d nProc = %d\n", rJob, nJob, nProc );
			free(( void * ) kidids ); free(( void * ) kidstatus ); free(( void * ) kidattempt );
			free( skip_job );
			free_matrix(( void ** ) rerundir, nProc ); if( type == 0 ) free_matrix(( void ** ) kidhost, nProc );
			return( -1 );
		}
		job_wait = 1;
		if( rJob > 0 && nKids < ( nHosts - nFailed ) )
		{
			if( p->cd->pardebug > 1 )
			{
				for( i = 0; i < nHosts; i++ )
					printf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
				printf( "Free processor to rerun: " );
			}
			for( j = 0; j < nHosts; j++ )
				if( kidattempt[j] <= 0 )
				{
					j++;
					if( p->cd->pardebug > 1 ) printf( "Processor %d\n", j );
					job_wait = 0;
				}
			if( p->cd->pardebug > 1 ) printf( "NONE! %d %d %d %d\n", nKids, nHosts, nFailed, rJob );
		}
		if( job_wait )
		{
			if( !done )
			{
				if( nKids >= nProc )
				{
					wait = 1;
					if( p->cd->pardebug > 1 ) printf( "Waiting ...\n" );
					sigsuspend( &act.sa_mask );
				}
			}
			else
			{
				if( wait != 2 && p->cd->pardebug > 1 ) printf( "All the jobs are started!\n" );
				if( nKids > 0 )
				{
					wait = 2;
					if( p->cd->pardebug > 1 ) printf( "Waiting ... (%d)\n", nKids );
					sigsuspend( &act.sa_mask );
				}
				else
					break;
			}
			for( j = 0; j < nHosts; ) if( kidstatus[j++] != 1 ) break; // find a kid with status != 1
		}
		if( j > nHosts )
		{
			if( p->cd->pardebug > 1 )
			{
				printf( "All the processors are busy! kids = %d\n", nKids );
				for( i = 0; i < nHosts; i++ )
					printf( "Processor %i [%d] : status = %d %d\n", i + 1, kidids[i], kidstatus[i], kidattempt[i] );
			}
			continue;
		}
		child1 = j; // available kid with status != 1
		child  = j - 1;
		atime = timestamp();
		if( p->cd->pardebug > 1 ) printf( "Processor %i [%s:%d] %s : ", child1, kidhost[child], kidids[child], atime );
		destroy = refresh = 0;
		if( rJob > 0 ) { rerun = 1; next = 0; }
		else           { rerun = 0; next = 1; }
		if( kidstatus[child] == 0 ) // kid suspended
		{
			if( done && !rerun ) { if( p->cd->pardebug > 1 ) printf( "suspended; will be killed! " ); kiddir[child][0] = 0; destroy = 1; }
			else { if( p->cd->pardebug > 1 ) printf( "suspended; new job " ); refork = 0; refresh = 1; }
		}
		else if( kidstatus[child] == -1 ) // kid finished
		{
			if( done && !rerun )
			{
				kiddir[child][0] = 0;
				kidhost[child][0] = 0;
				kidstatus[child] = 1;
				kidattempt[child] = 0;
				if( p->cd->pardebug > 1 ) printf( "finished!\n" );
				continue;
			}
			else
			{
				nKids++;
				refork = 1;  // refork; do not refresh
				refresh = 0;
				if( kidattempt[child] == -1 ) { if( p->cd->pardebug > 1 ) printf( "Initializing " ); }
				else                          { if( p->cd->pardebug > 1 ) printf( "finished; Starting " ); }
				kidattempt[child] = 1;
			}
		}
		else if( kidstatus[child] == -2 ) // kid killed
		{
			if( kiddir[child][0] != 0 )
			{
				if( p->cd->pardebug > 1 ) printf( "killed; " );
				if( kidattempt[child] >= MAXATTEMPTS ) // ignore the processor
				{
					nFailed++;
					nProc--;
					if( p->cd->pardebug > 1 ) printf( "(child %d on %s ignored)\n", child1, kidhost[child] );
					strcpy( rerundir[rJob], kiddir[child] );
					rJob++;
					kiddir[child][0] = 0;
					kidhost[child][0] = 0;
					kidstatus[child] = 1;
					printf( "WARNING: The number of currently used processors is decreased: %d of %d!\n", nProc, nHosts );
					continue;
				}
				else // try to restart
				{
					nKids++;
					kidattempt[child]++;
					if( p->cd->pardebug > 1 ) printf( "restart (attempt %d of %d)", kidattempt[child], MAXATTEMPTS );
					strcpy( dir, kiddir[child] );
					rerun = 0;
					next = 0;
					refork = 1;
					refresh = 0;
//					sleep( kidattempt[child] );
				}
			}
			else { if( p->cd->pardebug > 1 ) printf( "killed internally!\n" ); continue; }
		}
		if( destroy )
		{
			sleep( 1 );
			kill( kidids[child] * -1, SIGCONT );
			if( p->cd->pardebug > 1 ) printf( "\n" );
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
			if( p->cd->pardebug ) printf( "Job %d", cJob );
			if( skip_job[cJob - 1] == 1 )
			{
				if( p->cd->pardebug ) printf( " skipped!\n" );
				if( cJob >= nJob ) done = 1;
				continue;
			}
			else
			{
				sprintf( dir, "../%s_%08d", p->cd->mydir_hosts, ieval + cJob ); // Name of directory for parallel runs
				strcpy( kiddir[child], dir );
				if( cJob >= nJob ) done = 1;
			}
		}
		if( p->cd->pardebug > 1 ) printf( " : \'%s\' in \'%s\'\n", exec_name, dir );
		else if( p->cd->pardebug ) { if( nJob > 1 ) printf( " ...\n" ); else printf( " ... " ); }
		if( refork )
		{
			if( kidhost[child][0] == 0 ) continue;
			if(( return_fork = fork() ) == 0 )
			{
				pid = getpid();
				setpgid( pid, pid );
				sprintf( buf, "cd %s; %s", dir, exec_name );
				if( p->cd->pardebug > 3 ) printf( "Forked Process %i [%s:%d] : \'%s\' in \'%s\'\n", child1, kidhost[child], pid, exec_name, dir );
				if( type ) execlp( "bpsh", "bpsh", kidhost[child], "/bin/tcsh", "-f", "-c", buf, ( char * ) 0 );
				else       execlp( "/bin/tcsh", "-f", "-c", buf, ( char * ) 0 );
				exit( 7 );
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
				printf( "WARNING: fork failed!!!\n" );
			}
		}
		else if( refresh )
		{
			strcpy( kiddir[child], dir );
			printf( "Refreshed Process %i [%s:%d] : \'%s\' in \'%s\'\n", child1, kidhost[child], kidids[child], exec_name, kiddir[child] );
			sleep( 1 );
			kill( kidids[child] * -1, SIGCONT );
		}
	}
	printf( "Done.\n" );
	free(( void * ) kidids ); free(( void * ) kidstatus ); free(( void * ) kidattempt );
	free( skip_job );
	free_matrix(( void ** ) rerundir, nProc ); if( type == 0 ) free_matrix(( void ** ) kidhost, nProc );
	p->cd->neval += nJob;
	return( 1 );
}

static void handler( int sig )
{
	pid_t pid;
	int status, i, child, child1;
	while(( pid = waitpid(( pid_t ) - 1, &status, WNOHANG ) ) > 0 )
	{
		for( i = 0; i < nHosts; )
			if( kidids[i++] == pid )
			{
				child1 = i;
				break;
			}
		if( debug ) printf( "PROCESS %d [%d] : ", child1, pid );
		child = child1 - 1;
		if( WIFEXITED( status ) )
		{
			nKids--;
			if( debug ) printf( "finished (exit %d)\n", WEXITSTATUS( status ) );
			if( WEXITSTATUS( status ) == 1 ) kidstatus[child] = -2;
			else kidstatus[child] = -1;
		}
		else if( WIFSIGNALED( status ) )
		{
			nKids--;
			if( debug ) printf( "killed (signal %d)\n", WTERMSIG( status ) );
			kidstatus[child] = -2;
		}
		else { if( debug ) printf( "EXIT STATUS IS NOT UNDERSTOOD?!?\n" ); }
	}
}
