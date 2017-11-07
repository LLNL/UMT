/*----------------------------------------------------------*/
/* Simple elapsed-time timers for hand instrumenting codes. */
/*----------------------------------------------------------*/
/* link with libtimers.a                                    */
/* need locks if using inside parallel regions              */
/*----------------------------------------------------------*/
/* Fortran:                                                 */
/*    call timer_beg('label')  start timing                 */
/*    call timer_end('label')  stop  timing                 */
/*    call timer_reset()  reset all timers to zero          */
/*    call timer_print()  print timer values and labels     */
/*----------------------------------------------------------*/
/* C / C++:                                                 */
/*    Timer_Beg("label");  start timing                     */
/*    Timer_End("label");  stop  timing                     */
/*    Timer_Reset();  reset all timers to zero              */
/*    Timer_Print();  print timer values and labels         */
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <mpi.h>

#define FIFO_DEPTH 8
static int fifo[FIFO_DEPTH];

#pragma weak timer_beg_=timer_beg
#pragma weak timer_end_=timer_end
#pragma weak timer_print_=timer_print
#pragma weak timer_reset_=timer_reset
#pragma weak timer_enable_=timer_enable
#pragma weak timer_disable_=timer_disable

#pragma weak Timer_beg=Timer_Beg
#pragma weak Timer_end=Timer_End
#pragma weak Timer_print=Timer_Print


/*---------------------*/
/* Function prototypes */
/*---------------------*/
void Timer_Beg(char *);
void Timer_End(char *);
void Timer_Print(void);
void Timer_Reset(void);
void Timer_Enable(void);
void Timer_Disable(void);

void timer_beg(char *, int);
void timer_end(char *, int);
void timer_print(void);
void timer_reset(void);
void timer_enable(void);
void timer_disable(void);


static int index_from_label(char *);

static int printall = 0;


/*---------------------------*/
/* variables with file-scope */
/*---------------------------*/
static int initialized = 0;
static int code_block  = 0;
static int blocked = 0;

#define MAX_CODE_BLOCKS 50

static double timer_in[MAX_CODE_BLOCKS];
static double timer_sum[MAX_CODE_BLOCKS];
static char code_block_label[MAX_CODE_BLOCKS][80];
static int block_starts[MAX_CODE_BLOCKS];
static int block_stops[MAX_CODE_BLOCKS];



/*=======================================================*/
/* Fortran interface for tbeg: terminate the string      */
/* note: the length argument is hidden in Fortran source */
/*=======================================================*/
void timer_beg(char * f_label, int length)
{
   int i, j, rc;
   double xint, xfrac;
   char this_label[80];
   struct timeval tv;
   char * ptr;

   if (blocked) return;

   strncpy(this_label, f_label, length);
   this_label[length] = '\0';

   /*---------------------------------------------------------------*/
   /* if the timers are not initialized, then do the initialization */
   /*---------------------------------------------------------------*/
   if (!initialized)
   {
       ptr = getenv("PRINT_ALL");
       if (ptr != NULL) {
         if (strncasecmp(ptr, "yes", 3) == 0) printall = 1;
       }

       initialized = 1;

       memset(code_block_label, '\0', MAX_CODE_BLOCKS*80);

       /*--------------------------------------*/
       /* set the initial timer values to zero */
       /*--------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++) timer_sum[j] = 0.0;

       /*------------------------------------------*/
       /* zero-out the code block starts and stops */
       /*------------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++)
       {
           block_starts[j] = 0;
           block_stops[j]  = 0;
       }

   }

   j = index_from_label(this_label);

   block_starts[j] += 1;

   gettimeofday(&tv, NULL);
   timer_in[j] = ((double) (tv.tv_sec)) + 1.0e-6*((double) tv.tv_usec);

   return;
}


/*=======================================================*/
/* Fortran interface for tend: terminate the string      */
/* note: the length argument is hidden in Fortran source */
/*=======================================================*/
void timer_end(char * f_label, int length)
{
   int i, j, rc;
   double tnow;
   char this_label[80];
   struct timeval tv;

   if (blocked) return;

   strncpy(this_label, f_label, length);
   this_label[length] = '\0';

   gettimeofday(&tv, NULL);
   tnow = ((double) (tv.tv_sec)) + 1.0e-6*((double) tv.tv_usec);

   if (code_block >= MAX_CODE_BLOCKS) return;

   j = index_from_label(this_label);

   block_stops[j] += 1;

   timer_sum[j] += tnow - timer_in[j];

   return;
}

/*================================*/
/* Fortran interface for tprt     */
/*================================*/
void timer_print(void)
{
   Timer_Print();
}


/*=====================================================*/
/* Initialize and start timing.                        */
/*=====================================================*/
void Timer_Beg(char * this_label)
{
   int i, j, rc;
   double xint, xfrac;
   struct timeval tv;
   char * ptr;

   if (blocked) return;

   /*---------------------------------------------------------------*/
   /* if the timers are not initialized, then do the initialization */
   /*---------------------------------------------------------------*/
   if (!initialized)
   {
       ptr = getenv("PRINT_ALL");
       if (ptr != NULL) {
         if (strncasecmp(ptr, "yes", 3) == 0) printall = 1;
       }

       initialized = 1;

       memset(code_block_label, '\0', MAX_CODE_BLOCKS*80);

       /*--------------------------------------*/
       /* set the initial timer values to zero */
       /*--------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++)
              timer_sum[j] = 0.0;

       /*-------------------------------------------*/
       /* keep track of code block starts and stops */
       /*-------------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++)
       {
           block_starts[j] = 0;
           block_stops[j]  = 0;
       }
   }

   j = index_from_label(this_label);

   block_starts[j] += 1;

   gettimeofday(&tv, NULL);
   timer_in[j] = ((double) (tv.tv_sec)) + 1.0e-6*((double) tv.tv_usec);

   return;
}



/*================================================*/
/* stop timing, save the sum, and continue timing */
/*================================================*/
void Timer_End(char * this_label)
{
   int i, j, rc;
   struct timeval tv;
   double tnow;

   if (blocked) return;

   gettimeofday(&tv, NULL);
   tnow = ((double) (tv.tv_sec)) + 1.0e-6*((double) tv.tv_usec);

   if (code_block >= MAX_CODE_BLOCKS) return;

   j = index_from_label(this_label);

   block_stops[j] += 1;

   timer_sum[j] += tnow - timer_in[j];

   return;
}

/*==========================*/
/* reset all timers to zero */
/*==========================*/
void Timer_Reset(void)
{
   if (code_block > 0) {
      code_block = 0;
      initialized = 0;
   }
}


/*=============================================*/
/* reset all timers to zero, Fortran interface */
/*=============================================*/
void timer_reset(void)
{
   if (code_block > 0) {
      code_block = 0;
      initialized = 0;
   }
}

// service routines to temporarily enable/disable timers
void Timer_Enable(void) 
{
   blocked = 0;
}

void Timer_Disable(void) 
{
   blocked = 1;
}

void timer_enable(void) 
{
   blocked = 0;
}

void timer_disable(void) 
{
   blocked = 1;
}

/*====================================*/
/* print the timer values with labels */
/*====================================*/
void Timer_Print(void)
{
   int rc, i, j, nblocks, myrank, nranks, jobid;
   double elapsed_seconds, avg_seconds;
   int * all_starts;
   double * all_times;
   char filename[240];
   char label[8][80];
   FILE * fp;
   char * ptr = NULL;
   char * trace_dir = NULL;
   struct dataStruct {
                        double value;
                        int rank;
                     };
   struct dataStruct myVal, minVal, maxVal;

   PMPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   PMPI_Comm_size(MPI_COMM_WORLD, &nranks);

   if (myrank == 0) {
      ptr = getenv("LSB_JOBID");
//    ptr = getenv("SLURM_JOBID");
      if (ptr == NULL) jobid = getpid();
      else             jobid = atoi(ptr);

      if (trace_dir != NULL) {
         sprintf(filename, "%s/timing_summary.%d", trace_dir, jobid);
      }
      else {
         sprintf(filename, "timing_summary.%d", jobid);
      }

      fp = fopen(filename, "w");
      if (fp == NULL)
      {
         fprintf(stderr, "from Timer_Print: could not open %s\n", filename);
         fprintf(stderr, "timing summary printed to stderr\n");
         fp = stderr;
      }
   }

   if (myrank == 0) fprintf(fp, "\n");
   if (myrank == 0) fprintf(fp, "------------------------------------------------------------------------------------------------\n");
   if (myrank == 0) fprintf(fp, "Timing  summary:                 #calls       avg(sec)     min(sec) minRank     max(sec) maxRank\n");
   if (myrank == 0) fprintf(fp, "------------------------------------------------------------------------------------------------\n");
   if (code_block >= MAX_CODE_BLOCKS) nblocks = MAX_CODE_BLOCKS;
   else                               nblocks = code_block;

   for (j=0; j<nblocks; j++)
   { 
       if (block_starts[j] == block_stops[j])
       {
           elapsed_seconds = timer_sum[j];
           myVal.value = elapsed_seconds;
           myVal.rank  = myrank;
           PMPI_Allreduce(&myVal, &minVal, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
           PMPI_Allreduce(&myVal, &maxVal, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
           PMPI_Allreduce(&elapsed_seconds, &avg_seconds, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
           avg_seconds = avg_seconds / ((double) nranks);
           
           if (myrank == 0) fprintf(fp, "%-28s  %9d  %12.3f %12.3lf %7d %12.3lf %7d\n", 
                  code_block_label[j], block_starts[j], avg_seconds, minVal.value, minVal.rank, maxVal.value, maxVal.rank);
       }
       else
       {
           if (myrank == 0) fprintf(fp, "mismatch in starts/stops for code block '%s'\n", code_block_label[j]);
           if (myrank == 0) fprintf(fp, "  starts = %d\n", block_starts[j]);
           if (myrank == 0) fprintf(fp, "  stops  = %d\n", block_stops[j]);
       }
   }

   if (printall) {
      all_starts = (int *) malloc(nranks*sizeof(int));
      all_times  = (double *) malloc(nranks*sizeof(double));
      for (j=0; j<nblocks; j++)
      {
         PMPI_Gather(&block_starts[j], 1, MPI_INT, all_starts, 1, MPI_INT, 0, MPI_COMM_WORLD);
         PMPI_Gather(&timer_sum[j], 1, MPI_DOUBLE, all_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         if (myrank == 0) {
            if (block_starts[j] == block_stops[j])
            {
               fprintf(fp, "\nData for code block %s :\n", code_block_label[j]);
               fprintf(fp, "  rank      #calls     elapsed(sec)\n");
               for (i=0; i<nranks; i++) 
               {
                  fprintf(fp, "%6d  %9d  %12.3f\n", i, all_starts[i], all_times[i]); 
               }
            }
         }
      }
   }

   if (myrank == 0) fclose(fp);

   return;
}


/*===========================================*/
/* Find the code-block number from the label.*/
/*===========================================*/
int index_from_label(char * this_label)
{
   int i, k, match;
   char * ptr;

   match = 0;

   /*----------------------------*/
   /* first check the fifo queue */
   /*----------------------------*/
   if (code_block > FIFO_DEPTH)
   {
      for (k=0; k<FIFO_DEPTH; k++)
      {
          i = fifo[k];
          if (i >= 0)
          {
              if (0 == strcmp(code_block_label[i], this_label))
              {
                  match = 1;
                  break;
              }
          }
      }
   }

   if (match == 1) return i;

   /*-------------------------------------------------------*/
   /* not found in the fifo, so check all known code blocks */
   /*-------------------------------------------------------*/
   for (i=code_block-1; i>=0; i--)
   {
       if (0 == strcmp(code_block_label[i], this_label))
       {
           match = 1;
           break;
       }
   }

   /*------------------------------------------------*/
   /* if there is no match, this is a new code block */
   /*------------------------------------------------*/
   if (match == 0)
   {
       i = code_block;
       ptr = strcpy(code_block_label[i], this_label);
       if (ptr == NULL) code_block_label[i][0] = '\0';
       code_block ++;
   }

   /*-----------------------------------------------*/
   /* save the latest code block in the fifo        */
   /*-----------------------------------------------*/
   for (k=FIFO_DEPTH-2; k>=0; k--) fifo[k+1] = fifo[k];

   fifo[0] = i;


   return i;

}
