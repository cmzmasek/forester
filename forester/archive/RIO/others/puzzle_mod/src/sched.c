/*
 * sched.c
 *
 *
 * Part of TREE-PUZZLE 5.0 (June 2000)
 *
 * (c) 1999-2000 by Heiko A. Schmidt, Korbinian Strimmer,
 *                  M. Vingron, and Arndt von Haeseler
 * (c) 1995-1999 by Korbinian Strimmer and Arndt von Haeseler
 *
 * All parts of the source except where indicated are distributed under
 * the GNU public licence.  See http://www.opensource.org for details.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sched.h"
/* #include "ppuzzle.h" */

#define STDOUT     stdout
#ifndef PARALLEL             /* because printf() runs significantly faster */
                             /* than fprintf(stdout) on an Apple McIntosh  */
                             /* (HS) */
#       define FPRINTF    printf
#       define STDOUTFILE
#else
#       define FPRINTF    fprintf
#       define STDOUTFILE STDOUT,
#endif

int scinit;
int ssinit;
int fscinit;
int gssinit;
int tssinit;

int n, chunksize;
int p;

#ifdef SCHEDTEST
   schedtype testsched;
#endif

void printsched(schedtype sch)
{
   FPRINTF(STDOUTFILE "Current scheduling status:\n");
   FPRINTF(STDOUTFILE "  truetasks=%5ld - alltasks=%5ld - numtasks=%5ld - numprocs=%5d\n",
          sch.truetasks, sch.alltasks, sch.numtasks, sch.numprocs); 
   FPRINTF(STDOUTFILE "  delta    =%5d - overhead=%5d - rest    =%5d - inited  =%5d\n",
          sch.delta, sch.overhead, sch.rest, sch.inited);
   FPRINTF(STDOUTFILE "  nconst   =%5d - fconst  =%5f - lconst  =%5f - kconst  =%5f\n",
          sch.nconst, sch.fconst, sch.lconst, sch.kconst);
}

void initsched(schedtype *sch, uli tasks, int procs, uli minchunk)
{
   if (minchunk < 1) minchunk = 1;
   (*sch).minchunk  = minchunk;
   (*sch).truetasks = tasks;
   (*sch).rest      = (int)((*sch).truetasks % (*sch).minchunk);
   (*sch).alltasks  = (tasks - (*sch).rest);
   (*sch).numtasks  = (*sch).alltasks;
   (*sch).numprocs  = procs;
   (*sch).delta     = 0;
   (*sch).overhead  = 0;
   (*sch).nconst    = 0;
   (*sch).fconst    = 0;
   (*sch).lconst    = 0;
   (*sch).kconst    = 0;
   (*sch).inited    = 0;

#  ifdef PVERBOSE1
      printsched(*sch);
#  endif /* PVERBOSE1 */
}

/**************************************
*  Static Chunking
**************************************/
uli sc(schedtype *sch)
{  
  uli tmp;

  if ((*sch).inited == 0) {
    (*sch).overhead = (*sch).alltasks % (*sch).numprocs;
    (*sch).delta    = ((*sch).alltasks - (*sch).overhead) / (*sch).numprocs;
    (*sch).inited ++;
  }

  if (!(*sch).overhead) {
       if ((*sch).numtasks >= (*sch).delta)
          tmp = (uli)(*sch).delta;
       else
          tmp = 0;
  } else {
       if ((*sch).numtasks >= ((*sch).delta + 1)) {
          tmp = (uli)(*sch).delta + 1;
          (*sch).overhead--;
       } else 
          tmp = 0;
  }

  /* correction */
  if ((tmp % (*sch).minchunk) > 0) {
      tmp += (*sch).minchunk - (tmp % (*sch).minchunk);
  }

  (*sch).numtasks -= tmp;

  if ((*sch).numtasks == 0) {
     tmp += (uli)(*sch).rest;
     (*sch).rest = 0;
  }
  return tmp;
}  /* SC */


/**************************************
*  Self Scheduling
**************************************/
uli ss(schedtype *sch)
{  
  uli tmp;

  if ((*sch).inited == 0) {
     (*sch).inited ++;
  }

  if ((*sch).numtasks >= 1)
     tmp = 1;
  else
     tmp = (*sch).numtasks;

  /* correction */
  if ((tmp % (*sch).minchunk) > 0) {
      tmp += (*sch).minchunk - (tmp % (*sch).minchunk);
  }

  (*sch).numtasks -= tmp;

  if ((*sch).numtasks == 0) {
     tmp += (uli)(*sch).rest;
     (*sch).rest = 0;
  }

  return tmp;
}  /* SS */


/**************************************
*  fixed-size chunking
**************************************/
int fsc()
{  
  static int R ;
  static int delta ;
  static int overhead;

         int tmp;

  if (fscinit == 0) {
    R = n;
    overhead = n % p;
    delta    = (n - overhead) / p;
    fscinit ++;
  }

  if (!overhead) {
       if (R >= delta)
          tmp = delta;
       else
          tmp = 0;
  } else {
       if (R >= (delta + 1)) {
          tmp = delta + 1;
          overhead--;
       } else 
          tmp = 0;
  }

  R -= tmp;
  return tmp;
}  /* FSC */


/**************************************
*  Guided Self Scheduling
**************************************/
uli gss(schedtype *sch)
{  
  uli tmp;

  if ((*sch).inited == 0) {
    (*sch).inited ++;
  }

  if ((*sch).numtasks >= 1) {
     tmp = (uli)ceil((*sch).numtasks / (*sch).numprocs);
     if (tmp == 0) tmp = 1;
  } else
     tmp = 0;

  /* correction */
  if ((tmp % (*sch).minchunk) > 0) {
      tmp += (*sch).minchunk - (tmp % (*sch).minchunk);
  }

  (*sch).numtasks -= tmp;

  if ((*sch).numtasks == 0) {
     tmp += (uli)(*sch).rest;
     (*sch).rest = 0;
  }
  return tmp;
}  /* GSS */

/**************************************
*  Smooth Guided Self Scheduling
**************************************/
uli sgss(schedtype *sch)
{  
  uli tmp;

  if ((*sch).inited == 0) {
    (*sch).inited ++;
  }

  if ((*sch).numtasks >= 1) {
     tmp = (uli)ceil(((*sch).numtasks / (*sch).numprocs) / 2);
     if (tmp == 0) tmp = 1;
  } else
     tmp = 0;

  /* correction */
  if ((tmp % (*sch).minchunk) > 0) {
      tmp += (*sch).minchunk - (tmp % (*sch).minchunk);
  }

  (*sch).numtasks -= tmp;

  if ((*sch).numtasks == 0) {
     tmp += (uli)(*sch).rest;
     (*sch).rest = 0;
  }
  return tmp;
}  /* SGSS */


/**************************************
*  Trapezoid Self Scheduling
**************************************/
uli tss(schedtype *sch)
{
  uli tmp;

  if ((*sch).inited == 0) {
    (*sch).fconst = ceil((*sch).numtasks / (2*(*sch).numprocs));
    if ((*sch).fconst == 0) (*sch).fconst = 1;
    (*sch).lconst = 1;
    (*sch).nconst = ceil( (2*n) / ((*sch).fconst + (*sch).lconst) );
    (*sch).ddelta  = (((*sch).fconst - (*sch).lconst) / ((*sch).nconst - 1));
    (*sch).kconst = (*sch).fconst;
    FPRINTF(STDOUTFILE "f = n/2p = %.2f ; l = %.2f\n", (*sch).fconst, (*sch).lconst);
    FPRINTF(STDOUTFILE "N = 2n/(f+l) = %d ; delta = (f-l)/(N-1) = %.2f\n", (*sch).nconst, (*sch).ddelta);
    (*sch).inited ++;
  }

  if ((*sch).kconst <= (double) (*sch).numtasks) {
     tmp = (uli)ceil((*sch).kconst);
     (*sch).kconst -= (*sch).ddelta;
  } else {
     tmp = (uli)(*sch).numtasks;
     (*sch).kconst = 0.0;
  }

  /* correction */
  if ((tmp % (*sch).minchunk) > 0) {
      tmp += (*sch).minchunk - (tmp % (*sch).minchunk);
  }

  (*sch).numtasks -= tmp;

  if ((*sch).numtasks == 0) {
     tmp += (uli)(*sch).rest;
     (*sch).rest = 0;
  }
  return tmp;

} /* TSS */


/******************/


#ifdef SCHEDTEST
   uli numquarts(int maxspc)
   { 
      uli tmp;
      int a, b, c, d;
 
      if (maxspc < 4) 
         return (uli)0;
      else {
         maxspc--;
         a = maxspc-3;
         b = maxspc-2;
         c = maxspc-1;
         d = maxspc;

         tmp = (uli) 1 + a +
               (uli) b * (b-1) / 2 +
               (uli) c * (c-1) * (c-2) / 6 +
               (uli) d * (d-1) * (d-2) * (d-3) / 24;
         return (tmp);
      }
   }  /* numquarts */
#endif




/**************************************
*  main
**************************************/
#ifdef SCHEDTEST
int main(int argc, char *argv[])
{
  int tcount,
      count,
      lastsize,
      size;
  if ((argc > 4) || (argc < 3)) {
     FPRINTF(STDOUTFILE "\n\n   Usage: %s  <# species> <# processors> [<min chunk size>]\n\n", argv[0]);
     exit(1);
  }

  chunksize = 1;

  switch(argc) {
    case 4:
	  chunksize = atoi(argv[3]);
    case 3:
          n = numquarts(atoi(argv[1]));
          p = atoi(argv[2]);
  }

  FPRINTF(STDOUTFILE "proc=%6d\n", p);
  FPRINTF(STDOUTFILE "task=%6d\n", n);

  initsched(&testsched, n, p, chunksize);
  printsched(testsched);

  count=1; tcount = 0;
  FPRINTF(STDOUTFILE "\n\n---------------------------\n");
  FPRINTF(STDOUTFILE "SC(sched) - Static Chunking\n");
  FPRINTF(STDOUTFILE "---------------------------\n\n");
  do { size = sc(&testsched); 
       if (size > 0) {FPRINTF(STDOUTFILE "%6d. chunk = %6d %c\n", count++, size , (size%chunksize) ? '!' : ' ');
                      tcount+=size;}
       else FPRINTF(STDOUTFILE "%d tasks in %d chunks\n", tcount, (count-1));
     } while (size > 0);


  initsched(&testsched, n, p, chunksize);
  printsched(testsched);

  count=1; tcount = 0;
  FPRINTF(STDOUTFILE "\n\n---------------------------\n");
  FPRINTF(STDOUTFILE "SS(sched) - Self Scheduling\n");
  FPRINTF(STDOUTFILE "---------------------------\n\n");
  do { size = ss(&testsched); 
       if (size > 0) {if (count==1) FPRINTF(STDOUTFILE "%6d. chunk = %6d %c\n", count++, size , (size%chunksize) ? '!' : ' ');
                      count++;
                      tcount+=size;
                      lastsize = size;}
       else          {FPRINTF(STDOUTFILE "      ...\n");
                      FPRINTF(STDOUTFILE "%6d. chunk = %6d %c\n", count++, lastsize , (lastsize%chunksize) ? '!' : ' ');
                      FPRINTF(STDOUTFILE "%d tasks in %d chunks\n", tcount, (count-1));}
     } while (size > 0);


/**/
  count=1; tcount = 0;
  FPRINTF(STDOUTFILE "\n\n---------------------------\n");
  FPRINTF(STDOUTFILE "FSC() - Fixed-Size Chunking\n");
  FPRINTF(STDOUTFILE "---------------------------\n\n");
  do { size = fsc(); 
       if (size > 0) {FPRINTF(STDOUTFILE "%6d. chunk = %6d %c\n", count++, size , (size%chunksize) ? '!' : ' ');
                      tcount+=size;}
       else FPRINTF(STDOUTFILE "%d tasks in %d chunks\n", tcount, (count-1));
     } while (size > 0);
/**/

  initsched(&testsched, n, p, chunksize);
  printsched(testsched);

  count=1; tcount = 0;
  FPRINTF(STDOUTFILE "\n\n-----------------------------------\n");
  FPRINTF(STDOUTFILE "GSS(sched) - Guided Self Scheduling\n");
  FPRINTF(STDOUTFILE "-----------------------------------\n\n");
  do { size = gss(&testsched); 
       if (size > 0) {FPRINTF(STDOUTFILE "%6d. chunk = %6d %c\n", count++, size , (size%chunksize) ? '!' : ' ');
                      tcount+=size;}
       else FPRINTF(STDOUTFILE "%d tasks in %d chunks\n", tcount, (count-1));
     } while (size > 0);

  initsched(&testsched, n, p, chunksize);
  printsched(testsched);

  count=1; tcount = 0;
  FPRINTF(STDOUTFILE "\n\n--------------------------------------\n");
  FPRINTF(STDOUTFILE "TSS(sched) - Trapezoid Self Scheduling\n");
  FPRINTF(STDOUTFILE "--------------------------------------\n\n");
  do { size = tss(&testsched); 
       if (size > 0) {FPRINTF(STDOUTFILE "%6d. chunk = %6d %c\n", count++, size , (size%chunksize) ? '!' : ' ');
                      tcount+=size;}
       else FPRINTF(STDOUTFILE "%d tasks in %d chunks\n", tcount, (count-1));
     } while (size > 0);
  return (0);
}
#endif
