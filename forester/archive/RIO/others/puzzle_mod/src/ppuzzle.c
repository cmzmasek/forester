/*
 *   ppuzzle.c
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

#define EXTERN extern
 
#include <mpi.h>
#include <time.h>
#include "ppuzzle.h"
 

int PP_IamMaster;
int PP_IamSlave;
int PP_Myid;
int PP_MyMaster;
int PP_NumProcs;
MPI_Comm PP_Comm;

int *freeslaves;           /* Queue of free slaves */
int firstslave,            /* headpointer of queue */
    lastslave;             /* tailpointer of queue */

int *permutsent,
    *permutrecved,
    *quartsent,
    *quartrecved,
    *doquartsent,
    *doquartrecved,
    *splitsent,
    *splitrecved,
    *permutsentn,
    *permutrecvedn,
    *quartsentn,
    *quartrecvedn,
    *doquartsentn,
    *doquartrecvedn,
    *splitsentn,
    *splitrecvedn;

double *walltimes,
       *cputimes;
double *fullwalltimes,
       *fullcputimes;
double *altwalltimes,
       *altcputimes;

int PP_permutsent = 0;          /* # of  */
int PP_permutrecved = 0;        /* # of  */
int PP_quartsent = 0;           /* # of  */
int PP_quartrecved = 0;         /* # of  */
int PP_doquartsent = 0;         /* # of  */
int PP_doquartrecved = 0;       /* # of  */
int PP_splitsent = 0;           /* # of  */
int PP_splitrecved = 0;         /* # of  */
int PP_permutsentn = 0;          /* # of  */
int PP_permutrecvedn = 0;        /* # of  */
int PP_quartsentn = 0;           /* # of  */
int PP_quartrecvedn = 0;         /* # of  */
int PP_doquartsentn = 0;         /* # of  */
int PP_doquartrecvedn = 0;       /* # of  */
int PP_splitsentn = 0;           /* # of  */
int PP_splitrecvedn = 0;         /* # of  */

double PP_starttime     = 0,
 PP_stoptime      = 0,
 PP_inittime      = 0,
 PP_paramcomptime = 0,
 PP_paramsendtime = 0,
 PP_quartcomptime = 0,
 PP_quartsendtime = 0,
 PP_puzzletime    = 0,
 PP_treetime      = 0,
 PP_lasttime      = 0;

int PP_MaxSlave = 0;
 

/*********************************************************************
*  miscellaneous utilities                                           *
*********************************************************************/

int dcmp(const void *a, const void *b)
{
	if (*(double *)a > *(double *)b) return (-1);
	else if (*(double *)a < *(double *)b) return 1;
	else return 0;
}

/******************/

void PP_cmpd(int rank, double a, double b)
{
  if (a != b)
     FPRINTF(STDOUTFILE "(%2d) *** %.3f != %.3f\n", rank, a, b);
}

/******************/

void PP_cmpi(int rank, int a, int b)
{
  if (a != b)
     FPRINTF(STDOUTFILE "(%2d) *** %d != %d\n", rank, a, b);
}

/******************/

double PP_timer()
{
  double tmptime;
  if (PP_lasttime == 0) {
     PP_lasttime = MPI_Wtime();
     return(0);
     }
  else {
     tmptime = PP_lasttime;
     PP_lasttime = MPI_Wtime();
     return(PP_lasttime - tmptime);
  }
}

/******************/
 
void PP_Printerror(FILE *of, int id, int err)
{
        char  errstr[MPI_MAX_ERROR_STRING];
        int   errstrlen;

  if ((err > MPI_SUCCESS) && (err <= MPI_ERR_LASTCODE)) {
    MPI_Error_string(err, errstr, &errstrlen);
    fprintf(of, "(%2d) MPI ERROR %d : %s\n", id, err, errstr);
    }
  else {
    if (err == MPI_SUCCESS) 
      fprintf(of, "(%2d) MPI ERROR %d : No error\n", id, err);
    else
      fprintf(of, "(%2d) MPI ERROR %d : unknown error number\n", id, err);
  }
} /* PP_Printerror */

/******************/

void PP_Printbiparts(cmatrix biparts)
{ int n1, n2;
    for (n1=0; n1<(Maxspc-3); n1++) {
       if (n1==0) FPRINTF(STDOUTFILE "(%2d) bipartition : ", PP_Myid);
       else       FPRINTF(STDOUTFILE "(%2d)             : ", PP_Myid);
       for (n2=0; n2<Maxspc; n2++)
          FPRINTF(STDOUTFILE "%c", biparts[n1][n2]);
    FPRINTF(STDOUTFILE "\n");
    }
}  /* PP_Printbiparts */

/******************/

#if 0
void num2quart(uli qnum, int *a, int *b, int *c, int *d)
{
	double temp;
	uli aa, bb, cc, dd;
	uli lowval=0, highval=0;

	aa=0; bb=1; cc=2; dd=3;

	temp = (double)(24 * qnum);
	temp = sqrt(temp);
	temp = sqrt(temp);
	/* temp = pow(temp, (double)(1/4)); */
	dd = (uli) floor(temp) + 1;
	if (dd < 3) dd = 3; 
	lowval =  (uli) dd*(dd-1)*(dd-2)*(dd-3)/24;
	highval = (uli) (dd+1)*dd*(dd-1)*(dd-2)/24;
	if (lowval >= qnum)
	    while ((lowval > qnum)) {
	        dd -= 1; lowval = (uli) dd*(dd-1)*(dd-2)*(dd-3)/24; 
	    }
	else {
	    while (highval <= qnum) {
	        dd += 1; highval = (uli) (dd+1)*dd*(dd-1)*(dd-2)/24; 
	    }
	    lowval = (uli) dd*(dd-1)*(dd-2)*(dd-3)/24; 
	}
	qnum -= lowval;
	if (qnum > 0) {
	    temp = (double)(6 * qnum);
	    temp = pow(temp, (double)(1/3));
	    cc = (uli) floor(temp);
	    if (cc < 2) cc= 2;
	    lowval =  (uli) cc*(cc-1)*(cc-2)/6;
	    highval = (uli) (cc+1)*cc*(cc-1)/6;
	    if (lowval >= qnum)
	        while ((lowval > qnum)) {
	           cc -= 1; lowval = (uli) cc*(cc-1)*(cc-2)/6; 
	        }
	    else {
	        while (highval <= qnum) {
	           cc += 1; highval = (uli) (cc+1)*cc*(cc-1)/6; 
	        }
	        lowval = (uli) cc*(cc-1)*(cc-2)/6; 
	    }
	    qnum -= lowval;
	    if (qnum > 0) {
	        temp = (double)(2 * qnum);
	        temp = sqrt(temp);
	        bb = (uli) floor(temp);
	        if (bb < 1) bb= 1;
	        lowval =  (uli) bb*(bb-1)/2;
	        highval = (uli) (bb+1)*bb/2;
	        if (lowval >= qnum)
	            while ((lowval > qnum)) {
	                bb -= 1; lowval = (uli) bb*(bb-1)/2; 
	            }
	        else {
	            while (highval <= qnum) {
	               bb += 1; highval = (uli) (bb+1)*bb/2;
	            }
	            lowval = (uli) bb*(bb-1)/2; 
	        }
	        qnum -= lowval;
	        if (qnum > 0) {
	           aa = (uli) qnum;
	           if (aa < 0) aa= 0;
	        }
	    }
	}
	*d = (int)dd;
	*c = (int)cc;
	*b = (int)bb;
	*a = (int)aa;
}  /* num2quart */

/******************/

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

/******************/

uli quart2num (int a, int b, int c, int d)
{
      uli tmp;
      if ((a>b) || (b>c) || (c>d)) {
         fprintf(stderr, "Error PP5 not (%d <= %d <= %d <= %d) !!!\n", a, b, c, d);
         exit (1);
      }
      tmp = (uli) a +
            (uli) b * (b-1) / 2 +
            (uli) c * (c-1) * (c-2) / 6 +
            (uli) d * (d-1) * (d-2) * (d-3) / 24;
      return (tmp);
}  /* quart2num */
#endif
/******************/


/*********************************************************************
*  queue for storing the ranks of slaves waiting for work            *
*********************************************************************/

void PP_initslavequeue()
{
  int n;
  freeslaves = new_ivector(PP_NumProcs);
    firstslave = 0;
    PP_MaxSlave = PP_NumProcs-1;
    lastslave  = PP_MaxSlave-1;
    freeslaves[PP_MaxSlave] = PP_MaxSlave;
    for (n=0; n<PP_MaxSlave; n++){
      freeslaves[n] = n+1;
    }
}

/******************/

int PP_fullslave()
{
return (freeslaves[PP_MaxSlave] == PP_MaxSlave);
}

/******************/

int PP_emptyslave()
{
return (freeslaves[PP_MaxSlave] == 0);
}

/******************/

void PP_putslave(int sl)
{
  if (freeslaves[PP_MaxSlave] == PP_MaxSlave) {
    FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR PP1 TO DEVELOPERS\n\n\n");
    MPI_Finalize();
    exit(1);
  }
  (freeslaves[PP_MaxSlave])++;
  lastslave = (lastslave + 1) % PP_MaxSlave;
  freeslaves[lastslave] = sl;
}

/******************/

int PP_getslave()
{
  int tmp;
  if (freeslaves[PP_MaxSlave] == 0)
    return MPI_PROC_NULL;
  else {
    tmp = freeslaves[firstslave];
    firstslave = (firstslave + 1) % PP_MaxSlave;
    (freeslaves[PP_MaxSlave])--;
    return tmp;
  }
}

/*********************************************************************
*  procedures to parallelize the puzzling step                       *
*********************************************************************/

void PP_slave_do_puzzling(ivector trueID)
{
int    i, a, b, c;
uli    nq;
#if SEQUENTIAL
       double tc2, mintogo, minutes, hours;
#endif

     /* initialize tree */
     inittree();

     /* adding all other leafs */
     for (i = 3; i < Maxspc; i++) { 
 
         /* clear all edgeinfos */
         resetedgeinfo();
 
         /* clear counter of quartets */
         nq = 0;
 
         /*
           * core of quartet puzzling algorithm
          */
 
          for (a = 0; a < nextleaf - 2; a++)
              for (b = a + 1; b < nextleaf - 1; b++)
                  for (c = b + 1; c < nextleaf; c++) {
 
                      /*      check which two _leaves_ out of a, b, c
                          are closer related to each other than
                          to leaf i according to a least squares
                          fit of the continous Baysian weights to the
                          seven trivial "attractive regions". We assign 
                          a score of 1 to all edges between these two leaves
                          chooseA and chooseB */
 
                      checkquartet(a, b, c, i);
                      incrementedgeinfo(chooseA, chooseB);
 
                      nq++;
 
#                     if SEQUENTIAL
                         /* generate message every 15 minutes */
    
                         /* check timer */
                         time(&time2);
                         if ( (time2 - time1) > 900) {
                              /* every 900 seconds */
                              /* percentage of completed trees */
                              if (mflag == 0) {
                                      FPRINTF(STDOUTFILE "\n");
                                      mflag = 1;
                              }
                              tc2 = 100.0*Currtrial/Numtrial + 
                                    100.0*nq/Numquartets/Numtrial;
                              mintogo = (100.0-tc2) *
                                        (double) (time2-time0)/60.0/tc2;
                              hours = floor(mintogo/60.0);
                              minutes = mintogo - 60.0*hours;
                              FPRINTF(STDOUTFILE "%2.2f%%", tc2);
                              FPRINTF(STDOUTFILE " completed (remaining");
                              FPRINTF(STDOUTFILE " time: %.0f", hours);
                              FPRINTF(STDOUTFILE " hours %.0f", minutes);
                              FPRINTF(STDOUTFILE " minutes)\n");
                              time1 = time2;
                         }
#                     endif /* SEQUENTIAL */
              }

          /* find out which edge has the lowest edgeinfo */
          minimumedgeinfo();
 
          /* add the next leaf on minedge */
          addnextleaf(minedge);
     }
 
     /* compute bipartitions of current tree */
     computebiparts();

#if  PARALLEL
       if (PP_IamMaster) makenewsplitentries();
#    else
       makenewsplitentries();
#    endif 

	{
		int *ctree, startnode;
		char *trstr;
		ctree = initctree();
		copytree(ctree);
		startnode = sortctree(ctree);
		trstr=sprintfctree(ctree, psteptreestrlen);
		(void) addtree2list(&trstr, 1, &psteptreelist, &psteptreenum, &psteptreesum);
#		ifdef PVERBOSE2
			/* fprintf(STDOUT, "%s\n", trstr); */
			printfpstrees(psteptreelist);
#		endif
		freectree(&ctree);
	}


     /* free tree before building the next tree */
     freetree();

} /* PP_slave_do_puzzling */

/******************/

void PP_do_puzzling(ivector trueID)
{
int    dest;

#    if PARALLEL
       dest = PP_getslave();
       PP_SendPermut(dest, Maxspc, trueID);
#    endif

     /* initialize tree */
     inittree();

     PP_RecvSplits(Maxspc, biparts); 

#    ifdef PVERBOSE3 
       PP_Printbiparts(biparts);
#    endif /* PVERBOSE3 */

     makenewsplitentries();

     /* free tree before building the next tree */
     freetree();

} /* PP_do_puzzling */
 
/******************/


void PP_do_write_quart(int e,
                       int f,
                       int g,
                       int h,
                       double d1,
                       double d2,
                       double d3,
                       uli   *numbq,
                       uli   *bqarr)
{
	double lhs[3],
	       temp,
	       wlist[6],
	       plist[6];
        unsigned char qpbranching;
        int badquartet;

        lhs[0] = d1;
	lhs[1] = d2;
	lhs[2] = d3;

         badquartet = FALSE;

         /* compute Bayesian weights */
         temp = (lhs[0] + lhs[1] + lhs[2])/3.0;
         lhs[0] = exp(lhs[0] - temp);
         lhs[1] = exp(lhs[1] - temp);
         lhs[2] = exp(lhs[2] - temp);
         temp = lhs[0] + lhs[1] + lhs[2];
         wlist[0] = lhs[0] / temp;
         wlist[1] = 1.0;
         wlist[2] = lhs[1] / temp;
         wlist[3] = 2.0;
         wlist[4] = lhs[2] / temp;
         wlist[5] = 4.0;

         /* sort in descending order */
         qsort(wlist, 3, 2*sizeof(double), dcmp);

         /* check out the three possibilities */

         /* 100 distribution */
         plist[0] = (1.0 - wlist[0])*(1.0 - wlist[0]) +
                    (0.0 - wlist[2])*(0.0 - wlist[2]) +
                    (0.0 - wlist[4])*(0.0 - wlist[4]);
         plist[1] = wlist[1];

         /* 110 distribution */
         plist[2] = (0.5 - wlist[0])*(0.5 - wlist[0]) +
                    (0.5 - wlist[2])*(0.5 - wlist[2]) +
                    (0.0 - wlist[4])*(0.0 - wlist[4]);
         plist[3] = wlist[1] + wlist[3];

         /* 111 distribution */
         temp = 1.0/3.0;
         plist[4] = (temp - wlist[0])*(temp - wlist[0]) +
                    (temp - wlist[2])*(temp - wlist[2]) +
                    (temp - wlist[4])*(temp - wlist[4]);
         plist[5] = wlist[1] + wlist[3] + wlist[5];

         /* sort in descending order */
         qsort(plist, 3, 2*sizeof(double), dcmp);

         qpbranching = (unsigned char) plist[5];
         writequartet(e, f, g, h, qpbranching);

         /* a bad quartet is a quartet that shows
            equal weights for all three possible topologies */
         if (qpbranching == 7) badquartet = TRUE;

         if (badquartet) {
               bqarr[(*numbq)++] = quart2num(e, f, g, h);
#	       ifdef PVERBOSE3
	          FPRINTF(STDOUTFILE "(%2d) bad quartet: %d  %d  %d %d -> %ld\n", 
	                  PP_Myid, e, f, g, h, quart2num(e, f, g, h));
#	       endif /* PVERBOSE3 */
               badqs++;
               badtaxon[e]++;
               badtaxon[f]++;
               badtaxon[g]++;
               badtaxon[h]++;
         } /* if badquartet */
} /* PP_do_write_quart */

/*********************************************************************
*  sending/receiving the important sizes and parameter  (M->S)       *
*********************************************************************/
 
void PP_SendSizes(int    mspc, 
                  int    msite,
                  int    ncats,
                  int    nptrn,
                  int    rad,
                  int    outgr,
                  double frconst,
                  int    rseed)
{
# define NUMINT 7
# define NUMDBL 1
  int    ints[NUMINT];
  double doubles[NUMDBL];
  MPI_Datatype Dtypes[2] =    {MPI_INT, MPI_DOUBLE};
  int          Dtypelens[2] = {NUMINT , NUMDBL};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_Sizes;
  int          dest;
  int          error;

# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: Maxspc=%d  Maxsite=%d  numcats=%d\n", PP_Myid, mspc, msite, ncats);
    FPRINTF(STDOUTFILE "(%2d)          Numprtn=%d  tpmradix=%d  fracconst=%.3f\n", PP_Myid, nptrn, rad, frconst);
# endif /* PVERBOSE2 */

  ints[0] = mspc;
  ints[1] = msite;
  ints[2] = ncats;
  ints[3] = nptrn;
  ints[4] = rad;
  ints[5] = outgr;
  ints[6] = rseed;
  doubles[0] = frconst;

  MPI_Address(ints,     Dtypeaddr);
  MPI_Address(doubles, (Dtypeaddr+1));
  
  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_Sizes);
  MPI_Type_commit(&PP_Sizes);

  for (dest=1; dest<PP_NumProcs; dest++) {

    error = MPI_Ssend(MPI_BOTTOM, 1, PP_Sizes, dest, PP_SIZES, PP_Comm);
    if (error != MPI_SUCCESS)
      PP_Printerror(STDOUT, 600+PP_Myid, error);

#   ifdef PVERBOSE3
       FPRINTF(STDOUTFILE "(%2d) -> (%2d) Sent Sizes\n", PP_Myid, dest);
#   endif /* PVERBOSE3 */

  } /* for each slave */

  MPI_Type_free(&PP_Sizes);

# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent Sizes\n", PP_Myid);
# endif /* PVERBOSE3 */

# undef NUMINT
# undef NUMDBL
} /* PP_SendSizes */


/******************/

void PP_RecvSizes(int    *mspc, 
                  int    *msite,
                  int    *ncats,
                  int    *nptrn,
                  int    *rad,
                  int    *outgr,
                  double *frconst,
                  int    *rseed)
{
# define NUMINT 7
# define NUMDBL 1
  int    ints[NUMINT];
  double doubles[NUMDBL];
  MPI_Datatype Dtypes[2] =    {MPI_INT, MPI_DOUBLE};
  int          Dtypelens[2] = {NUMINT , NUMDBL};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_Sizes;
  MPI_Status   stat;
  int          error;
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) Receiving Sizes ...\n", PP_Myid);
# endif /* PVERBOSE3 */

  MPI_Address(ints,     Dtypeaddr);
  MPI_Address(doubles, (Dtypeaddr+1));
  
  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_Sizes);
  MPI_Type_commit(&PP_Sizes);
 
  error = MPI_Probe(PP_MyMaster, MPI_ANY_TAG, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 700+PP_Myid, error);
  if (stat.MPI_TAG != PP_SIZES) {
	if (stat.MPI_TAG == PP_DONE) {
		PP_RecvDone();
#		ifdef PVERBOSE1
			FPRINTF(STDOUTFILE "(%2d) Finishing...\n", PP_Myid);
#		endif /* PVERBOSE1 */
		MPI_Finalize();
		exit(1);
	} else {
		FPRINTF(STDOUTFILE "(%2d) Error: unexpected TAG received...\n", PP_Myid);
		MPI_Finalize();
		exit(1);
	}
  }

  error = MPI_Recv(MPI_BOTTOM, 1, PP_Sizes, PP_MyMaster, MPI_ANY_TAG, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 700+PP_Myid, error);
  if (stat.MPI_TAG != PP_SIZES) {
	FPRINTF(STDOUTFILE "(%2d) Error: unexpected TAG received...\n", PP_Myid);
	MPI_Finalize();
	exit(1);
  }

  *mspc    = ints[0];
  *msite   = ints[1];
  *ncats   = ints[2];
  *nptrn   = ints[3];
  *rad     = ints[4];
  *outgr   = ints[5];
  *rseed   = ints[6];
  *frconst = doubles[0];
 
  MPI_Type_free(&PP_Sizes);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) <- (%2d) Received: Maxspec=%d  Maxsite=%d  numcats=%d\n", PP_Myid, PP_MyMaster, *mspc, *msite, *ncats);
    FPRINTF(STDOUTFILE "(%2d)                   Numprtn=%d  tpmradix=%d  fracconst=%.3f\n", PP_Myid, *nptrn, *rad, *frconst);
# endif /* PVERBOSE2 */

# undef NUMINT
# undef NUMDBL
} /* PP_RecvSizes */
 


/*********************************************************************
*  sending/receiving the data matrizes  (M->S)                       *
*********************************************************************/

void PP_RecvData(
	cmatrix Seqpat,           /* cmatrix (Maxspc x Numptrn)             */
	ivector Alias,            /* ivector (Maxsite)                      */
	ivector Weight,           /* ivector (Numptrn)                      */
	ivector constpat,
	dvector Rates,            /* dvector (numcats)                      */
	dvector Eval,             /* dvector (tpmradix)                     */
	dvector Freqtpm,
	dmatrix Evec,             /* dmatrix (tpmradix x tpmradix)          */
	dmatrix Ievc,
	dmatrix iexp,
	dmatrix Distanmat,        /* dmatrix (Maxspc x Maxspc)              */
	dcube   ltprobr)          /* dcube (numcats x tpmradix x tpmradix)  */
{
  MPI_Datatype Dtypes[12];
  int          Dtypelens[12];
  MPI_Aint     Dtypeaddr[12];
  MPI_Datatype PP_Data;
  MPI_Status   stat;
  int          error;
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) Receiving Sizes ...\n", PP_Myid);
# endif /* PVERBOSE2 */

  Dtypes  [0] = MPI_CHAR; Dtypelens [0] = Maxspc * Numptrn;
  MPI_Address(&(Seqpat[0][0]), &(Dtypeaddr[0]));
  Dtypes  [1] = MPI_INT; Dtypelens  [1] = Maxsite ;
  MPI_Address(&(Alias[0]), &(Dtypeaddr[1]));
  Dtypes  [2] = MPI_INT; Dtypelens  [2] = Numptrn ;
  MPI_Address(&(Weight[0]), &(Dtypeaddr[2]));
  Dtypes  [3] = MPI_INT; Dtypelens  [3] = Numptrn ;
  MPI_Address(&(constpat[0]), &(Dtypeaddr[3]));
  Dtypes  [4] = MPI_DOUBLE; Dtypelens  [4] = numcats ;
  MPI_Address(&(Rates[0]), &(Dtypeaddr[4]));
  Dtypes  [5] = MPI_DOUBLE; Dtypelens  [5] = tpmradix ;
  MPI_Address(&(Eval[0]), &(Dtypeaddr[5]));
  Dtypes  [6] = MPI_DOUBLE; Dtypelens  [6] = tpmradix ;
  MPI_Address(&(Freqtpm[0]), &(Dtypeaddr[6]));
  Dtypes  [7] = MPI_DOUBLE; Dtypelens  [7] = tpmradix * tpmradix ;
  MPI_Address(&(Evec[0][0]), &(Dtypeaddr[7]));
  Dtypes  [8] = MPI_DOUBLE; Dtypelens  [8] = tpmradix * tpmradix ;
  MPI_Address(&(Ievc[0][0]), &(Dtypeaddr[8]));
  Dtypes  [9] = MPI_DOUBLE; Dtypelens [9] = tpmradix * tpmradix ;
  MPI_Address(&(iexp[0][0]), &(Dtypeaddr[9]));
  Dtypes [10] = MPI_DOUBLE; Dtypelens [10] = Maxspc * Maxspc ;
  MPI_Address(&(Distanmat[0][0]), &(Dtypeaddr[10]));
  Dtypes [11] = MPI_DOUBLE; Dtypelens [11] = numcats * tpmradix * tpmradix ;
  MPI_Address(&(ltprobr[0][0][0]), &(Dtypeaddr[11]));
 
  MPI_Type_struct(12, Dtypelens, Dtypeaddr, Dtypes, &PP_Data);
  MPI_Type_commit(&PP_Data);
 

  error = MPI_Probe(PP_MyMaster, MPI_ANY_TAG, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 700+PP_Myid, error);
  if (stat.MPI_TAG != PP_DATA) {
        if (stat.MPI_TAG == PP_DONE) {
                PP_RecvDone();
#               ifdef PVERBOSE1
                        FPRINTF(STDOUTFILE "(%2d) Finishing...\n", PP_Myid);
#               endif /* PVERBOSE1 */
                MPI_Finalize();
                exit(1);
        } else {
                FPRINTF(STDOUTFILE "(%2d) Error: unexpected TAG received...\n", PP_Myid);
                MPI_Finalize();
                exit(1);
        }
  }


  error = MPI_Recv(MPI_BOTTOM, 1, PP_Data, PP_MyMaster, PP_DATA, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 900+PP_Myid, error);

# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) <- (%2d) Received : Alias(0)=%d - Weight(0)=%d - constpat(0)=%d\n", PP_Myid, PP_MyMaster, Alias[0], Weight[0], constpat[0]);
    FPRINTF(STDOUTFILE "(%2d)                    Rates(0)=%.3f - Eval(0)=%.3f - Freqtpm(0)=%.3f\n", PP_Myid, Rates[0], Eval[0], Freqtpm[0]);
    FPRINTF(STDOUTFILE "(%2d)                    Evec(0,0)=%.3f - Ievc(0,0)=%.3f - iexp(0,0)=%.3f - Distanmat(0,1)=%.3f\n", PP_Myid, Evec[0][0], Ievc[0][0], iexp[0][0], Distanmat[0][1]);
    FPRINTF(STDOUTFILE "(%2d)                    Distanmat(0,1)=%.3f\n", PP_Myid, Distanmat[0][1]);
    FPRINTF(STDOUTFILE "(%2d)                    ltprobr(0,0,0)=%.3f\n", PP_Myid, ltprobr[0][0][0]);
# endif /* PVERBOSE2 */
 
  MPI_Type_free(&PP_Data);
 
} /* PP_RecvData */


/******************/

void PP_SendData(
        cmatrix Seqpat,           /* cmatrix (Maxspc x Numptrn)             */
        ivector Alias,            /* ivector (Maxsite)                      */
        ivector Weight,           /* ivector (Numptrn)                      */
        ivector constpat,
        dvector Rates,            /* dvector (numcats)                      */
        dvector Eval,             /* dvector (tpmradix)                     */
        dvector Freqtpm,
        dmatrix Evec,             /* dmatrix (tpmradix x tpmradix)          */
        dmatrix Ievc,
        dmatrix iexp,
        dmatrix Distanmat,        /* dmatrix (Maxspc x Maxspc)              */
        dcube   ltprobr)          /* dcube (numcats x tpmradix x tpmradix)  */
{
  MPI_Datatype Dtypes[12];
  int          Dtypelens[12];
  MPI_Aint     Dtypeaddr[12];
  MPI_Datatype PP_Data;
  int          dest;
  int          error;
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: Alias(0)=%d - Weight(0)=%d - constpat(0)=%d\n", PP_Myid, Alias[0], Weight[0], constpat[0]);
    FPRINTF(STDOUTFILE "(%2d)          Rates(0)=%.3f - Eval(0)=%.3f - Freqtpm(0)=%.3f\n", PP_Myid, Rates[0], Eval[0], Freqtpm[0]);
    FPRINTF(STDOUTFILE "(%2d)          Evec(0,0)=%.3f - Ievc(0,0)=%.3f - iexp(0,0)=%.3f - Distanmat(0,1)=%.3f\n", PP_Myid, Evec[0][0], Ievc[0][0], iexp[0][0], Distanmat[0][1]);
    FPRINTF(STDOUTFILE "(%2d)          ltprobr(0,0,0)=%.3f\n", PP_Myid, ltprobr[0][0][0]);
# endif /* PVERBOSE2 */
 
  Dtypes  [0] = MPI_CHAR; Dtypelens [0] = Maxspc * Numptrn;
  MPI_Address(&(Seqpat[0][0]), &(Dtypeaddr[0]));
  Dtypes  [1] = MPI_INT; Dtypelens  [1] = Maxsite ;
  MPI_Address(&(Alias[0]), &(Dtypeaddr[1]));
  Dtypes  [2] = MPI_INT; Dtypelens  [2] = Numptrn ;
  MPI_Address(&(Weight[0]), &(Dtypeaddr[2]));
  Dtypes  [3] = MPI_INT; Dtypelens  [3] = Numptrn ;
  MPI_Address(&(constpat[0]), &(Dtypeaddr[3]));
  Dtypes  [4] = MPI_DOUBLE; Dtypelens  [4] = numcats ;
  MPI_Address(&(Rates[0]), &(Dtypeaddr[4]));
  Dtypes  [5] = MPI_DOUBLE; Dtypelens  [5] = tpmradix ;
  MPI_Address(&(Eval[0]), &(Dtypeaddr[5]));
  Dtypes  [6] = MPI_DOUBLE; Dtypelens  [6] = tpmradix ;
  MPI_Address(&(Freqtpm[0]), &(Dtypeaddr[6]));
  Dtypes  [7] = MPI_DOUBLE; Dtypelens  [7] = tpmradix * tpmradix ;
  MPI_Address(&(Evec[0][0]), &(Dtypeaddr[7]));
  Dtypes  [8] = MPI_DOUBLE; Dtypelens  [8] = tpmradix * tpmradix ;
  MPI_Address(&(Ievc[0][0]), &(Dtypeaddr[8]));
  Dtypes  [9] = MPI_DOUBLE; Dtypelens  [9] = tpmradix * tpmradix ;
  MPI_Address(&(iexp[0][0]), &(Dtypeaddr [9]));
  Dtypes [10] = MPI_DOUBLE; Dtypelens [10] = Maxspc * Maxspc ;
  MPI_Address(&(Distanmat[0][0]), &(Dtypeaddr[10]));
  Dtypes [11] = MPI_DOUBLE; Dtypelens [11] = numcats * tpmradix * tpmradix ;
  MPI_Address(&(ltprobr[0][0][0]), &(Dtypeaddr[11]));
 
  MPI_Type_struct(12, Dtypelens, Dtypeaddr, Dtypes, &PP_Data);
  MPI_Type_commit(&PP_Data);
 
  for (dest=1; dest<PP_NumProcs; dest++) {
 
    error = MPI_Ssend(MPI_BOTTOM, 1, PP_Data, dest, PP_DATA, PP_Comm);
    if (error != MPI_SUCCESS)
      PP_Printerror(STDOUT, 1100+PP_Myid, error);
 
#     ifdef PVERBOSE3
         FPRINTF(STDOUTFILE "(%2d) -> (%2d) Sent Data\n", PP_Myid, dest);
#     endif /* PVERBOSE2 */
 
  } /* for each slave */
 
  MPI_Type_free(&PP_Data);

# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent Data\n", PP_Myid);
# endif /* PVERBOSE2 */
 
} /* PP_SendData */


/**************************************************************************
*  procedures to send the request to compute a single quartet  (M->S)     *
**************************************************************************/
 
void PP_SendDoQuart(int    dest,
                    int    a, 
                    int    b,
                    int    c,
                    int    d,
                    int    approx)
{
# define NUMINT 5
  int    ints[NUMINT];
  int    error;
 
  ints[0] = a;
  ints[1] = b;
  ints[2] = c;
  ints[3] = d;
  ints[4] = approx;

  PP_doquartsent++;
  PP_doquartsentn++;

# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending  -> (%2d): Quart(%d,%d,%d,%d)\n", PP_Myid, dest, a, b, c, d);
# endif /* PVERBOSE2 */
 
  error = MPI_Ssend(ints, NUMINT, MPI_INT, dest, PP_DOQUART, PP_Comm);
  if (error != MPI_SUCCESS)
     PP_Printerror(STDOUT, PP_Myid, error);
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent \n", PP_Myid);
# endif /* PVERBOSE3 */
# undef NUMINT
 
} /* PP_SendDoQuart */



/******************/

void PP_RecvDoQuart(int    *a, 
                    int    *b,
                    int    *c,
                    int    *d,
                    int    *approx)
{
# define NUMINT 5
  int    ints[NUMINT];
  int    error;
  MPI_Status   stat;
  PP_doquartrecved++;
  PP_doquartrecvedn++;

# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) Receiving: Quart\n", PP_Myid);
# endif /* PVERBOSE3 */
 
  error = MPI_Recv(ints, NUMINT, MPI_INT, PP_MyMaster, PP_DOQUART, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 200+PP_Myid, error);
 
  *a      = ints[0];
  *b      = ints[1];
  *c      = ints[2];
  *d      = ints[3];
  *approx = ints[4];
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Received: Quart(%d,%d,%d,%d,%c)\n", PP_Myid, *a, *b, *c, *d, (approx ? 'A' : 'E'));
# endif /* PVERBOSE2 */
# undef NUMINT
 
} /* PP_RecvDoQuart */
 
 
/**************************************************************************
*  procedures to send the result of a single quartet  (S->M)              *
**************************************************************************/
 
void PP_SendQuart(int    a, 
                  int    b,
                  int    c,
                  int    d,
                  double d1,
                  double d2,
                  double d3,
                  int    approx)
{
# define NUMINT 5
# define NUMDBL 3
  int    ints[NUMINT];
  double doubles[NUMDBL];
  MPI_Datatype Dtypes[2] =    {MPI_INT, MPI_DOUBLE};
  int          Dtypelens[2] = {NUMINT , NUMDBL};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_Quart;
  int          error;
 
  PP_quartsent++;
  PP_quartsentn++;
  ints[0] = a;
  ints[1] = b;
  ints[2] = c;
  ints[3] = d;
  ints[4] = approx;
  doubles[0] = d1;
  doubles[1] = d2;
  doubles[2] = d3;
 
  MPI_Address(ints,     Dtypeaddr);
  MPI_Address(doubles, (Dtypeaddr+1));
  
  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_Quart);
  MPI_Type_commit(&PP_Quart);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: Quart(%d,%d,%d,%d) = (%.3f, %.3f, %.3f)\n", PP_Myid, a, b, c, d, d1, d2, d3);
# endif /* PVERBOSE2 */
 
  error = MPI_Ssend(MPI_BOTTOM, 1, PP_Quart, PP_MyMaster, PP_QUART, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 400+PP_Myid, error);
 
  MPI_Type_free(&PP_Quart);

# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent \n", PP_Myid);
# endif /* PVERBOSE3 */
# undef NUMINT
# undef NUMDBL

} /* PP_SendQuart */
 
 
 
/******************/
 
void PP_RecvQuart(int    *a, 
                  int    *b,
                  int    *c,
                  int    *d,
                  double *d1,
                  double *d2,
                  double *d3,
                  int    *approx)
{
# define NUMINT 5
# define NUMDBL 3
  int          ints[NUMINT];
  double       doubles[NUMDBL];
  MPI_Datatype Dtypes[2] =    {MPI_INT, MPI_DOUBLE};
  int          Dtypelens[2] = {NUMINT , NUMDBL};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_Quart;
  int          error;
  MPI_Status   stat;
 
  PP_quartrecved++;
  PP_quartrecvedn++;
  MPI_Address(ints,     Dtypeaddr);
  MPI_Address(doubles, (Dtypeaddr+1));
  
  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_Quart);
  MPI_Type_commit(&PP_Quart);
 
  error = MPI_Recv(MPI_BOTTOM, 1, PP_Quart, MPI_ANY_SOURCE, PP_QUART, PP_Comm, &stat);

  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 500+PP_Myid, error);

  PP_putslave(stat.MPI_SOURCE);
 
  *a      = ints[0];
  *b      = ints[1];
  *c      = ints[2];
  *d      = ints[3];
  *d1     = doubles[0];
  *d2     = doubles[1];
  *d3     = doubles[2];
  *approx = ints[4];

  MPI_Type_free(&PP_Quart);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Received <- (%2d): Quart(%d,%d,%d,%d)=(%.3f, %.3f, %.3f)\n", PP_Myid, stat.MPI_SOURCE, *a, *b, *c, *d, *d1, *d2, *d3);
# endif /* PVERBOSE2 */
# undef NUMINT
# undef NUMDBL
 
} /* PP_RecvQuart */


 
/**************************************************************************
*  procedures to send the request to compute a block of quartets  (M->S)  *
**************************************************************************/

void PP_SendDoQuartBlock(int dest, uli firstq, uli amount, int approx)
{
#	define NUMULI 3
	uli           ulongs[NUMULI];
	int           error;
 
        PP_doquartsent += amount;
        PP_doquartsentn++;
#	ifdef PVERBOSE2
	   FPRINTF(STDOUTFILE "(%2d) Sending: DOQuartBlock Signal\n", PP_Myid);
#	endif /* PVERBOSE2 */
 
	ulongs[0] = firstq;
	ulongs[1] = amount;
	ulongs[2] = (uli)approx;

	error = MPI_Ssend(ulongs, NUMULI, MPI_UNSIGNED_LONG, dest, PP_DOQUARTBLOCK, PP_Comm);
	if (error != MPI_SUCCESS) PP_Printerror(STDOUT, 2100+PP_Myid, error);
 
#	ifdef PVERBOSE3
	   FPRINTF(STDOUTFILE "(%2d) ... Sent DOQuartBlock Signal (addr:%ld, num:%ld)\n", PP_Myid, firstq, amount);
#	endif /* PVERBOSE3 */
#       undef NUMULI
 
} /* PP_SendDoQuartBlock */

/******************/

void PP_RecvDoQuartBlock(uli *firstq, uli *amount, uli **bq, int *approx)
{
#	define NUMULI 3
	uli           ulongs[NUMULI];
	MPI_Status    stat;
	int           error;
 
#	ifdef PVERBOSE2
	  FPRINTF(STDOUTFILE "(%2d) Receiving: DOQuartBlock Signal\n", PP_Myid);
#	endif /* PVERBOSE2 */
 
	error = MPI_Recv(&ulongs, NUMULI, MPI_UNSIGNED_LONG, PP_MyMaster, PP_DOQUARTBLOCK, PP_Comm, &stat);
	if (error != MPI_SUCCESS)
	   PP_Printerror(STDOUT, 2100+PP_Myid, error);

	*firstq=ulongs[0];
	*amount=ulongs[1];
	*approx= (int)ulongs[2];

	*bq = malloc((unsigned)*amount * sizeof(uli));

        PP_doquartrecved += *amount;
        PP_doquartrecvedn++;

#	ifdef PVERBOSE3
	   FPRINTF(STDOUTFILE "(%2d) ... DOQuartBlock (addr:%ld, num:%ld)\n", 
	                      PP_Myid, *firstq, *amount);
#	endif /* PVERBOSE3 */
 
#	undef NUMULI
} /* PP_RecvDoQuartBlock */

/*********************************************************************
*  procedures to send the results of a block of quartets  (S->M)     *
*********************************************************************/

void PP_SendQuartBlock(uli startq,
                       uli numofq,
                       unsigned char *quartetinfo,
                       uli  numofbq,
                       uli *bq,
                       int approx)
{
# define NUMULI 3
# define NUMINT 1
  unsigned char *trueaddr;
  uli            truenum;
  int            error;
  int    ints[NUMINT];
  uli    ulis[NUMULI];
  MPI_Datatype Dtypes[2] =    {MPI_UNSIGNED_LONG, MPI_INT};
  int          Dtypelens[2] = {NUMULI,            NUMINT};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_QBlockSpecs;
  MPI_Datatype DtypesRes[2] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_LONG};
  int          DtypelensRes[2];
  MPI_Aint     DtypeaddrRes[2];
  MPI_Datatype PP_QBlockRes;

/*
  uli *bq;
  uli  numofbq;
*/

  PP_quartsent += numofq;
  PP_quartsentn++;

  truenum = (uli)((numofq+1)/2);
  trueaddr = (unsigned char *)(quartetinfo + (uli)(startq/2));

# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: startq=%lud numofq=%lud\n", PP_Myid, startq, numofq);
    FPRINTF(STDOUTFILE "(%2d)          approx=%c\n", PP_Myid, (approx ? 'A' : 'E'));
# endif /* PVERBOSE2 */

  ints[0] = approx;
  ulis[0] = startq;
  ulis[1] = numofq;
  ulis[2] = numofbq;

  MPI_Address(ulis,  Dtypeaddr);
  MPI_Address(ints, (Dtypeaddr+1));
  
  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_QBlockSpecs);
  MPI_Type_commit(&PP_QBlockSpecs);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: xxPP_QuartBlockSpecs(0,%lu)=%d,%d\n", PP_Myid, truenum-1, trueaddr[0], trueaddr[truenum-1]);
# endif /* PVERBOSE2 */
 

  error = MPI_Ssend(MPI_BOTTOM, 1, PP_QBlockSpecs, PP_MyMaster, PP_QUARTBLOCKSPECS, PP_Comm);
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent QuartBlockSpecs (%ld, %ld, %ld, %d)\n", PP_Myid, ulis[0], ulis[1], ulis[2], ints[0]);
# endif /* PVERBOSE3 */

  MPI_Address(trueaddr, DtypeaddrRes);
  DtypelensRes[0] = truenum;

  MPI_Address(bq, (DtypeaddrRes + 1));
  DtypelensRes[1] = numofbq;
  MPI_Type_struct(2, DtypelensRes, DtypeaddrRes, DtypesRes, &PP_QBlockRes);
  MPI_Type_commit(&PP_QBlockRes);

  error = MPI_Ssend(MPI_BOTTOM, 1, PP_QBlockRes, PP_MyMaster, PP_QUARTBLOCK, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, PP_Myid, error);
 
  MPI_Type_free(&PP_QBlockSpecs);
  MPI_Type_free(&PP_QBlockRes);
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent xxPP_QuartBlock(0,%lu)=%d,%d\n", PP_Myid, truenum-1, trueaddr[0], trueaddr[truenum-1]);
# endif /* PVERBOSE3 */
 
# undef NUMULI
# undef NUMINT
} /* PP_SendQuartBlock */
 
 

/******************/

void PP_RecvQuartBlock(int slave,
                       uli *startq,
                       uli *numofq,
                       unsigned char *quartetinfo,
                       int *approx)
{
#	define NUMULI 3
#	define NUMINT 1
	unsigned char *trueaddr;
	uli            truenum;
	int            error;
	int            dest;
	int    ints[NUMINT];
	uli    ulis[NUMULI];
	MPI_Datatype Dtypes[2] =    {MPI_UNSIGNED_LONG, MPI_INT};
	int          Dtypelens[2] = {NUMULI,             NUMINT};
	MPI_Aint     Dtypeaddr[2];
	MPI_Datatype PP_QBlockSpecs;
	MPI_Datatype DtypesRes[2] =    {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_LONG};
	int          DtypelensRes[2];
	MPI_Aint     DtypeaddrRes[2];
	MPI_Datatype PP_QBlockRes;
	MPI_Status   stat;
	uli          count;
uli num;
uli *numofbq;
uli *bq;
numofbq=&num;
#	ifdef PVERBOSE3
	    FPRINTF(STDOUTFILE "(%2d) Receiving QuartBlock ...\n", PP_Myid);
#	endif /* PVERBOSE3 */
	MPI_Address(ulis,  Dtypeaddr);
	MPI_Address(ints, (Dtypeaddr+1));
 
	MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_QBlockSpecs);
	MPI_Type_commit(&PP_QBlockSpecs);
 
	MPI_Probe(MPI_ANY_SOURCE, PP_QUARTBLOCKSPECS, PP_Comm, &stat);
	dest = stat.MPI_SOURCE;
	error = MPI_Recv(MPI_BOTTOM, 1, PP_QBlockSpecs, dest, PP_QUARTBLOCKSPECS, PP_Comm, &stat);
	if (error != MPI_SUCCESS)
	      PP_Printerror(STDOUT, PP_Myid, error);

	*approx  = ints[0];
	*startq  = ulis[0];
	*numofq  = ulis[1];
	*numofbq = ulis[2];

        PP_quartrecved += *numofq;
        PP_quartrecvedn++;
	truenum = (uli)((*numofq+1)/2);
	trueaddr = (unsigned char *)(quartetinfo + (uli)(*startq/2));
#       ifdef PVERBOSE3
           FPRINTF(STDOUTFILE "(%2d) ... Recv QuartBlockSpecs (%ld, %ld, %ld, %d)\n", PP_Myid, ulis[0], ulis[1], ulis[2], ints[0]);
#       endif /* PVERBOSE3 */

	DtypelensRes[0] =  truenum;
	MPI_Address(trueaddr,  DtypeaddrRes);

	bq = malloc((unsigned) *numofbq * sizeof(uli));

	DtypelensRes[1] = *numofbq;
	MPI_Address(bq, (DtypeaddrRes+1));
	MPI_Type_struct(2, DtypelensRes, DtypeaddrRes, DtypesRes, &PP_QBlockRes);
	MPI_Type_commit(&PP_QBlockRes);
 
	error = MPI_Recv(MPI_BOTTOM, 1, PP_QBlockRes, dest, PP_QUARTBLOCK, PP_Comm, &stat);
	if (error != MPI_SUCCESS)
	    PP_Printerror(STDOUT, PP_Myid, error);
#       ifdef PVERBOSE3
           FPRINTF(STDOUTFILE "(%2d) ... Recv QuartBlock \n", PP_Myid);
#       endif /* PVERBOSE3 */

        PP_putslave(dest);

	for(count = 0; count < *numofbq; count++){
	  int a, b, c, d;
	  num2quart(bq[count], &a, &b, &c, &d);
#	  ifdef PVERBOSE2
	     FPRINTF(STDOUTFILE "(%2d) %ld. bad quarted (%d, %d, %d, %d) = %ld\n", PP_Myid, count, a, b, c, d, bq[count]);
#	  endif /* PVERBOSE2 */

	  badqs++;
	  badtaxon[a]++;
	  badtaxon[b]++;
	  badtaxon[c]++;
	  badtaxon[d]++;
	  if (show_optn) {
	     fputid10(unresfp, a);
	     fprintf(unresfp, "  ");
	     fputid10(unresfp, b);
	     fprintf(unresfp, "  ");
	     fputid10(unresfp, c);
	     fprintf(unresfp, "  ");
	     fputid(unresfp, d);
	     fprintf(unresfp, "\n");
	  }
	}
	free(bq);
	MPI_Type_free(&PP_QBlockSpecs);
	MPI_Type_free(&PP_QBlockRes);
#	ifdef PVERBOSE2
           FPRINTF(STDOUTFILE "(%2d) <- (%2d) ... Recv xxPP_QuartBlock(0,%lu)=%d,%d\n", PP_Myid, dest, truenum-1, trueaddr[0], trueaddr[truenum-1]);
#	endif /* PVERBOSE2 */
 
#	undef NUMULI
#	undef NUMINT
} /* PP_RecvQuartBlock */
 

/*********************************************************************
*  send/receive array with all quartets  (M->S)                      *
*********************************************************************/

void PP_SendAllQuarts(unsigned long  Numquartets,
                      unsigned char *quartetinfo)
{
  MPI_Datatype  Dtypes[1] =    {MPI_UNSIGNED_CHAR};
  int           Dtypelens[1];
  MPI_Aint      Dtypeaddr[1];
  MPI_Datatype  PP_AllQuarts;
  int           dest;
  int           error;
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: PP_AllQuart(0)=%d\n", PP_Myid, quartetinfo[0]);
# endif /* PVERBOSE2 */
 
  /* compute number of quartets */
  if (Numquartets % 2 == 0) { /* even number */
     Dtypelens[0] = (Numquartets)/2;
  } else { /* odd number */
     Dtypelens[0] = (Numquartets + 1)/2;
  }

  MPI_Address(&(quartetinfo[0]), Dtypeaddr);
  MPI_Type_struct(1, Dtypelens, Dtypeaddr, Dtypes, &PP_AllQuarts);
  MPI_Type_commit(&PP_AllQuarts);
 
  for (dest=1; dest<PP_NumProcs; dest++) {
 
    error = MPI_Ssend(MPI_BOTTOM, 1, PP_AllQuarts, dest, PP_ALLQUARTS, PP_Comm);
    if (error != MPI_SUCCESS)
      PP_Printerror(STDOUT, 1200+PP_Myid, error);
 
#   ifdef PVERBOSE3
       FPRINTF(STDOUTFILE "(%2d) -> (%2d) ... Sent xxAllQuart(0,%d)=%d,%d (%luq -> %db)\n", 
                           PP_Myid, dest, Dtypelens[0]-1, quartetinfo[0], quartetinfo[Dtypelens[0]-1],
                           Numquartets, Dtypelens[0]-1);
#   endif /* PVERBOSE3 */
  } /* for each slave */
 
  MPI_Type_free(&PP_AllQuarts);
 
 
} /* PP_SendAllQuarts */
 
 

/******************/

void PP_RecvAllQuarts(int            taxa,
                      unsigned long *Numquartets,
                      unsigned char *quartetinfo)
{
  MPI_Datatype  Dtypes[1] =    {MPI_UNSIGNED_CHAR};
  int           Dtypelens[1];
  MPI_Aint      Dtypeaddr[1];
  MPI_Datatype  PP_AllQuarts;
  MPI_Status    stat;
  int           error;
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) Receiving AllQuarts ...\n", PP_Myid);
# endif /* PVERBOSE3 */
 
  /* compute number of quartets */
  *Numquartets = (uli) taxa*(taxa-1)*(taxa-2)*(taxa-3)/24;
  if (*Numquartets % 2 == 0) { /* even number */
    Dtypelens[0] = (*Numquartets)/2;
  } else { /* odd number */
    Dtypelens[0] = (*Numquartets + 1)/2;
  }
 
  MPI_Address(&(quartetinfo[0]), Dtypeaddr);
  MPI_Type_struct(1, Dtypelens, Dtypeaddr, Dtypes, &PP_AllQuarts);
  MPI_Type_commit(&PP_AllQuarts);
 
  error = MPI_Recv(MPI_BOTTOM, 1, PP_AllQuarts, PP_MyMaster, PP_ALLQUARTS, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1300+PP_Myid, error);

  MPI_Type_free(&PP_AllQuarts);
 
# ifdef PVERBOSE2
     FPRINTF(STDOUTFILE "(%2d) <- (%2d) ... Recv xxAllQuart(0,%d)=%d,%d (%luq -> %db)\n", 
                        PP_Myid, PP_MyMaster, Dtypelens[0]-1, quartetinfo[0], quartetinfo[Dtypelens[0]-1],
                        *Numquartets, Dtypelens[0]-1);
# endif /* PVERBOSE2 */
 
} /* PP_RecvAllQuarts */
 


/*********************************************************************
*  procedures to send request for a single puzzle tree               *
*********************************************************************/

void PP_SendPermut(int      dest,
                   int      taxa, 
                   ivector  permut)
{
  MPI_Datatype  Dtypes[1] =    {MPI_INT};
  int           Dtypelens[1];
  MPI_Aint      Dtypeaddr[1];
  MPI_Datatype  PP_Permut;
  int           error;
 
  PP_permutsent++;
  PP_permutsentn++;
  Dtypelens[0] = taxa;
 
  MPI_Address(&(permut[0]), Dtypeaddr);
  MPI_Type_struct(1, Dtypelens, Dtypeaddr, Dtypes, &PP_Permut);
  MPI_Type_commit(&PP_Permut);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending  -> (%2d): PP_Permut(0)=%d\n", PP_Myid, dest, permut[0]);
# endif /* PVERBOSE2 */
 
  error = MPI_Ssend(MPI_BOTTOM, 1, PP_Permut, dest, PP_DOPUZZLE, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1500+PP_Myid, error);
 
  MPI_Type_free(&PP_Permut);
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent PP_Permut\n", PP_Myid);
# endif /* PVERBOSE3 */
 
} /* PP_SendPermut */
 
/******************/
 
void PP_RecvPermut(int      taxa,
                   ivector  permut)
{
  MPI_Datatype  Dtypes[1] =    {MPI_INT};
  int           Dtypelens[1];
  MPI_Aint      Dtypeaddr[1];
  MPI_Datatype  PP_Permut;
  MPI_Status    stat;
  int           error;
  
  PP_permutrecved++;
  PP_permutrecvedn++;
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) Receiving: PP_Permut\n", PP_Myid);
# endif /* PVERBOSE3 */
 
  Dtypelens[0] = taxa;
 
  MPI_Address(&(permut[0]), Dtypeaddr);
  MPI_Type_struct(1, Dtypelens, Dtypeaddr, Dtypes, &PP_Permut);
  MPI_Type_commit(&PP_Permut);
 
  error = MPI_Recv(MPI_BOTTOM, 1, PP_Permut, PP_MyMaster, PP_DOPUZZLE, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1700+PP_Myid, error);
 
  MPI_Type_free(&PP_Permut);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Received: PP_Permut(0)=%d\n", PP_Myid, permut[0]);
# endif /* PVERBOSE2 */
 
} /* PP_RecvPermut */

/*********************************************************************
*  procedures to send the splits of a puzzle tree to the master      *
*********************************************************************/

void PP_SendSplitsBlock(int               taxa, 
                        uli               blocksize,
                        cmatrix          *biparts,
                        int               pstnum,
                        treelistitemtype *pstlist)
{
  MPI_Datatype     *Dtypes;
  int              *Dtypelens;
  MPI_Aint         *Dtypeaddr;
  MPI_Datatype      PP_Biparts;
  int               error;
  int               n;
  int               ints[3];
  int              *pstnumarr;
  treelistitemtype *pstptr;

  PP_splitsent+=blocksize;
  PP_splitsentn++;
 
  ints[0] = taxa;
  ints[1] = (int) blocksize;
  ints[2] = pstnum;
  error = MPI_Ssend(ints, 3, MPI_INT, PP_MyMaster, PP_PUZZLEBLOCKSPECS, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1800+PP_Myid, error);
 
  Dtypes    = malloc((blocksize + pstnum + 1) * sizeof(MPI_Datatype));
  Dtypelens = malloc((blocksize + pstnum + 1) * sizeof(int));
  Dtypeaddr = malloc((blocksize + pstnum + 1) * sizeof(MPI_Aint));
  pstnumarr = malloc(pstnum * sizeof(int));

# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: PP_bipartsblock(0..%lu,0,0)8=\"%c\"\n", PP_Myid, blocksize, biparts[0][0][0]);
# endif /* PVERBOSE2 */
 
  for (n=0; n<(int)blocksize; n++) {
    Dtypes[n]    = MPI_CHAR;
    Dtypelens[n] = (taxa - 3) * taxa;
    MPI_Address(&(biparts[n][0][0]), &(Dtypeaddr[n]));
  }
  pstptr = pstlist;
  for (n=0; n<pstnum; n++) {
    Dtypes[(int)blocksize + n]    = MPI_CHAR;
    Dtypelens[(int)blocksize + n] = psteptreestrlen;
    MPI_Address((*pstptr).tree, &(Dtypeaddr[(int)blocksize + n]));
    pstnumarr[n] = (*pstptr).count;
#   ifdef PVERBOSE3
       FPRINTF(STDOUTFILE "(%2d) Sent tree item ->%d: [%d/%d] #=%d  \"%s\"\n",
               PP_Myid, PP_MyMaster, n, pstnum, pstnumarr[n], (*pstptr).tree);
#   endif /* PVERBOSE3 */
    pstptr = (*pstptr).succ;
  }
  Dtypes[((int)blocksize + pstnum)]    = MPI_INT;
  Dtypelens[((int)blocksize + pstnum)] = pstnum;
  MPI_Address(&(pstnumarr[0]), &(Dtypeaddr[((int)blocksize + pstnum)]));

  MPI_Type_struct(((int)blocksize + pstnum + 1), Dtypelens, Dtypeaddr, Dtypes, &PP_Biparts);
  MPI_Type_commit(&PP_Biparts);
 
  error = MPI_Ssend(MPI_BOTTOM, 1, PP_Biparts, PP_MyMaster, PP_PUZZLEBLOCK, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1800+PP_Myid, error);
 
  MPI_Type_free(&PP_Biparts);
  free(Dtypes);
  free(Dtypelens);
  free(Dtypeaddr);
  free(pstnumarr);
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent PP_bipartsblock\n", PP_Myid);
# endif /* PVERBOSE3 */
 
} /* PP_SendSplitsBlock */
 
/******************/

void PP_RecvSplitsBlock(int               *taxa, 
                        uli               *blocksize,
                        cmatrix          **bip,
                        treelistitemtype **pstlist,
                        int               *pstnum,
                        int               *pstsum)
/* bp -> (*bip) */
{
  MPI_Datatype      *Dtypes;
  int               *Dtypelens;
  MPI_Aint          *Dtypeaddr;
  MPI_Datatype       PP_Biparts;
  MPI_Status         stat;
  int                error;
  int                n;
  int                dest;
  int                ints[3];
  int                pstlistnum;
  int                tmpnum;
  int                tmpsum;
  int               *pstnumarr;
  char             **pstarr;
  treelistitemtype  *treeitem;

  error = MPI_Recv(ints, 3, MPI_INT, MPI_ANY_SOURCE, PP_PUZZLEBLOCKSPECS, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1900+PP_Myid, error);

    dest       =       stat.MPI_SOURCE; 
   *taxa       =       ints[0];
   *blocksize  = (uli) ints[1];
    pstlistnum =       ints[2];
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Received<-%d: PP_bipartsblockspec(t=%d,b=%ld,p=%d)\n", PP_Myid, dest, *taxa, *blocksize, pstlistnum);
# endif /* PVERBOSE2 */

  PP_splitrecved += *blocksize;
  PP_splitrecvedn++;
 
  Dtypes    = malloc((*blocksize + pstlistnum + 1) * sizeof(MPI_Datatype));
  Dtypelens = malloc((*blocksize + pstlistnum + 1) * sizeof(int));
  Dtypeaddr = malloc((*blocksize + pstlistnum + 1) * sizeof(MPI_Aint));
  (*bip)    = (cmatrix *) malloc(*blocksize * sizeof(void *));
  pstnumarr = (int *) malloc(pstlistnum * sizeof(int));
  pstarr    = (char **) malloc(pstlistnum * sizeof(char *));
  
/*  pstarr[0] = (char *) malloc(psteptreestrlen * pstlistnum * sizeof(char)); 
    for(n=1; n<pstlistnum; n++)
    pstarr[n] = &(pstarr[0][n * psteptreestrlen]);
*/

  for (n=0; n<(int)*blocksize; n++) {
    (*bip)[n]        = new_cmatrix(*taxa - 3, *taxa);
    Dtypes[n]    = MPI_CHAR;
    Dtypelens[n] = (*taxa - 3) * *taxa;
    MPI_Address(&((*bip)[n][0][0]), &(Dtypeaddr[n]));
  }
  for (n=0; n<pstlistnum; n++) {
    pstarr[n]        = (char *)malloc(psteptreestrlen * sizeof(char)); 
    Dtypes[(int)*blocksize + n]    = MPI_CHAR;
    Dtypelens[(int)*blocksize + n] = psteptreestrlen;
    MPI_Address(&(pstarr[n][0]), &(Dtypeaddr[(int)*blocksize + n]));
  }
  
  Dtypes[(int)*blocksize + pstlistnum]    = MPI_INT;
  Dtypelens[(int)*blocksize + pstlistnum] = pstlistnum;
  MPI_Address(&(pstnumarr[0]), &(Dtypeaddr[(int)*blocksize + pstlistnum]));

  MPI_Type_struct(((int)*blocksize + pstlistnum + 1), Dtypelens, Dtypeaddr, Dtypes, &PP_Biparts);
  MPI_Type_commit(&PP_Biparts);
 
  error = MPI_Recv(MPI_BOTTOM, 1, PP_Biparts, dest, PP_PUZZLEBLOCK, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1910+PP_Myid, error);

  tmpsum = *pstsum;
  for (n=0; n<pstlistnum; n++) tmpsum += pstnumarr[n];
  tmpnum = *pstnum;

  for (n=0; n<pstlistnum; n++) {
#   ifdef PVERBOSE3
       FPRINTF(STDOUTFILE "(%2d) Received tree item <-%d: [%d/%d] #=%d  \"%s\"\n",
               PP_Myid, dest, n, pstlistnum, pstnumarr[n], pstarr[n]);
#   endif /* PVERBOSE3 */

    treeitem = addtree2list(&(pstarr[n]), pstnumarr[n], pstlist, pstnum, pstsum);
    if((listqptrees == PSTOUT_LIST) 
      || (listqptrees == PSTOUT_LISTORDER)) {
       /* print: order no/# topol per this id/tree id/sum of unique topologies/sum of trees so far */
       fprintf(qptlist, "%d.\t%d\t%d\t%d\t%d\t%d\n", 
              PP_splitrecvedn, pstnumarr[n], (*treeitem).count, (*treeitem).id, *pstnum, tmpsum);
    }


  }

  MPI_Type_free(&PP_Biparts);
  free(Dtypes);
  free(Dtypelens);
  free(Dtypeaddr);
  free(pstnumarr);
  free(pstarr);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Received<-%d: PP_bipartsblock(0..%lu,0,0):=%c\n", PP_Myid, dest, *blocksize, (*bip)[0][0][0]);
# endif /* PVERBOSE2 */
 
  PP_putslave(dest); 
 
  /* MPI_Type_free(&PP_Biparts); */

} /* PP_RecvSplitsBlock */
 
/******************/

void PP_SendSplits(int      taxa, 
                   cmatrix  biparts)
{
  MPI_Datatype  Dtypes[1] =    {MPI_CHAR};
  int           Dtypelens[1];
  MPI_Aint      Dtypeaddr[1];
  MPI_Datatype  PP_Biparts;
  int           error;

  PP_splitsent++;
  PP_splitsentn++;
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: PP_biparts(0,0)=%c\n", PP_Myid, biparts[0][0]);
# endif /* PVERBOSE2 */
 
  Dtypelens[0] = (taxa - 3) * taxa;

  MPI_Address(&(biparts[0][0]), Dtypeaddr);
  MPI_Type_struct(1, Dtypelens, Dtypeaddr, Dtypes, &PP_Biparts);
  MPI_Type_commit(&PP_Biparts);
 
  error = MPI_Ssend(MPI_BOTTOM, 1, PP_Biparts, PP_MyMaster, PP_PUZZLE, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1800+PP_Myid, error);
 
  MPI_Type_free(&PP_Biparts);
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent PP_biparts\n", PP_Myid);
# endif /* PVERBOSE3 */
 
} /* PP_SendSplits */
 
/******************/

void PP_RecvSplits(int      taxa,
                   cmatrix  biparts)
{
  MPI_Datatype  Dtypes[1] =    {MPI_CHAR};
  int           Dtypelens[1];
  MPI_Aint      Dtypeaddr[1];
  MPI_Datatype  PP_Biparts;
  MPI_Status    stat;
  int           error;

  PP_splitrecved++;
  PP_splitrecvedn++;
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) Receiving: PP_biparts\n", PP_Myid);
# endif /* PVERBOSE3 */
 
  Dtypelens[0] = (taxa - 3) * taxa;
 
  MPI_Address(&(biparts[0][0]), Dtypeaddr);
  MPI_Type_struct(1, Dtypelens, Dtypeaddr, Dtypes, &PP_Biparts);
  MPI_Type_commit(&PP_Biparts);
 
  error = MPI_Recv(MPI_BOTTOM, 1, PP_Biparts, MPI_ANY_SOURCE, PP_PUZZLE, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 1920+PP_Myid, error);

  PP_putslave(stat.MPI_SOURCE); 
 
  MPI_Type_free(&PP_Biparts);
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Received <- (%2d): PP_biparts(0,0)=%c\n", PP_Myid, stat.MPI_SOURCE, biparts[0][0]);
# endif /* PVERBOSE2 */
 
} /* PP_RecvSplits */
 
/*********************************************************************
*  procedures to send the request to compute puzzle tree             *
*********************************************************************/

void PP_SendDoPermutBlock(uli puzzlings)
{
  int           dest;
  int           error;
  uli           numpuzzlings;
  uli           puzzlingssent = 0;
  schedtype     sched; 
  int           bipnum;
  int           tx;
  uli           bs;
  cmatrix      *bp;


 
  initsched(&sched, puzzlings, PP_NumProcs-1, 4);
  /* numpuzzlings = sc(&sched); */
  numpuzzlings = sgss(&sched);
  while (numpuzzlings > 0) {
     if (PP_emptyslave()) {
        PP_RecvSplitsBlock(&tx, &bs, &bp, &psteptreelist, &psteptreenum, &psteptreesum);
        for (bipnum=0; bipnum<bs; bipnum++) {
               inittree();
               biparts = bp[bipnum];
               makenewsplitentries();
               freetree();
               free_cmatrix(bp[bipnum]);
        }
	free(bp);
        puzzlingssent -= bs;
     }

     dest = PP_getslave();

#    ifdef PVERBOSE2
        FPRINTF(STDOUTFILE "(%2d) Sending: DOPuzzleBlock Signal\n", PP_Myid);
#    endif /* PVERBOSE2 */

     error = MPI_Ssend(&numpuzzlings, 1, MPI_UNSIGNED_LONG, dest, PP_DOPUZZLEBLOCK, PP_Comm);
     if (error != MPI_SUCCESS)
        PP_Printerror(STDOUT, 2100+PP_Myid, error);

#    ifdef PVERBOSE3
        FPRINTF(STDOUTFILE "(%2d) ... Sent DOPuzzleBlock Signal\n", PP_Myid);
#    endif /* PVERBOSE3 */

     puzzlingssent += numpuzzlings;
     PP_permutsent += numpuzzlings;
     PP_permutsentn ++;
 
     numpuzzlings = sgss(&sched);

  } /* while */
 
  while (puzzlingssent > 0) {
        PP_RecvSplitsBlock(&tx, &bs, &bp, &psteptreelist, &psteptreenum, &psteptreesum);
        for (bipnum=0; bipnum<bs; bipnum++) {
               inittree();
               biparts = bp[bipnum];
               makenewsplitentries();
               freetree();
               free_cmatrix(bp[bipnum]);
        }
	free(bp);
        puzzlingssent -= bs;
  }
 
} /* PP_SendDoPermutBlock */

/******************/


void PP_RecvDoPermutBlock(uli *taxa)
{
  MPI_Status    stat;
  int           error;
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Receiving: DOPuzzleBlock Signal\n", PP_Myid);
# endif /* PVERBOSE2 */
 
  error = MPI_Recv(taxa, 1, MPI_UNSIGNED_LONG, PP_MyMaster, PP_DOPUZZLEBLOCK, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 2100+PP_Myid, error);

  PP_permutrecved += *taxa;
  PP_permutrecvedn ++;

# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... DOPuzzleBlock Signal\n", PP_Myid);
# endif /* PVERBOSE3 */
 
} /* PP_RecvDoPermutBlock */

/*********************************************************************
*  procedures to make the slaves to finalize                         *
*********************************************************************/

void PP_SendDone()
{
# define NUMINT 8
# define NUMDBL 6
  int           dest;
  int           error;
  MPI_Status    stat;
  int    ints[NUMINT];
  double doubles[NUMDBL];
  MPI_Datatype Dtypes[2] =    {MPI_INT, MPI_DOUBLE};
  int          Dtypelens[2] = {NUMINT , NUMDBL};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_Stats;
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Sending: DONE Signal\n", PP_Myid);
# endif /* PVERBOSE2 */
 
  for (dest=1; dest<PP_NumProcs; dest++) {
 
    error = MPI_Ssend(&dest, 0, MPI_INT, dest, PP_DONE, PP_Comm);
    if (error != MPI_SUCCESS)
      PP_Printerror(STDOUT, 2100+PP_Myid, error);
 
  } /* for each slave */
 
# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... Sent DONE Signal\n", PP_Myid);
# endif /* PVERBOSE3 */
 
  MPI_Address(ints,     Dtypeaddr);
  MPI_Address(doubles, (Dtypeaddr+1));

  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_Stats);
  MPI_Type_commit(&PP_Stats);

  doquartrecved[0]  = 0;
  doquartrecvedn[0] = 0;
  quartsent[0]      = 0;
  quartsentn[0]     = 0;
  permutrecved[0]   = 0;
  permutrecvedn[0]  = 0;
  splitsent[0]      = 0;
  splitsentn[0]     = 0;
  cputimes[0]       = 0.0;
  walltimes[0]      = 0.0;
  fullcputimes[0]   = 0.0;
  fullwalltimes[0]  = 0.0;
  altcputimes[0]    = 0.0;
  altwalltimes[0]   = 0.0;

  for (dest=1; dest<PP_NumProcs; dest++) {

     error = MPI_Recv(MPI_BOTTOM, 1, PP_Stats, dest, PP_STATS, PP_Comm, &stat);
     if (error != MPI_SUCCESS)
       PP_Printerror(STDOUT, 2100+PP_Myid, error);

     doquartrecved[dest]  = ints[0];
     doquartrecvedn[dest] = ints[1];
     quartsent[dest]      = ints[2];
     quartsentn[dest]     = ints[3];
     permutrecved[dest]   = ints[4];
     permutrecvedn[dest]  = ints[5];
     splitsent[dest]      = ints[6];
     splitsentn[dest]     = ints[7];
     cputimes[dest]       = doubles[0];
     walltimes[dest]      = doubles[1];
     fullcputimes[dest]   = doubles[2];
     fullwalltimes[dest]  = doubles[3];
     altcputimes[dest]    = doubles[4];
     altwalltimes[dest]   = doubles[5];

#    ifdef PVERBOSE1
        FPRINTF(STDOUTFILE "(%2d) ... Stats received (%d/%d, %d/%d, %d/%d, %d/%d, %f, %f, %f, %f, %f, %f)\n", 
        PP_Myid, doquartrecved[dest], doquartrecvedn[dest], quartsent[dest], quartsentn[dest],
                 permutrecved[dest], permutrecvedn[dest], splitsent[dest], splitsentn[dest],
                 cputimes[dest], walltimes[dest], fullcputimes[dest], fullwalltimes[dest], 
                 altcputimes[dest], altwalltimes[dest]);
#    endif /* PVERBOSE1 */
  } /* for each slave */
 
  MPI_Type_free(&PP_Stats);

# undef NUMINT
# undef NUMDBL
} /* PP_SendDone */

/******************/



void PP_RecvDone()
{
# define NUMINT 8
# define NUMDBL 6

  int           dummy;
  MPI_Status    stat;
  int           error;
  int    ints[NUMINT];
  double doubles[NUMDBL];
  MPI_Datatype Dtypes[2] =    {MPI_INT, MPI_DOUBLE};
  int          Dtypelens[2] = {NUMINT , NUMDBL};
  MPI_Aint     Dtypeaddr[2];
  MPI_Datatype PP_Stats;
 
# ifdef PVERBOSE2
    FPRINTF(STDOUTFILE "(%2d) Receiving: DONE Signal\n", PP_Myid);
# endif /* PVERBOSE2 */
 
  error = MPI_Recv(&dummy, 0, MPI_INT, PP_MyMaster, PP_DONE, PP_Comm, &stat);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 2100+PP_Myid, error);

# ifdef PVERBOSE3
    FPRINTF(STDOUTFILE "(%2d) ... DONE Signal\n", PP_Myid);
# endif /* PVERBOSE3 */
 
  addtimes(OVERALL, &tarr);
# ifdef TIMEDEBUG
  printtimearr(&tarr);
# endif /* TIMEDEBUG */

  time(&walltimestop);
  walltime = difftime(walltimestop, walltimestart);
  cputimestop = clock();
  cputime = ((double) (cputimestop - cputimestart) / CLOCKS_PER_SEC);

# ifdef PVERBOSE1
    FPRINTF(STDOUTFILE "(%2d) total wallclock time: %0.2f sec (= %0.2f min = %0.2f hours)\n",
            PP_Myid, walltime, walltime / 60, walltime / 3600);
    FPRINTF(STDOUTFILE "(%2d) total CPU time: %0.2f sec (= %0.2f min = %0.2f hours)\n",
            PP_Myid, cputime, cputime / 60, cputime / 3600);
# endif /* PVERBOSE1 */

  ints[0] = PP_doquartrecved;
  ints[1] = PP_doquartrecvedn;
  ints[2] = PP_quartsent;
  ints[3] = PP_quartsentn;
  ints[4] = PP_permutrecved;
  ints[5] = PP_permutrecvedn;
  ints[6] = PP_splitsent;
  ints[7] = PP_splitsentn;
  doubles[0] = cputime;
  doubles[1] = walltime;
  doubles[2] = tarr.fullcpu;
  doubles[3] = tarr.fulltime;
  doubles[4] = tarr.cpu;
  doubles[5] = tarr.time;

  MPI_Address(ints,     Dtypeaddr);
  MPI_Address(doubles, (Dtypeaddr+1));

  MPI_Type_struct(2, Dtypelens, Dtypeaddr, Dtypes, &PP_Stats);
  MPI_Type_commit(&PP_Stats);

  error = MPI_Ssend(MPI_BOTTOM, 1, PP_Stats, PP_MyMaster, PP_STATS, PP_Comm);
  if (error != MPI_SUCCESS)
    PP_Printerror(STDOUT, 2100+PP_Myid, error);

  MPI_Type_free(&PP_Stats);

# ifdef PVERBOSE3
        FPRINTF(STDOUTFILE "(%2d) ... Stats sent (%d/%d, %d/%d, %d/%d, %d/%d, %f, %f)\n", 
        PP_Myid, PP_doquartrecved, PP_doquartrecvedn, PP_quartsent, PP_quartsentn,
                 PP_permutrecved, PP_permutrecvedn, PP_splitsent, PP_splitsentn,
                 cputime, walltime, fullcputime, fullwalltime, altcputime, altwalltime);
# endif /* PVERBOSE3 */

# undef NUMINT
# undef NUMDBL
} /* PP_RecvDone */

/*********************************************************************
*  procedures to initialize and cleanup                              *
*********************************************************************/

void PP_Init(int *argc, char **argv[])
{
  MPI_Init(argc, argv);
  PP_Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(PP_Comm, &PP_Myid);
  MPI_Comm_size(PP_Comm, &PP_NumProcs);

  if (PP_NumProcs <= 1) {
	FPRINTF(STDOUTFILE "\nHalt: TREE-PUZZLE needs at least 2 processes for a parallel run!!!\n\n");
	MPI_Finalize();
	exit(1);
  }

  PP_MyMaster = 0;

  if (PP_Myid == PP_MyMaster) {
    PP_IamMaster=1;
    PP_IamSlave =0;
    PP_initslavequeue();

    permutsent     = new_ivector(PP_NumProcs);
    permutrecved   = new_ivector(PP_NumProcs);
    quartsent      = new_ivector(PP_NumProcs);
    quartrecved    = new_ivector(PP_NumProcs);
    doquartsent    = new_ivector(PP_NumProcs);
    doquartrecved  = new_ivector(PP_NumProcs);
    splitsent      = new_ivector(PP_NumProcs);
    splitrecved    = new_ivector(PP_NumProcs);
    permutsentn    = new_ivector(PP_NumProcs);
    permutrecvedn  = new_ivector(PP_NumProcs);
    quartsentn     = new_ivector(PP_NumProcs);
    quartrecvedn   = new_ivector(PP_NumProcs);
    doquartsentn   = new_ivector(PP_NumProcs);
    doquartrecvedn = new_ivector(PP_NumProcs);
    splitsentn     = new_ivector(PP_NumProcs);
    splitrecvedn   = new_ivector(PP_NumProcs);

    walltimes      = new_dvector(PP_NumProcs);
    cputimes       = new_dvector(PP_NumProcs);
    fullwalltimes  = new_dvector(PP_NumProcs);
    fullcputimes   = new_dvector(PP_NumProcs);
    altwalltimes   = new_dvector(PP_NumProcs);
    altcputimes    = new_dvector(PP_NumProcs);
    { int n;
      for (n = 0; n<PP_NumProcs; n++){
        permutsent     [n] = -1; 
	permutrecved   [n] = -1; 
	quartsent      [n] = -1; 
	quartrecved    [n] = -1; 
	doquartsent    [n] = -1; 
	doquartrecved  [n] = -1; 
	splitsent      [n] = -1; 
	splitrecved    [n] = -1; 
	permutsentn    [n] = -1; 
	permutrecvedn  [n] = -1; 
	quartsentn     [n] = -1; 
	quartrecvedn   [n] = -1; 
	doquartsentn   [n] = -1; 
	doquartrecvedn [n] = -1; 
	splitsentn     [n] = -1; 
	splitrecvedn   [n] = -1; 
	walltimes      [n] = -1.0; 
	cputimes       [n] = -1.0; 
	fullwalltimes  [n] = -1.0; 
	fullcputimes   [n] = -1.0; 
	altwalltimes   [n] = -1.0; 
	altcputimes    [n] = -1.0; 
      }
    }
  } else {
    PP_IamMaster=0;
    PP_IamSlave =1;
  }

# ifdef PVERBOSE1
    { char ma[7] = "master";
      char sl[7] = "slave ";
      char *st;
      char PP_ProcName[MPI_MAX_PROCESSOR_NAME] = "empty";
      int  flag;
      int  PP_ProcNameLen = 6;
      int  PP_Clocksync = 0;
      int  PP_Maxtag    = 0;
      int  PP_IO        = 0;
      int  PP_Hostexist = 0;
    if (PP_IamMaster)
      st = ma;
    else
      st = sl;
    FPRINTF(STDOUTFILE "(%2d) Init: ID %d, %d processes running \n", PP_Myid, PP_Myid, PP_NumProcs); 
    FPRINTF(STDOUTFILE "(%2d)       I am %s \n", PP_Myid, st); 
    MPI_Get_processor_name(PP_ProcName, &PP_ProcNameLen);
    FPRINTF(STDOUTFILE "(%2d)       I am running on \"%s\"\n", PP_Myid, PP_ProcName); 
    MPI_Attr_get(PP_Comm, MPI_IO, &PP_IO, &flag);
    if (flag) FPRINTF(STDOUTFILE "(%2d)       next IO Proc=%d - ", PP_Myid, PP_IO); 
    MPI_Attr_get(PP_Comm, MPI_TAG_UB, &PP_Maxtag, &flag);
    if (flag) FPRINTF(STDOUTFILE "MaxTag=%d - ", PP_Maxtag); 
    MPI_Attr_get(PP_Comm, MPI_HOST, &PP_Hostexist, &flag);
    if (flag) {
      if (PP_Hostexist == MPI_PROC_NULL)
        FPRINTF(STDOUTFILE "No HostProc - "); 
      else
        FPRINTF(STDOUTFILE "HostProc=%d - ", PP_Hostexist); 
    }
    MPI_Attr_get(PP_Comm, MPI_WTIME_IS_GLOBAL, &PP_Clocksync, &flag);
    if (PP_Clocksync)
      FPRINTF(STDOUTFILE "global time"); 
    else
      FPRINTF(STDOUTFILE "local time"); 
    FPRINTF(STDOUTFILE "\n");
    }
# endif /* PVERBOSE1 */
} /* PP_Init */

/******************/

void PP_Finalize()
{
  if (PP_IamMaster) {
    free_ivector(freeslaves);

    free_ivector(permutsent);
    free_ivector(permutrecved);
    free_ivector(quartsent);
    free_ivector(quartrecved);
    free_ivector(doquartsent);
    free_ivector(doquartrecved);
    free_ivector(splitsent);
    free_ivector(splitrecved);
    free_ivector(permutsentn);
    free_ivector(permutrecvedn);
    free_ivector(quartsentn);
    free_ivector(quartrecvedn);
    free_ivector(doquartsentn);
    free_ivector(doquartrecvedn);
    free_ivector(splitsentn);
    free_ivector(splitrecvedn);

    free_dvector(walltimes);
    free_dvector(cputimes);
    free_dvector(fullwalltimes);
    free_dvector(fullcputimes);
    free_dvector(altwalltimes);
    free_dvector(altcputimes);
  }

# ifdef PVERBOSE1
    FPRINTF(STDOUTFILE "(%2d) Finished ...\n", PP_Myid);

    FPRINTF(STDOUTFILE "(%2d) doqu(s%6d/%6d,%6d/%6d) qu(s%6d/%6d,%6d/%6d)\n",
            PP_Myid, PP_doquartsent, PP_doquartsentn, PP_doquartrecved, PP_doquartrecvedn, 
                     PP_quartsent, PP_quartsentn, PP_quartrecved, PP_quartrecvedn);
    FPRINTF(STDOUTFILE "(%2d) perm(s%6d/%6d,%6d/%6d) spli(s%6d/%6d,%6d/%6d)\n",
            PP_Myid, PP_permutsent, PP_permutsentn, PP_permutrecved, PP_permutrecvedn, 
                     PP_splitsent, PP_splitsentn, PP_splitrecved, PP_splitrecvedn);
    if (PP_IamMaster) {
       FPRINTF(STDOUTFILE "(%2d)  Init: %2.4f   Param(Comp: %2.4f  Send: %2.4f)\n", PP_Myid, PP_inittime, PP_paramcomptime, PP_paramsendtime);
       FPRINTF(STDOUTFILE "(%2d)  Quart(Comp: %2.4f  Send: %2.4f)   Puzzle: %2.4f   Tree: %2.4f\n", PP_Myid, PP_quartcomptime, PP_quartsendtime, PP_puzzletime, PP_treetime);
    }
# endif /* PVERBOSE1 */

  MPI_Finalize();
}

/***********************************************************
*  main part of the slave process                          *
***********************************************************/

int slave_main(int argc, char *argv[])
{	
	int i, a, b, c, d, approx;
	double d1, d2, d3;

	int        notdone;
	MPI_Status stat;
	int        PP_AllQuartsReceived = 0;
	
        PP_RecvSizes(&Maxspc, &Maxsite, &numcats, &Numptrn, &tpmradix, &outgroup, &fracconst, &randseed);
	psteptreestrlen = (Maxspc * (int)(1 + log10(Maxspc))) +
	                  (Maxspc * 3);


        Maxbrnch = 2*Maxspc - 3;

	initrandom(randseed);

	/* initialise ML (from PP_mlstart) */
	Seqpat    = new_cmatrix(Maxspc, Numptrn);
	Alias     = new_ivector(Maxsite);
	Weight    = new_ivector(Numptrn);
	constpat  = new_ivector(Numptrn);
	Rates     = new_dvector(numcats);
	Eval      = new_dvector(tpmradix);
	Freqtpm   = new_dvector(tpmradix);
	Evec      = new_dmatrix(tpmradix,tpmradix);
	Ievc      = new_dmatrix(tpmradix,tpmradix);
	iexp      = new_dmatrix(tpmradix,tpmradix);
	Distanmat = new_dmatrix(Maxspc, Maxspc);
	ltprobr   = new_dcube(numcats, tpmradix,tpmradix);
	PP_RecvData(Seqpat,                       /* cmatrix */
	            Alias, Weight, constpat,      /* ivector */
	            Rates, Eval, Freqtpm,         /* dvector */
	            Evec, Ievc, iexp, Distanmat,  /* dmatrix */
	            ltprobr);                     /* dcube   */

	Ctree = new_quartet(Numptrn, Seqpat);
	Numbrnch = NUMQBRNCH;
	Numibrnch = NUMQIBRNCH;
	Numspc = NUMQSPC;

	/* allocate variable used for randomizing input order */
	trueID = new_ivector(Maxspc);

	/* allocate and initialize badtaxonvector */
	badtaxon = new_ulivector(Maxspc);
	for (i = 0; i < Maxspc; i++) badtaxon[i] = 0;

	/* allocate memory for quartets */
	quartetinfo = mallocquartets(Maxspc);

        /* prepare for consensus tree analysis */
	initconsensus();

	MPI_Probe(PP_MyMaster, MPI_ANY_TAG, PP_Comm, &stat);
	
	notdone = (stat.MPI_TAG != PP_DONE);
        if (!notdone) {
	  PP_RecvDone();
#         ifdef TIMEDEBUG
             printtimearr(&tarr);
#         endif /* TIMEDEBUG */

	}

	while (notdone) {
	  switch (stat.MPI_TAG) {
	    case PP_ALLQUARTS: {
	  	PP_RecvAllQuarts(Maxspc, &Numquartets, quartetinfo);
		PP_AllQuartsReceived = 1;
		break;
		}
	    case PP_DOQUARTBLOCK:
	    case PP_DOQUARTBLOCKSPECS:
                {
		uli qtodo, qstart, qend, idx, nofbq, *bqarr;
		int a, b, c, i;
		double d1, d2, d3;

	        nofbq=0;
		PP_RecvDoQuartBlock(&qstart, &qtodo, &bqarr, &approx);
		qend   = (qstart + qtodo) - 1;
#		ifdef PVERBOSE1
			FPRINTF(STDOUTFILE "(%2d) Quartets %4ld->%4ld (%dx%ld)\n", PP_Myid, qstart, qend, PP_NumProcs-1, qtodo);
#		endif

		addtimes(GENERAL, &tarr);
		for (i = 3; i < Maxspc; i++)
		   for (c = 2; c < i; c++)
		      for (b = 1; b < c; b++)
		         for (a = 0; a < b; a++) {

		            idx = (uli) a + 
		                  (uli) b*(b-1)/2 +
		                  (uli) c*(c-1)*(c-2)/6 +
		                  (uli) i*(i-1)*(i-2)*(i-3)/24;
		            if ((idx >= qstart) && (idx <= qend)) {
#			    ifdef PVERBOSE4
			       FPRINTF(STDOUTFILE "(%2d)  %4ld <---> (%d,%d,%d,%d)\n",PP_Myid, idx, a,b,c,i);
#			    endif
			       compute_quartlklhds(a,b,c,i,&d1,&d2,&d3,approx);
		               PP_do_write_quart(a,b,c,i,d1,d2,d3,&nofbq,bqarr);
			       addtimes(QUARTETS, &tarr);
		            } /* if idx */
		         } /* for for for for */   
		PP_SendQuartBlock(qstart, qtodo, quartetinfo, nofbq, bqarr, approx);

		free(bqarr); bqarr=NULL;

		break;
		}

	    case PP_DOPUZZLEBLOCK: {
		if (PP_AllQuartsReceived){
		   uli Numtrial, ptodo;
		   cmatrix *bp;
		   int n;

		   PP_RecvDoPermutBlock(&Numtrial);
		   ptodo = Numtrial;

		   bp = (cmatrix *) malloc(Numtrial * sizeof(void *));
                   for(n=0; n<Numtrial; n++) {
		      bp[n] = new_cmatrix(Maxspc - 3, Maxspc);
		   }

		   addtimes(GENERAL, &tarr);
		   for (Currtrial = 0; Currtrial < ptodo; Currtrial++) {
		      /* randomize input order */
                      biparts=bp[Currtrial];
		      chooser(Maxspc, Maxspc, trueID);

		      PP_slave_do_puzzling(trueID);
		      addtimes(PUZZLING, &tarr);
		   }
		   PP_SendSplitsBlock(Maxspc, Numtrial, bp, psteptreenum, psteptreelist);
                   for (Currtrial = 0; Currtrial < ptodo; Currtrial++) {
                      free_cmatrix(bp[Currtrial]);
                   }
		   free(bp);
		} else {
		   FPRINTF(STDOUTFILE "ERROR: Requested Puzzle Tree without receiving Quartets!!!\n");
		   notdone = 0;
		}
		freetreelist(&psteptreelist, &psteptreenum, &psteptreesum);
		break;
		}
	    case PP_DOQUART: {
		PP_RecvDoQuart(&a,&b,&c,&d, &approx);
		addtimes(GENERAL, &tarr);
		compute_quartlklhds(a,b,c,d,&d1,&d2,&d3,approx);
		addtimes(QUARTETS, &tarr);
		PP_SendQuart(a,b,c,d,d1,d2,d3, approx);
		break;
		}
	    case PP_DOPUZZLE: {
		if (PP_AllQuartsReceived){
		   PP_RecvPermut(Maxspc, trueID);
		   addtimes(GENERAL, &tarr);
		   PP_slave_do_puzzling(trueID);
		   addtimes(PUZZLING, &tarr);
		   }
		else {
		   FPRINTF(STDOUTFILE "ERROR: Requested Puzzle Tree without receiving Quartets!!!\n");
		   notdone = 0;
		}
		break;
		}
	    default: {
		FPRINTF(STDOUTFILE "ERROR: Unknown Message Tag received\n");
		MPI_Abort(PP_Comm, 1);
		}
	  } /* switch stat.MPI_TAG */ 

	  MPI_Probe(PP_MyMaster, MPI_ANY_TAG, PP_Comm, &stat);
	  notdone = (stat.MPI_TAG != PP_DONE);
	  if (!notdone)
	    PP_RecvDone();

	} /* while notdone */

        return (0);

} /* slave_main */

