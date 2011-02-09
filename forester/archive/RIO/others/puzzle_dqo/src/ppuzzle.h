/*
 *   ppuzzle.h
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


#ifndef _PPUZZLE_
#define _PPUZZLE_

#include "puzzle.h"
#include "util.h"
#include "ml.h"
#include "sched.h"

extern int PP_IamSlave;
extern int PP_IamMaster;

#ifdef PARALLEL
#    ifdef SEQUENTIAL
#        undef SEQUENTIAL
#    endif
#    define SEQUENTIAL 0
#    undef  PARALLEL
#    define PARALLEL 1
#    include "mpi.h"
#else
#    ifdef SEQUENTIAL
#        undef SEQUENTIAL
#    endif
#    define SEQUENTIAL 1
#    define PARALLEL   0
#    undef PVERBOSE
#    undef PVERBOSE1
#    undef PVERBOSE2
#    undef PVERBOSE3
#endif

/* PVERBOSE3 includes PVERBOSE2 includes PVERBOSE1 */
/* PVERBOSE1 is default (PVERBOSE)                 */

#ifdef PVERBOSE
#    undef  PVERBOSE1
#    define PVERBOSE1
#endif
#ifdef PVERBOSE3
#    undef  PVERBOSE2
#    define PVERBOSE2
#endif
#ifdef PVERBOSE2
#    undef  PVERBOSE1
#    define PVERBOSE1
#endif

#if PARALLEL
#  define PP_DONE             0 /* Finished            M->S */
#  define PP_SIZES               1 /* Array sizes needed  M->S */
#  define PP_DATA                2 /* Data Arrays         M->S */

#  define PP_ALLQUARTS           3 /* All Quartets        M->S */

#  define PP_DOQUART             4 /* do 4Specs           M->S */
#  define PP_DOQUARTX2           5 /* do 4Specs + X^2     M->S */
#  define PP_QUART               6 /* quartet back        S->M */
#  define PP_QUARTX2             7 /* quartet + X^2 back  S->M */

#  define PP_DOQUARTBLOCKSPECS   8 /* do block Specs      M->S */
#  define PP_DOQUARTBLOCK        9 /* do block of Quarts  M->S */
#  define PP_QUARTBLOCKSPECS    10 /* block Specs         S->M */
#  define PP_QUARTBLOCK         11 /* block of Quarts     S->M */

#  define PP_DOPUZZLE           12 /* do Puzzling step    M->S */
#  define PP_PUZZLE             13 /* Puzzling tree back  S->M */
#  define PP_DOPUZZLEBLOCK      14 /* do Puzzling block   M->S */
#  define PP_DOPUZZLEBLOCKSPECS 15 /* do Puzzling block   M->S */
#  define PP_PUZZLEBLOCK        16 /* Puzzling block      S->M */
#  define PP_PUZZLEBLOCKSPECS   17 /* Puzzling block      S->M */

#  define PP_STATS              18 /* Slave Statistics    S->M */

#  define PP_WAIT               18 /* waiting for work    S->M */
#  define PP_TEST              100 /* testing                  */

#  define PERMUTQUEUESIZE 100
#  define QUARTQUEUESIZE 100

   extern int      PP_IamMaster;
   extern int      PP_IamSlave;
   extern int      PP_Myid;
   extern int      PP_MyMaster;
   extern int      PP_NumProcs;
   extern MPI_Comm PP_Comm;
#endif /* PARALLEL */

extern int *permutsent,
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
extern double *walltimes,
              *cputimes;
extern double *fullwalltimes,
              *fullcputimes;
extern double *altwalltimes,
              *altcputimes;

extern int PP_permutsent,
           PP_permutrecved,
           PP_quartsent,
           PP_quartrecved,
           PP_doquartsent,
           PP_doquartrecved,
           PP_splitsent,
           PP_splitrecved,
           PP_permutsentn,
           PP_permutrecvedn,
           PP_quartsentn,
           PP_quartrecvedn,
           PP_doquartsentn,
           PP_doquartrecvedn,
           PP_splitsentn,
           PP_splitrecvedn;

extern double PP_starttime,
       PP_stoptime,
       PP_inittime,
       PP_paramcomptime,
       PP_paramsendtime,
       PP_quartcomptime,
       PP_quartsendtime,
       PP_puzzletime,
       PP_treetime;

void num2quart(uli qnum, int *a, int *b, int *c, int *d);
uli numquarts(int maxspc);
uli quart2num (int a, int b, int c, int d);

int  slave_main(int argc, char *argv[]);
void PP_Init(int *argc, char **argv[]);
void PP_Finalize();
void PP_Printerror(FILE *of, int id, int err);
void PP_do_puzzling(ivector trueID);

void PP_RecvDoQuart(int  *a,
                  int    *b,
                  int    *c,
                  int    *d,
                  int    *approx);
void PP_SendDoQuart(int    dest,
                  int    a, 
                  int    b,
                  int    c,
                  int    d,
                  int    approx);
void PP_RecvQuart(int    *a,
                  int    *b,
                  int    *c,
                  int    *d,
                  double *d1,
                  double *d2,
                  double *d3,
                  int    *approx);
void PP_SendQuart(int    a,
                  int    b,
                  int    c,
                  int    d,
                  double d1,
                  double d2,
                  double d3,
                  int    approx);
void PP_SendSizes(int    mspc, 
                  int    msite,
                  int    ncats,
                  int    nptrn,
                  int    rad,
                  int    outgr,
                  double frconst,
                  int    rseed);
void PP_RecvSizes(int    *mspc, 
                  int    *msite,
                  int    *ncats,
                  int    *nptrn,
                  int    *rad,
                  int    *outgr,
                  double *frconst,
                  int    *rseed);
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
        dcube   ltprobr);          /* dcube (numcats x tpmradix x tpmradix)  */
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
        dcube   ltprobr);          /* dcube (numcats x tpmradix x tpmradix)  */
void PP_SendAllQuarts(unsigned long  Numquartets,
                      unsigned char *quartetinfo);
void PP_RecvAllQuarts(int            taxa,
                      unsigned long *Numquartets,
                      unsigned char *quartetinfo);

void PP_SendDoQuartBlock(int dest, uli firstq, uli amount, int approx);
void PP_RecvDoQuartBlock(uli *firstq, uli *amount, uli **bq, int *approx);
void PP_SendQuartBlock(uli startq,
                       uli numofq,
                       unsigned char *quartetinfo,
                       uli  numofbq,
                       uli *bq,
                       int approx);
void PP_RecvQuartBlock(int slave,
                       uli *startq,
                       uli *numofq,
                       unsigned char *quartetinfo,
                       int *approx);

void PP_SendPermut(int      dest,
                   int      taxa, 
                   ivector  permut);
void PP_RecvPermut(int      taxa,
                   ivector  permut);
void PP_SendDoPermutBlock(uli puzzlings);
void PP_RecvDoPermutBlock(uli *taxa);

void PP_SendSplits(int      taxa, 
                   cmatrix  biparts);
void PP_RecvSplits(int      taxa,
                   cmatrix  biparts);
void PP_SendDone();
void PP_RecvDone();

int PP_emptyslave();
void PP_putslave(int sl);
int PP_getslave();

void PP_cmpd(int rank, double a, double b);
void PP_cmpi(int rank, int a, int b);

#endif /* _PPUZZLE_ */
