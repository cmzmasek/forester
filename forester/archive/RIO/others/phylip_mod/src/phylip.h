/*Modified by Christian Zmasek. Use at your own risk.*/


#ifndef _PHYLIP_H_
#define _PHYLIP_H_

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   Mike Palczewski, Doug Buxton and Dan Fineman.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define VERSION "3.65"

/* machine-specific stuff:
   based on a number of factors in the library stdlib.h, we will try
   to determine what kind of machine/compiler this program is being
   built on.  However, it doesn't always succeed.  However, if you have
   ANSI conforming C, it will probably work.

   We will try to figure out machine type
   based on defines in stdio, and compiler-defined things as well.: */

#include <stdio.h>
#include <stdlib.h>
#ifdef WIN32
#include <windows.h>

void phyClearScreen(void);
void phySaveConsoleAttributes(void);
void phySetConsoleAttributes(void);
void phyRestoreConsoleAttributes(void);
void phyFillScreenColor(void);

#endif

#ifdef  GNUDOS
#define DJGPP
#define DOS
#endif

#ifdef THINK_C
#define MAC
#endif
#ifdef __MWERKS__
#ifndef WIN32
#define MAC
#endif
#endif

#ifdef __CMS_OPEN
#define CMS
#define EBCDIC true
#define INFILE "infile data"
#define OUTFILE "outfile data"
#define FONTFILE "fontfile data"
#define PLOTFILE "plotfile data"
#define INTREE "intree data"
#define INTREE2 "intree data 2"
#define OUTTREE "outtree data"
#define CATFILE "categories data"
#define WEIGHTFILE "weights data"
#define ANCFILE "ancestors data"
#define MIXFILE "mixture data"
#define FACTFILE "factors data"
#else
#define EBCDIC false
#define INFILE "infile"
#define OUTFILE "outfile"
#define FONTFILE "fontfile" /* on unix this might be /usr/local/lib/fontfile */
#define PLOTFILE "plotfile"
#define INTREE "intree"
#define INTREE2 "intree2"
#define OUTTREE "outtree"
#define CATFILE "categories"
#define WEIGHTFILE "weights"
#define ANCFILE "ancestors"
#define MIXFILE "mixture"
#define FACTFILE "factors"
#endif

#ifdef L_ctermid            /* try and detect for sysV or V7. */
#define SYSTEM_FIVE
#endif

#ifdef sequent
#define SYSTEM_FIVE
#endif

#ifndef SYSTEM_FIVE
#include <stdlib.h>
# if defined(_STDLIB_H_) || defined(_H_STDLIB) || defined(H_SCCSID) || defined(unix)
# define UNIX
# define MACHINE_TYPE "BSD Unix C"
# endif
#endif


#ifdef __STDIO_LOADED
#define VMS
#define MACHINE_TYPE "VAX/VMS C"
#endif

#ifdef __WATCOMC__
#define QUICKC
#define WATCOM
#define DOS
#include "graph.h"
#endif
/* watcom-c has graphics library calls that are almost identical to    *
 * quick-c, so the "QUICKC" symbol name stays.                         */


#ifdef _QC
#define MACHINE_TYPE "MS-DOS / Quick C"
#define QUICKC
#include "graph.h"
#define DOS
#endif

#ifdef _DOS_MODE
#define MACHINE_TYPE "MS-DOS /Microsoft C "
#define DOS           /* DOS is always defined if on a DOS machine */
#define MSC           /* MSC is defined for microsoft C              */
#endif

#ifdef __MSDOS__      /* TURBO c compiler, ONLY (no other DOS C compilers) */
#define DOS
#define TURBOC
#include <stdlib.h>
#include <graphics.h>
#endif

#ifdef DJGPP          /* DJ Delorie's original gnu  C/C++ port */
#include <graphics.h>
#endif

#ifndef MACHINE_TYPE
#define MACHINE_TYPE "ANSI C"
#endif

#ifdef DOS
#define MALLOCRETURN void 
#else
#define MALLOCRETURN void
#endif
#ifdef VMS
#define signed /* signed doesn't exist in VMS */
#endif

/* default screen types */
/*  if on a DOS but not a Windows system can use IBM PC screen controls */
#ifdef DOS
#ifndef WIN32
#define IBMCRT true
#define ANSICRT false
#endif
#endif
/*  if on a Mac cannot use screen controls */
#ifdef MAC
#define IBMCRT false 
#define ANSICRT false
#endif
/*  if on a Windows system can use IBM PC screen controls */
#ifdef WIN32
#define IBMCRT true 
#define ANSICRT false
#endif
/* otherwise, let's assume we are on a Linux or Unix system
   with ANSI terminal controls */
#ifndef MAC
#ifndef DOS
#ifndef WIN32
#define IBMCRT false 
#define ANSICRT true
#endif
#endif
#endif

#ifdef DJGPP
#undef MALLOCRETURN
#define MALLOCRETURN void
#endif


/* includes: */
#ifdef UNIX
#include <strings.h>
#else
#include <string.h>
#endif

#include <math.h>
#include <ctype.h>

#ifdef MAC
#ifdef DRAW
#include "interface.h"
#else
#include "macface.h"
#endif
#define getch gettch
#endif

/* directory delimiters */
#ifdef MAC
#define DELIMITER ':'
#else 
#ifdef WIN32
#define DELIMITER '\\'
#else 
#define DELIMITER '/'
#endif
#endif


#define FClose(file) if (file) fclose(file) ; file=NULL
#define Malloc(x) mymalloc((long)x)

typedef void *Anyptr;
#define Signed     signed
#define Const     const
#define Volatile  volatile
#define Char        char      /* Characters (not bytes) */
#define Static     static     /* Private global funcs and vars */
#define Local      static     /* Nested functions */

typedef unsigned char boolean;

#define true    1
#define false   0
#define SETBITS 31

MALLOCRETURN    *mymalloc(long);

#define FNMLNGTH        200  /* length of array to store a file name */
#define MAXNCH          26   /*changed from to 20 to 26 by CZ 2006-07-28 */
#define nmlngth         26   /*changed from to 10 to 26 by CZ 2006-07-28 */ /* number of characters in species name    */
#define maxcategs       9    /* maximum number of site types */
#define maxcategs2     11    /* maximum number of site types + 2 */
#define point           "."
#define pointe          '.'
#define down            2
#define MAXSHIMOTREES 100

#define smoothings      4    /* number of passes through smoothing algorithm */
#define iterations      4    /* number of iterates for each branch           */
#define epsilon         0.0001  /* small number used in makenewv */
#define EPSILON         0.00001  /* small number used in hermite root-finding */
#define initialv        0.1  /* starting branch length unless otherwise */
#define over            60   /* maximum width all branches of tree on screen */
#define SQRTPI 1.7724538509055160273
#define SQRT2  1.4142135623730950488

#define NLRSAVES 5 /* number of views that need to be saved during local  *
                    * rearrangement                                       */

typedef long *steptr;
typedef long longer[6];
typedef char naym[MAXNCH];
typedef long *bitptr;
typedef double raterootarray[maxcategs2][maxcategs2];

typedef struct bestelm {
  long *btree;
  boolean gloreange;
  boolean locreange;
  boolean collapse;
} bestelm;

extern FILE *infile, *outfile, *intree, *intree2, *outtree,
    *weightfile, *catfile, *ancfile, *mixfile, *factfile;
extern long spp, words, bits;
extern boolean ibmpc, ansi, tranvsp;
extern naym *nayme;                     /* names of species */


#define ebcdic          EBCDIC

typedef Char plotstring[MAXNCH];

/* Approx. 1GB, used to test for memory request errors */
#define TOO_MUCH_MEMORY 1000000000


/* The below pre-processor commands define the type used to store
   group arrays.  We can't use #elif for metrowerks, so we use
   cascaded if statements */
#include <limits.h>

/* minimum double we feel safe with, anything less will be considered 
   underflow */
#define MIN_DOUBLE 10e-100

/* K&R says that there should be a plus in front of the number, but no
   machine we've seen actually uses one; we'll include it just in
   case. */
#define MAX_32BITS        2147483647
#define MAX_32BITS_PLUS  +2147483647

/* If ints are 4 bytes, use them */
#if INT_MAX == MAX_32BITS
typedef int  group_type;

#else  
     #if INT_MAX == MAX_32BITS_PLUS
     typedef int  group_type;

     #else 
          /* Else, if longs are 4 bytes, use them */
          #if LONG_MAX == MAX_32BITS
          typedef long group_type;

          #else 
               #if LONG_MAX == MAX_32BITS_PLUS
                typedef long group_type;

               /* Default to longs */
               #else
                    typedef long group_type;
               #endif

          #endif
     #endif
#endif

/* for many programs */

#define maxuser         1000  /* maximum number of user-defined trees    */

typedef Char **sequence;

typedef enum {
  A, C, G, T, O
} bases;

typedef enum {
  alanine, arginine, asparagine, aspartic, cysteine, 
  glutamine, glutamic, glycine, histidine, isoleucine,
  leucine, lysine, methionine, phenylalanine, proline,
  serine, threonine, tryptophan, tyrosine, valine
} acids;

/* for Pars */

typedef enum {
  zero = 0, one, two, three, four, five, six, seven
} discbases;

/* for Protpars */

typedef enum {
  ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro,
  ser1, ser2, thr, trp, tyr, val, del, stop, asx, glx, ser, unk, quest
} aas;

typedef double sitelike[(long)T - (long)A + 1];   /* used in dnaml, dnadist */
typedef double psitelike[(long)valine - (long)alanine + 1];
                             /* used in proml                                    */        
     
typedef long *baseptr;       /* baseptr used in dnapars, dnacomp & dnapenny */
typedef long *baseptr2;      /* baseptr used in dnamove                     */
typedef unsigned char *discbaseptr;         /* discbaseptr used in pars     */
typedef sitelike *ratelike;                    /* used in dnaml ...            */
typedef psitelike *pratelike;                    /* used in proml                    */
typedef ratelike *phenotype;    /* phenotype used in dnaml, dnamlk, dnadist */
typedef pratelike *pphenotype;  /* phenotype used in proml                    */ 
typedef double *sitelike2;
typedef sitelike2 *phenotype2;              /* phenotype2 used in restml    */
typedef double *phenotype3;                 /* for continuous char programs */

typedef double *vector;                     /* used in distance programs    */

typedef long nucarray[(long)O - (long)A + 1];
typedef long discnucarray[(long)seven - (long)zero + 1];

typedef enum { nocollap, tocollap, undefined } collapstates;

typedef enum { bottom, nonbottom, hslength, tip, iter, length,
                 hsnolength, treewt, unittrwt } initops;


typedef double **transmatrix;
typedef transmatrix *transptr;                /* transptr used in restml */

typedef long sitearray[3];
typedef sitearray *seqptr;                    /* seqptr used in protpars */

typedef struct node {
  struct node *next, *back;
  plotstring nayme;
  long naymlength, tipsabove, index;
  double times_in_tree;            /* Previously known as cons_index */
  double xcoord, ycoord;
  long long_xcoord, long_ycoord;         /* for use in cons.               */
  double oldlen, length, r, theta, oldtheta, width, depth,
         tipdist, lefttheta, righttheta;
  group_type *nodeset;                   /* used by accumulate      -plc   */
  long ymin, ymax;                       /* used by printree        -plc   */
  boolean haslength;               /* haslength used in dnamlk             */
  boolean iter;                    /* iter used in dnaml, fitch & restml   */
  boolean initialized;             /* initialized used in dnamlk & restml  */
  long branchnum;                  /* branchnum used in restml             */
  phenotype x;                     /* x used in dnaml, dnamlk, dnadist     */
  phenotype2 x2;                   /* x2 used in restml                    */
  phenotype3 view;                 /* contml etc                           */
  pphenotype protx;                /* protx used in proml */ 
  aas *seq;                  /* the sequence used in protpars              */
  seqptr siteset;            /* temporary storage for aa's used in protpars*/
  double v, deltav, ssq;       /* ssq used only in contrast                */
  double bigv;                 /* bigv used in contml                      */
  double tyme, oldtyme;        /* used in dnamlk                           */
  double t;                    /* time in kitsch                           */
  boolean sametime;            /* bookkeeps scrunched nodes in kitsch      */
  double weight;               /* weight of node used by scrunch in kitsch */
  boolean processed;           /* used by evaluate in kitsch               */
  boolean deleted;        /* true if node is deleted (retree)              */
  boolean hasname;        /* true if tip has a name (retree)               */
  double beyond;          /* distance beyond this node to most distant tip */
                            /* (retree) */
  boolean deadend;          /* true if no undeleted nodes beyond this node */
                            /* (retree) */
  boolean onebranch;        /* true if there is one undeleted node beyond  */
                            /* this node (retree)                          */
  struct node *onebranchnode;
                            /* if there is, a pointer to that node (retree)*/
  double onebranchlength;   /* if there is, the distance from here to there*/
                                /* (retree)                                */
  boolean onebranchhaslength;   /* true if there is a valid combined length*/
                                 /* from here to there (retree)            */
  collapstates collapse;         /* used in dnapars & dnacomp              */
  boolean tip;
  boolean bottom;                /* used in dnapars & dnacomp, disc char   */
  boolean visited;               /* used in dnapars & dnacomp  disc char   */
  baseptr base;                  /* the sequence in dnapars/comp/penny     */
  discbaseptr discbase;          /* the sequence in pars                   */
  baseptr2 base2;                /* the sequence in dnamove                */
  baseptr oldbase;               /* record previous sequence               */
  discbaseptr olddiscbase;       /* record previous sequence               */
  long numdesc;                  /* number of immediate descendants        */
  nucarray *numnuc;              /* bookkeeps number of nucleotides        */
  discnucarray *discnumnuc;      /* bookkeeps number of nucleotides        */
  steptr numsteps;               /* bookkeeps steps                        */
  steptr oldnumsteps;            /* record previous steps                  */
  double sumsteps;               /* bookkeeps sum of steps                 */
  nucarray cumlengths;           /* bookkeeps cummulative minimum lengths  */
  discnucarray disccumlengths;   /* bookkeeps cummulative minimum lengths  */
  nucarray numreconst;           /* bookkeeps number of  reconstructions   */
  discnucarray discnumreconst;   /* bookkeeps number of  reconstructions   */
  vector d, w;                   /* for distance matrix programs           */
  double dist;                   /* dist used in fitch                     */
  bitptr stateone, statezero;    /* discrete char programs                 */
  long maxpos;                   /* maxpos used in Clique                  */
  Char state;                    /* state used in Dnamove, Dolmove & Move  */
  double* underflows;            /* used to record underflow               */
} node;

typedef node **pointarray;

typedef struct tree {
  pointarray nodep;
  double likelihood;
  transptr trans;                 /* restml */
  long *freetrans;                /* restml */
  long transindex;                /* restml */
  node *start;                    
  node *root;                    
} tree;

typedef void (*initptr)(node **, node **, node *, long, long,
                         long *, long *, initops, pointarray,
                         pointarray, Char *, Char *, FILE *);

#ifndef OLDC
/* function prototypes */
void   scan_eoln(FILE *);
boolean    eoff(FILE *);
boolean    eoln(FILE *);
int    filexists(char *);
const char*  get_command_name (const char *);
void   getstryng(char *);
void   openfile(FILE **,const char *,const char *,const char *,const char *,
                char *);
void   cleerhome(void);
void   loopcount(long *, long);
double randum(longer);
void   randumize(longer, long *);
double normrand(longer);
long   readlong(const char *);

void   uppercase(Char *);
void   initseed(long *, long *, longer);
void   initjumble(long *, long *, longer, long *);
void   initoutgroup(long *, long);
void   initthreshold(double *);
void   initcatn(long *);
void   initcategs(long, double *);
void   initprobcat(long, double *, double *);
double logfac (long);
double halfroot(double (*func)(long , double), long, double, double);
double hermite(long, double);
void initlaguerrecat(long, double, double *, double *);
void   root_hermite(long, double *);
void   hermite_weight(long, double *, double *);
void   inithermitcat(long, double, double *, double *);
void   lgr(long, double, raterootarray);
double glaguerre(long, double, double);
void   initgammacat(long, double, double *, double *);
void   inithowmany(long *, long);
void   inithowoften(long *);

void   initlambda(double *);
void   initfreqs(double *, double *, double *, double *);
void   initratio(double *);
void   initpower(double *);
void   initdatasets(long *);
void   justweights(long *);
void   initterminal(boolean *, boolean *);
void   initnumlines(long *);
void   initbestrees(bestelm *, long, boolean);
void   newline(FILE *, long, long, long);

void   inputnumbers(long *, long *, long *, long);
void   inputnumbersold(long *, long *, long *, long);
void   inputnumbers2(long *, long *, long n);
void   inputnumbers3(long *, long *);
void   samenumsp(long *, long);
void   samenumsp2(long);
void   readoptions(long *, const char *);
void   matchoptions(Char *, const char *);
void   inputweights(long, steptr, boolean *);
void   inputweightsold(long, steptr, boolean *);
void   inputweights2(long, long, long *, steptr, boolean *, const char *);
void   printweights(FILE *, long, long, steptr, const char *);

void   inputcategs(long, long, steptr, long, const char *);
void   printcategs(FILE *, long, steptr, const char *);
void   inputfactors(long, Char *, boolean *);
void   inputfactorsnew(long, Char *, boolean *);
void   printfactors(FILE *, long, Char *, const char *);
void   headings(long, const char *, const char *);
void   initname(long);
void   findtree(boolean *,long *,long,long *,bestelm *);
void   addtree(long,long *,boolean,long *,bestelm *);
long   findunrearranged(bestelm *, long, boolean);
boolean torearrange(bestelm *, long);

void   reducebestrees(bestelm *, long *);
void   shellsort(double *, long *, long);
void   getch(Char *, long *, FILE *);
void   getch2(Char *, long *);
void   findch(Char, Char *, long);
void   findch2(Char, long *, long *, Char *);
void   findch3(Char, Char *, long, long);
void   processlength(double *,double *,Char *,boolean *,FILE *,long *);
void   writename(long, long, long *);
void   memerror(void);

void   odd_malloc(long);

void   gnu(node **, node **);
void   chuck(node **, node *);
void   zeronumnuc(node *, long);
void   zerodiscnumnuc(node *, long);
void   allocnontip(node *, long *, long);
void   allocdiscnontip(node *, long *, unsigned char *, long );
void   allocnode(node **, long *, long);
void   allocdiscnode(node **, long *, unsigned char *, long );
void   gnutreenode(node **, node **, long, long, long *);
void   gnudisctreenode(node **, node **, long , long, long *,
                unsigned char *);

void   chucktreenode(node **, node *);
void   setupnode(node *, long);
long   count_sibs (node *);
void   inittrav (node *);
void   commentskipper(FILE ***, long *);
long   countcomma(FILE **, long *);
long   countsemic(FILE **);
void   hookup(node *, node *);
void   link_trees(long, long , long, pointarray);
void   allocate_nodep(pointarray *, FILE **, long  *);
  
void   malloc_pheno(node *, long, long);
void   malloc_ppheno(node *, long, long);
long   take_name_from_tree (Char *, Char *, FILE *);
void   match_names_to_data (Char *, pointarray, node **, long);
void   addelement(node **, node *, Char *, long *, FILE *, pointarray,
                boolean *, boolean *, pointarray, long *, long *, boolean *,
                node **, initptr,boolean,long);
void   treeread (FILE *, node **, pointarray, boolean *, boolean *,
                pointarray, long *, boolean *, node **, initptr,boolean,long);
void   addelement2(node *, Char *, long *, FILE *, pointarray, boolean,
                double *, boolean *, long *, long *, long, boolean *,boolean,
                long);
void   treeread2 (FILE *, node **, pointarray, boolean, double *,
                boolean *, boolean *, long *,boolean,long);
void   exxit (int);
void countup(long *loopcount, long maxcount);
char gettc(FILE* file);
void unroot_r(node* p,node ** nodep, long nonodes);
void unroot(tree* t,long nonodes);
void unroot_here(node* root, node** nodep, long nonodes);
void clear_connections(tree *t, long nonodes);
void init(int argc, char** argv);
#endif /* OLDC */
#endif /* _PHYLIP_H_ */
