/*
 * puzzle.h
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


#ifndef _PUZZLE_
#define _PUZZLE_

#ifndef PACKAGE
#  define PACKAGE    "tree-puzzle"
#endif
#ifndef VERSION
#  define VERSION    "5.0"
#endif
#define DATE       "October 2000"

/* prototypes */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "util.h"
#include "ml.h"
#ifdef PARALLEL
#  include "ppuzzle.h"
#endif

#define STDOUT stdout
#ifndef PARALLEL	/* because printf() runs significantly faster */
			/* than fprintf(stdout) on an Apple McIntosh  */
			/* (HS) */
#	define FPRINTF    printf
#	define STDOUTFILE
#else
#	define FPRINTF    fprintf
#	define STDOUTFILE STDOUT,
#endif

/* filenames */
#  define FILENAMELENTH 2048


#  define INFILEDEFAULT     "infile"
#  define OUTFILEDEFAULT    "outfile"
#  define TREEFILEDEFAULT   "outtree"
#  define INTREEDEFAULT     "intree"
#  define DISTANCESDEFAULT  "outdist"
#  define TRIANGLEDEFAULT   "outlm.eps"
#  define UNRESOLVEDDEFAULT "outqlist"
#  define ALLQUARTDEFAULT   "outallquart"
#  define ALLQUARTLHDEFAULT "outallquartlh"
#  define OUTPTLISTDEFAULT  "outpstep"
#  define OUTPTORDERDEFAULT "outptorder"

#  define INFILE     infilename
#  define OUTFILE    outfilename
#  define TREEFILE   outtreename
#  define INTREE     intreename
#  define DISTANCES  outdistname
#  define TRIANGLE   outlmname
#  define UNRESOLVED outqlistname
#  define ALLQUART   outallquartname
#  define ALLQUARTLH outallquartlhname
#  define OUTPTLIST  outpstepname
#  define OUTPTORDER outptordername

EXTERN char infilename        [FILENAMELENTH];
EXTERN char outfilename       [FILENAMELENTH];
EXTERN char outtreename       [FILENAMELENTH];
EXTERN char intreename        [FILENAMELENTH];
EXTERN char outdistname       [FILENAMELENTH];
EXTERN char outlmname         [FILENAMELENTH];
EXTERN char outqlistname      [FILENAMELENTH];
EXTERN char outallquartname   [FILENAMELENTH];
EXTERN char outallquartlhname [FILENAMELENTH];
EXTERN char outpstepname      [FILENAMELENTH];
EXTERN char outptordername    [FILENAMELENTH];

#define OUTFILEEXT    "puzzle"
#define TREEFILEEXT   "tree"
#define DISTANCESEXT  "dist"
#define TRIANGLEEXT   "eps"
#define UNRESOLVEDEXT "qlist"
#define ALLQUARTEXT   "allquart"
#define ALLQUARTLHEXT "allquartlh"
#define OUTPTLISTEXT  "pstep"
#define OUTPTORDEREXT "ptorder"

#ifndef PARALLEL             /* because printf() runs significantly faster */
                             /* than fprintf(stdout) on an Apple McIntosh  */
                             /* (HS) */
#	define FPRINTF    printf
#	define STDOUTFILE
#else
#	define FPRINTF    fprintf
#	define STDOUT     stdout
#	define STDOUTFILE STDOUT,
#endif


/* auto_aamodel/auto_datatype values  (xxx) */
#define AUTO_OFF      0
#define AUTO_GUESS    1
#define AUTO_DEFAULT  2


/* qptlist values  (xxx) */
#define PSTOUT_NONE      0
#define PSTOUT_ORDER     1
#define PSTOUT_LISTORDER 2
#define PSTOUT_LIST      3

/* dtat_optn values  (xxx) */
#define NUCLEOTIDE 0
#define AMINOACID  1
#define BINARY     2

/* typ_optn values  (xxx) */
#define LIKMAPING_OPTN 1
#define TREERECON_OPTN 0

/* puzzlemodes (xxx) */
#define QUARTPUZ 0
#define USERTREE 1
#define PAIRDIST 2

/* rhetmodes (xxx) Modes of rate heterogeneity */
#define UNIFORMRATE 0
#define GAMMARATE   1
#define TWORATE     2
#define MIXEDRATE   3

/* defines for types of quartet likelihood computation (xxx) */
#define EXACT  0
#define APPROX 1

/* tree structure */
typedef struct oneedge {
	/* pointer to other three edges */
	struct oneedge *up;
	struct oneedge *downleft;
	struct oneedge *downright;
	int numedge;    /* number of edge */
	uli edgeinfo;   /* value of this edge */
	int *edgemap;   /* pointer to the local edgemap */
} ONEEDGE;


/* variables */
EXTERN cmatrix biparts;      /* bipartitions of tree of current puzzling step */
EXTERN cmatrix consbiparts;  /* bipartitions of majority rule consensus tree */
EXTERN cmatrix seqchars;     /* characters contained in data set */
EXTERN cmatrix treepict;     /* picture of consensus tree */
EXTERN double minscore;      /* value of edgescore on minedge */
EXTERN double tstvf84;       /* F84 transition/transversion ratio */
EXTERN double tstvratio;     /* expected transition/transversion ratio */
EXTERN double yrtsratio;     /* expected pyrimidine/purine transition ratio */
EXTERN dvector ulkl;         /* log L of user trees */
EXTERN dmatrix allsites;     /* log L per sites of user trees */
EXTERN dvector ulklc;        /* log L of user trees (clock) */
EXTERN dmatrix allsitesc;    /* log L per sites of user trees (clock) */
EXTERN FILE *utfp;           /* pointer to user tree file */
EXTERN FILE *ofp;            /* pointer to output file */
EXTERN FILE *seqfp;          /* pointer to sequence input file */
EXTERN FILE *tfp;            /* pointer to tree file */
EXTERN FILE *dfp;            /* pointer to distance file */
EXTERN FILE *trifp;          /* pointer to triangle file */
EXTERN FILE *unresfp;        /* pointer to file with unresolved quartets */
EXTERN FILE *tmpfp;          /* pointer to temporary file */
EXTERN FILE *qptlist;        /* pointer to file with puzzling step trees */
EXTERN FILE *qptorder;       /* pointer to file with unique puzzling step trees */
EXTERN int SHcodon;          /* whether SH should be applied to 1st, 2nd codon positions */
EXTERN int utree_optn;       /* use first user tree for estimation */
EXTERN int listqptrees;      /* list puzzling step trees */
EXTERN int approxqp;         /* approximate QP quartets */
EXTERN int *edgeofleaf;      /* vector with edge number of all leaves */
EXTERN int codon_optn;       /* declares what positions in a codon should be used */
EXTERN int compclock;        /* computation of clocklike branch lengths */
EXTERN int chooseA;          /* leaf variable */
EXTERN int chooseB;          /* leaf variable */
EXTERN int clustA, clustB, clustC, clustD; /* number of members of LM clusters */
EXTERN int column;           /* used for breaking lines (writing tree to treefile) */
EXTERN int Frequ_optn;       /* use empirical base frequencies */
EXTERN int Maxbrnch;         /* 2*Maxspc - 3 */
EXTERN int Maxseqc;          /* number of sequence characters per taxum */
EXTERN int mflag;            /* flag used for correct printing of runtime messages */
EXTERN int minedge;          /* edge with minimum edgeinfo */
EXTERN int nextedge;         /* number of edges in the current tree */
EXTERN int nextleaf;         /* next leaf to add to tree */
EXTERN int numclust;         /* number of clusters in LM analysis */
EXTERN int outgroup;         /* outgroup */
EXTERN int puzzlemode;       /* computation of QP tree and/or ML distances */
EXTERN int rootsearch;       /* how location of root is found */
EXTERN int rhetmode;         /* model of rate heterogeneity */
EXTERN int splitlength;      /* length of one entry in splitpatterns */
EXTERN int *splitsizes;      /* size of all different splits of all trees */
EXTERN int usebestq_optn;    /* use only best quartet topology, no bayesian weights */
EXTERN int show_optn;        /* show unresolved quartets                        */
EXTERN int savequart_optn;   /* save memory block which quartets to file */
EXTERN int savequartlh_optn; /* save quartet likelihoods to file */
EXTERN int saveqlhbin_optn;  /* save quartet likelihoods binary */
EXTERN int readquart_optn;   /* read memory block which quartets from file */
EXTERN int sym_optn;         /* symmetrize doublet frequencies */
EXTERN int xsize;            /* depth of consensus tree picture */
EXTERN int ytaxcounter;      /* counter for establishing y-coordinates of all taxa */
EXTERN int numutrees;        /* number of users trees in input tree file */
EXTERN ivector clusterA, clusterB, clusterC, clusterD;  /* clusters for LM analysis */
EXTERN ivector consconfid;   /* confidence values of majority rule consensus tree */
EXTERN ivector conssizes;    /* partition sizes of majority rule consensus tree */
EXTERN ivector trueID;       /* leaf -> taxon on this leaf */
EXTERN ivector xcor;         /* x-coordinates of consensus tree nodes */
EXTERN ivector ycor;         /* y-coordinates of consensus tree nodes */
EXTERN ivector ycormax;      /* maximal y-coordinates of consensus tree nodes */
EXTERN ivector ycormin;      /* minimal y-coordinates of consensus tree nodes */
EXTERN ivector ycortax;      /* y-coordinates of all taxa */
EXTERN ONEEDGE *edge;        /* vector with all the edges of the tree */
EXTERN uli *splitcomp;       /* bipartition storage */
EXTERN uli *splitfreqs;      /* frequencies of all different splits of all trees */
EXTERN uli *splitpatterns;   /* all different splits of all trees */
EXTERN uli badqs;            /* number of bad quartets */
EXTERN uli consincluded;     /* number of included biparts in the consensus tree */
EXTERN uli Currtrial;        /* counter for puzzling steps */
EXTERN uli maxbiparts;       /* space is reserved for that many bipartitions */
EXTERN uli mininfo;          /* value of edgeinfo on minedge */
EXTERN uli numbiparts;       /* number of different bipartitions */
EXTERN uli Numquartets;      /* number of quartets */
EXTERN uli Numtrial;         /* number of puzzling steps */
EXTERN uli lmqts;            /* quartets investigated in LM analysis (0 = ALL) */

EXTERN int auto_datatype;       /* guess datatype ? */
EXTERN int guessdata_optn;      /* guessed datatype */

EXTERN int auto_aamodel;        /* guess amino acid modell ? */
EXTERN int guessauto_aamodel;   /* guessed amino acid modell ? */
EXTERN int guessDayhf_optn;     /* guessed Dayhoff model option */
EXTERN int guessJtt_optn;       /* guessed JTT model option */
EXTERN int guessblosum62_optn;  /* guessed BLOSUM 62 model option */
EXTERN int guessmtrev_optn;     /* guessed mtREV model option */
EXTERN int guesscprev_optn;     /* guessed cpREV model option */
EXTERN int guessvtmv_optn;      /* guessed VT model option */
EXTERN int guesswag_optn;       /* guessed WAG model option */

/* counter variables needed in likelihood mapping analysis */
EXTERN uli ar1, ar2, ar3;
EXTERN uli reg1, reg2, reg3, reg4, reg5, reg6, reg7;
EXTERN uli reg1l, reg1r, reg2u, reg2d, reg3u, reg3d,
 reg4u, reg4d, reg5l, reg5r, reg6u, reg6d;
EXTERN unsigned char *quartetinfo; /* place where quartets are stored */
EXTERN dvector qweight; /* for use in QP and LM analysis */
EXTERN dvector sqdiff;
EXTERN ivector qworder;
EXTERN ivector sqorder;

EXTERN int randseed;
EXTERN int psteptreestrlen;

typedef struct treelistitemtypedummy {
	struct treelistitemtypedummy *pred;
	struct treelistitemtypedummy *succ;
	struct treelistitemtypedummy *sortnext;
	struct treelistitemtypedummy *sortlast;
	char  *tree;
	int    count;
	int    id;
	int    idx;
} treelistitemtype;

EXTERN treelistitemtype *psteptreelist;
EXTERN treelistitemtype *psteptreesortlist;
EXTERN int               psteptreenum;
EXTERN int               psteptreesum;


/* prototypes */
void makeF84model(void);
void compnumqts(void);
void setoptions(void);
void openfiletoread(FILE **, char[], char[]);
void openfiletowrite(FILE **, char[], char[]);
void openfiletoappend(FILE **, char[], char[]);
void closefile(FILE *);
void symdoublets(void);
void computeexpectations(void);
void putdistance(FILE *);
void findidenticals(FILE *);
double averagedist(void);
void initps(FILE *);
void plotlmpoint(FILE *, double, double);
void finishps(FILE *);
void makelmpoint(FILE *, double, double, double);
void printtreestats(FILE *);
void timestamp(FILE *);
void writeoutputfile(FILE *, int);

/* definitions for writing output */
#define WRITEALL    0
#define WRITEPARAMS 1
#define WRITEREST   2

void writetimesstat(FILE *ofp);
void writecutree(FILE *, int);
void starttimer(void);
void checktimer(uli);
void estimateparametersnotree(void);
void estimateparameterstree(void);
int main(int, char *[]);
int ulicmp(const void *, const void *);
int intcmp(const void *, const void *);

void readid(FILE *, int);
char readnextcharacter(FILE *, int, int);
void skiprestofline(FILE *, int, int);
void skipcntrl(FILE *, int, int);
void getseqs(FILE *);
void initid(int);
void fputid10(FILE *, int);
int fputid(FILE *, int);
void getsizesites(FILE *);
void getdataset(FILE *);
int guessdatatype(void);
void translatedataset(void);
void estimatebasefreqs(void);
void guessmodel(void);
void inittree(void);
void addnextleaf(int);
void freetree(void);
void writeOTU(FILE *, int);
void writetree(FILE *);
int *initctree();
void copytree(int *ctree);
void freectree(int **snodes);
void printctree(int *ctree);
char *sprintfctree(int *ctree, int strlen);
void fprintffullpstree(FILE *outf, char *treestr);
int printfsortctree(int *ctree);
int sortctree(int *ctree);
int ct_1stedge(int node);
int ct_2ndedge(int node);
int ct_3rdedge(int node);

void printfpstrees(treelistitemtype *list);
void printfsortedpstrees(treelistitemtype *list);
void fprintfsortedpstrees(FILE *output, treelistitemtype *list, int itemnum, int itemsum, int comment, float cutoff);

void sortbynum(treelistitemtype *list, treelistitemtype **sortlist);
treelistitemtype *addtree2list(char             **tree,
                               int                numtrees,
                               treelistitemtype **list,
                               int               *numitems,
                               int               *numsum);
void freetreelist(treelistitemtype **list,
                  int               *numitems,
                  int               *numsum);
void resetedgeinfo(void);
void incrementedgeinfo(int, int);
void minimumedgeinfo(void);
void initconsensus(void);
void makepart(int, int);
void computebiparts(void);
void printsplit(FILE *, uli);
void makenewsplitentries(void);
void copysplit(uli, int);
void makeconsensus(void);
void writenode(FILE *, int);
void writeconsensustree(FILE *);
void nodecoordinates(int);
void drawnode(int, int);
void plotconsensustree(FILE *);
unsigned char *mallocquartets(int);
void freequartets(void);
unsigned char readquartet(int, int, int, int);
void writequartet(int, int, int, int, unsigned char);
void sort3doubles(dvector, ivector);
void computeallquartets(void);
void checkquartet(int, int, int, int);
void num2quart(uli qnum, int *a, int *b, int *c, int *d);
uli numquarts(int maxspc);
uli quart2num (int a, int b, int c, int d);

void writetpqfheader(int nspec, FILE *ofp, int flag);


/* extracted from main (xxx) */
void compute_quartlklhds(int a, int b, int c, int d, double *d1, double *d2, double *d3, int approx);


/* definitions for timing */

#define OVERALL   0
#define GENERAL   1
#define OPTIONS   2
#define PARAMEST  3
#define QUARTETS  4
#define PUZZLING  5
#define TREEEVAL  6

typedef struct {
	int      currentjob;
	clock_t  tempcpu;
	clock_t  tempfullcpu;
	clock_t  tempcpustart;
	time_t   temptime;
	time_t   tempfulltime;
	time_t   temptimestart;

        clock_t  maxcpu;
	clock_t  mincpu;
	time_t   maxtime;
	time_t   mintime;

	double   maxcpublock;
	double   mincpublock;
	double   mincputick;
	double   mincputicktime;
	double   maxtimeblock;
	double   mintimeblock;

	double   generalcpu;
	double   optionscpu;
	double   paramestcpu;
	double   quartcpu;
	double   quartblockcpu;
	double   quartmaxcpu;
	double   quartmincpu;
	double   puzzcpu;
	double   puzzblockcpu;
	double   puzzmaxcpu;
	double   puzzmincpu;
	double   treecpu;
	double   treeblockcpu;
	double   treemaxcpu;
	double   treemincpu;
	double   cpu;
	double   fullcpu;

	double   generaltime;
	double   optionstime;
	double   paramesttime;
	double   quarttime;
	double   quartblocktime;
	double   quartmaxtime;
	double   quartmintime;
	double   puzztime;
	double   puzzblocktime;
	double   puzzmaxtime;
	double   puzzmintime;
	double   treetime;
	double   treeblocktime;
	double   treemaxtime;
	double   treemintime;
	double   time;
	double   fulltime;
} timearray_t;

EXTERN double cputime, walltime;
EXTERN double fullcpu, fulltime;
EXTERN double fullcputime, fullwalltime;
EXTERN double altcputime, altwalltime;
EXTERN clock_t cputimestart,  cputimestop, cputimedummy;
EXTERN time_t  walltimestart, walltimestop, walltimedummy;
EXTERN clock_t Startcpu;     /* start cpu time */
EXTERN clock_t Stopcpu;      /* stop cpu time */
EXTERN time_t Starttime;     /* start time */
EXTERN time_t Stoptime;      /* stop time */
EXTERN time_t time0;         /* timer variable */
EXTERN time_t time1;         /* yet another timer */
EXTERN time_t time2;         /* yet another timer */
EXTERN timearray_t tarr;

void resetqblocktime(timearray_t *ta);
void resetpblocktime(timearray_t *ta);
void inittimearr(timearray_t *ta);
void addtimes(int jobtype, timearray_t *ta);
#ifdef TIMEDEBUG
  void printtimearr(timearray_t *ta);
#endif /* TIMEDEBUG */

#endif /* _PUZZLE_ */

