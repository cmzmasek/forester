/*
 * ml.h
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


#ifndef _ML_
#define _ML_

/* definitions */

#define MINTS 0.20  /* Ts/Tv parameter */
#define MAXTS 30.0
#define MINYR 0.10  /* Y/R Ts parameter */
#define MAXYR 6.00
#define MINFI 0.00  /* fraction invariable sites */
#define MAXFI 0.99  /* only for input */
#define MINGE 0.01  /* rate heterogeneity parameter */
#define MAXGE 0.99
#define MINCAT 4    /* discrete Gamma categories */
#define MAXCAT 16

#define RMHROOT    5.0     /* upper relative bound for height of root         */    
#define MAXARC     900.0   /* upper limit on branch length (PAM) = 6.0        */
#define MINARC     0.001   /* lower limit on branch length (PAM) = 0.00001    */
#define EPSILON    0.0001  /* error in branch length (PAM) = 0.000001         */
#define HEPSILON   0.0001  /* error in node and root heights                  */
#define MAXIT      100     /* maximum number of iterates of smoothing         */
#define MINFDIFF   0.00002 /* lower limit on base frequency differences       */
#define MINFREQ    0.0001  /* lower limit on base frequencies = 0.01%         */
#define NUMQBRNCH  5       /* number of branches in a quartet                 */
#define NUMQIBRNCH 1       /* number of internal branches in a quartet        */
#define NUMQSPC    4       /* number of sequences in a quartet                */

/* 2D minimisation */
#define PEPS1  0.01    /* epsilon substitution process estimation */
#define PEPS2  0.01    /* epsilon rate heterogeneity estimation  */

/* quartet series */
#define MINPERTAXUM 2
#define MAXPERTAXUM 6
#define TSDIFF 0.20
#define YRDIFF 0.10

/* type definitions */

typedef struct node
{
	struct node *isop;
	struct node *kinp;
	int descen;
	int number;
	double length;
	double lengthc;
	double varlen;
	double height;
	double varheight;
	ivector paths;
	cvector eprob;
	dcube partials; /* partial likelihoods */
	char *label;  /* internal labels */
} Node;

typedef struct tree
{
	Node *rootp;
	Node **ebrnchp;        /* list of pointers to external branches           */
	Node **ibrnchp;        /* list of pointers to internal branches           */
	double lklhd;          /* total log-likelihood                            */
	double lklhdc;         /* total log-likelihood clock                      */
	dmatrix condlkl;       /* likelihoods for each pattern and non-zero rate  */
	double rssleast;
} Tree;


/* global variables */

EXTERN Node *chep;         /* pointer to current height node                  */
EXTERN Node *rootbr;       /* pointer to root branch                          */
EXTERN Node **heights;     /* pointer to height nodes in unrooted tree        */
EXTERN int Numhts;         /* number of height nodes in unrooted tree         */
EXTERN double hroot;       /* height of root                                  */
EXTERN double varhroot;    /* variance of height of root                      */
EXTERN double maxhroot;    /* maximal height of root                          */
EXTERN int locroot;        /* location of root                                */
EXTERN int numbestroot;    /* number of best locations for root               */
EXTERN int clockmode;      /* clocklike vs. nonclocklike computation          */
EXTERN cmatrix Identif;    /* sequence names                                  */
EXTERN cmatrix Seqchar;    /* ML sequence data                                */
EXTERN cmatrix Seqpat;     /* ordered site patterns                           */
EXTERN ivector constpat;   /* indicates constant site patterns                */
EXTERN cvector seqchi;
EXTERN cvector seqchj;
EXTERN dcube partiali;
EXTERN dcube partialj;
EXTERN dcube ltprobr;      /* transition probabilites (for all non-zero rates */
EXTERN dmatrix Distanmat;  /* matrix with maximum likelihood distances        */
EXTERN dmatrix Evec;       /* Eigenvectors                                    */
EXTERN dmatrix Ievc;       /* Inverse eigenvectors                            */
EXTERN double TSparam;     /* Ts/Tv parameter                                 */
EXTERN double tsmean, yrmean;
EXTERN double YRparam;     /* Y/R Ts parameter                                */
EXTERN double geerr;       /* estimated error of rate heterogeneity           */
EXTERN double Geta;        /* rate heterogeneity parameter                    */
EXTERN double fracconst;   /* fraction of constant sites                      */
EXTERN double fracconstpat;/* fraction of constant patterns                   */
EXTERN double Proportion;  /* for tree drawing                                */
EXTERN double tserr;       /* estimated error of TSparam                      */
EXTERN double yrerr;       /* estimated error of YRparam                      */
EXTERN double fracinv;     /* fraction of invariable sites                    */
EXTERN double fierr;       /* estimated error of fracinv                      */
EXTERN dvector Brnlength;
EXTERN dvector Distanvec;
EXTERN dvector Eval;       /* Eigenvalues of 1 PAM rate matrix                */
EXTERN dvector Freqtpm;    /* base frequencies                                */
EXTERN dvector Rates;      /* rate of each of the categories                  */
EXTERN dmatrix iexp;
EXTERN imatrix Basecomp;   /* base composition of each taxon                  */
EXTERN ivector usedtaxa;   /* list needed in the input treefile procedure     */
EXTERN int numtc;          /* auxiliary variable for printing rooted tree     */
EXTERN int qcalg_optn;     /* use quartet subsampling algorithm               */
EXTERN int approxp_optn;   /* approximate parameter estimation                */
EXTERN int chi2fail;       /* flag for chi2 test                              */
EXTERN int Converg;        /* flag for ML convergence (no clock)              */
EXTERN int Convergc;       /* flag for ML convergence (clock)                 */
EXTERN int data_optn;      /* type of sequence input data                     */  
EXTERN int Dayhf_optn;     /* Dayhoff model                                   */
EXTERN int HKY_optn;       /* use HKY model                                   */
EXTERN int Jtt_optn;       /* JTT model                                       */
EXTERN int blosum62_optn;  /* BLOSUM 62 model                                 */
EXTERN int mtrev_optn;     /* mtREV model                                     */
EXTERN int cprev_optn;     /* cpREV model                                     */
EXTERN int vtmv_optn;      /* VT model                                        */
EXTERN int wag_optn;       /* WAG model                                       */
EXTERN int Maxsite;        /* number of ML characters per taxum               */
EXTERN int Maxspc;         /* number of sequences                             */
EXTERN int mlmode;         /* quartet ML or user defined tree ML              */
EXTERN int nuc_optn;       /* nucleotide (4x4) models                         */
EXTERN int Numbrnch;       /* number of branches of current tree              */
EXTERN int numcats;        /* number of rate categories                       */
EXTERN int Numconst;       /* number of constant sites                        */
EXTERN int Numconstpat;    /* number of constant patterns                     */
EXTERN int Numibrnch;      /* number of internal branches of current tree     */
EXTERN int Numitc;         /* number of ML iterations assumning clock         */
EXTERN int Numit;          /* number of ML iterations if there is convergence */
EXTERN int Numptrn;        /* number of site patterns                         */
EXTERN int Numspc;         /* number of sequences of current tree             */
EXTERN int optim_optn;     /* optimize model parameters                       */
EXTERN int grate_optim;    /* optimize Gamma rate heterogeneity parameter     */
EXTERN int SH_optn;        /* SH nucleotide (16x16) model                     */
EXTERN int TN_optn;        /* use TN model                                    */
EXTERN int tpmradix;       /* number of different states                      */
EXTERN int fracinv_optim;  /* optimize fraction of invariable sites           */
EXTERN int typ_optn;       /* type of PUZZLE analysis                         */
EXTERN ivector Weight;     /* weight of each site pattern                     */
EXTERN Tree *Ctree;        /* pointer to current tree                         */
EXTERN ulivector badtaxon; /* involment of each taxon in a bad quartet        */
EXTERN int qca, qcb, qcc, qcd; /* quartet currently optimized                 */
EXTERN ivector Alias;      /* link site -> corresponding site pattern         */
EXTERN ivector bestrate;   /* optimal assignment of rates to sequence sites   */

EXTERN int bestratefound;

/* function prototypes of all ml function */

void convfreq(dvector);
void radixsort(cmatrix, ivector, int, int, int *);
void condenceseq(cmatrix, ivector, cmatrix, ivector, int, int, int);
void countconstantsites(cmatrix, ivector, int, int, int *, int*);
void evaluateseqs(void);
void elmhes(dmatrix, ivector, int);
void eltran(dmatrix, dmatrix, ivector, int);
void mcdiv(double, double, double, double, double *, double *);
void hqr2(int, int, int, dmatrix, dmatrix, dvector, dvector);
void onepamratematrix(dmatrix);
void eigensystem(dvector, dmatrix);
void luinverse(dmatrix, dmatrix, int);
void checkevector(dmatrix, dmatrix, int);
void tranprobmat(void);
void tprobmtrx(double, dmatrix);
double comptotloglkl(dmatrix);
void allsitelkl(dmatrix, dvector);
double pairlkl(double);
double mldistance(int, int);
void initdistan(void);
void computedistan(void);
void productpartials(Node *);
void partialsinternal(Node *);
void partialsexternal(Node *);
void initpartials(Tree *);
double intlkl(double);
void optinternalbranch(Node *);
double extlkl(double);
void optexternalbranch(Node *);
void finishlkl(Node *);
double optlkl(Tree *);
double treelkl(Tree *);
void luequation(dmatrix, dvector, int);
void lslength(Tree *, dvector, int, int, dvector);

void getusertree(FILE *, cvector, int);
Node *internalnode(Tree *, char **, int *);
void constructtree(Tree *, cvector);
void removebasalbif(cvector);
void makeusertree(FILE *);
Tree *new_tree(int, int, cmatrix);
Tree *new_quartet(int, cmatrix);
void free_tree(Tree *, int);
void make_quartet(int, int, int, int);
void changedistan(dmatrix, dvector, int);
double quartet_lklhd(int, int, int, int);
double quartet_alklhd(int, int, int, int);
void readusertree(FILE *);
double usertree_lklhd(void);
double usertree_alklhd(void);
void mlstart(void);
void distupdate(int, int, int, int);
void mlfinish(void);
void prbranch(Node *, int, int, int, ivector, ivector, FILE *);
void getproportion(double *, dvector, int);
void prtopology(FILE *);
void fputphylogeny(FILE *);
void resulttree(FILE *);
void njtree(FILE *);
void njdistantree(Tree *);
void findbestratecombination(void);
void printbestratecombination(FILE *);
int checkedge(int);
void fputsubstree(FILE *, Node *);
void fputrooted(FILE *, int);
void findheights(Node *);
void initclock(int);
double clock_alklhd(int);
double heightlkl(double);
void optheight(void);
double rheightlkl(double);
void optrheight(void);
double clock_lklhd(int);
int findrootedge(void);
void resultheights(FILE *);

double homogentest(int);
void YangDiscreteGamma(double, int, double *);
void updaterates(void);
void computestat(double *, int, double *, double *);
double quartetml(int, int, int, int);
double opttsq(double);
double optyrq(double);
void optimseqevolparamsq(void);
double opttst(double);
double optyrt(double);
void optimseqevolparamst(void);
double optfi(double);
double optge(double);
void optimrateparams(void);

int gettpmradix(void);
void rtfdata(dmatrix, double *);
int code2int(cvector);
char *int2code(int);

void jttdata(dmatrix, double *);
void dyhfdata(dmatrix, double *);
void mtrevdata(dmatrix, double *);
void cprev45data(dmatrix, double *);
void blosum62data(dmatrix, double *);
void vtmvdata(dmatrix, double *);
void wagdata(dmatrix, double *);

#endif
