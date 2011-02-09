/*
 * util.h
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


#ifndef _UTIL_
#define _UTIL_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/*
 * general definitions
 */

#define TRUE 1
#define FALSE 0

#ifdef PARALLEL
	extern long int PP_randn;
	extern long int PP_rand;
#endif

/*
 * type definitions
 */

typedef unsigned long int uli;

typedef double *dvector, **dmatrix, ***dcube;
typedef char *cvector, **cmatrix, ***ccube;
typedef int *ivector, **imatrix, ***icube;
typedef uli *ulivector, **ulimatrix, ***ulicube;


/*
 * prototypes of functions defined in util.c
 */

void maerror(char *message);

dvector new_dvector(int n);
dmatrix new_dmatrix(int nrow, int ncol);
dcube new_dcube(int ntri, int nrow, int ncol);
void free_dvector(dvector v);
void free_dmatrix(dmatrix m);
void free_dcube(dcube c);

cvector new_cvector(int n);
cmatrix new_cmatrix(int nrow, int ncol);
ccube new_ccube(int ntri, int nrow, int ncol);
void free_cvector(cvector v);
void free_cmatrix(cmatrix m);
void free_ccube(ccube c);

ivector new_ivector(int n);
imatrix new_imatrix(int nrow, int ncol);
icube new_icube(int ntri, int nrow, int ncol);
void free_ivector(ivector v);
void free_imatrix(imatrix m);
void free_icube(icube c);

ulivector new_ulivector(int n);
ulimatrix new_ulimatrix(int nrow, int ncol);
ulicube new_ulicube(int ntri, int nrow, int ncol);
void free_ulivector(ulivector v);
void free_ulimatrix(ulimatrix m);
void free_ulicube(ulicube c);

double randomunitintervall(void);
int initrandom(int seed);
int randominteger(int n);
void chooser(int t, int s, ivector slist);
void *myrealloc(void *, size_t);
cvector mygets(void);

#define MAXITS 10      /* maximum number of iterations in twoedimenmin */
double onedimenmin(double, double, double, double (*f )(double ), double, double *, double *);
void twodimenmin(double, int, double, double *, double, double (*func1 )(double ), double *, int, double, double *, double, double (*func2 )(double ), double *);



#endif
