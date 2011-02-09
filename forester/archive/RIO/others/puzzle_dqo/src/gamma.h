/*
 * gamma.h
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

#ifndef _GAMMA_
#define _GAMMA_

double densityGamma (double, double);
double cdfGamma (double, double);
double icdfGamma (double, double);
double momentGamma (int, double);

double LnGamma (double);
double IncompleteGammaQ (double, double);

double chi2prob (int, double);
double chi2test (double *, int *, int , int *);


#endif /* _GAMMA_ */
