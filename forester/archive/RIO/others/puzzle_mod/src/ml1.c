/*
 * ml1.c
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


/******************************************************************************/
/* definitions and prototypes                                                 */
/******************************************************************************/

#define EXTERN extern

/* prototypes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "util.h"
#include "ml.h"

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


/******************************************************************************/
/* compacting sequence data information                                       */
/******************************************************************************/


/***************************** internal functions *****************************/


/* make all frequencies a little different */
void convfreq(dvector freqemp)
{
	int i, j, maxi=0;
	double freq, maxfreq, sum;


	sum = 0.0;
	maxfreq = 0.0;
	for (i = 0; i < tpmradix; i++) {
		freq = freqemp[i];
		if (freq < MINFREQ) freqemp[i] = MINFREQ;
		if (freq > maxfreq) {
			maxfreq = freq;
			maxi = i;
		}
		sum += freqemp[i];
	}
	freqemp[maxi] += 1.0 - sum;
	
	for (i = 0; i < tpmradix - 1; i++) {
		for (j = i + 1; j < tpmradix; j++) {
			if (freqemp[i] == freqemp[j]) {
				freqemp[i] += MINFDIFF/2.0;
				freqemp[j] -= MINFDIFF/2.0;
			}
		}
	}
}

/* sort site patters of original input data */
void radixsort(cmatrix seqchar, ivector ali, int maxspc, int maxsite,
	int *numptrn)
{
	int i, j, k, l, n, pass;
	int *awork;
	int *count;


	awork = new_ivector(maxsite);
	count = new_ivector(tpmradix+1);
	for (i = 0; i < maxsite; i++)
		ali[i] = i;
	for (pass = maxspc - 1; pass >= 0; pass--) {
		for (j = 0; j < tpmradix+1; j++)
			count[j] = 0;
		for (i = 0; i < maxsite; i++)
			count[(int) seqchar[pass][ali[i]]]++;
		for (j = 1; j < tpmradix+1; j++)
			count[j] += count[j-1];
		for (i = maxsite-1; i >= 0; i--)
			awork[ --count[(int) seqchar[pass][ali[i]]] ] = ali[i];
		for (i = 0; i < maxsite; i++)
			ali[i] = awork[i];
	}
	free_ivector(awork);
	free_ivector(count);
	n = 1;
	for (j = 1; j < maxsite; j++) {
		k = ali[j];
		l = ali[j-1];
		for (i = 0; i < maxspc; i++) {
			if (seqchar[i][l] != seqchar[i][k]) {
				n++;
				break;
			}
		}
	}
	*numptrn = n;
}


void condenceseq(cmatrix seqchar, ivector ali, cmatrix seqconint,
	ivector weight, int maxspc, int maxsite, int numptrn)
{
	int i, j, k, n;
	int agree_flag; /* boolean */


	n = 0;
	k = ali[n];
	for (i = 0; i < maxspc; i++) {
		seqconint[i][n] = seqchar[i][k];
	}
	weight[n] = 1;
	Alias[k] = 0;
	for (j = 1; j < maxsite; j++) {
		k = ali[j];
		agree_flag = TRUE;
		for (i = 0; i < maxspc; i++) {
			if (seqconint[i][n] != seqchar[i][k]) {
				agree_flag = FALSE;
				break;
			}
		}
		if (agree_flag == FALSE) {
			n++;
			for (i = 0; i < maxspc; i++) {
				seqconint[i][n] = seqchar[i][k];
			}
			weight[n] = 1;
			Alias[k] = n;
		} else {
			weight[n]++;
			Alias[k] = n;
		}
	}
	n++;
	if (numptrn != n) {
		/* Problem in condenceseq */
		FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR A TO DEVELOPERS\n\n\n");
		exit(1);
	}
}

void countconstantsites(cmatrix seqpat, ivector weight, int maxspc, int numptrn,
	int *numconst, int *numconstpat)
{
	int character, s, i, constflag;

	*numconst    = 0;
	*numconstpat = 0;
	for (s = 0; s < numptrn; s++) { /* check all patterns */
		constpat[s] = FALSE;
		constflag = TRUE;
		character = seqpat[0][s];
		for (i = 1; i < maxspc; i++) {
			if (seqpat[i][s] != character) {
				constflag = FALSE;
				break;
			}
		}
		if (character != tpmradix && constflag) {
			(*numconst) = (*numconst) + weight[s];
			(*numconstpat)++;
			constpat[s] = TRUE;
		}
	}
}

/***************************** exported functions *****************************/


void evaluateseqs()
{	
	ivector ali;

	convfreq(Freqtpm); /* make all frequencies slightly different */
	ali = new_ivector(Maxsite);
	radixsort(Seqchar, ali, Maxspc, Maxsite, &Numptrn);
	Seqpat = new_cmatrix(Maxspc, Numptrn);
	constpat = new_ivector(Numptrn);
	Weight = new_ivector(Numptrn);
	condenceseq(Seqchar, ali, Seqpat, Weight, Maxspc, Maxsite, Numptrn);
	free_ivector(ali);
	countconstantsites(Seqpat, Weight, Maxspc, Numptrn, &Numconst, &Numconstpat);
	fracconstpat = (double) Numconstpat / (double) Numptrn;	
	fracconst    = (double) Numconst / (double) Maxsite;	
}


/******************************************************************************/
/* computation of Pij(t)                                                      */
/******************************************************************************/


/***************************** internal functions *****************************/


void elmhes(dmatrix a, ivector ordr, int n)
{
	int m, j, i;
	double y, x;


	for (i = 0; i < n; i++)
		ordr[i] = 0;
	for (m = 2; m < n; m++) {
		x = 0.0;
		i = m;
		for (j = m; j <= n; j++) {
			if (fabs(a[j - 1][m - 2]) > fabs(x)) {
				x = a[j - 1][m - 2];
				i = j;
			}
		}
		ordr[m - 1] = i;      /* vector */
		if (i != m) {
			for (j = m - 2; j < n; j++) {
				y = a[i - 1][j];
				a[i - 1][j] = a[m - 1][j];
				a[m - 1][j] = y;
			}
			for (j = 0; j < n; j++) {
				y = a[j][i - 1];
				a[j][i - 1] = a[j][m - 1];
				a[j][m - 1] = y;
			}
		}
		if (x != 0.0) {
			for (i = m; i < n; i++) {
				y = a[i][m - 2];
				if (y != 0.0) {
					y /= x;
					a[i][m - 2] = y;
					for (j = m - 1; j < n; j++)
						a[i][j] -= y * a[m - 1][j];
					for (j = 0; j < n; j++)
						a[j][m - 1] += y * a[j][i];
				}
			}
		}
	}
}


void eltran(dmatrix a, dmatrix zz, ivector ordr, int n)
{
	int i, j, m;


	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			zz[i][j] = 0.0;
			zz[j][i] = 0.0;
		}
		zz[i][i] = 1.0;
	}
	if (n <= 2)
		return;
	for (m = n - 1; m >= 2; m--) {
		for (i = m; i < n; i++)
			zz[i][m - 1] = a[i][m - 2];
		i = ordr[m - 1];
		if (i != m) {
			for (j = m - 1; j < n; j++) {
				zz[m - 1][j] = zz[i - 1][j];
				zz[i - 1][j] = 0.0;
			}
			zz[i - 1][m - 1] = 1.0;
		}
	}
}


void mcdiv(double ar, double ai, double br, double bi,
	double *cr, double *ci)
{
	double s, ars, ais, brs, bis;


	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs * brs + bis * bis;
	*cr = (ars * brs + ais * bis) / s;
	*ci = (ais * brs - ars * bis) / s;
}


void hqr2(int n, int low, int hgh, dmatrix h,
	dmatrix zz, dvector wr, dvector wi)
{
	int i, j, k, l=0, m, en, na, itn, its;
	double p=0, q=0, r=0, s=0, t, w, x=0, y, ra, sa, vi, vr, z=0, norm, tst1, tst2;
	int notlas; /* boolean */


	norm = 0.0;
	k = 1;
	/* store isolated roots and compute matrix norm */
	for (i = 0; i < n; i++) {
		for (j = k - 1; j < n; j++)
			norm += fabs(h[i][j]);
		k = i + 1;
		if (i + 1 < low || i + 1 > hgh) {
			wr[i] = h[i][i];
			wi[i] = 0.0;
		}
	}
	en = hgh;
	t = 0.0;
	itn = n * 30;
	while (en >= low) {	/* search for next eigenvalues */
		its = 0;
		na = en - 1;
		while (en >= 1) {
			/* look for single small sub-diagonal element */
			for (l = en; l > low; l--) {
				s = fabs(h[l - 2][l - 2]) + fabs(h[l - 1][l - 1]);
				if (s == 0.0)
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l - 1][l - 2]);
				if (tst2 == tst1)
					goto L100;
			}
			l = low;
	L100:
			x = h[en - 1][en - 1];	/* form shift */
			if (l == en || l == na)
				break;
			if (itn == 0) {
				/* all eigenvalues have not converged */
				FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR B TO DEVELOPERS\n\n\n");
				exit(1);
			}
			y = h[na - 1][na - 1];
			w = h[en - 1][na - 1] * h[na - 1][en - 1];
			/* form exceptional shift */
			if (its == 10 || its == 20) {
				t += x;
				for (i = low - 1; i < en; i++)
					h[i][i] -= x;
				s = fabs(h[en - 1][na - 1]) + fabs(h[na - 1][en - 3]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
			}
			its++;
			itn--;
			/* look for two consecutive small sub-diagonal elements */
			for (m = en - 2; m >= l; m--) {
				z = h[m - 1][m - 1];
				r = x - z;
				s = y - z;
				p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
				q = h[m][m] - z - r - s;
				r = h[m + 1][m];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break;
				tst1 = fabs(p) *
						(fabs(h[m - 2][m - 2]) + fabs(z) + fabs(h[m][m]));
				tst2 = tst1 + fabs(h[m - 1][m - 2]) * (fabs(q) + fabs(r));
				if (tst2 == tst1)
					break;
			}
			for (i = m + 2; i <= en; i++) {
				h[i - 1][i - 3] = 0.0;
				if (i != m + 2)
					h[i - 1][i - 4] = 0.0;
			}
			for (k = m; k <= na; k++) {
				notlas = (k != na);
				if (k != m) {
					p = h[k - 1][k - 2];
					q = h[k][k - 2];
					r = 0.0;
					if (notlas)
						r = h[k + 1][k - 2];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x != 0.0) {
						p /= x;
						q /= x;
						r /= x;
					}
				}
				if (x != 0.0) {
					if (p < 0.0) /* sign */
						s = - sqrt(p * p + q * q + r * r);
					else
						s = sqrt(p * p + q * q + r * r);
					if (k != m)
						h[k - 1][k - 2] = -s * x;
					else {
						if (l != m)
							h[k - 1][k - 2] = -h[k - 1][k - 2];
					}
					p += s;
					x = p / s;
					y = q / s;
					z = r / s;
					q /= p;
					r /= p;
					if (!notlas) {
						for (j = k - 1; j < n; j++) {	/* row modification */
							p = h[k - 1][j] + q * h[k][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {	/* column modification */
							p = x * h[i][k - 1] + y * h[i][k];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
						}
					} else {
						for (j = k - 1; j < n; j++) {	/* row modification */
							p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
							h[k + 1][j] -= p * z;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {	/* column modification */
							p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
							h[i][k + 1] -= p * r;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k] +
								z * zz[i][k + 1];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
							zz[i][k + 1] -= p * r;
						}
					}
				}
			}	       /* for k */
		}		       /* while infinite loop */
		if (l == en) {	       /* one root found */
			h[en - 1][en - 1] = x + t;
			wr[en - 1] = h[en - 1][en - 1];
			wi[en - 1] = 0.0;
			en = na;
			continue;
		}
		y = h[na - 1][na - 1];
		w = h[en - 1][na - 1] * h[na - 1][en - 1];
		p = (y - x) / 2.0;
		q = p * p + w;
		z = sqrt(fabs(q));
		h[en - 1][en - 1] = x + t;
		x = h[en - 1][en - 1];
		h[na - 1][na - 1] = y + t;
		if (q >= 0.0) {	       /* real pair */
			if (p < 0.0) /* sign */
				z = p - fabs(z);
			else
				z = p + fabs(z);
			wr[na - 1] = x + z;
			wr[en - 1] = wr[na - 1];
			if (z != 0.0)
				wr[en - 1] = x - w / z;
			wi[na - 1] = 0.0;
			wi[en - 1] = 0.0;
			x = h[en - 1][na - 1];
			s = fabs(x) + fabs(z);
			p = x / s;
			q = z / s;
			r = sqrt(p * p + q * q);
			p /= r;
			q /= r;
			for (j = na - 1; j < n; j++) {	/* row modification */
				z = h[na - 1][j];
				h[na - 1][j] = q * z + p * h[en - 1][j];
				h[en - 1][j] = q * h[en - 1][j] - p * z;
			}
			for (i = 0; i < en; i++) {	/* column modification */
				z = h[i][na - 1];
				h[i][na - 1] = q * z + p * h[i][en - 1];
				h[i][en - 1] = q * h[i][en - 1] - p * z;
			}
			/* accumulate transformations */
			for (i = low - 1; i < hgh; i++) {
				z = zz[i][na - 1];
				zz[i][na - 1] = q * z + p * zz[i][en - 1];
				zz[i][en - 1] = q * zz[i][en - 1] - p * z;
			}
		} else {	       /* complex pair */
			wr[na - 1] = x + p;
			wr[en - 1] = x + p;
			wi[na - 1] = z;
			wi[en - 1] = -z;
		}
		en -= 2;
	}			       /* while en >= low */
	/* backsubstitute to find vectors of upper triangular form */
	if (norm != 0.0) {
		for (en = n; en >= 1; en--) {
			p = wr[en - 1];
			q = wi[en - 1];
			na = en - 1;
			if (q == 0.0) {/* real vector */
				m = en;
				h[en - 1][en - 1] = 1.0;
				if (na != 0) {
					for (i = en - 2; i >= 0; i--) {
						w = h[i][i] - p;
						r = 0.0;
						for (j = m - 1; j < en; j++)
							r += h[i][j] * h[j][en - 1];
						if (wi[i] < 0.0) {
							z = w;
							s = r;
						} else {
							m = i + 1;
							if (wi[i] == 0.0) {
								t = w;
								if (t == 0.0) {
									tst1 = norm;
									t = tst1;
									do {
										t = 0.01 * t;
										tst2 = norm + t;
									} while (tst2 > tst1);
								}
								h[i][en - 1] = -(r / t);
							} else {	/* solve real equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
								t = (x * s - z * r) / q;
								h[i][en - 1] = t;
								if (fabs(x) > fabs(z))
									h[i + 1][en - 1] = (-r - w * t) / x;
								else
									h[i + 1][en - 1] = (-s - y * t) / z;
							}
							/* overflow control */
							t = fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++)
										h[j][en - 1] /= t;
								}
							}
						}
					}
				}
			} else if (q > 0.0) {
				m = na;
				if (fabs(h[en - 1][na - 1]) > fabs(h[na - 1][en - 1])) {
					h[na - 1][na - 1] = q / h[en - 1][na - 1];
					h[na - 1][en - 1] = (p - h[en - 1][en - 1]) /
														h[en - 1][na - 1];
				} else
					mcdiv(0.0, -h[na - 1][en - 1], h[na - 1][na - 1] - p, q,
						&h[na - 1][na - 1], &h[na - 1][en - 1]);
				h[en - 1][na - 1] = 0.0;
				h[en - 1][en - 1] = 1.0;
				if (en != 2) {
					for (i = en - 3; i >= 0; i--) {
						w = h[i][i] - p;
						ra = 0.0;
						sa = 0.0;
						for (j = m - 1; j < en; j++) {
							ra += h[i][j] * h[j][na - 1];
							sa += h[i][j] * h[j][en - 1];
						}
						if (wi[i] < 0.0) {
							z = w;
							r = ra;
							s = sa;
						} else {
							m = i + 1;
							if (wi[i] == 0.0)
								mcdiv(-ra, -sa, w, q, &h[i][na - 1],
									&h[i][en - 1]);
							else {	/* solve complex equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								vr = (wr[i] - p) * (wr[i] - p);
								vr = vr + wi[i] * wi[i] - q * q;	
								vi = (wr[i] - p) * 2.0 * q;
								if (vr == 0.0 && vi == 0.0) {
									tst1 = norm * (fabs(w) + fabs(q) + fabs(x) +
												   fabs(y) + fabs(z));
									vr = tst1;
									do {
										vr = 0.01 * vr;
										tst2 = tst1 + vr;
									} while (tst2 > tst1);
								}
								mcdiv(x * r - z * ra + q * sa,
									 x * s - z * sa - q * ra, vr, vi,
								     &h[i][na - 1], &h[i][en - 1]);
								if (fabs(x) > fabs(z) + fabs(q)) {
									h[i + 1]
										[na - 1] = (q * h[i][en - 1] -
													w * h[i][na - 1] - ra) / x;
									h[i + 1][en - 1] = (-sa - w * h[i][en - 1] -
														q * h[i][na - 1]) / x;
								} else
									mcdiv(-r - y * h[i][na - 1],
										 -s - y * h[i][en - 1], z, q,
									     &h[i + 1][na - 1], &h[i + 1][en - 1]);
							}
							/* overflow control */
							t = (fabs(h[i][na - 1]) > fabs(h[i][en - 1])) ?
								 fabs(h[i][na - 1]) : fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++) {
										h[j][na - 1] /= t;
										h[j][en - 1] /= t;
									}
								}
							}
						}
					}
				}
			}
		}
		/* end back substitution. vectors of isolated roots */
		for (i = 0; i < n; i++) {
			if (i + 1 < low || i + 1 > hgh) {
				for (j = i; j < n; j++)
					zz[i][j] = h[i][j];
			}
		}
		/* multiply by transformation matrix to give vectors of
		 * original full matrix. */
		for (j = n - 1; j >= low - 1; j--) {
			m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */
			for (i = low - 1; i < hgh; i++) {
				z = 0.0;
				for (k = low - 1; k < m; k++)
					z += zz[i][k] * h[k][j];
				zz[i][j] = z;
			}
		}
	}
	return;
}


/* make rate matrix with 0.01 expected substitutions per unit time */
void onepamratematrix(dmatrix a)
{
	int i, j;
	double delta, temp, sum;
	dvector m;

	for (i = 0; i < tpmradix; i++)
	{
		for (j = 0; j < tpmradix; j++)
		{
			a[i][j] = Freqtpm[j]*a[i][j];
		}
	}

	m = new_dvector(tpmradix);
	for (i = 0, sum = 0.0; i < tpmradix; i++)
	{
		for (j = 0, temp = 0.0; j < tpmradix; j++)
			temp += a[i][j];
		m[i] = temp; /* row sum */
		sum += temp*Freqtpm[i]; /* exp. rate */
	}
	delta = 0.01 / sum; /* 0.01 subst. per unit time */
	for (i = 0; i < tpmradix; i++) {
		for (j = 0; j < tpmradix; j++) {
			if (i != j)
				a[i][j] = delta * a[i][j];
			else
				a[i][j] = delta * (-m[i]);
		}
	}
	free_dvector(m);
}


void eigensystem(dvector eval, dmatrix evec)
{
	dvector evali, forg;
	dmatrix a, b;
	ivector ordr;
	int i, j, k, error;
	double zero;


	ordr = new_ivector(tpmradix);
	evali = new_dvector(tpmradix);
	forg = new_dvector(tpmradix);
	a = new_dmatrix(tpmradix,tpmradix);
	b = new_dmatrix(tpmradix,tpmradix);

	rtfdata(a, forg); /* get relative transition matrix and frequencies */
	
	onepamratematrix(a); /* make 1 PAM rate matrix */
	
	/* copy a to b */
	for (i = 0; i < tpmradix; i++)
		for (j = 0; j < tpmradix; j++)
			b[i][j] = a[i][j];

	elmhes(a, ordr, tpmradix); /* compute eigenvalues and eigenvectors */
	eltran(a, evec, ordr, tpmradix);
	hqr2(tpmradix, 1, tpmradix, a, evec, eval, evali);
	
	/* check eigenvalue equation */
	error = FALSE;
	for (j = 0; j < tpmradix; j++) {
		for (i = 0, zero = 0.0; i < tpmradix; i++) {
			for (k = 0; k < tpmradix; k++) zero += b[i][k] * evec[k][j];
			zero -= eval[j] * evec[i][j];
			if (fabs(zero) > 1.0e-5)
				error = TRUE;	
		}
	}
	if (error)
		FPRINTF(STDOUTFILE "\nWARNING: Eigensystem doesn't satisfy eigenvalue equation!\n");

	free_ivector(ordr);
	free_dvector(evali);
	free_dvector(forg);
	free_dmatrix(a);
	free_dmatrix(b);
}


void luinverse(dmatrix inmat, dmatrix imtrx, int size)
{
    double eps = 1.0e-20; /* ! */
	int i, j, k, l, maxi=0, idx, ix, jx;
	double sum, tmp, maxb, aw;
	ivector index;
	double *wk;
	dmatrix omtrx;

	
	index = new_ivector(tpmradix);
	omtrx = new_dmatrix(tpmradix,tpmradix);
	
	/* copy inmat to omtrx */
	for (i = 0; i < tpmradix; i++)
		for (j = 0; j < tpmradix; j++)
			omtrx[i][j] = inmat[i][j];
	
	wk = (double *) malloc((unsigned)size * sizeof(double));
	aw = 1.0;
	for (i = 0; i < size; i++) {
		maxb = 0.0;
		for (j = 0; j < size; j++) {
			if (fabs(omtrx[i][j]) > maxb)
				maxb = fabs(omtrx[i][j]);
		}
		if (maxb == 0.0) {
			/* Singular matrix */
			FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR C TO DEVELOPERS\n\n\n");
			exit(1);
		}
		wk[i] = 1.0 / maxb;
	}
	for (j = 0; j < size; j++) {
		for (i = 0; i < j; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < i; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
		}
		maxb = 0.0;
		for (i = j; i < size; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < j; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
			tmp = wk[i] * fabs(sum);
			if (tmp >= maxb) {
				maxb = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < size; k++) {
				tmp = omtrx[maxi][k];
				omtrx[maxi][k] = omtrx[j][k];
				omtrx[j][k] = tmp;
			}
			aw = -aw;
			wk[maxi] = wk[j];
		}
		index[j] = maxi;
		if (omtrx[j][j] == 0.0)
			omtrx[j][j] = eps;
		if (j != size - 1) {
			tmp = 1.0 / omtrx[j][j];
			for (i = j + 1; i < size; i++)
				omtrx[i][j] *= tmp;
		}
	}
	for (jx = 0; jx < size; jx++) {
		for (ix = 0; ix < size; ix++)
			wk[ix] = 0.0;
		wk[jx] = 1.0;
		l = -1;
		for (i = 0; i < size; i++) {
			idx = index[i];
			sum = wk[idx];
			wk[idx] = wk[i];
			if (l != -1) {
				for (j = l; j < i; j++)
					sum -= omtrx[i][j] * wk[j];
			} else if (sum != 0.0)
				l = i;
			wk[i] = sum;
		}
		for (i = size - 1; i >= 0; i--) {
			sum = wk[i];
			for (j = i + 1; j < size; j++)
				sum -= omtrx[i][j] * wk[j];
			wk[i] = sum / omtrx[i][i];
		}
		for (ix = 0; ix < size; ix++)
			imtrx[ix][jx] = wk[ix];
	}
	free((char *)wk);
	wk = NULL;
	free_ivector(index);
	free_dmatrix(omtrx);
}


void checkevector(dmatrix evec, dmatrix ivec, int nn)
{
	int i, j, ia, ib, ic, error;
	dmatrix matx;
	double sum;


	matx = new_dmatrix(nn, nn);
	/* multiply matrix of eigenvectors and its inverse */
	for (ia = 0; ia < nn; ia++) {
		for (ic = 0; ic < nn; ic++) {
			sum = 0.0;
			for (ib = 0; ib < nn; ib++) sum += evec[ia][ib] * ivec[ib][ic];
			matx[ia][ic] = sum;
		}
	}
	/* check whether the unitary matrix is obtained */
	error = FALSE;
	for (i = 0; i < nn; i++) {
		for (j = 0; j < nn; j++) {
			if (i == j) {
				if (fabs(matx[i][j] - 1.0) > 1.0e-5)
					error = TRUE;
			} else {
				if (fabs(matx[i][j]) > 1.0e-5)
					error = TRUE;
			}
		}
	}
	if (error) {
		FPRINTF(STDOUTFILE "\nWARNING: Inversion of eigenvector matrix not perfect!\n");
	}
	free_dmatrix(matx);
}


/***************************** exported functions *****************************/


/* compute 1 PAM rate matrix, its eigensystem, and the inverse matrix thereof */
void tranprobmat()
{
	eigensystem(Eval, Evec); /* eigensystem of 1 PAM rate matrix */
	luinverse(Evec, Ievc, tpmradix); /* inverse eigenvectors are in Ievc */
	checkevector(Evec, Ievc, tpmradix); /* check whether inversion was OK */
}


/* compute P(t) */
void tprobmtrx(double arc, dmatrix tpr)
{
	register int i, j, k;
	register double temp;


	for (k = 0; k < tpmradix; k++) {
		temp = exp(arc * Eval[k]);
		for (j = 0; j < tpmradix; j++)
			iexp[k][j] = Ievc[k][j] * temp;
	}
	for (i = 0; i < tpmradix; i++) {
		for (j = 0; j < tpmradix; j++) {
			temp = 0.0;
			for (k = 0; k < tpmradix; k++)
				temp += Evec[i][k] * iexp[k][j];
			tpr[i][j] = fabs(temp);
		}
	}
}


/******************************************************************************/
/* estimation of maximum likelihood distances                                 */
/******************************************************************************/

/* compute total log-likelihood
   input: likelihoods for each site and non-zero rate
   output: total log-likelihood (incl. zero rate category) */
double comptotloglkl(dmatrix cdl)
{
	int k, r;
	double loglkl, fv, fv2, sitelkl;

	loglkl = 0.0;
	fv = 1.0-fracinv;
	fv2 = (1.0-fracinv)/(double) numcats;
	
	if (numcats == 1) {

		for (k = 0; k < Numptrn; k++) {	
		
			/* compute likelihood for pattern k */
			sitelkl = cdl[0][k]*fv;
			if (constpat[k] == TRUE)
				sitelkl += fracinv*Freqtpm[(int) Seqpat[0][k]];
		
			/* total log-likelihood */
			loglkl += log(sitelkl)*Weight[k];
		
		}

	} else {
	
		for (k = 0; k < Numptrn; k++) {
			
			/* this general routine works always but it's better 
               to run it only when it's really necessary */
		
			/* compute likelihood for pattern k */
			sitelkl = 0.0;
			for (r = 0; r < numcats; r++)
				sitelkl += cdl[r][k];
			sitelkl = fv2*sitelkl;
			if (constpat[k] == TRUE)
				sitelkl += fracinv*Freqtpm[(int) Seqpat[0][k]];
		
			/* total log-likelihood */
			loglkl += log(sitelkl)*Weight[k];
		
		}

	}
	
	return loglkl;
}


/* computes the site log-likelihoods 
   input: likelihoods for each site and non-zero rate
   output: log-likelihood for each site */
void allsitelkl(dmatrix cdl, dvector aslkl)
{
	int k, r;
	double fv, fv2, sitelkl;

	fv = 1.0-fracinv;
	fv2 = (1.0-fracinv)/(double) numcats;
	
	if (numcats == 1) {

		for (k = 0; k < Numptrn; k++) {	
		
			/* compute likelihood for pattern k */
			sitelkl = cdl[0][k]*fv;
			if (constpat[k] == TRUE)
				sitelkl += fracinv*Freqtpm[(int) Seqpat[0][k]];
		
			/* site log-likelihood */
			aslkl[k] = log(sitelkl);
		}

	} else {
	
		for (k = 0; k < Numptrn; k++) {
			
			/* this general routine works always but it's better 
               to run it only when it's really necessary */
		
			/* compute likelihood for pattern k */
			sitelkl = 0.0;
			for (r = 0; r < numcats; r++)
				sitelkl += cdl[r][k];
			sitelkl = fv2*sitelkl;
			if (constpat[k] == TRUE)
				sitelkl += fracinv*Freqtpm[(int) Seqpat[0][k]];
		
			/* total log-likelihood */
			aslkl[k] = log(sitelkl);
		
		}
	}
}


/***************************** internal functions *****************************/

/* compute negative log-likelihood of distance arc between sequences seqchi/j */
double pairlkl(double arc)
{
	int k, r, ci, cj;
	double loglkl, fv, sitelkl;

	
	/* compute tpms */
	for (r = 0; r < numcats; r++)
		/* compute tpm for rate category r */
		tprobmtrx(arc*Rates[r], ltprobr[r]);

	loglkl = 0.0;
	fv = 1.0-fracinv;
	
	if (numcats == 1) {

		for (k = 0; k < Numptrn; k++) {	
		
			/* compute likelihood for site k */
			ci = seqchi[k];
			cj = seqchj[k];
			if (ci != tpmradix && cj != tpmradix)
				sitelkl = ltprobr[0][ci][cj]*fv;
			else
				sitelkl = fv;
			if (ci == cj && ci != tpmradix)
				sitelkl += fracinv*Freqtpm[ci];
		
			/* total log-likelihood */
			loglkl += log(sitelkl)*Weight[k];
		
		}

	} else {
	
		for (k = 0; k < Numptrn; k++) {
			
			/* this general routine works always but it's better 
               to run it only when it's really necessary */
		
			/* compute likelihood for site k */
			ci = seqchi[k];
			cj = seqchj[k];
			if (ci != tpmradix && cj != tpmradix) {
				sitelkl = 0.0;
				for (r = 0; r < numcats; r++)
					sitelkl += ltprobr[r][ci][cj];
				sitelkl = fv*sitelkl/(double) numcats;	
			} else
				sitelkl = fv;
			if (ci == cj && ci != tpmradix)
				sitelkl += fracinv*Freqtpm[ci];
		
			/* total log-likelihood */
			loglkl += log(sitelkl)*Weight[k];
		
		}

	}
	
	/* return negative log-likelihood as we use a minimizing procedure */
	return -loglkl;
}


/***************************** exported functions *****************************/


/* maximum likelihood distance between sequence i and j */
double mldistance(int i, int j)
{
	double dist, fx, f2x;
	
	if (i == j) return 0.0;
	
	/* use old distance as start value */	
	dist = Distanmat[i][j];

	if (dist == 0.0) return 0.0;

	seqchi = Seqpat[i];
	seqchj = Seqpat[j];

	if (dist <= MINARC) dist = MINARC+1.0;
	if (dist >= MAXARC) dist = MAXARC-1.0;
	
 	dist = onedimenmin(MINARC, dist, MAXARC, pairlkl, EPSILON, &fx, &f2x);
  	
	return dist;
}


/* initialize distance matrix */
void initdistan()
{
	int i, j, k, diff, x, y;
	double obs, temp;
	
	for (i = 0; i < Maxspc; i++) {
		Distanmat[i][i] = 0.0;
		for (j = i + 1; j < Maxspc; j++) {
			seqchi = Seqpat[i];
			seqchj = Seqpat[j];
			
			/* count observed differences */
			diff = 0;
			for (k = 0; k < Numptrn; k++) {
				x = seqchi[k];
				y = seqchj[k];
				if (x != y &&
					x != tpmradix &&
					y != tpmradix)
					diff += Weight[k];
			}
			if (diff == 0)
				Distanmat[i][j] = 0.0;
			else {
				/* use generalized JC correction to get first estimate
				   (for the SH model the observed distance is used) */
				/* observed distance */
				obs = (double) diff / (double) Maxsite;
				temp = 1.0 - (double) obs*tpmradix/(tpmradix-1.0);
				if (temp > 0.0 && !(data_optn == 0 && SH_optn))
					/* use JC corrected distance */
					Distanmat[i][j] = -100.0*(tpmradix-1.0)/tpmradix * log(temp);
				else
					/* use observed distance */
					Distanmat[i][j] = obs * 100.0;
				if (Distanmat[i][j] < MINARC) Distanmat[i][j] = MINARC;
				if (Distanmat[i][j] > MAXARC) Distanmat[i][j] = MAXARC;
			}
			Distanmat[j][i] = Distanmat[i][j];
		}
	}
}

/* compute distance matrix */
void computedistan()
{
	int i, j;

	for (i = 0; i < Maxspc; i++)
		for (j = i; j < Maxspc; j++) {
			Distanmat[i][j] = mldistance(i, j);
			Distanmat[j][i] = Distanmat[i][j];
		}
}


/******************************************************************************/
/* computation of maximum likelihood edge lengths for a given tree            */
/******************************************************************************/


/***************************** internal functions *****************************/


/* multiply partial likelihoods */
void productpartials(Node *op)
{
	Node *cp;
	int i, j, r;
	dcube opc, cpc;

	cp = op;
	opc = op->partials;
	while (cp->isop->isop != op) {
		cp = cp->isop;
		cpc = cp->partials;
		for (r = 0; r < numcats; r++)
			for (i = 0; i < Numptrn; i++)
				for (j = 0; j < tpmradix; j++)
					opc[r][i][j] *= cpc[r][i][j];
	}
}


/* compute internal partial likelihoods */
void partialsinternal(Node *op)
{
	int i, j, k, r;
	double sum;
	dcube oprob, cprob;

	if (clockmode == 1) { /* clocklike branch lengths */
		for (r = 0; r < numcats; r++) {
			tprobmtrx((op->lengthc)*Rates[r], ltprobr[r]);
		}
	} else {  /* non-clocklike branch lengths */
		for (r = 0; r < numcats; r++) {
			tprobmtrx((op->length)*Rates[r], ltprobr[r]);
		}
	}
	
	oprob = op->partials;
	cprob = op->kinp->isop->partials;
	for (r = 0; r < numcats; r++) {
		for (k = 0; k < Numptrn; k++) {
			for (i = 0; i < tpmradix; i++) {
				sum = 0.0;
				for (j = 0; j < tpmradix; j++)
					sum += ltprobr[r][i][j] * cprob[r][k][j];
				oprob[r][k][i] = sum;
			}
		}
	}
}


/* compute external partial likelihoods */
void partialsexternal(Node *op)
{
	int i, j, k, r;
	dcube oprob;
	cvector dseqi;

	if (clockmode == 1) { /* clocklike branch lengths */
		for (r = 0; r < numcats; r++) {
			tprobmtrx((op->lengthc)*Rates[r], ltprobr[r]);
		}
	} else {  /* nonclocklike branch lengths */
		for (r = 0; r < numcats; r++) {
			tprobmtrx((op->length)*Rates[r], ltprobr[r]);
		}
	}

	oprob = op->partials;
	dseqi = op->kinp->eprob;
	for (r = 0; r < numcats; r++) {
		for (k = 0; k < Numptrn; k++) {
			if ((j = dseqi[k]) == tpmradix) {
				for (i = 0; i < tpmradix; i++)
					oprob[r][k][i] = 1.0;
			} else {
				for (i = 0; i < tpmradix; i++)
					oprob[r][k][i] = ltprobr[r][i][j];
			}
		}
	}
}


/* compute all partial likelihoods */
void initpartials(Tree *tr)
{
	Node *cp, *rp;

	cp = rp = tr->rootp;
	do {
		cp = cp->isop->kinp;
		if (cp->isop == NULL) { /* external node */			
			cp = cp->kinp; /* not descen */
			partialsexternal(cp);
		} else { /* internal node */	
			if (!cp->descen) {
				productpartials(cp->kinp->isop);
				partialsinternal(cp);
			}
		}
	} while (cp != rp);
}


/* compute log-likelihood given internal branch with length arc
   between partials partiali and partials partialj */
double intlkl(double arc)
{
	double sumlk, slk;
	int r, s, i, j;
	dmatrix cdl;

	cdl = Ctree->condlkl;
	for (r = 0; r < numcats; r++) {
		tprobmtrx(arc*Rates[r], ltprobr[r]);
	}
	for (r = 0; r < numcats; r++) {
		for (s = 0; s < Numptrn; s++) {
			sumlk = 0.0;
			for (i = 0; i < tpmradix; i++) {			
				slk = 0.0;
				for (j = 0; j < tpmradix; j++) 
					slk += partialj[r][s][j] * ltprobr[r][i][j];		
				sumlk += Freqtpm[i] * partiali[r][s][i] * slk;	
			}
			cdl[r][s] = sumlk;
		}
	}

	/* compute total log-likelihood for current tree */
	Ctree->lklhd = comptotloglkl(cdl);

	return -(Ctree->lklhd); /* we use a minimizing procedure */
}


/* optimize internal branch */
void optinternalbranch(Node *op)
{
	double arc, fx, f2x;

	partiali = op->isop->partials;
	partialj = op->kinp->isop->partials;
	arc = op->length; /* nonclocklike branch lengths */
	if (arc <= MINARC) arc = MINARC+1.0;
	if (arc >= MAXARC) arc = MAXARC-1.0;
	arc = onedimenmin(MINARC, arc, MAXARC, intlkl, EPSILON, &fx, &f2x);
	op->kinp->length = arc;
	op->length = arc;
	
	/* variance of branch length */
	f2x = fabs(f2x);
	if (1.0/(MAXARC*MAXARC) < f2x)
		op->varlen = 1.0/f2x;
	else
		op->varlen = MAXARC*MAXARC;
}


/* compute log-likelihood given external branch with length arc
   between partials partiali and sequence seqchi */
double extlkl(double arc)
{
	double sumlk;
	int r, s, i, j;
	dvector opb;
	dmatrix cdl;
	
	cdl = Ctree->condlkl;
	for (r = 0; r < numcats; r++) {
		tprobmtrx(arc*Rates[r], ltprobr[r]);
	}
	for (r = 0; r < numcats; r++) {
		for (s = 0; s < Numptrn; s++) {
			opb = partiali[r][s];
			sumlk = 0.0;
			if ((j = seqchi[s]) != tpmradix) {
				for (i = 0; i < tpmradix; i++)
					sumlk += (Freqtpm[i] * (opb[i] * ltprobr[r][i][j]));
			} else {
				for (i = 0; i < tpmradix; i++)
					sumlk += Freqtpm[i] * opb[i];
			}
			cdl[r][s] = sumlk;
		}
	}
	
	/* compute total log-likelihood for current tree */
	Ctree->lklhd = comptotloglkl(cdl);

	return -(Ctree->lklhd); /* we use a minimizing procedure */
}

/* optimize external branch */
void optexternalbranch(Node *op)
{
	double arc, fx, f2x;

	partiali = op->isop->partials;
	seqchi = op->kinp->eprob;
	arc = op->length; /* nonclocklike branch lengths */
	if (arc <= MINARC) arc = MINARC+1.0;
	if (arc >= MAXARC) arc = MAXARC-1.0;
	arc = onedimenmin(MINARC, arc, MAXARC, extlkl, EPSILON, &fx, &f2x);
	op->kinp->length = arc;
	op->length = arc;
	
	 /* variance of branch length */
	f2x = fabs(f2x);
	if (1.0/(MAXARC*MAXARC) < f2x)
		op->varlen = 1.0/f2x;
	else
		op->varlen = MAXARC*MAXARC;
}


/* finish likelihoods for each rate and site */
void finishlkl(Node *op)
{
	int r, k, i, j;
	double arc, sumlk, slk;
	dmatrix cdl;

	partiali = op->isop->partials;
	partialj = op->kinp->isop->partials;
	cdl = Ctree->condlkl;
	arc = op->length; /* nonclocklike branch lengths */
	for (r = 0; r < numcats; r++) {
		tprobmtrx(arc*Rates[r], ltprobr[r]);
	}
	for (r = 0; r < numcats; r++) {
		for (k = 0; k < Numptrn; k++) {
			sumlk = 0.0;
			for (i = 0; i < tpmradix; i++) {
				slk = 0.0;
				for (j = 0; j < tpmradix; j++)
					slk += partialj[r][k][j] * ltprobr[r][i][j];
				sumlk += Freqtpm[i] * partiali[r][k][i] * slk;
			}			
			cdl[r][k] = sumlk;
		}
	}
}


/***************************** exported functions *****************************/


/* optimize branch lengths to get maximum likelihood (nonclocklike branchs) */
double optlkl(Tree *tr)
{
	Node *cp, *rp;
	int nconv;
	double lendiff;
	
	clockmode = 0; /* nonclocklike branch lengths */
	nconv = 0;
	Converg = FALSE;
	initpartials(tr);
	for (Numit = 1; (Numit <= MAXIT) && (!Converg); Numit++) {
		
		cp = rp = tr->rootp;
		do {
			cp = cp->isop->kinp;
			productpartials(cp->kinp->isop);
			if (cp->isop == NULL) { /* external node */	
				cp = cp->kinp; /* not descen */
				
				lendiff = cp->length;
				optexternalbranch(cp);
				lendiff = fabs(lendiff - cp->length);
				if (lendiff < EPSILON) nconv++;
				else nconv = 0;				
				
				partialsexternal(cp);
			} else { /* internal node */
				if (cp->descen) {
					partialsinternal(cp);
				} else {
					
					lendiff = cp->length;
					optinternalbranch(cp);
					lendiff = fabs(lendiff - cp->length);
					if (lendiff < EPSILON) nconv++;
					else nconv = 0;
					
					/* eventually compute likelihoods for each site */
					if ((cp->number == Numibrnch-1 && lendiff < EPSILON) ||
					Numit == MAXIT-1) finishlkl(cp);

					partialsinternal(cp);
				}
			}			
			if (nconv >= Numbrnch) { /* convergence */
				Converg = TRUE;
				cp = rp; /* get out of here */
			}
		} while (cp != rp);
	}

	/* compute total log-likelihood for current tree */
	return comptotloglkl(tr->condlkl);
}


/* compute likelihood of tree for given branch lengths */
double treelkl(Tree *tr)
{
	int i, k, r;
	Node *cp;
	dmatrix cdl;
	dcube prob1, prob2;
	double sumlk;

	/* compute for each site and rate log-likelihoods */
	initpartials(tr);
	cp = tr->rootp;
	productpartials(cp->isop);
	prob1 = cp->partials;
	prob2 = cp->isop->partials;
	cdl = tr->condlkl;
	for (r = 0; r < numcats; r++) {
		for (k = 0; k < Numptrn; k++) {
			sumlk = 0.0;
			for (i = 0; i < tpmradix; i++)
				sumlk += Freqtpm[i] * (prob1[r][k][i] * prob2[r][k][i]);
			cdl[r][k] = sumlk;
		}
	}	
	
	/* return total log-likelihood for current tree */
	return comptotloglkl(cdl);
}


/******************************************************************************/
/* least-squares estimate of branch lengths                                   */
/******************************************************************************/


/***************************** internal functions *****************************/


void luequation(dmatrix amat, dvector yvec, int size)
{
    double eps = 1.0e-20; /* ! */
	int i, j, k, l, maxi=0, idx;
	double sum, tmp, maxb, aw;
	dvector wk;
	ivector index;


	wk = new_dvector(size);
	index = new_ivector(size);
	aw = 1.0;
	for (i = 0; i < size; i++) {
		maxb = 0.0;
		for (j = 0; j < size; j++) {
			if (fabs(amat[i][j]) > maxb)
				maxb = fabs(amat[i][j]);
		}
		if (maxb == 0.0) {
			/* Singular matrix */
			FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR D TO DEVELOPERS\n\n\n");
			exit(1);
		}
		wk[i] = 1.0 / maxb;
	}
	for (j = 0; j < size; j++) {
		for (i = 0; i < j; i++) {
			sum = amat[i][j];
			for (k = 0; k < i; k++)
				sum -= amat[i][k] * amat[k][j];
			amat[i][j] = sum;
		}
		maxb = 0.0;
		for (i = j; i < size; i++) {
			sum = amat[i][j];
			for (k = 0; k < j; k++)
				sum -= amat[i][k] * amat[k][j];
			amat[i][j] = sum;
			tmp = wk[i] * fabs(sum);
			if (tmp >= maxb) {
				maxb = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < size; k++) {
				tmp = amat[maxi][k];
				amat[maxi][k] = amat[j][k];
				amat[j][k] = tmp;
			}
			aw = -aw;
			wk[maxi] = wk[j];
		}
		index[j] = maxi;
		if (amat[j][j] == 0.0)
			amat[j][j] = eps;
		if (j != size - 1) {
			tmp = 1.0 / amat[j][j];
			for (i = j + 1; i < size; i++)
				amat[i][j] *= tmp;
		}
	}
	l = -1;
	for (i = 0; i < size; i++) {
		idx = index[i];
		sum = yvec[idx];
		yvec[idx] = yvec[i];
		if (l != -1) {
			for (j = l; j < i; j++)
				sum -= amat[i][j] * yvec[j];
		} else if (sum != 0.0)
			l = i;
		yvec[i] = sum;
	}
	for (i = size - 1; i >= 0; i--) {
		sum = yvec[i];
		for (j = i + 1; j < size; j++)
			sum -= amat[i][j] * yvec[j];
		yvec[i] = sum / amat[i][i];
	}
	free_ivector(index);
	free_dvector(wk);
}


/* least square estimation of branch lengths
   used for the approximate ML and as starting point 
   in the calculation of the exact value of the ML */
void lslength(Tree *tr, dvector distanvec, int numspc, int numibrnch, dvector Brnlength)
{
	int i, i1, j, j1, j2, k, numbrnch, numpair;
	double sum, leng, alllen, rss;
	ivector pths;
	dmatrix atmt, atamt;
	Node **ebp, **ibp;

	numbrnch = numspc + numibrnch;
	numpair = (numspc * (numspc - 1)) / 2;
	atmt = new_dmatrix(numbrnch, numpair);
	atamt = new_dmatrix(numbrnch, numbrnch);
	ebp = tr->ebrnchp;
	ibp = tr->ibrnchp;
	for (i = 0; i < numspc; i++) {
		for (j1 = 1, j = 0; j1 < numspc; j1++) {
			if (j1 == i) {
				for (j2 = 0; j2 < j1; j2++, j++) {
					atmt[i][j] = 1.0;
				}
			} else {
				for (j2 = 0; j2 < j1; j2++, j++) {
					if (j2 == i)
						atmt[i][j] = 1.0;
					else
						atmt[i][j] = 0.0;
				}
			}
		}
	}
	for (i1 = 0, i = numspc; i1 < numibrnch; i1++, i++) {
		pths = ibp[i1]->paths;
		for (j1 = 1, j = 0; j1 < numspc; j1++) {
			for (j2 = 0; j2 < j1; j2++, j++) {
				if (pths[j1] != pths[j2])
					atmt[i][j] = 1.0;
				else
					atmt[i][j] = 0.0;
			}
		}
	}
	for (i = 0; i < numbrnch; i++) {
		for (j = 0; j <= i; j++) {
			for (k = 0, sum = 0.0; k < numpair; k++)
				sum += atmt[i][k] * atmt[j][k];
			atamt[i][j] = sum;
			atamt[j][i] = sum;
		}
	}
	for (i = 0; i < numbrnch; i++) {
		for (k = 0, sum = 0.0; k < numpair; k++)
			sum += atmt[i][k] * distanvec[k];
		Brnlength[i] = sum;
	}
	luequation(atamt, Brnlength, numbrnch);
	for (i = 0, rss = 0.0; i < numpair; i++) {
		sum = distanvec[i];
		for (j = 0; j < numbrnch; j++) {
			if (atmt[j][i] == 1.0 && Brnlength[j] > 0.0)
				sum -= Brnlength[j];
		}
		rss += sum * sum;
	}
	tr->rssleast = sqrt(rss);
	alllen = 0.0;
	for (i = 0; i < numspc; i++) {
		leng = Brnlength[i];
		alllen += leng;
		if (leng < MINARC) leng = MINARC;
		if (leng > MAXARC) leng = MAXARC;
		if (clockmode) { /* clock */
			ebp[i]->lengthc = leng;
			ebp[i]->kinp->lengthc = leng;
		} else { /* no clock */
			ebp[i]->length = leng;
			ebp[i]->kinp->length = leng;
		}	
		Brnlength[i] = leng;
	}
	for (i = 0, j = numspc; i < numibrnch; i++, j++) {
		leng = Brnlength[j];
		alllen += leng;
		if (leng < MINARC) leng = MINARC;
		if (leng > MAXARC) leng = MAXARC;
		if (clockmode) { /* clock */
			ibp[i]->lengthc = leng;
			ibp[i]->kinp->lengthc = leng;
		} else { /* no clock */
			ibp[i]->length = leng;
			ibp[i]->kinp->length = leng;
		}
		Brnlength[j] = leng;
	}
	free_dmatrix(atmt);
	free_dmatrix(atamt);
}
