/*
 * ml3.c
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


/* prototypes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "ml.h"
#include "gamma.h"



/******************************************************************************/
/* discrete Gamma-distribution and related stuff                                              */
/******************************************************************************/

/* compare general base frequencies with frequencies of taxon i with chi square */
double homogentest(int taxon)
{
	return chi2test(Freqtpm, Basecomp[taxon], gettpmradix(), &chi2fail);
}


/* discrete Gamma according to Yang 1994 (JME 39:306-314) */
void YangDiscreteGamma (double shape, int c, dvector x)
{
	double twoc, mu;
	int i;
	
	twoc = 2.0*c;
	mu = 0.0;
	for (i = 0; i < c; i++)
	{
		/* corresponding rates */
		x[i] = icdfGamma ( (2.0*i+1.0)/twoc, shape);
		mu += x[i];
	}
	mu = mu/c;	

	/* rescale for avarage rate of 1.0 */
	for (i = 0; i < c; i++)
	{
		x[i] /= mu;
	}
}

/* compute rates of each category when rates are Gamma-distributed */
void updaterates()
{
	int i;
	double alpha;
	
	if (numcats == 1)
	{
		Rates[0] = 1.0;
		return;
	}
	if (Geta == 0.0)
	{
		for (i = 0; i < numcats; i++)
			Rates[i] = 1.0;
		return;
	}
	alpha = (1.0 - Geta)/Geta;

	YangDiscreteGamma (alpha, numcats, Rates);
	
	/* if invariable sites are present */
	for (i = 0; i < numcats; i++)
		Rates[i] = Rates[i]/(1.0-fracinv);
	
	/* check for very small rates */
	for (i = 0; i < numcats; i++)
		if (Rates[i] < 0.000001) Rates[i] = 0.000001;
}



/******************************************************************************/
/* parameter estimation                                                       */
/******************************************************************************/

/* compute sample mean and standard deviation of sample mean */
void computestat(double *data, int n, double *mean, double *err)
{	
	int i;
	double sum;
	
	sum = 0;
	for (i = 0; i < n; i++) sum += data[i];
	(*mean) = sum/(double) n;
	
	sum = 0;
	for (i = 0; i < n; i++) sum += (data[i] - (*mean))*(data[i] - (*mean));
	if (n != 1)
		(*err) = sqrt(sum)/sqrt((double)(n-1)*n); /* unbiased estimator */
	else
		(*err) = 0.0; /* if n == 1 */
}

/* compute ML value of quartet (a,b,c,d) */
double quartetml(int a, int b, int c, int d)
{
	double d1, d2, d3;

	/* compute ML for all topologies */
	if (approxp_optn) { /* approximate parameter mode */
		d1 = quartet_alklhd(a,b,c,d); /* (a,b)-(c,d) */
		d2 = quartet_alklhd(a,c,b,d); /* (a,c)-(b,d) */
		d3 = quartet_alklhd(a,d,b,c); /* (a,d)-(b,c) */
	} else {
		d1 = quartet_lklhd(a,b,c,d); /* (a,b)-(c,d) */
		d2 = quartet_lklhd(a,c,b,d); /* (a,c)-(b,d) */
		d3 = quartet_lklhd(a,d,b,c); /* (a,d)-(b,c) */
	}
	
	/* looking for max(d1, d2, d3) */						
	if (d1 < d2) { /* d2 > d1 */
		if (d2 < d3) { /* d3 > d2 > d1 */
			/* d3 maximum */ 
			return d3;
		} else { /* d2 >= d3 > d1 */
			/* d2 maximum */
			return d2;
		}
	} else { /* d1 >= d2 */
		if (d1 < d3) { /* d3 > d1 >= d2 */
			/* d3 maximum */
			return d3;
		} else { /* d1 >= d2 && d1 >= d3 */
			/* d1 maximum */
			return d1;
		}
	}
}

/* optimization function TSparam - quartets */
double opttsq(double x)
{
	if (x < MINTS) TSparam = MINTS;
	else if (x > MAXTS) TSparam = MAXTS;
	else TSparam = x;
	tranprobmat();
	distupdate(qca, qcb, qcc, qcd);	
	return (-quartetml(qca, qcb, qcc, qcd));
}

/* optimization function YRparam - quartets */
double optyrq(double x)
{
	if (x < MINYR) YRparam = MINYR;
	else if (x > MAXYR) YRparam = MAXYR;
	else YRparam = x;
	tranprobmat();
	distupdate(qca, qcb, qcc, qcd);
	return (-quartetml(qca, qcb, qcc, qcd));
}

/* estimate substitution process parameters - random quartets */
void optimseqevolparamsq()
{
	double tsmeanold, yrmeanold;
	dvector tslist, yrlist;
	int fin;
	ivector taxon;
	uli minqts, maxqts, n;


	taxon = new_ivector(4);

	/* number of quartets to be investigated */
	minqts = (uli) floor(0.25 * MINPERTAXUM * Maxspc) + 1;
	maxqts = (uli) floor(0.25 * MAXPERTAXUM * Maxspc) + 1;
	if (Maxspc == 4) {
		minqts = (uli) 1;
		maxqts = (uli) 1;
	}
	
	tslist = new_dvector(maxqts);
	yrlist = new_dvector(maxqts);
	
	/* initialize averages */
	tsmean = TSparam;
	yrmean = YRparam;

	fin = FALSE;
		
	/* investigate maxqts random quartets */
	for (n = 0; n < maxqts; n++) {
	
		/* choose random quartet */
		chooser(Maxspc, 4, taxon);

		/*
		 * optimize parameters on this quartet
		 */

		qca = taxon[0];
		qcb = taxon[1];
		qcc = taxon[2];
		qcd = taxon[3];
		
		/* initialize start values with average value */
		if ((SH_optn || nuc_optn) && optim_optn && (data_optn == 0)) TSparam = tsmean;
		if ((nuc_optn && TN_optn) && optim_optn && (data_optn == 0)) YRparam = yrmean;
		
		/* estimation */
		twodimenmin(PEPS1,
			(SH_optn || nuc_optn) && optim_optn && (data_optn == 0),
				MINTS, &TSparam, MAXTS, opttsq, &tserr,
			(nuc_optn && TN_optn) && optim_optn && (data_optn == 0),
				MINYR, &YRparam, MAXYR, optyrq, &yrerr);


		tsmeanold = tsmean;
		yrmeanold = yrmean;
		tslist[n] = TSparam;
		yrlist[n] = YRparam;
		computestat(tslist, n+1 , &tsmean, &tserr);
		computestat(yrlist, n+1 , &yrmean, &yrerr);

		/* check whether the means are converging */
		if (n > minqts-2) {
			if ((fabs(tsmean-tsmeanold) < TSDIFF) &&
				(fabs(yrmean-yrmeanold) < YRDIFF))
					fin = TRUE;
		}

		/* investigate at least minqts quartets */
		if (n > minqts-2 && (fin || n > maxqts-2)) break;
	}

	/* round estimated numbers to 2 digits after the decimal point */
	if (tserr != 0.0) tsmean = floor(100.0*tsmean+0.5)/100.0;
	if (yrerr != 0.0) yrmean = floor(100.0*yrmean+0.5)/100.0;

	/* update ML engine */
	TSparam = tsmean;
	YRparam = yrmean;
	tranprobmat();

	free_ivector(taxon);
}

/* optimization function TSparam - tree */
double opttst(double x)
{
	double result;

	if (x < MINTS) TSparam = MINTS;
	else if (x > MAXTS) TSparam = MAXTS;
	else TSparam = x;
	tranprobmat();
	computedistan();
	if (approxp_optn) result = usertree_alklhd();
	else result = usertree_lklhd();

	return (-result);
}

/* optimization function YRparam - tree */
double optyrt(double x)
{
	double result;
	
	if (x < MINYR) YRparam = MINYR;
	else if (x > MAXYR) YRparam = MAXYR;
	else YRparam = x;
	tranprobmat();
	computedistan();
	if (approxp_optn) result = usertree_alklhd();
	else result = usertree_lklhd();

	return (-result);
}


/* optimize substitution process parameters - tree */
void optimseqevolparamst()
{
	twodimenmin(PEPS1,
		(SH_optn || nuc_optn) && optim_optn && (data_optn == 0),
			MINTS, &TSparam, MAXTS, opttst, &tserr,
		(nuc_optn && TN_optn) && optim_optn && (data_optn == 0),
			MINYR, &YRparam, MAXYR, optyrt, &yrerr);
}


/* optimization function fracinv */
double optfi(double x)
{
	double result;

	if (x < MINFI) fracinv = MINFI;
	else if (x > MAXFI) fracinv = MAXFI;
	else fracinv = x;
	
	computedistan();
	if (approxp_optn) result = usertree_alklhd();
	else result = usertree_lklhd();

	return (-result);	
}


/* optimization function Geta */
double optge(double x)
{
	double result;

	if (x < MINGE) Geta = MINGE;
	else if (x > MAXGE) Geta = MAXGE;
	else Geta = x;
	
	updaterates();
	
	computedistan();
	if (approxp_optn) result = usertree_alklhd();
	else result = usertree_lklhd();

	return (-result);	
}


/* optimize rate heterogeneity parameters */
void optimrateparams()
{	
	twodimenmin(PEPS2,
		fracinv_optim,
			MINFI, &fracinv, fracconst, optfi, &fierr,
		grate_optim,
			MINGE, &Geta, MAXGE, optge, &geerr);

}
