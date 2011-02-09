/*
 * util.c
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


#include "util.h"

#define STDOUT     stdout
#ifndef PARALLEL             /* because printf() runs significantly faster */
                             /* than fprintf(stdout) on an Apple McIntosh  */
                             /* (HS) */
#       define FPRINTF    printf
#       define STDOUTFILE
#else
#       define FPRINTF    fprintf
#       define STDOUTFILE STDOUT,
	extern int PP_NumProcs;
	extern int PP_Myid;
        long int PP_randn;
        long int PP_rand;
#endif


/*
 * memory allocation error handler
 */

void maerror(char *message)
{
	FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (lack of memory: %s)\n\n", message);
	FPRINTF(STDOUTFILE "Hint for Macintosh users:\n");
	FPRINTF(STDOUTFILE "Use the <Get Info> command of the Finder to increase the memory partition!\n\n");
	exit(1);
}


/*
 * memory allocate double vectors, matrices, and cubes
 */

dvector new_dvector(int n)
{
	dvector v;

	v = (dvector) malloc((unsigned) (n * sizeof(double)));
	if (v == NULL) maerror("step 1 in new_dvector");

	return v;
}

dmatrix new_dmatrix(int nrow, int ncol)
{
	int i;
	dmatrix m;

	m = (dmatrix) malloc((unsigned) (nrow * sizeof(dvector)));
	if (m == NULL) maerror("step 1 in in new_dmatrix");

	*m = (dvector) malloc((unsigned) (nrow * ncol * sizeof(double)));
	if (*m == NULL) maerror("step 2 in in new_dmatrix");

	for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

	return m;
}

dcube new_dcube(int ntri, int nrow, int ncol)
{
	int i, j;
	dcube c;

	c = (dcube) malloc((unsigned) (ntri * sizeof(dmatrix)));
	if (c == NULL) maerror("step 1 in in new_dcube");

	*c = (dmatrix) malloc((unsigned) (ntri * nrow * sizeof(dvector)));
	if (*c == NULL) maerror("step 2 in in new_dcube");

	**c = (dvector) malloc((unsigned) (ntri * nrow * ncol * sizeof(double)));
	if (**c == NULL) maerror("step 3 in in new_dcube");

	for (j = 1; j < nrow; j++) c[0][j] = c[0][j-1] + ncol;

	for (i = 1; i < ntri; i++) {
		c[i] = c[i-1] + nrow;
		c[i][0] = c[i-1][0] + nrow * ncol;
		for (j = 1; j < nrow; j++) c[i][j] = c[i][j-1] + ncol;
	}

	return c;
}

void free_dvector(dvector v)
{
	free((double *) v);
}

void free_dmatrix(dmatrix m)
{
	free((double *) *m);
	free((double *) m);
}

void free_dcube(dcube c)
{
	free((double *) **c);
	free((double *) *c);
	free((double *) c);
}


/*
 * memory allocate char vectors, matrices, and cubes
 */

cvector new_cvector(int n)
{
	cvector v;

	v = (cvector) malloc((unsigned)n * sizeof(char));
	if (v == NULL) maerror("step1 in new_cvector");

	return v;
}

cmatrix new_cmatrix(int nrow, int ncol)
{
	int i;
	cmatrix m;

	m = (cmatrix) malloc((unsigned) (nrow * sizeof(cvector)));
	if (m == NULL) maerror("step 1 in new_cmatrix");

	*m = (cvector) malloc((unsigned) (nrow * ncol * sizeof(char)));
	if (*m == NULL) maerror("step 2 in new_cmatrix");

	for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

	return m;
}

ccube new_ccube(int ntri, int nrow, int ncol)
{
	int i, j;
	ccube c;

	c = (ccube) malloc((unsigned) (ntri * sizeof(cmatrix)));
	if (c == NULL) maerror("step 1 in new_ccube");

	*c = (cmatrix) malloc((unsigned) (ntri * nrow * sizeof(cvector)));
	if (*c == NULL) maerror("step 2 in new_ccube");

	**c = (cvector) malloc((unsigned) (ntri * nrow * ncol * sizeof(char)));
	if (**c == NULL) maerror("step 3 in new_ccube");

	for (j = 1; j < nrow; j++) c[0][j] = c[0][j-1] + ncol;

	for (i = 1; i < ntri; i++) {
		c[i] = c[i-1] + nrow;
		c[i][0] = c[i-1][0] + nrow * ncol;
		for (j = 1; j < nrow; j++) c[i][j] = c[i][j-1] + ncol;
	}

	return c;
}

void free_cvector(cvector v)
{
	free((char *) v);
}

void free_cmatrix(cmatrix m)
{
	free((char *) *m);
	free((char *) m);
}

void free_ccube(ccube c)
{
	free((char *) **c);
	free((char *) *c);
	free((char *) c);
}


/*
 * memory allocate int vectors, matrices, and cubes
 */

ivector new_ivector(int n)
{
	ivector v;

	v = (ivector) malloc((unsigned) (n * sizeof(int)));
	if (v == NULL) maerror("step 1 in new_ivector");

	return v;
}

imatrix new_imatrix(int nrow, int ncol)
{
	int i;
	imatrix m;

	m = (imatrix) malloc((unsigned) (nrow * sizeof(ivector)));
	if (m == NULL) maerror("step 1 in new_imatrix");

	*m = (ivector) malloc((unsigned) (nrow * ncol * sizeof(int)));
	if (*m == NULL) maerror("step 2 in new_imatrix");

	for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

	return m;
}

icube new_icube(int ntri, int nrow, int ncol)
{
	int i, j;
	icube c;

	c = (icube) malloc((unsigned) (ntri * sizeof(imatrix)));
	if (c == NULL) maerror("step 1 in new_icube");

	*c = (imatrix) malloc((unsigned) (ntri * nrow * sizeof(ivector)));
	if (*c == NULL) maerror("step 2 in new_icube");

	**c = (ivector) malloc((unsigned) (ntri * nrow * ncol * sizeof(int)));
	if (**c == NULL) maerror("step 3 in new_icube");

	for (j = 1; j < nrow; j++) c[0][j] = c[0][j-1] + ncol;

	for (i = 1; i < ntri; i++) {
		c[i] = c[i-1] + nrow;
		c[i][0] = c[i-1][0] + nrow * ncol;
		for (j = 1; j < nrow; j++) c[i][j] = c[i][j-1] + ncol;
	}
	
	return c;
}

void free_ivector(ivector v)
{
	free((int *) v);
}

void free_imatrix(imatrix m)
{
	free((int *) *m);
	free((int *) m);
}

void free_icube(icube c)
{
	free((int *) **c);
	free((int *) *c);
	free((int *) c);
}


/*
 * memory allocate uli vectors, matrices, and cubes
 */

ulivector new_ulivector(int n)
{
	ulivector v;

	v = (ulivector) malloc((unsigned) (n * sizeof(uli)));
	if (v == NULL) maerror("step 1 in new_ulivector");

	return v;
}

ulimatrix new_ulimatrix(int nrow, int ncol)
{
	int i;
	ulimatrix m;

	m = (ulimatrix) malloc((unsigned) (nrow * sizeof(ulivector)));
	if (m == NULL) maerror("step 1 in new_ulimatrix");

	*m = (ulivector) malloc((unsigned) (nrow * ncol * sizeof(uli)));
	if (*m == NULL) maerror("step 2 in new_ulimatrix");

	for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

	return m;
}

ulicube new_ulicube(int ntri, int nrow, int ncol)
{
	int i, j;
	ulicube c;

	c = (ulicube) malloc((unsigned) (ntri * sizeof(ulimatrix)));
	if (c == NULL) maerror("step 1 in new_ulicube");

	*c = (ulimatrix) malloc((unsigned) (ntri * nrow * sizeof(ulivector)));
	if (*c == NULL) maerror("step 2 in new_ulicube");

	**c = (ulivector) malloc((unsigned) (ntri * nrow * ncol * sizeof(uli)));
	if (**c == NULL) maerror("step 3 in new_ulicube");

	for (j = 1; j < nrow; j++) c[0][j] = c[0][j-1] + ncol;

	for (i = 1; i < ntri; i++) {
		c[i] = c[i-1] + nrow;
		c[i][0] = c[i-1][0] + nrow * ncol;
		for (j = 1; j < nrow; j++) c[i][j] = c[i][j-1] + ncol;
	}
	
	return c;
}

void free_ulivector(ulivector v)
{
	free((uli *) v);
}

void free_ulimatrix(ulimatrix m)
{
	free((uli *) *m);
	free((uli *) m);
}

void free_ulicube(ulicube c)
{
	free((uli *) **c);
	free((uli *) *c);
	free((uli *) c);
}


/******************************************************************************/
/* random numbers generator  (Numerical recipes)                              */
/******************************************************************************/

/* definitions */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* variable */
long idum;

double randomunitintervall()
/* Long period (> 2e18) random number generator. Returns a uniform random
   deviate between 0.0 and 1.0 (exclusive of endpoint values).

   Source:
   Press et al., "Numerical recipes in C", Cambridge University Press, 1992
   (chapter 7 "Random numbers", ran2 random number generator) */
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (idum <= 0) {
		if (-(idum) < 1)
			idum=1;
		else
			idum=-(idum);
		idum2=(idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0)
				idum += IM1;
			if (j < NTAB)
				iv[j] = idum;
		}
		iy=iv[0];
	}
	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0)
		idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0)
		idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1)
		iy += IMM1;
	if ((temp=AM*iy) > RNMX)
		return RNMX;
	else
		return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

int initrandom(int seed)
{
   srand((unsigned) time(NULL));
   if (seed < 0) 
	seed = rand();
   idum=-(long) seed;
#  ifdef PARALLEL
	{
	int n;
	for (n=0; n<PP_Myid; n++)
		(void) randomunitintervall();
#       ifdef PVERBOSE1
  	   fprintf(stderr, "(%2d) !!! random seed set to %d, %dx drawn !!!\n", PP_Myid, seed, n);
#       endif
	}
#  else
#       ifdef PVERBOSE1
	   fprintf(stderr, "!!! random seed set to %d !!!\n", seed);
#       endif
#  endif
   return (seed);
}  /* initrandom */ 


/* returns a random integer in the range [0; n - 1] */
int randominteger(int n)
{
	int t;
#  ifndef FIXEDINTRAND
#	ifndef PARALLEL
		t = (int) floor(randomunitintervall()*n);
		return t;
#	else
		int m;
		for (m=1; m<PP_NumProcs; m++)
			(void) randomunitintervall();
		PP_randn+=(m-1); PP_rand++;
		return (int) floor(randomunitintervall()*n);
#	endif
#  else
	fprintf(stderr, "!!! fixed \"random\" integers for testing purposes !!!\n");
	return (int)0;
#  endif /* FIXEDINTRAND */
}

/* Draw s numbers from the set 0,1,2,..,t-1 and put them
   into slist (every number can be drawn only one time) */
void chooser(int t, int s, ivector slist)
{
	int i, j, k, l;
	ivector isfree;

	isfree = new_ivector(t);
	for (i = 0; i < t; i++) isfree[i] = TRUE;
	for (i = 0; i < s; i++) {
		/* random integer in the range [0, t-1-i] */
		j = randominteger(t-i);
		k = -1;
		l = -1;
		do {
			k++;
			if (isfree[k] == TRUE) l++; 
		} while ( l != j);
		slist[i] = k;
		isfree[k] = FALSE;
	}
	free_ivector(isfree);
}

/* a realloc function that works also on non-ANSI compilers */ 
void *myrealloc(void *p, size_t s)
{
	if (p == NULL) return malloc(s);
	else return realloc(p, s);
}

/* safer variant of gets */
/* Reads characters from stdin until a newline character or EOF
   is received.  The newline is not made part of the string.
   If an error occurs a null string \0 is returned */
cvector mygets()
{
	int c, n;
	cvector str;

	str = new_cvector(100);
	
	n = 0;
	c = getchar();
	while (c != '\n' && c != '\r' && n < 99 && c != EOF && !ferror(stdin))
	{
		str[n] = (char) c;
		c = getchar();
		n++;
	}
	if (c == EOF || ferror(stdin))
	{
		str[0] = '\0';
	}
	else
	{
		str[n] = '\0';
	}
	
	return str;
}



/******************************************************************************/
/* minimization of a function by Brents method (Numerical Recipes)            */
/******************************************************************************/

double brent(double, double, double, double (*f )(double ), double, double *, double *, double, double, double);


#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Brents method in one dimension */
double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	double *foptx, double *f2optx, double fax, double fbx, double fcx)
{
	int iter;
	double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if (fax < fcx) {
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else {
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;	
	}
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
				(v*v*xw + x*x*wv + w*w*vx);
			return x;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
		(v*v*xw + x*x*wv + w*w*vx);
	return x;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN
#undef GOLD
#undef GLIMIT
#undef TINY

/* one-dimensional minimization - as input a lower and an upper limit and a trial
  value for the minimum is needed: xmin < xguess < xmax
  the function and a fractional tolerance has to be specified
  onedimenmin returns the optimal x value and the value of the function
  and its second derivative at this point
  */
double onedimenmin(double xmin, double xguess, double xmax, double (*f)(double),
	double tol, double *fx, double *f2x)
{
	double eps, optx, ax, bx, cx, fa, fb, fc;
		
	/* first attempt to bracketize minimum */
	eps = xguess*tol*50.0;
	ax = xguess - eps;
	if (ax < xmin) ax = xmin;
	bx = xguess;
	cx = xguess + eps;
	if (cx > xmax) cx = xmax;
	
	/* check if this works */
	fa = (*f)(ax);
	fb = (*f)(bx);
	fc = (*f)(cx);

	/* if it works use these borders else be conservative */
	if ((fa < fb) || (fc < fb)) {
		if (ax != xmin) fa = (*f)(xmin);
		if (cx != xmax) fc = (*f)(xmax);
		optx = brent(xmin, xguess, xmax, f, tol, fx, f2x, fa, fb, fc);
	} else
		optx = brent(ax, bx, cx, f, tol, fx, f2x, fa, fb, fc);

	return optx; /* return optimal x */
}

/* two-dimensional minimization with borders and calculations of standard errors */
/* we optimize along basis vectors - not very optimal but it seems to work well */
void twodimenmin(double tol,
	int active1, double min1, double *x1, double max1, double (*func1)(double), double *err1,
	int active2, double min2, double *x2, double max2, double (*func2)(double), double *err2)
{
	int it, nump, change;
	double x1old, x2old;
	double fx, f2x;

	it = 0;
	nump = 0;

	/* count number of parameters */
	if (active1) nump++;
	if (active2) nump++;

	do { /* repeat until nothing changes any more */
		it++;
		change = FALSE;

		/* optimize first variable */
		if (active1) {

			if ((*x1) <= min1) (*x1) = min1 + 0.2*(max1-min1);
			if ((*x1) >= max1) (*x1) = max1 - 0.2*(max1-min1);
			x1old = (*x1);
			(*x1) = onedimenmin(min1, (*x1), max1, func1, tol, &fx, &f2x);
			if ((*x1) < min1) (*x1) = min1;
			if ((*x1) > max1) (*x1) = max1;
			/* same tolerance as 1D minimization */
			if (fabs((*x1) - x1old) > 3.3*tol) change = TRUE;
			
			/* standard error */
			f2x = fabs(f2x);
			if (1.0/(max1*max1) < f2x) (*err1) = sqrt(1.0/f2x);
			else (*err1) = max1;

		}

		/* optimize second variable */
		if (active2) {

			if ((*x2) <= min2) (*x2) = min2 + 0.2*(max2-min2);
			if ((*x2) >= max2) (*x2) = max2 - 0.2*(max2-min2);
			x2old = (*x2);
			(*x2) = onedimenmin(min2, (*x2), max2, func2, tol, &fx, &f2x);
			if ((*x2) < min2) (*x2) = min2;
			if ((*x2) > max2) (*x2) = max2;
			/* same tolerance as 1D minimization */
			if (fabs((*x2) - x2old) > 3.3*tol) change = TRUE;
			
			/* standard error */
			f2x = fabs(f2x);
			if (1.0/(max2*max2) < f2x) (*err2) = sqrt(1.0/f2x);
			else (*err2) = max2;

		}

		if (nump == 1) return;
		
	} while (it != MAXITS && change);

	return;
}

