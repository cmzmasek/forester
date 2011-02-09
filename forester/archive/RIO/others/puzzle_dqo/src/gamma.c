/*
 * gamma.c
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

#include <math.h>
#include "util.h"
#include "gamma.h"

/* private prototypes */
static double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
static double PointNormal (double prob);
static double PointChi2 (double prob, double v);

/* Gamma density function */
double densityGamma (double x, double shape)
{
	return pow (shape, shape) * pow (x, shape-1) /
	 	exp (shape*x + LnGamma(shape));
}

/* Gamma cdf */
double cdfGamma (double x, double shape)
{
	double result;
	
	result = IncompleteGamma (shape*x, shape, LnGamma(shape));
	
	return result;
}

/* Gamma inverse cdf */
double icdfGamma (double y, double shape)
{
	double result;
	
	result = PointChi2 (y, 2.0*shape)/(2.0*shape);
	
	/* to avoid -1.0 */
	if (result < 0.0)
	{
		result = 0.0;
	}
	
	return result;
}

/* Gamma n-th moment */
double momentGamma (int n, double shape)
{
	int i;
	double tmp = 1.0;
	
	for (i = 1; i < n; i++)
	{
		tmp *= (shape + i)/shape;
	}
	
	return tmp;
}

/* The following code comes from tools.c in Yang's PAML package */

double LnGamma (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673 
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;  
}

static double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];
   
   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}


/* functions concerning the CDF and percentage points of the gamma and
   Chi2 distribution
*/
static double PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}


static double PointChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3: 
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:

   do
   {
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      return (-1);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   }
   while (fabs(q/ch-1) > e);

   return (ch);
}


/* Incomplete Gamma function Q(a,x)
   - this is a cleanroom implementation of NRs gammq(a,x)
*/
double IncompleteGammaQ (double a, double x)
{
	return 1.0-IncompleteGamma (x, a, LnGamma(a));
}


/* probability that the observed chi-square
   exceeds chi2 even if model is correct */
double chi2prob (int deg, double chi2)
{
	return IncompleteGammaQ (0.5*deg, 0.5*chi2);
}



/* chi square test
   ef      expected frequencies (sum up to 1 !!)
   of      observed frequencies (sum up to the number of samples)
   numcat  number of categories
   returns critical significance level */
double chi2test(double *ef, int *of, int numcat, int *chi2fail)
{	
	double chi2, criticals, efn;
	int i, below1, below5, reducedcat;
	int samples;
	
	*chi2fail = FALSE;
	reducedcat = numcat;
	below1 = 0;
	below5 = 0;
	
	/* compute number of samples */
	samples = 0;
	for (i = 0; i < numcat; i++)
		samples = samples + of[i];
	
	/* compute chi square */
	chi2 = 0;
	for (i = 0; i < numcat; i++) {
		efn = ef[i]*((double) samples);
		if (efn < 1.0) below1++;
		if (efn < 5.0) below5++;
		if (efn == 0.0) {
			reducedcat--;
			fprintf(stdout, "FPE error: samples=%d, ef[%d]=%f, of[%d]=%d, efn=%f, nc=%d, rc=%d\n",
					samples, i, ef[i], i, of[i], efn, numcat, reducedcat);
			fprintf(stdout, "PLEASE REPORT THIS ERROR TO DEVELOPERS !!!\n");
			fflush(stdout);
		} else chi2 = chi2 + ((double) of[i]-efn)*((double) of[i]-efn)/efn;
	}

	/* compute significance */
	criticals = chi2prob (numcat-1, chi2);
	
	/* no expected frequency category (sum up to # samples) below 1.0 */
	if (below1 > 0) *chi2fail = TRUE;
	/* no more than 1/5 of the frequency categories below 5.0 */
	if (below5 > (int) floor(samples/5.0)) *chi2fail = TRUE;

	return criticals;
}


/* chi square test
   ef      expected frequencies (sum up to 1 !!)
   of      observed frequencies (sum up to the number of samples)
   numcat  number of categories
   returns critical significance level */
double altchi2test(double *ef, int *of, int numcat, int *chi2fail)
{	
	double chi2, criticals, efn;
	int i, below1, below5;
	int samples;
	
	*chi2fail = FALSE;
	below1 = 0;
	below5 = 0;
	
	/* compute number of samples */
	samples = 0;
	for (i = 0; i < numcat; i++)
		samples = samples + of[i];
	
	/* compute chi square */
	chi2 = 0;
	for (i = 0; i < numcat; i++) {
		efn = ef[i]*((double) samples);
		if (efn < 1.0) below1++;
		if (efn < 5.0) below5++;
		chi2 = chi2 + ((double) of[i]-efn)*((double) of[i]-efn)/efn;
	}

	/* compute significance */
	criticals = chi2prob (numcat-1, chi2);
	
	/* no expected frequency category (sum up to # samples) below 1.0 */
	if (below1 > 0) *chi2fail = TRUE;
	/* no more than 1/5 of the frequency categories below 5.0 */
	if (below5 > (int) floor(samples/5.0)) *chi2fail = TRUE;

	return criticals;
}
