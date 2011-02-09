/*
 * puzzle2.c
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


/* Modified by Christian Zmasek to:
   - Allow 8000 seqs (for pairwise dist. calc.).
   - Names of 26 chars.
   

   !WARNING: Use ONLY together with FORESTER/RIO!
   !For all other puposes download the excellent original!
   
   last modification: 05/19/01


   void getsizesites(FILE *ifp):

   257 -> 8000



   void readid(FILE *infp, int t):

   for (i = 0; i < 10; i++) {                   -> for (i = 0; i < 26; i++) {

   for (i = 9; i > -1; i--) {                   -> for (i = 25; i > -1; i--) {

   for (j = 0; (j < 10) && (flag == TRUE); j++) -> for (j = 0; (j < 26) && (flag == TRUE); j++)



   void initid(int t):

   Identif = new_cmatrix(t, 10);                -> Identif = new_cmatrix(t, 26);

   for (j = 0; j < 10; j++)                     -> for (j = 0; j < 26; j++)



   fputid10(FILE *ofp, int t):
 
   for (i = 0; i < 10; i++)                     -> for (i = 0; i < 26; i++)



   int fputid(FILE *ofp, int t):

   while (Identif[t][i] != ' ' && i < 10) {     -> while (Identif[t][i] != ' ' && i < 26) {




*/

#define EXTERN extern

#include "puzzle.h"
#include <string.h>

#if PARALLEL
#	include "sched.h"
#endif /* PARALLEL */


/******************************************************************************/
/* sequences                                                                  */
/******************************************************************************/

/* read ten characters of current line as identifier */
void readid(FILE *infp, int t)
{
	int i, j, flag, ci;

	for (i = 0; i < 26; i++) { /* CZ 05/19/01 */
		ci = fgetc(infp);
		if (ci == EOF || !isprint(ci)) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (no name for sequence %d)\n\n\n", t+1);
			exit(1);
		}
		Identif[t][i] = (char) ci;
	}	
	/* convert leading blanks in taxon name to underscores */
	flag = FALSE;
	for (i = 25; i > -1; i--) { /* CZ 05/19/01 */
		if (flag == FALSE) {
			if (Identif[t][i] != ' ') flag = TRUE; 
		} else {
			if (Identif[t][i] == ' ') Identif[t][i] = '_';
		}
	}
	/* check whether this name is already used */
	for (i = 0; i < t; i++) { /* compare with all other taxa */
		flag = TRUE; /* assume identity */
		for (j = 0; (j < 26) && (flag == TRUE); j++) /* CZ 05/19/01 */
			if (Identif[t][j] != Identif[i][j])
				flag = FALSE;
		if (flag) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (multiple occurence of sequence name '");
			fputid(STDOUT, t);
			FPRINTF(STDOUTFILE "')\n\n\n");
			exit(1);
		}
	}
}

/* read next allowed character */
char readnextcharacter(FILE *ifp, int notu, int nsite)
{
	char c;

	/* ignore blanks and control characters except newline */
	do {
		if (fscanf(ifp, "%c", &c) != 1) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (missing character at position %d in sequence '", nsite + 1);
			fputid(STDOUT, notu);
			FPRINTF(STDOUTFILE "')\n\n\n");
			exit(1);
		}
	} while (c == ' ' || (iscntrl((int) c) && c != '\n'));
	return c;
}

/* skip rest of the line */
void skiprestofline(FILE* ifp, int notu, int nsite)
{
	int ci;

	/* read chars until the first newline */
	do{
		ci = fgetc(ifp);
		if (ci == EOF) {
			FPRINTF(STDOUTFILE "Unable to proceed (missing newline at position %d in sequence '", nsite + 1);
			fputid(STDOUT, notu);
			FPRINTF(STDOUTFILE "')\n\n\n");
			exit(1);
		}
	} while ((char) ci != '\n');
}

/* skip control characters and blanks */
void skipcntrl(FILE *ifp, int notu, int nsite)
{
	int ci;

	/* read over all control characters and blanks */
	do {
		ci = fgetc(ifp);
		if (ci == EOF) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (missing character at position %d in sequence '", nsite + 1);
			fputid(STDOUT, notu);
			FPRINTF(STDOUTFILE "')\n\n\n");
			exit(1);
		}
	} while (iscntrl(ci) || (char) ci == ' ');
	/* go one character back */
	if (ungetc(ci, ifp) == EOF) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (positioning error at position %d in sequence '", nsite + 1);
		fputid(STDOUT, notu);
		FPRINTF(STDOUTFILE "')\n\n\n");
		exit(1);
	}
}

/* read sequences of one data set */
void getseqs(FILE *ifp)
{
	int notu, nsite, endofline, linelength, i;
	char c;
	
	seqchars = new_cmatrix(Maxspc, Maxseqc);
	/* read all characters */
	nsite = 0; /* next site to be read */
	while (nsite < Maxseqc) {
		/* read first taxon */
		notu = 0;
		/* go to next true line */
		skiprestofline(ifp, notu, nsite);
		skipcntrl(ifp, notu, nsite);
		if (nsite == 0) readid(ifp, notu);
		endofline = FALSE;
		linelength = 0;		
		do {
			c = readnextcharacter(ifp, notu, nsite + linelength);
			if (c == '\n') endofline = TRUE;
			else if (c == '.') {
				FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (invalid character '.' at position ");
				FPRINTF(STDOUTFILE "%d in first sequence)\n\n\n", nsite + linelength + 1);
				exit(1);
			} else if (nsite + linelength < Maxseqc) {
				/* change to upper case */
				seqchars[notu][nsite + linelength] = (char) toupper((int) c);
				linelength++;
			} else {
				endofline = TRUE;
				skiprestofline(ifp, notu, nsite + linelength);
			}
		} while (!endofline);	
		if (linelength == 0) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (line with length 0 at position %d in sequence '", nsite + 1);
			fputid(STDOUT, notu);
			FPRINTF(STDOUTFILE "')\n\n\n");
			exit(1);
		}
		/* read other taxa */
		for (notu = 1; notu < Maxspc; notu++) {
			/* go to next true line */
			if (notu != 1) skiprestofline(ifp, notu, nsite);
			skipcntrl(ifp, notu, nsite);
			if (nsite == 0) readid(ifp, notu);
			for (i = nsite; i < nsite + linelength; i++) {
				c = readnextcharacter(ifp, notu, i);
				if (c == '\n') { /* too short */
					FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (line to short at position %d in sequence '", i + 1);
					fputid(STDOUT, notu);
					FPRINTF(STDOUTFILE "')\n\n\n");
					exit(1);
				} else if (c == '.') {
					seqchars[notu][i] = seqchars[0][i];
				} else {
					/* change to upper case */
					seqchars[notu][i] = (char) toupper((int) c);
				}
			}
		}
		nsite = nsite + linelength;
	}
}

/* initialize identifer array */
void initid(int t)
{
	int i, j;

	Identif = new_cmatrix(t, 26); /* CZ 05/19/01 */
	for (i = 0; i < t; i++)
		for (j = 0; j < 26; j++) /* CZ 05/19/01 */
			Identif[i][j] = ' ';
}

/* print identifier of specified taxon in full 10 char length */
void fputid10(FILE *ofp, int t)
{	
	int i;

	for (i = 0; i < 26; i++) fputc(Identif[t][i], ofp); /* CZ 05/19/01 */
}

/* print identifier of specified taxon up to first space */
int fputid(FILE *ofp, int t)
{	
	int i;
	
	i = 0;
	while (Identif[t][i] != ' ' && i < 26) { /* CZ 05/19/01 */
		fputc(Identif[t][i], ofp);
		i++;
	}
	return i;
}

/* read first line of sequence data set */
void getsizesites(FILE *ifp)
{
	if (fscanf(ifp, "%d", &Maxspc) != 1) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (missing number of sequences)\n\n\n");
			exit(1);
	}
	if (fscanf(ifp, "%d", &Maxseqc) != 1) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (missing number of sites)\n\n\n");
			exit(1);
	}
	
	if (Maxspc < 4) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (less than 4 sequences)\n\n\n");
		exit(1);
	}
	if (Maxspc > 8000) { /* CZ 05/19/01 */
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (more than 8000 sequences)\n\n\n");
		exit(1);
	}
	if (Maxseqc < 1) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (no sequence sites)\n\n\n");
		exit(1);
	}
	Maxbrnch = 2*Maxspc - 3;
}

/* read one data set - PHYLIP interleaved */
void getdataset(FILE *ifp)
{
	initid(Maxspc);
	getseqs(ifp);
}

/* guess data type */
int guessdatatype()
{
	uli numnucs, numchars, numbins;
	int notu, nsite;
	char c;
	
	/* count A, C, G, T, U, N */
	numnucs = 0;
	numchars = 0;
	numbins = 0;
	for (notu = 0; notu < Maxspc; notu++)
		for (nsite = 0; nsite < Maxseqc; nsite++) {
			c = seqchars[notu][nsite];
			if (c == 'A' || c == 'C' || c == 'G' ||
			    c == 'T' || c == 'U' || c == 'N') numnucs++;
			if (c != '-' && c != '?') numchars++;
			if (c == '0' || c == '1') numbins++;
		}
	if (numchars == 0) numchars = 1;
	/* more than 85 % frequency means nucleotide data */
	if ((double) numnucs / (double) numchars > 0.85) return 0;
	else if ((double) numbins / (double) numchars > 0.2) return 2; 
	else return 1;
}

/* translate characters into format used by ML engine */
void translatedataset()
{	
	int notu, sn, co;
	char c;
	cvector code;
	

	/* determine Maxsite - number of ML sites per taxon */
	if (data_optn == 0 && SH_optn) {
		if (SHcodon)
			Maxsite = Maxseqc / 3;
		else
			Maxsite = Maxseqc / 2; /* assume doublets */
		
	} else
		Maxsite = Maxseqc;
	if (data_optn == 0 && (Maxsite % 3) == 0  && !SH_optn) {	
		if (codon_optn == 1 || codon_optn == 2 || codon_optn == 3)
			Maxsite = Maxsite / 3; /* only one of the three codon positions */
		if (codon_optn == 4)
			Maxsite = 2*(Maxsite / 3); /* 1st + 2nd codon positions */
	}
	
	/* reserve memory */
	if (Seqchar != NULL) free_cmatrix(Seqchar);
	Seqchar = new_cmatrix(Maxspc, Maxsite);

	/* code length */
	if (data_optn == 0 && SH_optn)
		code = new_cvector(2);
	else
		code = new_cvector(1);
	
	/* decode characters */
	if (data_optn == 0 && SH_optn) { /* SH doublets */
		
		for (notu = 0; notu < Maxspc; notu++) {
			for (sn = 0; sn < Maxsite; sn++) {
					for (co = 0; co < 2; co++) {
						if (SHcodon)
							c = seqchars[notu][sn*3 + co];
						else
							c = seqchars[notu][sn*2 + co];
						code[co] = c;
					}
				Seqchar[notu][sn] = code2int(code);
			}
		}
		
	} else if (!(data_optn == 0 && (Maxseqc % 3) == 0)) { /* use all */

		for (notu = 0; notu < Maxspc; notu++) {
			for (sn = 0; sn < Maxsite; sn++) {
				code[0] = seqchars[notu][sn];
				Seqchar[notu][sn] = code2int(code);
			}
		}

	} else { /* codons */
		
		for (notu = 0; notu < Maxspc; notu++) {
			for (sn = 0; sn < Maxsite; sn++) {
				if (codon_optn == 1 || codon_optn == 2 || codon_optn == 3)
					code[0] = seqchars[notu][sn*3+codon_optn-1];
				else if (codon_optn == 4) {
					if ((sn % 2) == 0)
						code[0] = seqchars[notu][(sn/2)*3];
					else
						code[0] = seqchars[notu][((sn-1)/2)*3+1];
				} else
					code[0] = seqchars[notu][sn];
				Seqchar[notu][sn] = code2int(code);
			}
		}
	
	}
	free_cvector(code);
}

/* estimate mean base frequencies from translated data set */
void estimatebasefreqs()
{
	int tpmradix, i, j;
	uli all, *gene;
	
	tpmradix = gettpmradix();
	
	if (Freqtpm != NULL) free_dvector(Freqtpm);
	Freqtpm = new_dvector(tpmradix);
	
	if (Basecomp != NULL) free_imatrix(Basecomp);
	Basecomp = new_imatrix(Maxspc, tpmradix);
	
	gene = (uli *) malloc((unsigned) ((tpmradix + 1) * sizeof(uli)));
	if (gene == NULL) maerror("gene in estimatebasefreqs");
	
	for (i = 0; i < tpmradix + 1; i++) gene[i] = 0;
	for (i = 0; i < Maxspc; i++)
		for (j = 0; j < tpmradix; j++) Basecomp[i][j] = 0;
	for (i = 0; i < Maxspc; i++)
		for (j = 0; j < Maxsite; j++) {
			gene[(int) Seqchar[i][j]]++;
			if (Seqchar[i][j] != tpmradix) Basecomp[i][(int) Seqchar[i][j]]++;
			}

	all = Maxspc * Maxsite - gene[tpmradix];
	if (all != 0) { /* normal case */
		for (i = 0; i < tpmradix; i++)
			Freqtpm[i] = (double) gene[i] / (double) all;
	} else { /* pathological case with no unique character in data set */
		for (i = 0; i < tpmradix; i++)
			Freqtpm[i] = 1.0 / (double) tpmradix;
	}
	
	free(gene);
	
	Frequ_optn = TRUE;
}

/* guess model of substitution */
void guessmodel()
{
	double c1, c2, c3, c4, c5, c6;
	dvector f;
	dmatrix a;
	int i;

	Dayhf_optn = FALSE;
	Jtt_optn = TRUE;
	mtrev_optn = FALSE;
	cprev_optn = FALSE;
	blosum62_optn = FALSE;
	vtmv_optn = FALSE;
	wag_optn = FALSE;
	TSparam = 2.0;
	YRparam = 1.0;
	optim_optn = TRUE;
	HKY_optn = TRUE;
	TN_optn = FALSE;
	
	if (data_optn == 1) { /* amino acids */
		
		/* chi2 fit to amino acid frequencies */
		
		f = new_dvector(20);
		a = new_dmatrix(20,20);
		/* chi2 distance Dayhoff */
		dyhfdata(a, f);
		c1 = 0;
		for (i = 0; i < 20; i++)
			c1 = c1 + (Freqtpm[i]-f[i])*(Freqtpm[i]-f[i]);
		/* chi2 distance JTT */
		jttdata(a, f);
		c2 = 0;
		for (i = 0; i < 20; i++)
			c2 = c2 + (Freqtpm[i]-f[i])*(Freqtpm[i]-f[i]);
		/* chi2 distance mtREV */
		mtrevdata(a, f);
		c3 = 0;
		for (i = 0; i < 20; i++)
			c3 = c3 + (Freqtpm[i]-f[i])*(Freqtpm[i]-f[i]);
		/* chi2 distance VT */
		vtmvdata(a, f);
		c4 = 0;
		for (i = 0; i < 20; i++)
			c4 = c4 + (Freqtpm[i]-f[i])*(Freqtpm[i]-f[i]);
		/* chi2 distance WAG */
		wagdata(a, f);
		c5 = 0;
		for (i = 0; i < 20; i++)
			c5 = c5 + (Freqtpm[i]-f[i])*(Freqtpm[i]-f[i]);
		/* chi2 distance cpREV */
		cprev45data(a, f);
		c6 = 0;
		for (i = 0; i < 20; i++)
			c6 = c6 + (Freqtpm[i]-f[i])*(Freqtpm[i]-f[i]);

		free_dvector(f);
		free_dmatrix(a);

#ifndef CPREV
		if ((c1 < c2) && (c1 < c3) && (c1 < c4) && (c1 < c5)) {
		        /* c1 -> Dayhoff */
			Dayhf_optn = TRUE;
			Jtt_optn = FALSE;
			mtrev_optn = FALSE;
			cprev_optn = FALSE;
			vtmv_optn = FALSE;
			wag_optn = FALSE;
			FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
		} else {
			if ((c2 < c3) && (c2 < c4) && (c2 < c5)) {
				/* c2 -> JTT */
				Dayhf_optn = FALSE;
				Jtt_optn = TRUE;
				mtrev_optn = FALSE;
				cprev_optn = FALSE;
				vtmv_optn = FALSE;
				wag_optn = FALSE;
				FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
			} else {
				if ((c3 < c4) && (c3 < c5)) {
					/* c3 -> mtREV */
					Dayhf_optn = FALSE;
					Jtt_optn = FALSE;
					mtrev_optn = TRUE;
					cprev_optn = FALSE;
					vtmv_optn = FALSE;
					wag_optn = FALSE;
					FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on mtDNA)\n");
				} else {
					if ((c4 < c5)) {
						/* c4 -> VT */
						Dayhf_optn = FALSE;
						Jtt_optn = FALSE;
						mtrev_optn = FALSE;
						cprev_optn = FALSE;
						vtmv_optn = TRUE;
						wag_optn = FALSE;
						FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
					} else {
							/* c5 -> WAG */
							Dayhf_optn = FALSE;
							Jtt_optn = FALSE;
							mtrev_optn = FALSE;
							cprev_optn = FALSE;
							vtmv_optn = FALSE;
							wag_optn = TRUE;
							FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
					} /* if c4 else c5 */
				} /* if c3 else c4 */
			} /* if c2 */
		} /* if c1 */

#else /* CPREV */

		if ((c1 < c2) && (c1 < c3) && (c1 < c4) && (c1 < c5) && (c1 < c6)) {
		        /* c1 -> Dayhoff */
			Dayhf_optn = TRUE;
			Jtt_optn = FALSE;
			mtrev_optn = FALSE;
			cprev_optn = FALSE;
			vtmv_optn = FALSE;
			wag_optn = FALSE;
			FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
		} else {
			if ((c2 < c3) && (c2 < c4) && (c2 < c5) && (c2 < c6)) {
				/* c2 -> JTT */
				Dayhf_optn = FALSE;
				Jtt_optn = TRUE;
				mtrev_optn = FALSE;
				cprev_optn = FALSE;
				vtmv_optn = FALSE;
				wag_optn = FALSE;
				FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
			} else {
				if ((c3 < c4) && (c3 < c5) && (c3 < c6)) {
					/* c3 -> mtREV */
					Dayhf_optn = FALSE;
					Jtt_optn = FALSE;
					mtrev_optn = TRUE;
					cprev_optn = FALSE;
					vtmv_optn = FALSE;
					wag_optn = FALSE;
					FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on mtDNA)\n");
				} else {
					if ((c4 < c5) && (c4 < c6)) {
						/* c4 -> VT */
						Dayhf_optn = FALSE;
						Jtt_optn = FALSE;
						mtrev_optn = FALSE;
						cprev_optn = FALSE;
						vtmv_optn = TRUE;
						wag_optn = FALSE;
						FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
					} else {
						if (c5 < c6) {
							/* c5 -> WAG */
							Dayhf_optn = FALSE;
							Jtt_optn = FALSE;
							mtrev_optn = FALSE;
							cprev_optn = FALSE;
							vtmv_optn = FALSE;
							wag_optn = TRUE;
							FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on nuclear DNA)\n");
						} else {
							/* if (c6) */
							/* c6 -> cpREV */
							Dayhf_optn = FALSE;
							Jtt_optn = FALSE;
							mtrev_optn = FALSE;
							cprev_optn = TRUE;
							vtmv_optn = FALSE;
							wag_optn = FALSE;
							FPRINTF(STDOUTFILE "(consists very likely of amino acids encoded on cpDNA)\n");
						} /* if c5 else c6 */
					} /* if c4 else c5 */
				} /* if c3 else c4 */
			} /* if c2 */
		} /* if c1 */
#endif /* CPREV */

	} else if (data_optn == 0) {
		FPRINTF(STDOUTFILE "(consists very likely of nucleotides)\n");
	} else {
		FPRINTF(STDOUTFILE "(consists very likely of binary state data)\n");
	}
} /* guessmodel */


/******************************************************************************/
/* functions for representing and building puzzling step trees                */
/******************************************************************************/

/* initialize tree with the following starting configuration

                 2
          0  +------- C(=2)
  A(=0) -----+
             +------- B(=1)
                 1
                                */
void inittree()
{
	int i;

	/* allocate the memory for the whole tree */
	
	/* allocate memory for vector with all the edges of the tree */
	edge = (ONEEDGE *) calloc(Maxbrnch, sizeof(ONEEDGE) );
	if (edge == NULL) maerror("edge in inittree");

	/* allocate memory for vector with edge numbers of leaves */ 
	edgeofleaf = (int *) calloc(Maxspc, sizeof(int) );
	if (edgeofleaf == NULL) maerror("edgeofleaf in inittree");
	
	/* allocate memory for all the edges the edge map */
	for (i = 0; i < Maxbrnch; i++) {
		edge[i].edgemap = (int *) calloc(Maxbrnch, sizeof(int) );
		if (edge[i].edgemap == NULL) maerror("edgemap in inittree");
	}
		
	/* number all edges */
	for (i = 0; i < Maxbrnch; i++) edge[i].numedge = i;
	
	/* initialize tree */
	
	nextedge = 3;
	nextleaf = 3;
	
	/* edge maps */
	(edge[0].edgemap)[0] = 0; /* you are on the right edge */
	(edge[0].edgemap)[1] = 4; /* go down left for leaf 1 */
	(edge[0].edgemap)[2] = 5; /* go down right for leaf 2 */
	(edge[1].edgemap)[0] = 1; /* go up for leaf 0 */
	(edge[1].edgemap)[1] = 0; /* you are on the right edge */
	(edge[1].edgemap)[2] = 3; /* go up/down right for leaf 2 */
	(edge[2].edgemap)[0] = 1; /* go up for leaf 0 */
	(edge[2].edgemap)[1] = 2; /* go up/down left for leaf 1 */
	(edge[2].edgemap)[2] = 0; /* you are on the right edge */
	
	/* interconnection */
	edge[0].up = NULL;
	edge[0].downleft = &edge[1];
	edge[0].downright = &edge[2];
	edge[1].up = &edge[0];
	edge[1].downleft = NULL;
	edge[1].downright = NULL;
	edge[2].up = &edge[0];
	edge[2].downleft = NULL;
	edge[2].downright = NULL;
	
	/* edges of leaves */
	edgeofleaf[0] = 0;
	edgeofleaf[1] = 1;
	edgeofleaf[2] = 2;
} /* inittree */

/* add next leaf on the specified edge */
void addnextleaf(int dockedge)
{	
	int i;

	if (dockedge >= nextedge) {
		/* Trying to add leaf nextleaf to nonexisting edge dockedge */
		FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR F TO DEVELOPERS\n\n\n");
		exit(1);
	}
	
	if (nextleaf >= Maxspc) {
		/* Trying to add leaf nextleaf to a tree with Maxspc leaves */
		FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR G TO DEVELOPERS\n\n\n");
		exit(1);
	}
		
	/* necessary change in edgeofleaf if dockedge == edgeofleaf[0] */
	if (edgeofleaf[0] == dockedge) edgeofleaf[0] = nextedge;
	
	/* adding nextedge to the tree */
	edge[nextedge].up = edge[dockedge].up;
	edge[nextedge].downleft = &edge[dockedge];
	edge[nextedge].downright = &edge[nextedge+1];
	edge[dockedge].up = &edge[nextedge];

	if (edge[nextedge].up != NULL) {
		if ( ((edge[nextedge].up)->downleft) == &edge[dockedge] ) 
			(edge[nextedge].up)->downleft  = &edge[nextedge];
		else   
			(edge[nextedge].up)->downright = &edge[nextedge];
	}

	/* adding nextedge + 1 to the tree */
	edge[nextedge+1].up = &edge[nextedge];
	edge[nextedge+1].downleft = NULL;
	edge[nextedge+1].downright = NULL;
	edgeofleaf[nextleaf] = nextedge+1;
	
	/* the two new edges get info about the old edges */
	/* nextedge */
	for (i = 0; i < nextedge; i++) {
		switch ( (edge[dockedge].edgemap)[i] ) {
			
			/* down right changes to down left */
			case 5:		(edge[nextedge].edgemap)[i] = 4;
						break;
						
			/* null changes to down left */
			case 0:		(edge[nextedge].edgemap)[i] = 4;
						break;
			
			default: 	(edge[nextedge].edgemap)[i] =
							(edge[dockedge].edgemap)[i];
						break;
		}
	}

	/* nextedge + 1 */
	for (i = 0; i < nextedge; i++) {
		switch ( (edge[dockedge].edgemap)[i] ) {
			
			/* up/down left changes to up */
			case 2:		(edge[nextedge+1].edgemap)[i] = 1;
						break;
			
			/* up/down right changes to up */		
			case 3:		(edge[nextedge+1].edgemap)[i] = 1;
						break;
						
			/* down left changes to up/down left */			
			case 4:		(edge[nextedge+1].edgemap)[i] = 2;
						break;

			/* down right changes to up/down left */	
			case 5:		(edge[nextedge+1].edgemap)[i] = 2;
						break;
			
			/* null changes to up/down left */
			case 0:		(edge[nextedge+1].edgemap)[i] = 2;
						break;
			
			/* up stays up */
			default: 	(edge[nextedge+1].edgemap)[i] =
							(edge[dockedge].edgemap)[i];
						break;
		}
	}

	/* dockedge */
	for (i = 0; i < nextedge; i++) {
		switch ( (edge[dockedge].edgemap)[i] ) {
			
			/* up/down right changes to up */
			case 3:		(edge[dockedge].edgemap)[i] = 1;
						break;
						
			/* up/down left changes to up */
			case 2:		(edge[dockedge].edgemap)[i] = 1;
						break;
			
			default: 	break;
		}
	}

	/* all edgemaps are updated for the two new edges */
	/* nextedge */
	(edge[nextedge].edgemap)[nextedge] = 0;
	(edge[nextedge].edgemap)[nextedge+1] = 5; /* down right */
	
	/* nextedge + 1 */
	(edge[nextedge+1].edgemap)[nextedge] = 1; /* up */
	(edge[nextedge+1].edgemap)[nextedge+1] = 0;
	
	/* all other edges */
	for (i = 0; i < nextedge; i++) {
		(edge[i].edgemap)[nextedge] = (edge[i].edgemap)[dockedge];
		(edge[i].edgemap)[nextedge+1] = (edge[i].edgemap)[dockedge];
	}
	
	/* an extra for dockedge */
	(edge[dockedge].edgemap)[nextedge] = 1; /* up */
	(edge[dockedge].edgemap)[nextedge+1] = 3; /* up/down right */

	nextleaf++;
	nextedge = nextedge + 2;
} /* addnextleaf */


/* free memory (to be called after inittree) */
void freetree()
{
	int i;

	for (i = 0; i < 2 * Maxspc - 3; i++) free(edge[i].edgemap);
	free(edge);
	free(edgeofleaf);	
} /* freetree */

/* writes OTU sitting on edge ed */
void writeOTU(FILE *outfp, int ed)
{	
	int i;

	/* test whether we are on a leaf */
	if (edge[ed].downright == NULL && edge[ed].downleft == NULL) {
		for (i = 1; i < nextleaf; i++) {
			if (edgeofleaf[i] == ed) { /* i is the leaf of ed */
				column += fputid(outfp, trueID[i]);
				return;
			}
		}
	}

	/* we are NOT on a leaf */
	fprintf(outfp, "(");
	column++;
	writeOTU(outfp, edge[ed].downleft->numedge);
	fprintf(outfp, ",");
	column++;
	column++;
	if (column > 55) {
		column = 2;
		fprintf(outfp, "\n  ");
	}
	writeOTU(outfp, edge[ed].downright->numedge);	
	fprintf(outfp, ")");
	column++;
} /* writeOTU */

/* write tree */
void writetree(FILE *outfp)
{	
	column = 1;
	fprintf(outfp, "(");
	column += fputid(outfp, trueID[0]) + 3;
	fprintf(outfp, ",");
	writeOTU(outfp, edge[edgeofleaf[0]].downleft->numedge);
	column++;
	column++;
	fprintf(outfp, ",");
	writeOTU(outfp, edge[edgeofleaf[0]].downright->numedge);	
	fprintf(outfp, ");\n");	
} /* writetree */


/* clear all edgeinfos */
void resetedgeinfo()
{
	int i;
	
	for (i = 0; i < nextedge; i++)
		edge[i].edgeinfo = 0;
} /* resetedgeinfo */

/* increment all edgeinfo between leaf A and B */
void incrementedgeinfo(int A, int B)
{	
	int curredge, finaledge, nextstep;
	
	if (A == B) return;
	
	finaledge = edgeofleaf[B];
	
	curredge = edgeofleaf[A];
	edge[curredge].edgeinfo = edge[curredge].edgeinfo + 1;
	
	while (curredge != finaledge) {
		nextstep = (edge[curredge].edgemap)[finaledge];
		switch (nextstep) {

			/* up */
			case 1:	curredge = (edge[curredge].up)->numedge;
					break;
				
			/* up/down left */
			case 2:	curredge = ((edge[curredge].up)->downleft)->numedge;
					break;

			/* up/down right */
			case 3:	curredge = ((edge[curredge].up)->downright)->numedge;
					break;
				
			/* down left */
			case 4:	curredge = (edge[curredge].downleft)->numedge;
					break;
				
			/* down right */
			case 5:	curredge = (edge[curredge].downright)->numedge;
					break;
				
		}
		edge[curredge].edgeinfo = edge[curredge].edgeinfo + 1;
	}	
} /* incrementedgeinfo */

/* checks which edge has the lowest edgeinfo
   if there are several edges with the same lowest edgeinfo,
   one of them will be selected randomly */
void minimumedgeinfo()
{
	int i, k, howmany, randomnum;

	howmany = 1;
	minedge = 0;
	mininfo = edge[0].edgeinfo;
	for (i = 1; i < nextedge; i++)
		if (edge[i].edgeinfo <= mininfo) {
			if (edge[i].edgeinfo == mininfo) {
				howmany++;
			} else {
				minedge = i;
				mininfo = edge[i].edgeinfo;
				howmany = 1;
			}
		}
	
	if (howmany > 1) { /* draw random edge */
		randomnum = randominteger(howmany) + 1; /* 1 to howmany */
		i = -1;
		for (k = 0; k < randomnum; k++) {
			do {
				i++;
			} while (edge[i].edgeinfo != mininfo);
			minedge = i;
		}
	}
} /* minimumedgeinfo */




/*******************************************/
/* tree sorting                            */
/*******************************************/

/* compute address of the 4 int (sort key) in the 4 int node */
int ct_sortkeyaddr(int addr)
{
  int a, res;
  a = addr % 4;
  res = addr - a + 3;
  return res;
}


/**********/

/* compute address of the next edge pointer in a 4 int node (0->1->2->0) */
int ct_nextedgeaddr(int addr)
{
  int a, res;
  a = addr % 4;
  if ( a == 2 ) { res = addr - 2; }
  else          { res = addr + 1; }
  return res;
}


/**********/

/* compute address of 1st edge of a 4 int node from node number */
int ct_1stedge(int node)
{
  int res;
  res = 4 * node;
  return res;
}


/**********/

/* compute address of 2nd edge of a 4 int node from node number */
int ct_2ndedge(int node)
{
  int res;
  res = 4 * node +1;
  return res;
}


/**********/

/* compute address of 3rd edge of a 4 int node from node number */
int ct_3rdedge(int node)
{
  int res;
  res = 4 * node +2;
  return res;
}


/**********/

/* check whether node 'node' is a leaf (2nd/3rd edge pointer = -1) */
int ct_isleaf(int node, int *ctree)
{
  return (ctree[ct_3rdedge(node)] < 0);
}


/**********/

/* compute node number of 4 int node from an edge addr. */
int ct_addr2node(int addr)
{
  int a, res;
  a = addr % 4;
  res = (int) ((addr - a) / 4);
  return res;
}


/**********/

/* print graph pointers for checking */
void printctree(int *ctree)
{
	int n;
	for (n=0; n < 2*Maxspc; n++) {
		printf("n[%3d] = (%3d.%2d, %3d.%2d, %3d.%2d | %3d)\n", n,
		(int) ctree[ct_1stedge(n)]/4,
		(int) ctree[ct_1stedge(n)]%4,
		(int) ctree[ct_2ndedge(n)]/4,
		(int) ctree[ct_2ndedge(n)]%4,
		(int) ctree[ct_3rdedge(n)]/4,
		(int) ctree[ct_3rdedge(n)]%4,
		ctree[ct_3rdedge(n)+1]);
	}
        printf("\n");
} /* printctree */


/**********/

/* allocate memory for ctree 3 ints pointer plus 1 check byte */
int *initctree()
{
  int *snodes;
  int n;

  snodes = (int *) malloc(4 * 2 * Maxspc * sizeof(int));
  if (snodes == NULL) maerror("snodes in copytree");

  for (n=0; n<(4 * 2 * Maxspc); n++) {
      snodes[n]=-1;
  }
  return snodes;
}


/**********/

/* free memory of a tree for sorting */
void freectree(int **snodes)
{
	free(*snodes);
	*snodes = NULL;
}


/**********/

/* copy subtree recursively */
void copyOTU(int *ctree,                  /* tree array struct            */
             int *ct_nextnode,            /* next free node               */
             int ct_curredge,             /* currende edge to add subtree */
             int *ct_nextleaf,            /* next free leaf (0-maxspc)    */
             int ed)                      /* edge in puzzling step tree   */
{
        int i, nextcurredge;

        /* test whether we are on a leaf */
        if (edge[ed].downright == NULL && edge[ed].downleft == NULL) {
                for (i = 1; i < nextleaf; i++) {
                        if (edgeofleaf[i] == ed) { /* i is the leaf of ed */
				nextcurredge          = ct_1stedge(*ct_nextleaf);
				ctree[ct_curredge]    = nextcurredge;
				ctree[nextcurredge]   = ct_curredge;
                                ctree[ct_sortkeyaddr(nextcurredge)] = trueID[i];
				(*ct_nextleaf)++;
                                return;
                        }
                }
        }

        /* we are NOT on a leaf */
	nextcurredge        = ct_1stedge(*ct_nextnode);
        ctree[ct_curredge]     = nextcurredge;
	ctree[nextcurredge] = ct_curredge;
        (*ct_nextnode)++;
	nextcurredge = ct_nextedgeaddr(nextcurredge);
        copyOTU(ctree, ct_nextnode, nextcurredge, 
                ct_nextleaf, edge[ed].downleft->numedge);

	nextcurredge = ct_nextedgeaddr(nextcurredge);
        copyOTU(ctree, ct_nextnode, nextcurredge, 
                ct_nextleaf, edge[ed].downright->numedge);
}


/**********/

/* copy treestructure to sorting structure */
void copytree(int *ctree)
{
        int ct_curredge;
        int ct_nextleaf;
        int ct_nextnode;

        ct_nextnode = Maxspc;
        ct_curredge = ct_1stedge(ct_nextnode);
        ct_nextleaf = 1;

        ctree[ct_1stedge(0)] = ct_curredge;
        ctree[ct_curredge]   = ct_1stedge(0);
        ctree[ct_sortkeyaddr(0)] = trueID[0];

        ct_nextnode++;
        
        ct_curredge = ct_nextedgeaddr(ct_curredge);
        copyOTU(ctree, &ct_nextnode, ct_curredge, 
                &ct_nextleaf, edge[edgeofleaf[0]].downleft->numedge);

        ct_curredge = ct_nextedgeaddr(ct_curredge);
        copyOTU(ctree, &ct_nextnode, ct_curredge, 
                &ct_nextleaf, edge[edgeofleaf[0]].downright->numedge);
}


/**********/

/* sort subtree from edge recursively by indices */
int sortOTU(int edge, int *ctree)
{
	int key1, key2;
	int edge1, edge2;
	int tempedge;

	if (ctree[ct_2ndedge((int) (edge / 4))] < 0)
		return ctree[ct_sortkeyaddr(edge)];

	edge1 = ctree[ct_nextedgeaddr(edge)];
	edge2 = ctree[ct_nextedgeaddr(ct_nextedgeaddr(edge))];

        /* printf ("visiting [%5d] -> [%5d], [%5d]\n", edge, edge1, edge2); */
        /* printf ("visiting [%2d.%2d] -> [%2d.%2d], [%2d.%2d]\n", 
           (int)(edge/4), edge%4, (int)(edge1/4), edge1%4, 
           (int)(edge2/4), edge2%4); */

	key1  = sortOTU(edge1, ctree); 
	key2  = sortOTU(edge2, ctree); 
	
	if (key2 < key1) {
		tempedge            = ctree[ctree[edge1]];
		ctree[ctree[edge1]] = ctree[ctree[edge2]];
		ctree[ctree[edge2]] = tempedge;
		tempedge            = ctree[edge1];
		ctree[edge1]        = ctree[edge2];
		ctree[edge2]        = tempedge;
	  	ctree[ct_sortkeyaddr(edge)] = key2;
		
	} else {
	  ctree[ct_sortkeyaddr(edge)] = key1;
	}
	return ctree[ct_sortkeyaddr(edge)];
}


/**********/

/* sort ctree recursively by indices */
int sortctree(int *ctree)
{
	int n, startnode=-1;
	for(n=0; n<Maxspc; n++) {
		if (ctree[ct_sortkeyaddr(n*4)] == 0)
			startnode = n;
	}
	sortOTU(ctree[startnode * 4], ctree);
	return startnode;
}


/**********/

/* print recursively subtree of edge of sorted tree ctree */
void printfsortOTU(int edge, int *ctree)
{
        int edge1, edge2;

        if (ctree[ct_2ndedge((int) (edge / 4))] < 0) {
                printf("%d", ctree[ct_sortkeyaddr(edge)]);
                return;
        }

        edge1 = ctree[ct_nextedgeaddr(edge)];
        edge2 = ctree[ct_nextedgeaddr(ct_nextedgeaddr(edge))];

        printf("(");
        printfsortOTU(edge1, ctree); 
        printf(",");
        printfsortOTU(edge2, ctree); 
        printf(")");

}


/**********/

/* print recursively sorted tree ctree */
int printfsortctree(int *ctree)
{
        int n, startnode=-1;
        for(n=0; n<Maxspc; n++) {
                if (ctree[ct_sortkeyaddr(n*4)] == 0)
                        startnode = n;
        }
        printf ("(%d,", ctree[ct_sortkeyaddr(startnode*4)]);
        printfsortOTU(ctree[startnode * 4], ctree);
        printf (");\n");
        return startnode;
}


/**********/

/* print recursively subtree of edge of sorted tree ctree to string */
void sprintfOTU(char *str, int *len, int edge, int *ctree)
{
        int edge1, edge2;

        if (ctree[ct_2ndedge((int) (edge / 4))] < 0) {
                *len+=sprintf(&(str[*len]), "%d", ctree[ct_sortkeyaddr(edge)]);
		return;
	}

        edge1 = ctree[ct_nextedgeaddr(edge)];
        edge2 = ctree[ct_nextedgeaddr(ct_nextedgeaddr(edge))];

	sprintf(&(str[*len]), "(");
	(*len)++;
        sprintfOTU(str, len, edge1, ctree); 
	sprintf(&(str[*len]), ",");
	(*len)++;
        sprintfOTU(str, len, edge2, ctree); 
	sprintf(&(str[*len]), ")");
	(*len)++;
}

/**********/

/* print recursively sorted tree ctree to string */
char *sprintfctree(int *ctree, int strglen)
{
	char *treestr,
	     *tmpptr;
        int n,
	    len=0,
	    startnode=-1;
	treestr = (char *) malloc(strglen * sizeof(char));
	tmpptr  = treestr;
        for(n=0; n<Maxspc; n++) {
                if (ctree[ct_sortkeyaddr(n*4)] == 0)
                        startnode = n;
        }
	len+=sprintf (&(tmpptr[len]), "(%d,", ctree[ct_sortkeyaddr(startnode*4)]);
        sprintfOTU(tmpptr, &len, ctree[startnode * 4], ctree);
	len+=sprintf (&(tmpptr[len]), ");");
        return treestr;
}


/**********/


/***********************************************/
/* establish and handle a list of sorted trees */
/***********************************************/

int itemcount;

/* initialize structure */
treelistitemtype *inittreelist(int *treenum)
{
	*treenum = 0;
	return    NULL;
}


/**********/

/* malloc new tree list item */
treelistitemtype *gettreelistitem()
{
	treelistitemtype *tmpptr;
	tmpptr = (treelistitemtype *)malloc(sizeof(treelistitemtype));
	if (tmpptr == NULL) maerror("item of intermediate tree stuctures");
	(*tmpptr).pred = NULL;
	(*tmpptr).succ = NULL;
	(*tmpptr).tree = NULL;
	(*tmpptr).count = 0;
	(*tmpptr).idx = itemcount++;
	return tmpptr;
}

/**********/

/* free whole tree list */
void freetreelist(treelistitemtype **list,
                  int               *numitems,
                  int               *numsum)
{
	treelistitemtype *current; 
	treelistitemtype *next;
	current = *list;
	while (current != NULL) {
		next = (*current).succ;
		if ((*current).tree != NULL) {
			free ((*current).tree);
			(*current).tree = NULL;
		}
		free(current);
		current = next;
	}
	*list = NULL;
	*numitems = 0;
	*numsum = 0;
} /* freetreelist */


/**********/

/* add tree to the tree list */
treelistitemtype *addtree2list(char             **tree,         /* sorted tree string */
                               int                numtrees,     /* how many occurred, e.g. in parallel */
                               treelistitemtype **list,         /* addr. of tree list */
                               int               *numitems,     
                               int               *numsum)
{
	treelistitemtype *tmpptr = NULL;
	treelistitemtype *newptr = NULL;
	int               result;
	int               done = 0;

	if ((*list == NULL) || (numitems == 0)) {
		newptr = gettreelistitem();
		(*newptr).tree = *tree; 
		*tree = NULL;
		(*newptr).id    = *numitems;
		(*newptr).count = numtrees;
		*numitems = 1;
		*numsum   = numtrees;
		*list = newptr;
	} else {
		tmpptr = *list;
		while(done == 0) {
			result = strcmp( (*tmpptr).tree, *tree);
			if (result==0) {
				free(*tree); *tree = NULL;
				(*tmpptr).count += numtrees;
				*numsum += numtrees;
				done = 1;
				newptr = tmpptr;
			} else { if (result < 0) {
					if ((*tmpptr).succ != NULL)
						tmpptr = (*tmpptr).succ;
					else {
						newptr = gettreelistitem();
						(*newptr).tree = *tree; 
						*tree = NULL;
						(*newptr).id    = *numitems;
						(*newptr).count = numtrees;
						(*newptr).pred  = tmpptr;
						(*tmpptr).succ  = newptr;
						(*numitems)++;
						*numsum += numtrees;
						done = 1;
					}
			} else { /* result < 0 */
				newptr = gettreelistitem();
				(*newptr).tree = *tree; 
				*tree = NULL;
				(*newptr).id    = *numitems;
				(*newptr).count = numtrees;
				(*newptr).succ  = tmpptr;
				(*newptr).pred  = (*tmpptr).pred;
				(*tmpptr).pred  = newptr;
				*numsum += numtrees;

				if ((*newptr).pred != NULL) {
				   (*(*newptr).pred).succ = newptr;
				} else {
				   *list = newptr;
				}
				(*numitems)++;
				done = 1;
			} /* end if result < 0 */
			} /* end if result != 0 */
		} /* while  searching in list */
	} /* if list empty, else */
	return (newptr);
} /* addtree2list */


/**********/

/* resort list of trees by number of occurences for output */
void sortbynum(treelistitemtype *list, treelistitemtype **sortlist)
{
	treelistitemtype *tmpptr = NULL;
	treelistitemtype *curr = NULL;
	treelistitemtype *next = NULL;
	int xchange = 1;

	if (list == NULL) fprintf(stderr, "Grrrrrrrrr>>>>\n");
	tmpptr = list;
	*sortlist = list;
	while (tmpptr != NULL) {
		(*tmpptr).sortnext = (*tmpptr).succ;
		(*tmpptr).sortlast = (*tmpptr).pred;
		tmpptr = (*tmpptr).succ;
	}

	while (xchange > 0) {
		curr = *sortlist;
		xchange = 0;
		if (curr == NULL) fprintf(stderr, "Grrrrrrrrr>>>>\n");
		while((*curr).sortnext != NULL) {
			next = (*curr).sortnext;
			if ((*curr).count >= (*next).count)
				curr = (*curr).sortnext;
			else {
				if ((*curr).sortlast != NULL)
					(*((*curr).sortlast)).sortnext = next;
				if (*sortlist == curr)
					*sortlist = next;
				(*next).sortlast = (*curr).sortlast;

				if ((*next).sortnext != NULL)
					(*((*next).sortnext)).sortlast = curr;
				(*curr).sortnext = (*next).sortnext;

				(*curr).sortlast = next;
				(*next).sortnext = curr;

				xchange++;
			}
		}
	}
}  /* sortbynum */


/**********/

/* print puzzling step tree stuctures for checking */
void printfpstrees(treelistitemtype *list)
{
	char ch;
	treelistitemtype *tmpptr = NULL;
	tmpptr = list;
        ch = '-';
	while (tmpptr != NULL) {
		printf ("%c[%2d]  %5d     %s\n", ch, (*tmpptr).idx, (*tmpptr).count, (*tmpptr).tree);
		tmpptr = (*tmpptr).succ;
		ch = ' ';
	}
}

/**********/

/* print sorted puzzling step tree stucture with names */
void fprintffullpstree(FILE *outf, char *treestr)
{
	int count = 0;
	int idnum = 0;
	int n;
	for(n=0; treestr[n] != '\0'; n++){
		while(isdigit((int)treestr[n])){
			idnum = (10 * idnum) + ((int)treestr[n]-48);
			n++;
			count++;
		}
		if (count > 0){
#			ifdef USEQUOTES
				fprintf(outf, "'");
#			endif
			(void)fputid(outf, idnum);
#			ifdef USEQUOTES
				fprintf(outf, "'");
#			endif
			count = 0;
			idnum = 0;
		}
		fprintf(outf, "%c", treestr[n]);
	}
}


/**********/

/* print sorted puzzling step tree stuctures with names */
void fprintfsortedpstrees(FILE *output, 
                          treelistitemtype *list,  /* tree list */
                          int itemnum,             /* order number */
                          int itemsum,             /* number of trees */
                          int comment,             /* with statistics, or puzzle report ? */
                          float cutoff)            /* cutoff percentage */
{
	treelistitemtype *tmpptr = NULL;
	treelistitemtype *slist = NULL;
	int num = 1;
        float percent;

	if (list == NULL) fprintf(stderr, "Grrrrrrrrr>>>>\n");
	sortbynum(list, &slist); 

	tmpptr = slist;
	while (tmpptr != NULL) {
		percent = (float)(100.0 * (*tmpptr).count / itemsum);
		if ((cutoff == 0.0) || (cutoff <= percent)) {
			if (comment)
				fprintf (output, "[ %d. %d %.2f %d %d %d ]", num++, (*tmpptr).count, percent, (*tmpptr).id, itemnum, itemsum);
			else {
				if (num == 1){
					fprintf (output, "\n");
					fprintf (output, "The following tree(s) occured in more than %.2f%% of the %d puzzling steps.\n", cutoff, itemsum);
					fprintf (output, "The trees are orderd descending by the number of occurences.\n");
					fprintf (output, "\n");
					fprintf (output, "\n       occurences    ID  Phylip tree\n");
				}
				fprintf (output, "%2d. %5d %6.2f%% %5d  ", num++, (*tmpptr).count, percent, (*tmpptr).id);
			}
			fprintffullpstree(output, (*tmpptr).tree);
			fprintf (output, "\n");
		}
		tmpptr = (*tmpptr).sortnext;
	}

	if (!comment) {
		fprintf (output, "\n");
		switch(num) {
			case 1: fprintf (output, "There were no tree topologies (out of %d) occuring with a percentage >= %.2f%% of the %d puzzling steps.\n", itemnum, cutoff, itemsum); break;
			case 2: fprintf (output, "There was one tree topology (out of %d) occuring with a percentage >= %.2f%%.\n", itemnum, cutoff); break;
			default: fprintf (output, "There were %d tree topologies (out of %d) occuring with a percentage >= %.2f%%.\n", num-1, itemnum, cutoff); break;
		}
		fprintf (output, "\n");
		fprintf (output, "\n");
	}
	
}  /* fprintfsortedpstrees */

/**********/

/* print sorted tree topologies for checking */
void printfsortedpstrees(treelistitemtype *list)
{
	treelistitemtype *tmpptr = NULL;
	treelistitemtype *slist = NULL;

	sortbynum(list, &slist); 

	tmpptr = slist;
	while (tmpptr != NULL) {
		printf ("[%2d]  %5d     %s\n", (*tmpptr).idx, (*tmpptr).count, (*tmpptr).tree);
		tmpptr = (*tmpptr).sortnext;
	}
}  /* printfsortedpstrees */


/*******************************************/
/* end of tree sorting                     */
/*******************************************/



/******************************************************************************/
/* functions for computing the consensus tree                                 */
/******************************************************************************/

/* prepare for consensus tree analysis */
void initconsensus()
{
#	if ! PARALLEL
		biparts = new_cmatrix(Maxspc-3, Maxspc);
#	endif /* PARALLEL */

	if (Maxspc % 32 == 0)
		splitlength = Maxspc/32;
	else splitlength = (Maxspc + 32 - (Maxspc % 32))/32;
	numbiparts = 0; /* no pattern stored so far */
	maxbiparts = 0; /* no memory reserved so far */
	splitfreqs = NULL;
	splitpatterns = NULL;
	splitsizes = NULL;
	splitcomp = (uli *) malloc(splitlength * sizeof(uli) );
	if (splitcomp == NULL) maerror("splitcomp in initconsensus");
}

/* prototype needed for recursive function */
void makepart(int i, int curribrnch);

/* recursive function to get bipartitions */
void makepart(int i, int curribrnch)
{	
	int j;
	
	if ( edge[i].downright == NULL ||
		  edge[i].downleft == NULL) { /* if i is leaf */
			
			/* check out what leaf j sits on this edge i */		
			for (j = 1; j < Maxspc; j++) {
				if (edgeofleaf[j] == i) {
					biparts[curribrnch][trueID[j]] = '*';	
					return;
				}	
			}
	} else { /* still on inner branch */
		makepart(edge[i].downleft->numedge, curribrnch);
		makepart(edge[i].downright->numedge, curribrnch);
	}
}

/* compute bipartitions of tree of current puzzling step */
void computebiparts()
{
	int i, j, curribrnch;
	
	curribrnch = -1;
	
	for (i = 0; i < Maxspc - 3; i++)
		for (j = 0; j < Maxspc; j++)
			biparts[i][j] = '.';

	for (i = 0; i < Maxbrnch; i++) {
		if (!(     edgeofleaf[0] == i    ||
		       edge[i].downright == NULL ||
		        edge[i].downleft == NULL) ) { /* check all inner branches */
			curribrnch++;
			makepart(i, curribrnch);
			
			/* make sure that the root is always a '*' */
			if (biparts[curribrnch][outgroup] == '.') {
				for (j = 0; j < Maxspc; j++) {
					if (biparts[curribrnch][j] == '.')
						biparts[curribrnch][j] = '*';
					else
						biparts[curribrnch][j] = '.';
				}
			}
		}
	}
}

/* print out the bipartition n of all different splitpatterns */
void printsplit(FILE *fp, uli n)
{
	int i, j, col;
	uli z;
	
	col = 0;
	for (i = 0; i < splitlength; i++) {
		z = splitpatterns[n*splitlength + i];
		for (j = 0; j < 32 && col < Maxspc; j++) {
			if (col % 10 == 0 && col != 0) fprintf(fp, " ");
			if (z & 1) fprintf(fp, ".");
			else fprintf(fp, "*");
			z = (z >> 1);
			col++;
		}
	}
}

/* make new entries for new different bipartitions and count frequencies */
void makenewsplitentries()
{
	int i, j, bpc, identical, idflag, bpsize;
	uli nextentry, obpc;

	/* where the next entry would be in splitpatterns */
	nextentry = numbiparts;
	
	for (bpc = 0; bpc < Maxspc - 3; bpc++) { /* for every new bipartition */	
		/* convert bipartition into a more compact format */
		bpsize = 0;
		for (i = 0; i < splitlength; i++) {
			splitcomp[i] = 0;	
			for (j = 0; j < 32; j++) {
				splitcomp[i] = splitcomp[i] >> 1;
				if (i*32 + j < Maxspc)
					if (biparts[bpc][i*32 + j] == '.') {
						/* set highest bit */
						splitcomp[i] = (splitcomp[i] | 2147483648UL);
						bpsize++; /* count the '.' */
					}
			}
		}		
		/* compare to the *old* patterns */
		identical = FALSE;
		for (obpc = 0; (obpc < numbiparts) && (!identical); obpc++) {
			/* compare first partition size */
			if (splitsizes[obpc] == bpsize) idflag = TRUE;
			else idflag = FALSE;
			/* if size is identical compare whole partition */
			for (i = 0; (i < splitlength) && idflag; i++)
				if (splitcomp[i] != splitpatterns[obpc*splitlength + i])
					idflag = FALSE;
			if (idflag) identical = TRUE;
		}
		if (identical) { /* if identical increase frequency */
			splitfreqs[2*(obpc-1)]++;
		} else { /* create new entry */
			if (nextentry == maxbiparts) { /* reserve more memory */
				maxbiparts = maxbiparts + 2*Maxspc;
				splitfreqs = (uli *) myrealloc(splitfreqs,
					2*maxbiparts * sizeof(uli) );
				/* 2x: splitfreqs contains also an index (sorting!) */
				if (splitfreqs == NULL) maerror("splitfreqs in makenewsplitentries");
				splitpatterns = (uli *) myrealloc(splitpatterns,
					splitlength*maxbiparts * sizeof(uli) );
				if (splitpatterns == NULL) maerror("splitpatterns in makenewsplitentries");
				splitsizes = (int *) myrealloc(splitsizes,
					maxbiparts * sizeof(int) );
				if (splitsizes == NULL) maerror("splitsizes in makenewsplitentries");
			}
			splitfreqs[2*nextentry] = 1; /* frequency */
			splitfreqs[2*nextentry+1] = nextentry; /* index for sorting */
			for (i = 0; i < splitlength; i++)
				splitpatterns[nextentry*splitlength + i] = splitcomp[i];
			splitsizes[nextentry] = bpsize;
			nextentry++;
		}
	}
	numbiparts = nextentry;
}

/* general remarks:

   - every entry in consbiparts is one node of the consensus tree
   - for each node one has to know which taxa and which other nodes
     are *directly* descending from it
   - for every taxon/node number there is a flag that shows
     whether it descends from the node or not
   - '0' means that neither a taxon nor another node with the
         corresponding number decends from the node
     '1' means that the corresponding taxon descends from the node
     '2' means that the corresponding node descends from the node
     '3' means that the corresponding taxon and node descends from the node
*/

/* copy bipartition n of all different splitpatterns to consbiparts[k] */
void copysplit(uli n, int k)
{
	int i, j, col;
	uli z;
	
	col = 0;
	for (i = 0; i < splitlength; i++) {
		z = splitpatterns[n*splitlength + i];
		for (j = 0; j < 32 && col < Maxspc; j++) {
			if (z & 1) consbiparts[k][col] = '1';
			else consbiparts[k][col] = '0';
			z = (z >> 1);
			col++;
		}
	}
}

/* compute majority rule consensus tree */
void makeconsensus()
{
	int i, j, k, size, subnode;
	char chari, charj;

	/* sort bipartition frequencies */
	qsort(splitfreqs, numbiparts, 2*sizeof(uli), ulicmp);
	/* how many bipartitions are included in the consensus tree */
	consincluded = 0;
	for (i = 0; i < numbiparts && i == consincluded; i++) {
		if (2*splitfreqs[2*i] > Numtrial) consincluded = i + 1;
	}

	/* collect all info about majority rule consensus tree */
	/* the +1 is due to the edge with the root */
	consconfid = new_ivector(consincluded + 1);
	conssizes = new_ivector(2*consincluded + 2);
	consbiparts = new_cmatrix(consincluded + 1, Maxspc);
	
	for (i = 0; i < consincluded; i++) {
		/* copy partition to consbiparts */
		copysplit(splitfreqs[2*i+1], i);
		/* frequency in percent (rounded to integer) */
		consconfid[i] = (int) floor(100.0*splitfreqs[2*i]/Numtrial + 0.5);
		/* size of partition */
		conssizes[2*i] = splitsizes[splitfreqs[2*i+1]];
		conssizes[2*i+1] = i;
	}
	for (i = 0; i < Maxspc; i++) consbiparts[consincluded][i] = '1';
	consbiparts[consincluded][outgroup] = '0';
	consconfid[consincluded] = 100;
	conssizes[2*consincluded] = Maxspc - 1;
	conssizes[2*consincluded + 1] = consincluded;

	/* sort bipartitions according to cluster size */
	qsort(conssizes, consincluded + 1, 2*sizeof(int), intcmp);

	/* reconstruct consensus tree */
	for (i = 0; i < consincluded; i++) { /* try every node */
		size = conssizes[2*i]; /* size of current node */
		for (j = i + 1; j < consincluded + 1; j++) {
		
			/* compare only with nodes with more descendants */
			if (size == conssizes[2*j]) continue;
			
			/* check whether node i is a subnode of j */
			subnode = FALSE;
			for (k = 0; k < Maxspc && !subnode; k++) {
				chari = consbiparts[ conssizes[2*i+1] ][k];
				if (chari != '0') {
					charj = consbiparts[ conssizes[2*j+1] ][k];
					if (chari == charj || charj == '3') subnode = TRUE;
				}
			}
			
			/* if i is a subnode of j change j accordingly */
			if (subnode) {
				/* remove subnode i from j */
				for (k = 0; k < Maxspc; k++) {
					chari = consbiparts[ conssizes[2*i+1] ][k];
					if (chari != '0') {
						charj = consbiparts[ conssizes[2*j+1] ][k];
						if (chari == charj)
							consbiparts[ conssizes[2*j+1] ][k] = '0';
						else if (charj == '3') {
							if (chari == '1')
								consbiparts[ conssizes[2*j+1] ][k] = '2';
							else if (chari == '2')
								consbiparts[ conssizes[2*j+1] ][k] = '1';
							else {
								/* Consensus tree [1] */
								FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR H TO DEVELOPERS\n\n\n");
								exit(1);
							}	
						} else {
							/* Consensus tree [2] */
							FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR I TO DEVELOPERS\n\n\n");
							exit(1);
						}
					}
				}
				/* add link to subnode i in node j */
				charj = consbiparts[ conssizes[2*j+1] ][ conssizes[2*i+1] ];
				if (charj == '0')
					consbiparts[ conssizes[2*j+1] ][ conssizes[2*i+1] ] = '2';
				else if (charj == '1')
					consbiparts[ conssizes[2*j+1] ][ conssizes[2*i+1] ] = '3';
				else {
					/* Consensus tree [3] */
					FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR J TO DEVELOPERS\n\n\n");
					exit(1);
				}
			}
		}
	}
}

/* prototype for recursion */
void writenode(FILE *treefile, int node);

/* write node (writeconsensustree) */
void writenode(FILE *treefile, int node)
{
	int i, first;
	
	fprintf(treefile, "(");
	column++;
	/* write descending nodes */
	first = TRUE;
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '2' ||
		consbiparts[node][i] == '3') {
			if (first) first = FALSE;
			else {
				fprintf(treefile, ",");
				column++;
			}
			if (column > 60) {
				column = 2;
				fprintf(treefile, "\n");
			}
			/* write node i */
			writenode(treefile, i);

			/* reliability value as internal label */
			fprintf(treefile, "%d", consconfid[i]);

			column = column + 3;
		}
	}
	/* write descending taxa */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '1' ||
		consbiparts[node][i] == '3') {
			if (first) first = FALSE;
			else {
				fprintf(treefile, ",");
				column++;
			}
			if (column > 60) {
				column = 2;
				fprintf(treefile, "\n");
			}
			column += fputid(treefile, i);
		}
	}
	fprintf(treefile, ")");
	column++;
}

/* write consensus tree */
void writeconsensustree(FILE *treefile)
{
	int i, first;
	
	column = 1;
	fprintf(treefile, "(");
	column += fputid(treefile, outgroup) + 2;
	fprintf(treefile, ",");
	/* write descending nodes */
	first = TRUE;
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '2' ||
		consbiparts[consincluded][i] == '3') {
			if (first) first = FALSE;
			else {
				fprintf(treefile, ",");
				column++;
			}
			if (column > 60) {
				column = 2;
				fprintf(treefile, "\n");
			}
			/* write node i */
			writenode(treefile, i);

			/* reliability value as internal label */
			fprintf(treefile, "%d", consconfid[i]);

			column = column + 3;
		}
	}
	/* write descending taxa */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '1' ||
		consbiparts[consincluded][i] == '3') {
			if (first) first = FALSE;
			else {
				fprintf(treefile, ",");
				column++;
			}
			if (column > 60) {
				column = 2;
				fprintf(treefile, "\n");
			}
			column += fputid(treefile, i);
		}
	}
	fprintf(treefile, ");\n");
}

/* prototype for recursion */
void nodecoordinates(int node);

/* establish node coordinates (plotconsensustree) */
void nodecoordinates(int node)
{
	int i, ymin, ymax, xcoordinate;

	/* first establish coordinates of descending nodes */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '2' ||
		consbiparts[node][i] == '3') 
			nodecoordinates(i);
	}
	
	/* then establish coordinates of descending taxa */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '1' ||
		consbiparts[node][i] == '3') {
			/* y-coordinate of taxon i */
			ycortax[i] = ytaxcounter;
			ytaxcounter = ytaxcounter - 2;
		}
	}
	
	/* then establish coordinates of this node */
	ymin = 2*Maxspc - 2;
	ymax = 0;
	xcoordinate = 0;
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '2' ||
		consbiparts[node][i] == '3') {
			if (ycor[i] > ymax) ymax = ycor[i];
			if (ycor[i] < ymin) ymin = ycor[i];
			if (xcor[i] > xcoordinate) xcoordinate = xcor[i];
		}
	}
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '1' ||
		consbiparts[node][i] == '3') {
			if (ycortax[i] > ymax) ymax = ycortax[i];
			if (ycortax[i] < ymin) ymin = ycortax[i];
		}
	}
	ycormax[node] = ymax;
	ycormin[node] = ymin;
	ycor[node] = (int) floor(0.5*(ymax + ymin) + 0.5);
	if (xcoordinate == 0) xcoordinate = 9;
	xcor[node] = xcoordinate + 4;
}

/* prototype for recursion */
void drawnode(int node, int xold);

/* drawnode  (plotconsensustree) */
void drawnode(int node, int xold)
{
	int i, j;
	char buf[4];
	
	/* first draw vertical line */
	for (i = ycormin[node] + 1; i < ycormax[node]; i++)
		treepict[xcor[node]][i] = ':';
		
	/* then draw descending nodes */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '2' ||
		consbiparts[node][i] == '3') 
			drawnode(i, xcor[node]);
	}
	
	/* then draw descending taxa */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[node][i] == '1' ||
		consbiparts[node][i] == '3') {
			treepict[xcor[node]][ycortax[i]] = ':';
			for (j = xcor[node] + 1; j < xsize-10; j++)
				treepict[j][ycortax[i]] = '-';
			for (j = 0; j < 10; j++)
				treepict[xsize-10+j][ycortax[i]] = Identif[i][j];	
		}
	}
	
	/* then draw internal edge with consensus value */
	treepict[xold][ycor[node]] = ':';
	treepict[xcor[node]][ycor[node]] = ':';
	for (i = xold + 1; i < xcor[node]-3; i++)
		treepict[i][ycor[node]] = '-';
	sprintf(buf, "%d", consconfid[node]);
	if (consconfid[node] == 100) {
		treepict[xcor[node]-3][ycor[node]] = buf[0];
		treepict[xcor[node]-2][ycor[node]] = buf[1];
		treepict[xcor[node]-1][ycor[node]] = buf[2];	
	} else {
		treepict[xcor[node]-3][ycor[node]] = '-';
		treepict[xcor[node]-2][ycor[node]] = buf[0];
		treepict[xcor[node]-1][ycor[node]] = buf[1];
	}
}

/* plot consensus tree */
void plotconsensustree(FILE *plotfp)
{
	int i, j, yroot, startree;

	/* star tree or no star tree */
	if (consincluded == 0) {
		startree = TRUE;
		consincluded = 1; /* avoids problems with malloc */
	} else
		startree = FALSE;
	
	/* memory for x-y-coordinates of each bipartition */
	xcor = new_ivector(consincluded);
	ycor = new_ivector(consincluded);
	ycormax = new_ivector(consincluded);
	ycormin = new_ivector(consincluded);
	if (startree) consincluded = 0; /* avoids problems with malloc */

	/* y-coordinates of each taxon */
	ycortax = new_ivector(Maxspc);
	ycortax[outgroup] = 0;
	
	/* establish coordinates */
	ytaxcounter = 2*Maxspc - 2;
	
	/* first establish coordinates of descending nodes */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '2' ||
		consbiparts[consincluded][i] == '3') 
			nodecoordinates(i);
	}
	
	/* then establish coordinates of descending taxa */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '1' ||
		consbiparts[consincluded][i] == '3') {
			/* y-coordinate of taxon i */
			ycortax[i] = ytaxcounter;
			ytaxcounter = ytaxcounter - 2;
		}
	}

	/* then establish length of root edge and size of whole tree */
	yroot = 0;
	xsize = 0;
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '2' ||
		consbiparts[consincluded][i] == '3') {
			if (ycor[i] > yroot) yroot = ycor[i];
			if (xcor[i] > xsize) xsize = xcor[i];
		}
	}
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '1' ||
		consbiparts[consincluded][i] == '3') {
			if (ycortax[i] > yroot) yroot = ycortax[i];
		}
	}
	if (xsize == 0) xsize = 9;
	/* size in x direction inclusive one blank on the left */
	xsize = xsize + 6; 
	
	/* change all x-labels so that (0,0) is down-left */
	for (i = 0; i < consincluded; i++)
		xcor[i] = xsize-1-xcor[i];
	
	/* draw tree */
	treepict = new_cmatrix(xsize, 2*Maxspc-1);
	for (i = 0; i < xsize; i++)
		for (j = 0; j < 2*Maxspc-1; j++)
			treepict[i][j] = ' ';
	
	/* draw root */
	for (i = 1; i < yroot; i++)
		treepict[1][i] = ':';
	treepict[1][0] = ':';
	for (i = 2; i < xsize - 10; i++)
		treepict[i][0] = '-';
	for (i = 0; i < 10; i++)
		treepict[xsize-10+i][0] = Identif[outgroup][i];
	
	/* then draw descending nodes */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '2' ||
		consbiparts[consincluded][i] == '3') 
			drawnode(i, 1);
	}
	
	/* then draw descending taxa */
	for (i = 0; i < Maxspc; i++) {
		if (consbiparts[consincluded][i] == '1' ||
		consbiparts[consincluded][i] == '3') {
			treepict[1][ycortax[i]] = ':';
			for (j = 2; j < xsize-10; j++)
				treepict[j][ycortax[i]] = '-';
			for (j = 0; j < 10; j++)
				treepict[xsize-10+j][ycortax[i]] = Identif[i][j];	
		}
	}
	
	/* plot tree */
	for (i = 2*Maxspc-2; i > -1; i--) {
		for (j = 0; j < xsize; j++)
			fputc(treepict[j][i], plotfp);
		fputc('\n', plotfp);
	}	
	
	free_ivector(xcor);
	free_ivector(ycor);
	free_ivector(ycormax);
	free_ivector(ycormin);
	free_ivector(ycortax);
	free_cmatrix(treepict);
}



/******************************************************************************/
/* storing and evaluating quartet branching information                       */
/******************************************************************************/

/* general remarks:

	for a quartet with the taxa a, b, c, d there are
	three possible binary trees:
	
		1)  (a,b)-(c,d)
		2)  (a,c)-(b,d)
		3)  (a,d)-(b,c)
	
	For every quartet information about its branching structure is
	stored. With the functions  readquartet  and  writequartet
	this information can be accessed. For every quartet (a,b,c,d)
	with a < b < c < d (taxa) the branching information is encoded
	using 4 bits:
	
	value          8             4             2             1
	        +-------------+-------------+-------------+-------------+
	        |  not used   |   tree 3    |   tree 2    |   tree 1    |
	        +-------------+-------------+-------------+-------------+

	If the branching structure of the taxa corresponds to one of the
	three trees the corresponding bit is set. If the branching structure
	is unclear because two of the three trees have the same maximum
	likelihood value the corresponding two bits are set. If the branching
	structure is completely unknown all the bits are set (the highest
	bit is always cleared because it is not used).

*/

/* allocate memory for quartets */
unsigned char *mallocquartets(int taxa)
{
	uli nc, numch;
	unsigned char *qinfo;
	
	/* compute number of quartets */
	Numquartets = (uli) taxa*(taxa-1)*(taxa-2)*(taxa-3)/24;
	if (Numquartets % 2 == 0) { /* even number */
		numch = Numquartets/2;
	} else { /* odd number */
		numch = (Numquartets + 1)/2;
	}
	/* allocate memory */
	qinfo = (unsigned char *) malloc(numch * sizeof(unsigned char) );
	if (qinfo == NULL) maerror("quartetinfo in mallocquartets");
	for (nc = 0; nc < numch; nc++) qinfo[nc] = 0;
	return(qinfo);
}

/* free quartet memory */
void freequartets()
{	
	free(quartetinfo);
}

/* read quartet info - a < b < c < d */
unsigned char readquartet(int a, int b, int c, int d)
{
	uli qnum;

	qnum = (uli) a
			+ (uli) b*(b-1)/2
			+ (uli) c*(c-1)*(c-2)/6
			+ (uli) d*(d-1)*(d-2)*(d-3)/24;
	if (qnum % 2 == 0) { /* even number */
		/* bits 0 to 3 */
		return (quartetinfo[qnum/2] & (unsigned char) 15);
	} else { /* odd number */
		/* bits 4 to 7 */
		return ((quartetinfo[(qnum-1)/2] & (unsigned char) 240)>>4);
	}
}

/* write quartet info - a < b < c < d, 0 <= info <= 15 */
void writequartet(int a, int b, int c, int d, unsigned char info)
{
	uli qnum;

	qnum = (uli) a
			+ (uli) b*(b-1)/2
			+ (uli) c*(c-1)*(c-2)/6
			+ (uli) d*(d-1)*(d-2)*(d-3)/24;
	if (qnum % 2 == 0) { /* even number */
		/* bits 0 to 3 */
		quartetinfo[qnum/2] =
			((quartetinfo[qnum/2] & (unsigned char) 240) |
			(info & (unsigned char) 15));
	} else { /* odd number */
		/* bits 4 to 7 */
		quartetinfo[(qnum-1)/2] =
			((quartetinfo[(qnum-1)/2] & (unsigned char) 15) |
			((info & (unsigned char) 15)<<4));
	}
}

/* prototypes */
void openfiletowrite(FILE **, char[], char[]);
void closefile(FILE *);

/* sorts three doubles in descending order */
void sort3doubles(dvector num, ivector order)
{
	if (num[0] > num[1]) {
		if(num[2] > num[0]) {
			order[0] = 2;
			order[1] = 0;
			order[2] = 1;		
		} else if (num[2] < num[1]) {
			order[0] = 0;
			order[1] = 1;
			order[2] = 2;		
		} else {
			order[0] = 0;
			order[1] = 2;
			order[2] = 1;		
		}
	} else {
		if(num[2] > num[1]) {
			order[0] = 2;
			order[1] = 1;
			order[2] = 0;		
		} else if (num[2] < num[0]) {
			order[0] = 1;
			order[1] = 0;
			order[2] = 2;		
		} else {
			order[0] = 1;
			order[1] = 2;
			order[2] = 0;		
		}
	}
}

/* checks out all possible quartets */
void computeallquartets()
{	
	double onethird;
	uli nq;
	unsigned char treebits[3];
	FILE *lhfp;
#	if ! PARALLEL
		int a, b, c, i;
		double qc2, mintogo, minutes, hours, temp;
		double temp1, temp2, temp3;
		unsigned char discreteweight[3];
#	endif
	
	onethird = 1.0/3.0;
	treebits[0] = (unsigned char) 1;
	treebits[1] = (unsigned char) 2;
	treebits[2] = (unsigned char) 4;
	
	if (show_optn) { /* list all unresolved quartets */
		openfiletowrite(&unresfp, UNRESOLVED, "unresolved quartet trees");
		fprintf(unresfp, "List of all completely unresolved quartets:\n\n");
	}

	nq = 0;
	badqs = 0;
	
	/* start timer - percentage of completed quartets */
	time(&time0);
	time1 = time0;
	mflag = 0;
	
#	if PARALLEL
	{
           schedtype sched; 
           int flag;
           MPI_Status stat;
	   int dest = 1;
	   uli qaddr  =0;
	   uli qamount=0;
	   int qblocksent = 0;
	   int apr;
	   uli sq, noq;
	   initsched(&sched, numquarts(Maxspc), PP_NumProcs-1, 4);
	   qamount=sgss(&sched);
           while (qamount > 0) {
              if (PP_emptyslave()) {
	         PP_RecvQuartBlock(0, &sq, &noq, quartetinfo, &apr);
                 qblocksent -= noq;
              }
              dest = PP_getslave();
	      PP_SendDoQuartBlock(dest, qaddr, qamount, (approxqp ? APPROX : EXACT));
              qblocksent += qamount;
              qaddr += qamount;
	      qamount=sgss(&sched);
                 
              MPI_Iprobe(MPI_ANY_SOURCE, PP_QUARTBLOCKSPECS, PP_Comm, &flag, &stat);
              while (flag) {
	         PP_RecvQuartBlock(0, &sq, &noq, quartetinfo, &apr);
                 qblocksent -= noq;
                 MPI_Iprobe(MPI_ANY_SOURCE, PP_QUARTBLOCKSPECS, PP_Comm, &flag, &stat);
              }
           }
           while (qblocksent > 0) {
	      PP_RecvQuartBlock(0, &sq, &noq, quartetinfo, &apr);
              qblocksent -= noq;
           }
	}
#	else /* PARALLEL */

	addtimes(GENERAL, &tarr);
	if (savequartlh_optn) {
		openfiletowrite(&lhfp, ALLQUARTLH, "all quartet likelihoods");
		if (saveqlhbin_optn) writetpqfheader(Maxspc, lhfp, 3);
		else                 writetpqfheader(Maxspc, lhfp, 4); 
	}

	for (i = 3; i < Maxspc; i++) 
		for (c = 2; c < i; c++) 
			for (b = 1; b < c; b++)
				for (a = 0; a < b; a++) {
							nq++;

							/* generate message every 15 minutes */
							/* check timer */
							time(&time2);
							if ( (time2 - time1) > 900) {
								/* every 900 seconds */
								/* percentage of completed quartets */
								if (mflag == 0) {
									FPRINTF(STDOUTFILE "\n");
									mflag = 1;
								}
								qc2 = 100.*nq/Numquartets;
								mintogo = (100.0-qc2) *
									(double) (time2-time0)/60.0/qc2;
								hours = floor(mintogo/60.0);
								minutes = mintogo - 60.0*hours;
								FPRINTF(STDOUTFILE "%.2f%%", qc2);
								FPRINTF(STDOUTFILE " completed (remaining");
								FPRINTF(STDOUTFILE " time: %.0f", hours);
								FPRINTF(STDOUTFILE " hours %.0f", minutes);
								FPRINTF(STDOUTFILE " minutes)\n");
								fflush(STDOUT);
								time1 = time2;
							}
							
							/* maximum likelihood values */
							   
							/* exact or approximate maximum likelihood values */
							compute_quartlklhds(a,b,c,i,&qweight[0],&qweight[1],&qweight[2], (approxqp ? APPROX : EXACT));
							
							if (savequartlh_optn) {
								if (saveqlhbin_optn)
									fwrite(qweight, sizeof(double), 3, lhfp);
								else
									fprintf(lhfp, "(%d,%d,%d,%d)\t%f\t%f\t%f\n", a, b, c, i, 
									qweight[0], qweight[1], qweight[2]); 
							}

							/* sort in descending order */
							sort3doubles(qweight, qworder);

							if (usebestq_optn) {
								sqorder[2] = 2;
								discreteweight[sqorder[2]] = treebits[qworder[0]];
								if (qweight[qworder[0]] == qweight[qworder[1]]) {
								   discreteweight[sqorder[2]] = discreteweight[sqorder[2]] || treebits[qworder[1]];
								   if (qweight[qworder[1]] == qweight[qworder[2]]) {
								      discreteweight[sqorder[2]] = discreteweight[sqorder[2]] || treebits[qworder[2]];
								      discreteweight[sqorder[2]] = 7;
								   } 
								}
							} else {

								/* compute Bayesian weights */
								qweight[qworder[1]] = exp(qweight[qworder[1]]-qweight[qworder[0]]);
								qweight[qworder[2]] = exp(qweight[qworder[2]]-qweight[qworder[0]]);
								qweight[qworder[0]] = 1.0;
								temp = qweight[0] + qweight[1] + qweight[2];
								qweight[0] = qweight[0]/temp;
								qweight[1] = qweight[1]/temp;
								qweight[2] = qweight[2]/temp;
								
								/* square deviations */
								temp1 = 1.0 - qweight[qworder[0]];
								sqdiff[0] = temp1 * temp1 +
								           qweight[qworder[1]] * qweight[qworder[1]] +
								           qweight[qworder[2]] * qweight[qworder[2]];
								discreteweight[0] = treebits[qworder[0]];
     
								temp1 = 0.5 - qweight[qworder[0]];
								temp2 = 0.5 - qweight[qworder[1]];
								sqdiff[1] = temp1 * temp1 + temp2 * temp2 +
								           qweight[qworder[2]] * qweight[qworder[2]];           
								discreteweight[1] = treebits[qworder[0]] + treebits[qworder[1]];
							           
								temp1 = onethird - qweight[qworder[0]];
								temp2 = onethird - qweight[qworder[1]];
								temp3 = onethird - qweight[qworder[2]];
								sqdiff[2] = temp1 * temp1 + temp2 * temp2 + temp3 * temp3;
								discreteweight[2] = (unsigned char) 7;						           
							
								/* sort in descending order */
								sort3doubles(sqdiff, sqorder);
							}
							
							/* determine best discrete weight */
							writequartet(a, b, c, i, discreteweight[sqorder[2]]);
						
							/* counting completely unresolved quartets */
							if (discreteweight[sqorder[2]] == 7) {
								badqs++;
								badtaxon[a]++;
								badtaxon[b]++;
								badtaxon[c]++;
								badtaxon[i]++;
								if (show_optn) {
									fputid10(unresfp, a);
									fprintf(unresfp, "  ");
									fputid10(unresfp, b);
									fprintf(unresfp, "  ");
									fputid10(unresfp, c);
									fprintf(unresfp, "  ");
									fputid(unresfp, i);
									fprintf(unresfp, "\n");
								}
							}
							addtimes(QUARTETS, &tarr);
						}
	if (savequartlh_optn) {
		closefile(lhfp);
	}
	if (show_optn)
		closefile(unresfp);
	if (mflag == 1)
		FPRINTF(STDOUTFILE "\n");
#	endif /* PARALLEL */

}
							
/* check the branching structure between the leaves (not the taxa!)
   A, B, C, and I (A, B, C, I don't need to be ordered). As a result,
   the two leaves that are closer related to each other than to leaf I
   are found in chooseA and chooseB. If the branching structure is
   not uniquely defined, ChooseA and ChooseB are chosen randomly
   from the possible taxa */
void checkquartet(int A, int B, int C, int I)
{
	int i, j, a, b, taxon[5], leaf[5], ipos;
	unsigned char qresult;
	int notunique = FALSE;

	/* The relationship between leaves and taxa is defined by trueID */
	taxon[1] = trueID[A]; /* taxon number */
	leaf[1] = A;          /* leaf number  */
	taxon[2] = trueID[B];
	leaf[2] = B;
	taxon[3] = trueID[C];
	leaf[3] = C;
	taxon[4] = trueID[I];
	leaf[4] = I;

	/* sort for taxa */
	/* Source: Numerical Recipes (PIKSR2.C) */
	for (j = 2; j <= 4; j++) {
		a = taxon[j];
		b = leaf[j];
		i = j-1;
		while (i > 0 && taxon[i] > a) {
			taxon[i+1] = taxon[i];
			leaf[i+1] = leaf[i];
			i--;
		}
		taxon[i+1] = a;
		leaf[i+1] = b;
	}

	/* where is leaf I ? */
	ipos = 1;
	while (leaf[ipos] != I) ipos++;

	/* look at sequence quartet */
	qresult = readquartet(taxon[1], taxon[2], taxon[3], taxon[4]);

	/* chooseA and chooseB */
	do {	
		switch (qresult) {
		
			/* one single branching structure */
		
			/* 001 */
			case 1:		if (ipos == 1 || ipos == 2) {
							chooseA = leaf[3];
							chooseB = leaf[4];
						} else {
							chooseA = leaf[1];
							chooseB = leaf[2];
						}
						notunique = FALSE;
						break;

			/* 010 */
			case 2:		if (ipos == 1 || ipos == 3) {
							chooseA = leaf[2];
							chooseB = leaf[4];
						} else {
							chooseA = leaf[1];
							chooseB = leaf[3];
						}
						notunique = FALSE;
						break;

			/* 100 */
			case 4:		if (ipos == 1 || ipos == 4) {
							chooseA = leaf[2];
							chooseB = leaf[3];
						} else {
							chooseA = leaf[1];
							chooseB = leaf[4];
						}
						notunique = FALSE;
						break;

			/* two possible branching structures */		

			/* 011 */
			case 3:		if (randominteger(2)) qresult = 1;
						else qresult = 2;
						notunique = TRUE;
						break;

			/* 101 */
			case 5:		if (randominteger(2)) qresult = 1;
						else qresult = 4;
						notunique = TRUE;
						break;

			/* 110 */
			case 6:		if (randominteger(2)) qresult = 2;
						else qresult = 4;
						notunique = TRUE;
						break;

			/* three possible branching structures */

			/* 111 */
			case 7:		qresult = (1 << randominteger(3)); /* 1, 2, or 4 */
						notunique = TRUE;
						break;

			default:	/* Program error [checkquartet] */
#if PARALLEL
						FPRINTF(STDOUTFILE "\n\n\n(%2d)HALT: PLEASE REPORT ERROR K-PARALLEL TO DEVELOPERS (%d,%d,%d,%d) = %ld\n\n\n", 
						PP_Myid, taxon[1], taxon[2], taxon[3], taxon[4],
						quart2num(taxon[1], taxon[2], taxon[3], taxon[4]));
#else
						FPRINTF(STDOUTFILE "\n\n\nHALT: PLEASE REPORT ERROR K TO DEVELOPERS\n\n\n");
#endif
						
		}
	} while (notunique);

	return;
}

