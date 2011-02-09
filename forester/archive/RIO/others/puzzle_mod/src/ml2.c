/*
 * ml2.c
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
   - Names of 26 chars.

   !WARNING: Use ONLY together with FORESTER/RIO!
   !For all other puposes download the excellent original!
   
   last modification: 05/22/01


   Node *internalnode(Tree *tr, char **chpp, int *ninode):

   char ident[100], idcomp[11];           -> char ident[100], idcomp[27];

   idcomp[10] = '\0';                     -> idcomp[26] = '\0';

   } while (!stop && (ff != 10));         -> } while (!stop && (ff != 26)); 



*/



#define EXTERN extern

/* prototypes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
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

/* prototypes for two functions of puzzle2.c */
void fputid10(FILE *, int);
int fputid(FILE *, int);


/******************************************************************************/
/* user tree input                                                            */
/******************************************************************************/

/* read user tree, drop all blanks, tabs, and newlines.
   Drop edgelengths (after :) but keep internal
   labels. Check whether all pairs of brackets match. */
void getusertree(FILE *itfp, cvector tr, int maxlen)
{
	int n, brac, ci;
	int comment = 0;

	/* look for opening bracket */
	n = 0;
	brac = 0;
	do {
		ci = fgetc(itfp);
		if (ci == EOF) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (missing start bracket in tree)\n\n\n");
			exit(1);
		}
		if (ci == '[') comment = 1;
		if ((ci == ']') && comment) {
			comment = 0;
			ci = fgetc(itfp);
		}
	} while (comment || ((char) ci != '('));
	tr[n] = (char) ci;
	brac++;
	
	do {
		/* get next character (skip blanks, newlines, and tabs) */
		do {
			ci = fgetc(itfp);
			if (ci == EOF) {
				FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (no more characters in tree)\n\n\n");
				exit(1);
			}
			if (ci == '[') comment = 1;
			if ((ci == ']') && comment) {
				comment = 0;
				ci = fgetc(itfp);
			}
		} while (comment || (char) ci == ' ' || (char) ci == '\n' || (char) ci == '\t');
	
		if ((char) ci == ':') { /* skip characters until a ,) appears  */
			do {
				ci = fgetc(itfp);
				if (ci == EOF) {
					FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (missing ';' or ',' in tree)\n\n\n");
					exit(1);
				}
				if (ci == '[') comment = 1;
				if ((ci == ']') && comment) {
					comment = 0;
					ci = fgetc(itfp);
				}
			} while (comment || ((char) ci != ',' && (char) ci != ')') );
		}
		
		if ((char) ci == '(') {
			brac++;
		}
		if ((char) ci == ')') {
			brac--;
		}

		n++;
		tr[n] = (char) ci;
	
	} while (((char) ci != ';') && (n != maxlen-2));
	
	if (n == maxlen-2) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (tree description too long)\n\n\n");
		exit(1);
	}
	
	if (brac != 0) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (brackets don't match in tree)\n\n\n");
		exit(1);
	}
	
	n++;
	tr[n] = '\0';
}


Node *internalnode(Tree *tr, char **chpp, int *ninode)
{
	Node *xp, *np, *rp;
	int i, j, dvg, ff, stop, numc;
	char ident[100], idcomp[27];  /* CZ 05/22/01 */
	char *idp;

	(*chpp)++;
	if (**chpp == '(') { /* process subgroup */
		
		xp = internalnode(tr, chpp, ninode);
		xp->isop = xp;
		dvg = 1;
		while (**chpp != ')') {
			if (**chpp == '\0') {
				FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (unexpected end of tree)\n\n\n");
				exit(1);
			}
			dvg++;
			/* insert edges around node */
			np = internalnode(tr, chpp, ninode);
			np->isop = xp->isop;
			xp->isop = np; 
			xp = np;
		}
		/* closing bracket reached */
		
		(*chpp)++;
		if (dvg < 2) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (only one OTU inside pair of brackets)\n\n\n");
			exit(1);
		}
		
		if ((*ninode) >= Maxspc-3) { /* all internal nodes already used */
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (no unrooted tree)\n\n\n");
			exit(1);
		}

		rp = tr->ibrnchp[*ninode];
		rp->isop = xp->isop;
		xp->isop = rp;

		for (j = 0; j < Numspc; j++)
			rp->paths[j] = 0;
		xp = rp->isop;
		while (xp != rp) {
			for (j = 0; j < Numspc; j++) {
				if (xp->paths[j] == 1)
					rp->paths[j] = 1;
			}
			xp = xp->isop;
		}
		(*ninode)++;

		if ((**chpp) == ',' || (**chpp) == ')') return rp->kinp;
		if ((**chpp) == '\0')  {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (unexpected end of tree)\n\n\n");
			exit(1);
		}

		/* read internal label into rp->label (max. 20 characters) */
		rp->label = new_cvector(21);
		(rp->label)[0] = **chpp;
		(rp->label)[1] = '\0';
		for (numc = 1; numc < 20; numc++) {
			(*chpp)++;
			if ((**chpp) == ',' || (**chpp) == ')') return rp->kinp;
			if ((**chpp) == '\0')  {
				FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (unexpected end of tree)\n\n\n");
				exit(1);
			}
			(rp->label)[numc] = **chpp;
			(rp->label)[numc+1] = '\0';
		}	
		do { /* skip the rest of the internal label */
			(*chpp)++;
			if ((**chpp) == '\0')  {
				FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (unexpected end of tree)\n\n\n");
				exit(1);
			}
		} while (((**chpp) != ',' && (**chpp) != ')'));
			
		return rp->kinp;
		
	} else { /* process species names */
		/* read species name */
		for (idp = ident; **chpp != ',' &&
			**chpp != ')' && **chpp != '\0'; (*chpp)++) {
			*idp++ = **chpp;	
		}
		*idp = '\0';
		/* look for internal number */
		idcomp[26] = '\0'; /* CZ 05/22/01 */
		
		for (i = 0; i < Maxspc; i++) {
			ff = 0;
			stop = FALSE;
			do {
				idcomp[ff] = Identif[i][ff];
				ff++;
				if (idcomp[ff-1] == ' ') stop = TRUE;
			} while (!stop && (ff != 26)); /* CZ 05/22/01 */
			if (stop) idcomp[ff-1] = '\0';
			
			if (!strcmp(ident, idcomp)) {
				if (usedtaxa[i]) {
					FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (multiple occurence of sequence '");
					FPRINTF(STDOUTFILE "%s' in tree)\n\n\n", ident);
					exit(1);
				}
				usedtaxa[i] = TRUE;
				return tr->ebrnchp[i]->kinp;
			}
		}
		
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (unknown sequence '%s' in tree)\n\n\n", ident);
		exit(1);
	}
	return NULL; /* never returned but without some compilers complain */
}

/* make tree structure, the tree description may contain internal
    labels but no edge lengths */
void constructtree(Tree *tr, cvector strtree)
{
	char *chp;
	int ninode, i;
	int dvg, numc;
	Node *xp, *np;

	ninode = 0;
	chp = strtree;
	usedtaxa = new_ivector(Maxspc);
	for (i = 0; i < Maxspc; i++) usedtaxa[i] = FALSE;

	xp = internalnode(tr, &chp, &ninode);
	xp->isop = xp;
	dvg = 1;
	while (*chp != ')') { /* look for closing bracket */
		if (*chp == '\0') {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (unexpected end of tree)\n\n\n");
			exit(1);
		}
		dvg++;
		/* insert edges around node */
		np = internalnode(tr, &chp, &ninode);
		np->isop = xp->isop;
		xp->isop = np; 
		xp = np;
	}
	
	for (i = 0; i < Maxspc; i++)
		if (usedtaxa[i] == FALSE) {
			FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (sequences missing in tree)\n\n\n");
			exit(1);	
		}
	
	/* closing bracket reached */
	if (dvg < 3) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (no unrooted tree)\n\n\n");
		exit(1);
	}
	tr->rootp = xp;
	Numibrnch = ninode;
	Numbrnch = Numspc + ninode;

	chp++;
	if (*chp == ';' || *chp == '\0') {
		free_ivector(usedtaxa);
		return;
	}

	/* copy last internal label (max. 20 characters) */
	xp->label = new_cvector(21);
	(xp->label)[0] = *chp;
	(xp->label)[1] = '\0';		
	for (numc = 1; numc < 20; numc++) {
		chp++;
		if (*chp == ';' || *chp == '\0') {
			free_ivector(usedtaxa);
			return;
		} else {
			(xp->label)[numc] = *chp;
			(xp->label)[numc+1] = '\0';
		}
	}
	free_ivector(usedtaxa);	
	return;
}


/* remove possible basal bifurcation */
void removebasalbif(cvector strtree)
{
	int n, c, brak, cutflag, h;
	
	/* check how many OTUs on basal level */
	n = 0;
	c = 0;
	brak = 0;
	do {
		if (strtree[n] == '(') brak++;
		if (strtree[n] == ')') brak--;
		
		if (strtree[n] == ',' && brak == 1) c++; /* number of commas in outer bracket */
		
		n++;
	} while (strtree[n] != '\0');
	
	/* if only 1 OTU inside outer bracket stop now */
	if (c == 0) {
		FPRINTF(STDOUTFILE "\n\n\nUnable to proceed (Only 1 OTU inside outer bracket in tree)\n\n\n");
		exit(1);
	}	
	
	/* if only 2 OTUs inside outer bracket delete second pair of
	   brackets from the right to remove basal bifurcation */

	if (c == 1) {
			
		n = 0;
		brak = 0;
		cutflag = 0; /* not yet cutted */
		h = 0;
		do {
			if (strtree[n] == '(') brak++;
			if (strtree[n] == ')') brak--;
			
			if (brak == 2 && cutflag == 0) cutflag = 1; /* cutting */
			if (brak == 1 && cutflag == 1) {
				cutflag = 2; /* cutted */
				/* leave out internal label */
				do {
					h++;	
				} while (strtree[n+h] != ')' && strtree[n+h] != ',');

			}
			
			if (cutflag == 1) strtree[n] = strtree[n+1];
			if (cutflag == 2) strtree[n-1] = strtree[n+h];
		
			n++;
		} while (strtree[n] != '\0');
	}
}


void makeusertree(FILE *itfp)
{
	cvector strtree;
	
	strtree = new_cvector(23*Maxspc); /* for treefile */ 
	getusertree(itfp, strtree, 23*Maxspc);
	removebasalbif(strtree);
	constructtree(Ctree, strtree);
	free_cvector(strtree);
}


/******************************************************************************/
/* memory organisation for maximum likelihood tree                            */
/******************************************************************************/

/* initialise new tree */
Tree *new_tree(int maxspc, int numptrn, cmatrix seqconint)
{
	int n, i, maxibrnch;
	Tree *tr;
	Node *dp, *up;
	
	maxibrnch = maxspc - 3;
	heights = (Node **) malloc((unsigned)(maxspc-2) * sizeof(Node *));
	if (heights == NULL) maerror("heights in new_tree");
	tr = (Tree *) malloc(sizeof(Tree));
	if (tr == NULL) maerror("tr in new_tree");
	tr->ebrnchp = (Node **) malloc((unsigned)maxspc * sizeof(Node *));
	if (tr->ebrnchp == NULL) maerror("ebrnchp in new_tree");
	tr->ibrnchp = (Node **) malloc((unsigned)maxibrnch * sizeof(Node *));
	if (tr->ibrnchp == NULL) maerror("ibrnchp in new_tree");
	tr->condlkl = new_dmatrix(numcats, numptrn);	
	for (n = 0; n < maxspc; n++) {
		dp = (Node *) malloc(sizeof(Node));
		if (dp == NULL) maerror("dp in new_tree");
		up = (Node *) malloc(sizeof(Node));
		if (up == NULL) maerror("up in new_tree");
		dp->isop = NULL;
		up->isop = NULL;
		dp->kinp = up;
		up->kinp = dp;
		dp->descen = TRUE;
		up->descen = FALSE;
		dp->number = n;
		up->number = n;
		dp->length = 0.0;
		up->length = 0.0;
		dp->lengthc = 0.0;
		up->lengthc = 0.0;
		dp->varlen = 0.0;
		up->varlen = 0.0;
		dp->paths = new_ivector(maxspc);
		up->paths = dp->paths;
		for (i = 0; i < maxspc; i++) dp->paths[i] = 0;
		dp->paths[n] = 1;
		dp->eprob = seqconint[n];
		up->eprob = NULL;
		dp->partials = NULL;
		up->partials = new_dcube(numcats, numptrn, tpmradix);
		tr->ebrnchp[n] = dp;
		up->label = NULL;
		dp->label = NULL;
	}
	for (n = 0; n < maxibrnch; n++) {
		dp = (Node *) malloc(sizeof(Node));
		if (dp == NULL) maerror("dp in new_tree");
		up = (Node *) malloc(sizeof(Node));
		if (up == NULL) maerror("up in new_tree");
		dp->isop = NULL;
		up->isop = NULL;
		dp->kinp = up;
		up->kinp = dp;
		dp->descen = TRUE;
		up->descen = FALSE;
		dp->number = n;
		up->number = n;
		dp->length = 0.0;
		up->length = 0.0;
		dp->lengthc = 0.0;
		up->lengthc = 0.0;
		dp->varlen = 0.0;
		up->varlen = 0.0;
		dp->paths = new_ivector(maxspc);
		up->paths = dp->paths;
		for (i = 0; i < maxspc; i++) dp->paths[i] = 0;
		dp->eprob = NULL;
		up->eprob = NULL;
		dp->partials = new_dcube(numcats, numptrn, tpmradix);
		up->partials = new_dcube(numcats, numptrn, tpmradix);
		tr->ibrnchp[n] = dp;
		up->label = NULL;
		dp->label = NULL;
	}
	tr->rootp = NULL;
	
	/*
	 * reserve memory for lengths of the tree branches
	 * and for the distance matrix as a vector
	 * (needed for LS estimation of tree branch lengths)
	 */ 
	 
	Brnlength = new_dvector(2 * maxspc - 3); 
	Distanvec = new_dvector((maxspc * (maxspc - 1)) / 2);

	return tr;
}


/* initialise quartet tree */
Tree *new_quartet(int numptrn, cmatrix seqconint)
{
	int n, i;
	Tree *tr;
	Node *dp, *up;

	heights = (Node **) malloc((unsigned)2 * sizeof(Node *));
	if (heights == NULL) maerror("heights in new_quartet");
	/* reserve memory for tree */
	tr = (Tree *) malloc(sizeof(Tree));
	if (tr == NULL) maerror("tr in new_quartet");
	tr->ebrnchp = (Node **) malloc((unsigned) 4 * sizeof(Node *));
	if (tr->ebrnchp == NULL) maerror("ebrnchp in new_quartet");
	tr->ibrnchp = (Node **) malloc((unsigned) sizeof(Node *));
	if (tr->ibrnchp == NULL) maerror("ibrnchp in new_quartet");
	tr->condlkl = new_dmatrix(numcats, numptrn);
	/* reserve memory for nodes */
	for (n = 0; n < 4; n++) {
		dp = (Node *) malloc(sizeof(Node));
		if (dp == NULL) maerror("dp in new_quartet");
		up = (Node *) malloc(sizeof(Node));
		if (up == NULL) maerror("dp in new_quartet");
		dp->isop = NULL;
		dp->kinp = up;
		up->kinp = dp;
		dp->descen = TRUE;
		up->descen = FALSE;
		dp->number = n;
		up->number = n;
		dp->length = 0.0;
		up->length = 0.0;
		dp->lengthc = 0.0;
		up->lengthc = 0.0;
		dp->varlen = 0.0;
		up->varlen = 0.0;
		dp->paths = new_ivector(4);
		up->paths = dp->paths;
		for (i = 0; i < 4; i++) dp->paths[i] = 0;
		dp->paths[n] = 1;
		dp->eprob = seqconint[n]; /* make quartet (0,1)-(2,3) as default */
		up->eprob = NULL;		
		dp->partials = NULL;
		up->partials = new_dcube(numcats, numptrn, tpmradix);
		tr->ebrnchp[n] = dp;
	}

	/* reserve memory for internal branch */	
	dp = (Node *) malloc(sizeof(Node));
	if (dp == NULL) maerror("dp in new_quartet");	
	up = (Node *) malloc(sizeof(Node));
	if (up == NULL) maerror("dp in new_quartet");	
	dp->isop = tr->ebrnchp[3]->kinp; /* connect internal branch */
	up->isop = tr->ebrnchp[0]->kinp;
	dp->kinp = up;
	up->kinp = dp;
	dp->descen = TRUE;
	up->descen = FALSE;
	dp->number = 0;
	up->number = 0;
	dp->length = 0.0;
	up->length = 0.0;
	dp->lengthc = 0.0;
	up->lengthc = 0.0;
	dp->varlen = 0.0;
	up->varlen = 0.0;
	dp->paths = new_ivector(4);
	up->paths = dp->paths;
	up->paths[0] = 0;
	up->paths[1] = 0;
	up->paths[2] = 1;
	up->paths[3] = 1;
	dp->eprob = NULL;
	up->eprob = NULL;	
	dp->partials = new_dcube(numcats, numptrn, tpmradix);
	up->partials = new_dcube(numcats, numptrn, tpmradix);	
	tr->ibrnchp[0] = dp;
	
	/* place root */
	tr->rootp = up;

	/* connect external branches */ 
	tr->ebrnchp[0]->kinp->isop = tr->ebrnchp[1]->kinp;
	tr->ebrnchp[1]->kinp->isop = tr->rootp;
	tr->ebrnchp[3]->kinp->isop = tr->ebrnchp[2]->kinp;
	tr->ebrnchp[2]->kinp->isop = tr->rootp->kinp;
	
	/*
	 * reserve memory for lengths of the five branches
	 * of a quartet and for the six possible distances
	 * (needed for LS estimation of branch lengths)
	 */
	Brnlength = new_dvector(NUMQBRNCH); 
	Distanvec = new_dvector(NUMQSPC*(NUMQSPC-1)/2);

	return tr;
}


/* free tree memory */
void free_tree(Tree *tr, int taxa)
{	
	int n;
	Node *dp, *up;

	free(heights);
	free_dmatrix(tr->condlkl);
	for (n = 0; n < taxa; n++) {
		dp = tr->ebrnchp[n];
		up = dp->kinp;
		free_ivector(dp->paths);		
		free_dcube(up->partials);		
		free(dp);
		free(up);
	}
	free(tr->ebrnchp);
	for (n = 0; n < (taxa-3); n++) {
		dp = tr->ibrnchp[n];
		up = dp->kinp;
		free_dcube(dp->partials);
		free_dcube(up->partials);
		free_ivector(dp->paths);
		free(dp);
		free(up);
	}
	free(tr->ibrnchp);
	free(tr);
	free_dvector(Brnlength); /* branch lengths (for LS estimation) */
	free_dvector(Distanvec); /* distances (for LS estimation) */
}


/* make (a,b)-(c,d) quartet

	a ---+     +--- c
	     +-----+
	b ---+     +--- d

	species numbers range from 0 to Maxspc - 1  */

void make_quartet(int a, int b, int c, int d)
{
	/* place sequences */
	Ctree->ebrnchp[0]->eprob = Seqpat[a];
	Ctree->ebrnchp[1]->eprob = Seqpat[b];
	Ctree->ebrnchp[2]->eprob = Seqpat[c];
	Ctree->ebrnchp[3]->eprob = Seqpat[d];
	
	/* make distance vector */
	Distanvec[0] = Distanmat[b][a];
	Distanvec[1] = Distanmat[c][a];
	Distanvec[2] = Distanmat[c][b];
	Distanvec[3] = Distanmat[d][a];
	Distanvec[4] = Distanmat[d][b];
	Distanvec[5] = Distanmat[d][c];
}

/* write distance matrix as vector */
void changedistan(dmatrix distanmat, dvector distanvec, int numspc)
{
	int i, j, k;

	for (k = 0, i = 1; i < numspc; i++) {
		for (j = 0; j < i; j++, k++)
			distanvec[k] = distanmat[i][j];
	}
}


/******************************************************************************/
/* computation of maximum likelihood tree                                     */
/******************************************************************************/


/* compute the likelihood for (a,b)-(c,d) quartet */
double quartet_lklhd(int a, int b, int c, int d)
{
	/* reserve memory for quartet if necessary */
	if (mlmode != 1) { /* no quartet tree */
		if (Ctree != NULL)
			free_tree(Ctree, Numspc);
		Ctree = new_quartet(Numptrn, Seqpat);
		Numbrnch = NUMQBRNCH;
		Numibrnch = NUMQIBRNCH;
		Numspc = NUMQSPC;
		mlmode = 1;
	}
	
	/* make (a,b)-(c,d) quartet */
	make_quartet(a,b,c,d);

	clockmode = 0; /* nonclocklike branch lengths */
	
	/* least square estimate for branch length */	
	lslength(Ctree, Distanvec, Numspc, Numibrnch, Brnlength);

	/* compute likelihood */
	Ctree->lklhd = optlkl(Ctree);

	return Ctree->lklhd;
}


/* compute the approximate likelihood for (a,b)-(c,d) quartet */
double quartet_alklhd(int a, int b, int c, int d)
{
	/* reserve memory for quartet if necessary */
	if (mlmode != 1) { /* no quartet tree */
		if (Ctree != NULL)
			free_tree(Ctree, Numspc);
		Ctree = new_quartet(Numptrn, Seqpat);
		Numbrnch = NUMQBRNCH;
		Numibrnch = NUMQIBRNCH;
		Numspc = NUMQSPC;
		mlmode = 1; 
	}
	
	/* make (a,b)-(c,d) quartet */
	make_quartet(a,b,c,d);

	clockmode = 0; /* nonclocklike branch lengths */
	
	/* least square estimate for branch length */	
	lslength(Ctree, Distanvec, Numspc, Numibrnch, Brnlength);

	/* compute likelihood */
	Ctree->lklhd = treelkl(Ctree);

	return Ctree->lklhd;
}


/* read usertree from file to memory */
void readusertree(FILE *ifp)
{
	/* reserve memory for tree if necessary */
	if (mlmode != 2) { /* no tree */
		if (Ctree != NULL)
			free_tree(Ctree, Numspc);
		Ctree = new_tree(Maxspc, Numptrn, Seqpat);
		Numbrnch = 2*Maxspc-3;
		Numibrnch = Maxspc-3;
		Numspc = Maxspc;
		mlmode = 2;
	}

	/* read tree */
	makeusertree(ifp);
}


/* compute the likelihood of a usertree */
double usertree_lklhd()
{
	/* be sure to have a usertree in memory and
	   to have pairwise distances computed */

	clockmode = 0; /* nonclocklike branch lengths */

	/* least square estimate for branch length */
	changedistan(Distanmat, Distanvec, Numspc);	
	lslength(Ctree, Distanvec, Numspc, Numibrnch, Brnlength);

	/* compute likelihood */
	Ctree->lklhd = optlkl(Ctree);

	return Ctree->lklhd;
}


/* compute the approximate likelihood of a usertree */
double usertree_alklhd()
{
	/* be sure to have a usertree in memory and
	   to have pairwise distances computed */

	clockmode = 0; /* nonclocklike branch lengths */

	/* least square estimate for branch length */
	changedistan(Distanmat, Distanvec, Numspc);
	lslength(Ctree, Distanvec, Numspc, Numibrnch, Brnlength);

	/* compute likelihood */
	Ctree->lklhd = treelkl(Ctree);

	return Ctree->lklhd;
}


/* preparation for ML analysis */
void mlstart()
{
	/* number of states and code length */
	tpmradix = gettpmradix();
	
	/* declare variables */
	Eval = new_dvector(tpmradix);
	Evec = new_dmatrix(tpmradix,tpmradix);
	Ievc = new_dmatrix(tpmradix,tpmradix);
	iexp = new_dmatrix(tpmradix,tpmradix);
	Alias = new_ivector(Maxsite);
	
	/* process sequence information */
	evaluateseqs();
	bestrate = new_ivector(Numptrn);

	/* compute transition probability matrix */	
	tranprobmat();

	/* non-zero rate categories */
	Rates = new_dvector(numcats);
	updaterates();
	ltprobr = new_dcube(numcats, tpmradix,tpmradix);

	/* compute distance matrix */
	Distanmat = new_dmatrix(Maxspc, Maxspc);
	initdistan();
		
	/* initialize tree pointer for quartet tree */
	mlmode = 1;
	Ctree = new_quartet(Numptrn, Seqpat);
	Numbrnch = NUMQBRNCH;
	Numibrnch = NUMQIBRNCH;
	Numspc = NUMQSPC;

	/* computing ML distances */
	computedistan();
}


/* recompute ml distances for quartet only */
void distupdate(int a, int b, int c, int d)
{
	/* update distance matrix */
	/* consider only entries relevant to quartet */
	Distanmat[a][b] = mldistance(a, b);
	Distanmat[b][a] = Distanmat[a][b];
	Distanmat[a][c] = mldistance(a, c);
	Distanmat[c][a] = Distanmat[a][c];
	Distanmat[a][d] = mldistance(a, d);
	Distanmat[d][a] = Distanmat[a][d];
	Distanmat[b][c] = mldistance(b, c);
	Distanmat[c][b] = Distanmat[b][c];
	Distanmat[b][d] = mldistance(b, d);
	Distanmat[d][b] = Distanmat[b][d];
	Distanmat[c][d] = mldistance(c, d);
	Distanmat[d][c] = Distanmat[c][d];
}


/* cleanup after ML analysis */
void mlfinish()
{
	if (Ctree != NULL)
		free_tree(Ctree, Numspc);
	free_ivector(bestrate);
	free_ivector(Alias);
	free_cmatrix(Seqpat);
	free_ivector(constpat);
	free_ivector(Weight);
	free_dmatrix(Distanmat);
	free_dvector(Eval);
	free_dmatrix(Evec);
	free_dmatrix(Ievc);
	free_dvector(Rates);
	free_dcube(ltprobr);
	free_dmatrix(iexp);
}


/******************************************************************************/
/* tree output                                                                */
/******************************************************************************/


#define MAXOVER    50
#define MAXLENG    30
#define MAXCOLUMN  80


void prbranch(Node *up, int depth, int m, int maxm,
	ivector umbrella, ivector column, FILE *outfp)
{
	int i, num, n, maxn, lim;
	Node *cp;
	char bch;

	if ((int)((clockmode ? up->lengthc : up->length) * Proportion) >= MAXOVER) {
		column[depth] = MAXLENG;
		bch = '+';
	} else {
		column[depth] = (int)((clockmode ? up->lengthc : up->length) * Proportion) + 3;
		bch = '-';
	}

	if (up->isop == NULL) { /* external branch */
		num = up->number + 1; /* offset */
		if (m == 1) umbrella[depth - 1] = TRUE;
		for (i = 0; i < depth; i++) {
			if (umbrella[i])
				fprintf(outfp, "%*c", column[i], ':');
			else
				fprintf(outfp, "%*c", column[i], ' ');
		}
		if (m == maxm)
			umbrella[depth - 1] = FALSE;
		for (i = 0, lim = column[depth] - 3; i < lim; i++)
			fputc(bch, outfp);
		fprintf(outfp, "-%d ", num);
				
		fputid(outfp, up->number);
		
		
		fputc('\n', outfp);
		fputc(' ', outfp);
		return;
	}

	num = up->number + 1 + Numspc; /* offset, internal branch */
	for (cp = up->isop, maxn = 0; cp != up; cp = cp->isop, maxn++)
		;
	for (cp = up->isop, n = 1; cp != up; cp = cp->isop, n++) {
		prbranch(cp->kinp, depth + 1, n, maxn, umbrella, column, outfp);
		if (m == 1 && n == maxn / 2) umbrella[depth - 1] = TRUE;
		if (n != maxn) {
			for (i = 0; i < depth; i++) {
				if (umbrella[i])
					fprintf(outfp, "%*c", column[i], ':');
				else
					fprintf(outfp, "%*c", column[i], ' ');
			}
			if (n == maxn / 2) { /* internal branch */
				for (i = 0, lim = column[depth] - 3; i < lim; i++)
					fputc(bch, outfp);
				if (num < 10)
					fprintf(outfp, "--%d", num);
				else if (num < 100)
					fprintf(outfp, "-%2d", num);
				else
					fprintf(outfp, "%3d", num);
			} else {
				if (umbrella[depth])
					fprintf(outfp, "%*c", column[depth], ':');
				else
					fprintf(outfp, "%*c", column[depth], ' ');
			}
			fputc('\n', outfp);
			fputc(' ', outfp);
		}
		if (m == maxm) umbrella[depth - 1] = FALSE;
	}
	return;
}


void getproportion(double *proportion, dvector distanvec, int numspc)
{
	int i, maxpair;
	double maxdis;
	
	maxpair = (numspc*(numspc-1))/2;

	maxdis = 0.0;
	for (i = 0; i < maxpair; i++) {
		if (distanvec[i] > maxdis) {
			maxdis = distanvec[i];
		}
	}
	*proportion = (double) MAXCOLUMN / (maxdis * 3.0);
	if (*proportion > 1.0) *proportion = 1.0;
}


void prtopology(FILE *outfp)
{
	int n, maxn, depth;
	ivector umbrella;
	ivector column;
	Node *cp, *rp;
	
	getproportion(&Proportion, Distanvec, Numspc);

	umbrella = new_ivector(Numspc);
	column = new_ivector(Numspc);

	for (n = 0; n < Numspc; n++) {
		umbrella[n] = FALSE;
		column[n] = 3;
	}
	column[0] = 1;
	
	fputc(' ', outfp);

	/* original code: rp = Ctree->rootp */
	/* but we want to print the first group in the 
	   trichotomy as outgroup at the bottom! */
	rp = Ctree->rootp->isop;
	
	for (maxn = 1, cp = rp->isop; cp != rp; cp = cp->isop, maxn++)
		;
	depth = 1;
	n = 0;

	cp = rp;
	do {
		cp = cp->isop;
		n++;
		prbranch(cp->kinp, depth, n, maxn, umbrella, column, outfp);
		if (cp != rp) fprintf(outfp, "%*c\n ", column[0], ':');
	} while (cp != rp);

	free_ivector(umbrella);
	free_ivector(column);
}


/* print unrooted tree file with branch lengths */
void fputphylogeny(FILE *fp)
{
	Node *cp, *rp;
	int n;

	cp = rp = Ctree->rootp;
	putc('(', fp);
	n = 1;
	do {
		cp = cp->isop->kinp;
		if (cp->isop == NULL) { /* external node */
			if (n > 60) {
				fprintf(fp, "\n");
				n = 2;
			}
			n += fputid(fp, cp->number);
			fprintf(fp, ":%.5f", ((clockmode ? cp->lengthc : cp->length))*0.01);
			n += 7;
			cp = cp->kinp;
		} else { /* internal node */
			if (cp->descen) {
				if (n > 60) {
					fprintf(fp, "\n");
					n = 1;
				}
				putc('(', fp);
				n++;
			} else {
				putc(')', fp);
				n++;
				if (n > 60) {
					fprintf(fp, "\n");
					n = 1;
				}
				/* internal label */
				if (cp->kinp->label != NULL) {
					fprintf(fp, "%s", cp->kinp->label);
					n += strlen(cp->kinp->label);
				}
				fprintf(fp, ":%.5f", ((clockmode ? cp->lengthc : cp->length))*0.01);
				n += 7;
			}
		}
		if (!cp->descen && !cp->isop->descen && cp != rp) {
			putc(',', fp); /* not last subtree */
			n++;
		}
	} while (cp != rp);
	fprintf(fp, ")");
	/* internal label */
	if (cp->label != NULL)
		fprintf(fp, "%s", cp->label);
	fprintf(fp, ";\n");
}


void resulttree(FILE *outfp)
{
	int n, ne, closeflag;
	Node *ep, *ip;
	double blen;
	
	closeflag = FALSE;
	
	if (clockmode) {
		fprintf(outfp, "\n         branch  length     nc/c");
		fprintf(outfp, "   branch  length     nc/c (= non-clock/clock)\n");
	} else {
		fprintf(outfp, "\n         branch  length     S.E.");
		fprintf(outfp, "   branch  length     S.E.\n");
	}
	for (n = 0; n < Numspc; n++) {
		ep = Ctree->ebrnchp[n];
		ne = ep->number;
		fputid10(outfp, ne);
		fputs("  ", outfp);
		fprintf(outfp, "%3d", ne + 1);
		blen = (clockmode ? ep->lengthc : ep->length);
		fprintf(outfp, "%9.5f", blen*0.01);
		if (blen < 5.0*MINARC || blen > 0.95*MAXARC) closeflag = TRUE;
		if (clockmode)
			fprintf(outfp, "%9.3f", (ep->length)/(ep->lengthc));
		else
			fprintf(outfp, "%9.5f", 0.01*sqrt(ep->kinp->varlen));	
		if (n < Numibrnch) {
			ip = Ctree->ibrnchp[n];
			fprintf(outfp, "%8d", n + 1 + Numspc);
			blen = (clockmode ? ip->lengthc : ip->length);
			fprintf(outfp, "%9.5f", blen*0.01);
			if (blen < 5.0*MINARC || blen > 0.95*MAXARC) closeflag = TRUE;
			if (clockmode)
				fprintf(outfp, "%9.3f", (ip->length)/(ip->lengthc));
			else
				fprintf(outfp, "%9.5f", 0.01*sqrt(ip->kinp->varlen));	
			fputc('\n', outfp);
		} else {
			if (n == Numspc - 3) {
				fputc('\n', outfp);
			} else if (n == Numspc - 2) {
				if (clockmode) {
					if (!Convergc) 
						fprintf(outfp, "     No convergence after %d iterations!\n", Numitc);
					else
						fprintf(outfp, "     %d iterations until convergence\n", Numitc);
				} else {
					if (!Converg) 
						fprintf(outfp, "     No convergence after %d iterations!\n", Numit);
					else
						fprintf(outfp, "     %d iterations until convergence\n", Numit);				
				}		
			} else if (n == Numspc - 1) {
				fprintf(outfp, "     log L: %.2f\n", (clockmode ? Ctree->lklhdc : Ctree->lklhd));
			} else {
				fputc('\n', outfp);
			}
		}
	}
	if(closeflag)
		fprintf(outfp, "\nWARNING --- at least one branch length is close to an internal boundary!\n");
}


/******************************************************************************/
/* Neighbor-joining tree                                                      */
/******************************************************************************/


/* compute NJ tree and write to file */
void njtree(FILE *fp)
{
	/* reserve memory for tree if necessary */
	if (mlmode != 3) { /* no tree */
		if (Ctree != NULL)
			free_tree(Ctree, Numspc);
		Ctree = new_tree(Maxspc, Numptrn, Seqpat);
		Numbrnch = 2*Maxspc-3;
		Numibrnch = Maxspc-3;
		Numspc = Maxspc;
		mlmode = 3;
	}

	/* construct NJ tree from distance matrix */
	njdistantree(Ctree);
	
	fputphylogeny(fp);
}


/* construct  NJ tree from distance matrix */
void njdistantree(Tree *tr)
{
	int i, j, otui=0, otuj=0, otuk, nsp2, cinode, step, restsp, k;
	double dij, bix, bjx, bkx, sij, smax, dnsp2;
	dvector r;
	dmatrix distan;
	Node **psotu, *cp, *ip, *jp, *kp;

	distan = new_dmatrix(Maxspc,Maxspc);
	for (i = 0; i < Maxspc; i++)
		for (j = 0; j < Maxspc; j++)
			distan[i][j] = Distanmat[i][j];

	nsp2 = Maxspc - 2;
	dnsp2 = 1.0 / nsp2;
	
	r = new_dvector(Maxspc);
	
	psotu = (Node **) malloc((unsigned)Maxspc * sizeof(Node *));
	if (psotu == NULL) maerror("psotu in njdistantree");

	/* external branches are start OTUs */
	for (i = 0; i < Maxspc; i++)
		psotu[i] = tr->ebrnchp[i]->kinp;
	
	restsp = Maxspc;
	cinode = 0; /* counter for internal nodes */
	
	for (step = 0; restsp > 3; step++) { /* NJ clustering steps */
	
		for (i = 0; i < Maxspc; i++) {
			if (psotu[i] != NULL) {
				for (j = 0, r[i] = 0.0; j < Maxspc; j++)
					if (psotu[j] != NULL)
						r[i] += distan[i][j];
			}
		}
		
		smax = -1.0;
		for (i = 0; i < Maxspc-1; i++) {
			if (psotu[i] != NULL) {
				
				for (j = i+1; j < Maxspc; j++) {
					if (psotu[j] != NULL)
					{
						sij = ( r[i] + r[j] ) * dnsp2 - distan[i][j];
			
						if (sij > smax) {
							smax = sij;
							otui = i;
							otuj = j;
						}
					}
				}
			}
		}

		/* new pair: otui and otuj */

		dij = distan[otui][otuj];
		bix = (dij + r[otui]/nsp2 - r[otuj]/nsp2) * 0.5;
		bjx = dij - bix;
		
		cp = tr->ibrnchp[cinode];
		
		ip = psotu[otui];
		jp = psotu[otuj];
		cp->isop = ip;
		ip->isop = jp;
		jp->isop = cp;
		ip->length = bix;
		jp->length = bjx;
		ip->kinp->length = ip->length;
		jp->kinp->length = jp->length;
		
		cp = cp->kinp;
		
		for (k = 0; k < Maxspc; k++)
		{
			if (psotu[k] != NULL && k != otui && k != otuj)
			{
				dij = (distan[otui][k] + distan[otuj][k] - distan[otui][otuj]) * 0.5;
				distan[otui][k] = dij;
				distan[k][otui] = dij;
			}
		}
		distan[otui][otui] = 0.0;

		psotu[otui] = cp;
		psotu[otuj] = NULL;
		
		cinode++;
		
		restsp--;
		nsp2--;
		dnsp2 = 1.0 / nsp2;
	}
 
	otui = otuj = otuk = -1;
	for (i = 0; i < Maxspc; i++)
	{
		if (psotu[i] != NULL) {
			if (otui == -1) otui = i;
			else if (otuj == -1) otuj = i;
			else otuk = i;
		}
	}
	bix = (distan[otui][otuj] + distan[otui][otuk] - distan[otuj][otuk]) * 0.5;
	bjx = distan[otui][otuj] - bix;
	bkx = distan[otui][otuk] - bix;
	ip = psotu[otui];
	jp = psotu[otuj];
	kp = psotu[otuk];
	ip->isop = jp;
	jp->isop = kp;
	kp->isop = ip;
	ip->length = bix;
	jp->length = bjx;
	kp->length = bkx;
	ip->kinp->length = ip->length;
	jp->kinp->length = jp->length;
	kp->kinp->length = kp->length;

	tr->rootp = kp;

	free_dvector(r);
	free_dmatrix(distan);
	free((Node *) psotu);
}

/******************************************************************************/
/* find best assignment of rate categories                                    */
/******************************************************************************/

/* find best assignment of rate categories */
void findbestratecombination()
{
	int k, u;
	double bestvalue, fv2;
	dvector catprob;
	dmatrix cdl;

	cdl = Ctree->condlkl;
	catprob = new_dvector(numcats+1);
	fv2 = (1.0-fracinv)/(double) numcats;

	for (k = 0; k < Numptrn; k++) {
		/* zero rate */
		if (constpat[k] == TRUE)
			catprob[0] = fracinv*Freqtpm[(int) Seqpat[0][k]];
		else 
			catprob[0] = 0.0;
		/* non-zero-rates */
		for (u = 1; u < numcats+1; u++)
			catprob[u] = fv2*cdl[u-1][k];
		/* find best */
		bestvalue = catprob[0];
		bestrate[k] = 0;
		for (u = 1; u < numcats+1; u++)
			if (catprob[u] >= bestvalue) {
				bestvalue = catprob[u];
				bestrate[k] = u;
			}
	}
	free_dvector(catprob);
	bestratefound = 1;
}

/* print best assignment of rate categories */
void printbestratecombination(FILE *fp)
{
	int s, k;

	for (s = 0; s < Maxsite; s++) {
		k = Alias[s];
		fprintf(fp, "%2d", bestrate[k]);
		if ((s+1) % 30 == 0)
			fprintf(fp, "\n");
		else if ((s+1) % 10 == 0)
			fprintf(fp, "  ");
	}
	if (s % 70 != 0)
		fprintf(fp, "\n");
}


/******************************************************************************/
/* computation of clocklike branch lengths                                    */
/******************************************************************************/

/* checks wether e is a valid edge specification */
int checkedge(int e)
{
	/* there are Numspc external branches:
	     0 - Numspc-1
	   there are Numibrnch internal branches:
	     Numspc - Numspc+Numibrnch-1
	*/
	
	if (e < 0) return FALSE;
	if (e < Numspc+Numibrnch) return TRUE;
	else return FALSE;  
}

/* print topology of subtree */
void fputsubstree(FILE *fp, Node *ip)
{
	Node *cp;
	
	if (ip->isop == NULL) { /* terminal nodes */
		numtc += fputid(fp, ip->number);
	} else {	
		cp = ip;
		fprintf(fp, "(");
		numtc += 1;
		do {		
			cp = cp->isop->kinp;
			if (cp->isop == NULL) { /* external node */
				numtc += fputid(fp, cp->number);
				fprintf(fp, ":%.5f", (cp->lengthc)*0.01);
				numtc += 7;
				cp = cp->kinp;
			} else { /* internal node */
				if (cp->height > 0.0) {
					fprintf(fp, "(");
					numtc += 1;
				} else if (cp->height < 0.0) {
					fprintf(fp, ")");
					numtc += 1;
					if (numtc > 60) {
						fprintf(fp, "\n");
						numtc = 1;
					}
					/* internal label */
					if (cp->kinp->label != NULL) {
						fprintf(fp, "%s", cp->kinp->label);
						numtc += strlen(cp->kinp->label);
					}
					if (numtc > 60) {
						fprintf(fp, "\n");
						numtc = 1;
					}
					fprintf(fp, ":%.5f", (cp->lengthc)*0.01);
					numtc += 6;
					if (numtc > 60) {
						fprintf(fp, "\n");
						numtc = 1;
					}
				}
			}
			if (cp->height <= 0.0 && cp->isop->height <= 0.0 &&
				cp->isop != ip) {
				putc(',', fp); /* not last subtree */
				numtc += 1;
				if (numtc > 60) {
					fprintf(fp, "\n");
					numtc = 1;
				}
			}
		} while (cp->isop != ip);
		fprintf(fp, ")");
		numtc += 1;
	}
	if (numtc > 60) {
		fprintf(fp, "\n");
		numtc = 1;
	}

}

/* print rooted tree file  */
void fputrooted(FILE *fp, int e)
{
	Node *rootbr;
	
	/* to be called only after clocklike branch
	   lengths have been computed */
	 
	 /* pointer to root branch */
	if (e < Numspc) rootbr = Ctree->ebrnchp[e];
	else rootbr = Ctree->ibrnchp[e - Numspc];

	fprintf(fp, "(");
	numtc = 2;
	fputsubstree(fp, rootbr);
	/* internal label */
	if (rootbr->label != NULL) {
		fprintf(fp, "%s", rootbr->label);
		numtc += strlen(rootbr->label);
	}
	if (numtc > 60) {
		fprintf(fp, "\n");
		numtc = 1;
	}
	fprintf(fp, ":%.5f,", (hroot - rootbr->height)*0.01);
	numtc += 7;
	if (numtc > 60) {
		fprintf(fp, "\n");
		numtc = 1;
	}
	fputsubstree(fp, rootbr->kinp);
	/* internal label */
	if (rootbr->kinp->label != NULL) {
		fprintf(fp, "%s", rootbr->kinp->label);
		numtc += strlen(rootbr->kinp->label);
	}
	if (numtc > 60) {
		fprintf(fp, "\n");
		numtc = 1;
	}
	fprintf(fp, ":%.5f);\n", (hroot - rootbr->kinp->height)*0.01);
}

/* finds heights in subtree */
void findheights(Node *ip)
{
	Node *cp, *rp;
	
	if (ip->isop != NULL) { /* forget terminal nodes */
		
		cp = ip;
		
		/* initialise node  */
		cp->height = 1.0; /* up */
		rp = cp;
		while (rp->isop != cp) {
			rp = rp->isop;
			rp->height = -1.0; /* down */
		}

		do {		
			cp = cp->isop->kinp;
			if (cp->isop == NULL) { /* external node */
				cp = cp->kinp;
			} else { /* internal node */
				if (cp->height == 0.0) { /* node not yet visited */
					cp->height = 1.0; /* up */
					rp = cp;
					while (rp->isop != cp) {
						rp = rp->isop;
						rp->height = -1.0; /* down */
					}
				} else if (cp->kinp->height == 1.0) {
					/* cp->kinp is next height pointer  */
					heights[Numhts] = cp->kinp;
					Numhts++;
				}
			}
		} while (cp->isop != ip);
		/* ip is last height pointer */
		heights[Numhts] = ip;
		Numhts++;
	}	
}


/* initialise clocklike branch lengths (with root on edge e) */
void initclock(int e)
{
	int n, h, count;
	Node *cp, *rp;
	double sum, minh, aveh, len;

	/* be sure to have a Ctree in memory and
	   to have pairwise distances computed */

	clockmode = 1; /* clocklike branch lengths */

	/* least square estimate for branch length */
	changedistan(Distanmat, Distanvec, Numspc);
	lslength(Ctree, Distanvec, Numspc, Numibrnch, Brnlength);

	/* pointer to root branch */
	if (e < Numspc) rootbr = Ctree->ebrnchp[e];
	else rootbr = Ctree->ibrnchp[e - Numspc];
	   
	/* clear all heights */
	for (n = 0; n < Numspc; n++) {
		Ctree->ebrnchp[n]->height = 0.0;
		Ctree->ebrnchp[n]->kinp->height = 0.0;
		Ctree->ebrnchp[n]->varheight = 0.0;
		Ctree->ebrnchp[n]->kinp->varheight = 0.0;
		if (n < Numibrnch) {
			Ctree->ibrnchp[n]->height = 0.0;
			Ctree->ibrnchp[n]->kinp->height = 0.0;
			Ctree->ibrnchp[n]->varheight = 0.0;
			Ctree->ibrnchp[n]->kinp->varheight = 0.0;
		}
	}
	
	/* collect pointers to height nodes */
	Numhts = 0;
	findheights(rootbr); /* one side */
	findheights(rootbr->kinp); /* other side */

	/* assign preliminary approximate heights and
	   corresponding branch lengths */ 
	for (h = 0; h < Numhts; h++) {
		
		cp = rp = heights[h];
		sum = 0;
		count = 0;
		minh = 0.0;
		while (rp->isop != cp) {
			count++;
			rp = rp->isop;
			sum += rp->lengthc + rp->kinp->height;
			if (rp->kinp->height > minh) minh = rp->kinp->height;
		}
		aveh = sum / (double) count;
		if (aveh < minh + MINARC) aveh = minh + MINARC;
		cp->height = aveh;
		rp = cp;
		while (rp->isop != cp) {
			rp = rp->isop;
			len = aveh - rp->kinp->height;
			rp->kinp->lengthc = len;
			rp->lengthc = len;
		}
		
	}
	if (rootbr->height > rootbr->kinp->height) minh = rootbr->height;
	else minh = rootbr->kinp->height;
	aveh = 0.5*(rootbr->lengthc + rootbr->height + rootbr->kinp->height);
	if (aveh < minh + MINARC) aveh = minh + MINARC;	
	hroot = aveh;
	maxhroot = RMHROOT*hroot; /* maximal possible hroot */
	len = (hroot - rootbr->height) + (hroot - rootbr->kinp->height);
	rootbr->lengthc = len;
	rootbr->kinp->lengthc = len;
}

/* approximate likelihood under the constaining assumption of
   clocklike branch lengths (with root on edge e) */
double clock_alklhd(int e)
{
	initclock(e);
	Ctree->lklhdc = treelkl(Ctree);

	return Ctree->lklhdc;	
}

/* log-likelihood given height ht at node pointed to by chep */
double heightlkl(double ht)
{
	Node *rp;
	double len;	
	
	/* adjust branch lengths */
	chep->height = ht;
	/* descendent branches */
	rp = chep;
	while (rp->isop != chep) {
		rp = rp->isop;
		len = chep->height - rp->kinp->height;
		rp->kinp->lengthc = len;
		rp->lengthc = len;
	}
	/* upward branch */
	if (chep == rootbr || chep->kinp == rootbr) {
		len = (hroot - chep->height) + (hroot - chep->kinp->height);
		chep->lengthc = len;
		chep->kinp->lengthc = len;
	} else {
		rp = chep->kinp;
		while (rp->isop->height <= 0.0)
			rp = rp->isop;
		chep->lengthc = rp->isop->height - chep->height;
		chep->kinp->lengthc = rp->isop->height - chep->height;
	}

	/* compute likelihood */
	Ctree->lklhdc = treelkl(Ctree);

	return -(Ctree->lklhdc); /* we use a minimizing procedure */
}

/* optimize current height */
void optheight(void)
{
	double he, fx, f2x, minh, maxh, len;
	Node *rp;

	/* current height */
	he = chep->height;
	
	/* minimum */
	minh = 0.0;
	rp = chep;
	while (rp->isop != chep) {
		rp = rp->isop;
		if (rp->kinp->height > minh)
			minh = rp->kinp->height;
	}
	minh += MINARC;
	
	/* maximum */
	if (chep == rootbr || chep->kinp == rootbr) {
		maxh = hroot;
	} else {
		rp = chep->kinp;
		while (rp->isop->height <= 0.0)
			rp = rp->isop;
		maxh = rp->isop->height;
	}
	maxh -= MINARC;

	/* check borders for height */
	if (he < minh) he = minh;
	if (he > maxh) he = maxh;

	/* optimization */
        if (!(he == minh && he == maxh))
		he = onedimenmin(minh, he, maxh, heightlkl, HEPSILON, &fx, &f2x);
	
	/* variance of height */
	f2x = fabs(f2x);
	if (1.0/(maxhroot*maxhroot) < f2x)
		chep->varheight = 1.0/f2x;
	else
		chep->varheight = maxhroot*maxhroot;
	
	/* adjust branch lengths */
	chep->height = he;
	/* descendent branches */
	rp = chep;
	while (rp->isop != chep) {
		rp = rp->isop;
		len = chep->height - rp->kinp->height;
		rp->kinp->lengthc = len;
		rp->lengthc = len;
	}
	/* upward branch */
	if (chep == rootbr || chep->kinp == rootbr) {
		len = (hroot - chep->height) + (hroot - chep->kinp->height);
		chep->lengthc = len;
		chep->kinp->lengthc = len;
	} else {
		rp = chep->kinp;
		while (rp->isop->height <= 0.0)
			rp = rp->isop;
		chep->lengthc = rp->isop->height - chep->height;
		chep->kinp->lengthc = rp->isop->height - chep->height;
	}
}

/* log-likelihood given height ht at root */
double rheightlkl(double ht)
{
	double len;	
	
	/* adjust branch lengths */
	hroot = ht;
	len = (hroot - rootbr->height) + (hroot - rootbr->kinp->height);
	rootbr->lengthc = len;
	rootbr->kinp->lengthc = len;

	/* compute likelihood */
	Ctree->lklhdc = treelkl(Ctree);

	return -(Ctree->lklhdc); /* we use a minimizing procedure */
}

/* optimize height of root */
void optrheight(void)
{
	double he, fx, f2x, minh, len;

	/* current height */
	he = hroot;
	
	/* minimum */
	if (rootbr->height > rootbr->kinp->height)
		minh = rootbr->height;
	else
		minh = rootbr->kinp->height;
	minh += MINARC;
		
	/* check borders for height */
	if (he < minh) he = minh;
	if (he > maxhroot) he = maxhroot;

	/* optimization */
	he = onedimenmin(minh, he, maxhroot, rheightlkl, HEPSILON, &fx, &f2x);
	
	/* variance of height of root */
	f2x = fabs(f2x);
	if (1.0/(maxhroot*maxhroot) < f2x)
		varhroot = 1.0/f2x;
	else
		varhroot = maxhroot*maxhroot;

	/* adjust branch lengths */
	hroot = he;
	len = (hroot - rootbr->height) + (hroot - rootbr->kinp->height);
	rootbr->lengthc = len;
	rootbr->kinp->lengthc = len;
}

/* exact likelihood under the constaining assumption of
   clocklike branch lengths (with root on edge e) */
double clock_lklhd(int e)
{
	int h, nconv;
	double old;
	
	Numitc = 0;
	Convergc = FALSE;
	
	initclock(e);
	
	do {
	
		Numitc++;
		nconv = 0;

		/* optimize height of root */
		old = hroot;
		optrheight();
		if (fabs(old - hroot) < HEPSILON) nconv++;

		/* optimize height of nodes */
		for (h = Numhts-1; h >= 0; h--) {
		
			/* pointer chep to current height node */
			chep = heights[h];
			
			/* store old value */
			old = chep->height;

			/* find better height */
			optheight();
			
			/* converged ? */
			if (fabs(old - chep->height) < HEPSILON) nconv++;
		}

		if (nconv == Numhts+1) Convergc = TRUE;
			
	} while (Numitc < MAXIT && !Convergc);
		
	/* compute final likelihood */
	Ctree->lklhdc = treelkl(Ctree);

	return Ctree->lklhdc;
}

/* find out the edge containing the root */
int findrootedge()
{
	int e, ebest;
	double logbest, logtest;
	
	/* compute the likelihood for all edges and take the edge with
	   best likelihood (using approximate ML) */

	ebest = 0;
	logbest = clock_alklhd(0);
	numbestroot = 1;
	for (e = 1; e < Numspc+Numibrnch; e++) {
		logtest = clock_alklhd(e);
		if (logtest > logbest) {
			ebest = e;
			logbest = logtest;
			numbestroot = 1;
		} else if (logtest == logbest) {
			numbestroot++;
		}
	}

	return ebest;
}

/* show heights and corresponding standard errors */
void resultheights(FILE *fp)
{
	int h, num;
	Node *cp;
	
	fprintf(fp, " height    S.E.    of node common to branches\n");
	for (h = 0; h < Numhts; h++) {
		fprintf(fp, "%.5f  %.5f    ", (heights[h]->height)*0.01,
			sqrt(heights[h]->varheight)*0.01);
		cp = heights[h];
		do {
			num = (cp->number) + 1;
			if (cp->kinp->isop != NULL) num += Numspc; /* internal branch */			
			fprintf(fp, "%d  ", num);
			cp = cp->isop;
		} while (cp != heights[h]);	
		fprintf(fp, "\n");
		
	}
	fprintf(fp, "%.5f  %.5f   of root at branch %d\n",
		hroot*0.01, sqrt(varhroot)*0.01, locroot+1);
}

