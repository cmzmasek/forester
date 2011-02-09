/*
 * puzzle1.c
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


#define EXTERN

#include "puzzle.h"
#include "gamma.h"

void num2quart(uli qnum, int *a, int *b, int *c, int *d)
{
        double temp;
        uli aa, bb, cc, dd;
        uli lowval=0, highval=0;

        aa=0; bb=1; cc=2; dd=3;

        temp = (double)(24 * qnum);
        temp = sqrt(temp);
        temp = sqrt(temp);
        /* temp = pow(temp, (double)(1/4)); */
        dd = (uli) floor(temp) + 1;
        if (dd < 3) dd = 3;
        lowval =  (uli) dd*(dd-1)*(dd-2)*(dd-3)/24;
        highval = (uli) (dd+1)*dd*(dd-1)*(dd-2)/24;
        if (lowval >= qnum)
            while ((lowval > qnum)) {
                dd -= 1; lowval = (uli) dd*(dd-1)*(dd-2)*(dd-3)/24;
            }
        else {
            while (highval <= qnum) {
                dd += 1; highval = (uli) (dd+1)*dd*(dd-1)*(dd-2)/24;
            }
            lowval = (uli) dd*(dd-1)*(dd-2)*(dd-3)/24;
        }
        qnum -= lowval;
        if (qnum > 0) {
            temp = (double)(6 * qnum);
            temp = pow(temp, (double)(1/3));
            cc = (uli) floor(temp);
            if (cc < 2) cc= 2;
            lowval =  (uli) cc*(cc-1)*(cc-2)/6;
            highval = (uli) (cc+1)*cc*(cc-1)/6;
            if (lowval >= qnum)
                while ((lowval > qnum)) {
                   cc -= 1; lowval = (uli) cc*(cc-1)*(cc-2)/6;
                }
            else {
                while (highval <= qnum) {
                   cc += 1; highval = (uli) (cc+1)*cc*(cc-1)/6;
                }
                lowval = (uli) cc*(cc-1)*(cc-2)/6;
            }
            qnum -= lowval;
            if (qnum > 0) {
                temp = (double)(2 * qnum);
                temp = sqrt(temp);
                bb = (uli) floor(temp);
                if (bb < 1) bb= 1;
                lowval =  (uli) bb*(bb-1)/2;
                highval = (uli) (bb+1)*bb/2;
                if (lowval >= qnum)
                    while ((lowval > qnum)) {
                        bb -= 1; lowval = (uli) bb*(bb-1)/2;
                    }
                else {
                    while (highval <= qnum) {
                       bb += 1; highval = (uli) (bb+1)*bb/2;
                    }
                    lowval = (uli) bb*(bb-1)/2;
                }
                qnum -= lowval;
                if (qnum > 0) {
                   aa = (uli) qnum;
                   if (aa < 0) aa= 0;
                }
            }
        }
        *d = (int)dd;
        *c = (int)cc;
        *b = (int)bb;
        *a = (int)aa;
}  /* num2quart */

/******************/

uli numquarts(int maxspc)
{
   uli tmp;
   int a, b, c, d;

   if (maxspc < 4)
     return (uli)0;
   else {
      maxspc--;
      a = maxspc-3;
      b = maxspc-2;
      c = maxspc-1;
      d = maxspc;

      tmp = (uli) 1 + a +
            (uli) b * (b-1) / 2 +
            (uli) c * (c-1) * (c-2) / 6 +
            (uli) d * (d-1) * (d-2) * (d-3) / 24;
      return (tmp);
   }
}  /* numquarts */

/******************/

uli quart2num (int a, int b, int c, int d)
{
      uli tmp;
      if ((a>b) || (b>c) || (c>d)) {
         fprintf(stderr, "Error PP5 not (%d <= %d <= %d <= %d) !!!\n", a, b, c,
d);
         exit (1);
      }
      tmp = (uli) a +
            (uli) b * (b-1) / 2 +
            (uli) c * (c-1) * (c-2) / 6 +
            (uli) d * (d-1) * (d-2) * (d-3) / 24;
      return (tmp);
}  /* quart2num */

/******************/



/*  flag=0	old allquart binary  */
/*  flag=1	allquart binary  */
/*  flag=2	allquart ACSII   */
/*  flag=3	quartlh  binary  */
/*  flag=4	quartlh  ASCII   */

void writetpqfheader(int            nspec,
                     FILE          *ofp,
                     int            flag)
{ int            currspec;

  if (flag == 0) {
     unsigned long  nquart;
     unsigned long  blocklen;

     nquart = numquarts(nspec);
     /* compute number of bytes */
     if (nquart % 2 == 0) { /* even number */
        blocklen = (nquart)/2;
     } else { /* odd number */
        blocklen = (nquart + 1)/2;
     }
     /* FPRINTF(STDOUTFILE "Writing quartet file: %s\n", filename); */
     fprintf(ofp, "TREE-PUZZLE\n%s\n\n", VERSION);
     fprintf(ofp, "species: %d\n",       nspec);
     fprintf(ofp, "quartets: %lu\n",     nquart);
     fprintf(ofp, "bytes: %lu\n\n",      blocklen);


     /* fwrite(&(quartetinfo[0]), sizeof(char), blocklen, ofp); */
  }

  if (flag == 1) fprintf(ofp, "##TPQF-BB (TREE-PUZZLE %s)\n%d\n", VERSION, nspec);
  if (flag == 2) fprintf(ofp, "##TPQF-BA (TREE-PUZZLE %s)\n%d\n", VERSION, nspec);
  if (flag == 3) fprintf(ofp, "##TPQF-LB (TREE-PUZZLE %s)\n%d\n", VERSION, nspec);
  if (flag == 4) fprintf(ofp, "##TPQF-LA (TREE-PUZZLE %s)\n%d\n", VERSION, nspec);

  for (currspec=0; currspec<nspec; currspec++) {
     fputid(ofp, currspec);
     fprintf(ofp, "\n");
  } /* for each species */
  fprintf(ofp, "\n");

} /* writetpqfheader */



void writeallquarts(int            nspec,
                    char          *filename,
                    unsigned char *quartetinfo)
{ unsigned long  nquart;
  unsigned long  blocklen;
  FILE          *ofp;

  nquart = numquarts(nspec);

  /* compute number of bytes */
  if (nquart % 2 == 0) { /* even number */
     blocklen = (nquart)/2;
  } else { /* odd number */
     blocklen = (nquart + 1)/2;
  }

  FPRINTF(STDOUTFILE "Writing quartet file: %s\n", filename);

  openfiletowrite(&ofp, filename, "all quartets");

  writetpqfheader(nspec, ofp, 0);

  fwrite(&(quartetinfo[0]), sizeof(char), blocklen, ofp);
  fclose(ofp);

} /* writeallquart */



void readallquarts(int            nspec,
                   char          *filename,
                   unsigned char *quartetinfo)
{ unsigned long  nquart;
  unsigned long  blocklen;
  int            currspec;
  unsigned long  dummynquart;
  unsigned long  dummyblocklen;
  int            dummynspec;
  char           dummyversion[128];
  char           dummyname[128];
  FILE          *ifp;

  nquart = numquarts(nspec);

  /* compute number of bytes */
  if (nquart % 2 == 0) { /* even number */
     blocklen = (nquart)/2;
  } else { /* odd number */
     blocklen = (nquart + 1)/2;
  }

/*  &(quartetinfo[0] */

  openfiletoread(&ifp, filename, "all quartets");

  fscanf(ifp, "TREE-PUZZLE\n");
  fscanf(ifp, "%s\n\n", dummyversion);
  
  fscanf(ifp, "species: %d\n",   &dummynspec);
  fscanf(ifp, "quartets: %lu\n", &dummynquart);
  fscanf(ifp, "bytes: %lu\n\n",  &dummyblocklen);

  if ((nspec != dummynspec) ||
      (nquart != dummynquart) ||
      (blocklen != dummyblocklen)) {
     FPRINTF(STDOUTFILE "ERROR: Parameters in quartet file %s do not match!!!\n",
                     filename);
#    if PARALLEL
         PP_SendDone();
         MPI_Finalize();
#    endif /* PARALLEL */
     exit(1);
  }

  FPRINTF(STDOUTFILE "Reading quartet file: %s\n", filename);
  FPRINTF(STDOUTFILE "   generated by: TREE-PUZZLE %s\n", dummyversion);
  FPRINTF(STDOUTFILE "   current:      TREE-PUZZLE %s\n", VERSION);

  FPRINTF(STDOUTFILE "   %d species, %lu quartets, %lu bytes\n", 
                   dummynspec, dummynquart, dummyblocklen);

  for (currspec=0; currspec<nspec; currspec++) {

     fscanf(ifp, "%s\n", dummyname);
     FPRINTF(STDOUTFILE "   %3d. %s (", currspec+1, dummyname);
     fputid(STDOUT, currspec);
     FPRINTF(STDOUTFILE ")\n");

  } /* for each species */
  FPRINTF(STDOUTFILE "\n");

  fread(&(quartetinfo[0]), sizeof(char), blocklen, ifp);
  fclose(ifp);

} /* writeallquart */




/******************************************************************************/
/* options - file I/O - output                                                */
/******************************************************************************/

/* compute TN parameters according to F84 Ts/Tv ratio */
void makeF84model()
{
	double rho, piA, piC, piG, piT, piR, piY, ts, yr;
	
	piA = Freqtpm[0];
	piC = Freqtpm[1];
	piG = Freqtpm[2];
	piT = Freqtpm[3];
	piR = piA + piG;
	piY = piC + piT;
	if (piC*piT*piR + piA*piG*piY == 0.0) {		
		FPRINTF(STDOUTFILE "\n\n\nF84 model not possible ");
		FPRINTF(STDOUTFILE "(bad combination of base frequencies)\n");
		tstvf84 = 0.0;
		return;
	}
	rho = (piR*piY*(piR*piY*tstvf84 - (piA*piG + piC*piT)))/
		(piC*piT*piR + piA*piG*piY);
	
	if(piR == 0.0 || piY == 0.0 || (piR + rho) == 0.0) {
		FPRINTF(STDOUTFILE "\n\n\nF84 model not possible ");
		FPRINTF(STDOUTFILE "(bad combination of base frequencies)\n");
		tstvf84 = 0.0;
		return;
	}
	ts = 0.5 + 0.25*rho*(1.0/piR + 1.0/piY);
	yr = (piY + rho)/piY * piR/(piR + rho);
	if (ts < MINTS || ts > MAXTS) {
		FPRINTF(STDOUTFILE "\n\n\nF84 model not possible ");
		FPRINTF(STDOUTFILE "(bad Ts/Tv parameter)\n");
		tstvf84 = 0.0;
		return;
	}
	if (yr < MINYR || yr > MAXYR) {
		FPRINTF(STDOUTFILE "\n\n\nF84 model not possible ");
		FPRINTF(STDOUTFILE "(bad Y/R transition parameter)\n");
		tstvf84 = 0.0;
		return;
	}
	TSparam = ts;
	YRparam = yr;
	optim_optn = FALSE;
}

/* compute number of quartets used in LM analysis */
void compnumqts()
{
	if (lmqts == 0) {
		if (numclust == 4)
			Numquartets = (uli) clustA*clustB*clustC*clustD;
		if (numclust == 3)
			Numquartets = (uli) clustA*clustB*clustC*(clustC-1)/2;
		if (numclust == 2)
			Numquartets = (uli) clustA*(clustA-1)/2 * clustB*(clustB-1)/2;
		if (numclust == 1)	
			Numquartets = (uli) Maxspc*(Maxspc-1)*(Maxspc-2)*(Maxspc-3)/24;
	} else {
		Numquartets = lmqts;
	}
}

/* set options interactively */
void setoptions()
{	
	int i, valid;
	double sumfreq;
	char ch;

    puzzlemode = PAIRDIST; /*Only do pairwise dist. CZ, 05/16/01*/

	/* defaults */
	rhetmode = UNIFORMRATE;          /* assume rate homogeneity               */
	numcats = 1;
	Geta = 0.05;
	grate_optim = FALSE;
	fracinv = 0.0;
	fracinv_optim = FALSE;

	compclock = FALSE;           /* compute clocklike branch lengths      */
	locroot = -1;                /* search for optimal place of root      */ 
	qcalg_optn = FALSE;          /* don't use sampling of quartets        */
	approxp_optn = TRUE;         /* approximate parameter estimates       */
	listqptrees = PSTOUT_NONE;   /* list puzzling step trees              */
	
	/* approximate QP quartets? */
	if (Maxspc <= 6) approxqp = FALSE;
	else approxqp = TRUE;
	
	codon_optn = 0;        /* use all positions in a codon          */
	
	/* number of puzzling steps */
	if (Maxspc <= 25) Numtrial = 1000;
	else if (Maxspc <= 50) Numtrial = 10000;
	else if (Maxspc <= 75) Numtrial = 25000;
	else Numtrial = 50000;
	    
	utree_optn = TRUE;     /* use first user tree for estimation    */
	outgroup = 0;          /* use first taxon as outgroup           */
	sym_optn = FALSE;      /* symmetrize doublet frequencies        */
	tstvf84 = 0.0;         /* disable F84 model                     */
	show_optn = FALSE;     /* show unresolved quartets              */
	typ_optn = TREERECON_OPTN;          /* tree reconstruction                   */
	numclust = 1;          /* one clusters in LM analysis           */
	lmqts = 0;             /* all quartets in LM analysis           */
	compnumqts();
	if (Numquartets > 10000) {
		lmqts = 10000;     /* 10000 quartets in LM analysis         */
		compnumqts();
	}
	
 	do {
		FPRINTF(STDOUTFILE "\n\n\nGENERAL OPTIONS\n");
		FPRINTF(STDOUTFILE " b                     Type of analysis?  ");
		if (typ_optn == TREERECON_OPTN) FPRINTF(STDOUTFILE "Tree reconstruction\n");
		if (typ_optn == LIKMAPING_OPTN) FPRINTF(STDOUTFILE "Likelihood mapping\n");
		if (typ_optn == TREERECON_OPTN) {
			FPRINTF(STDOUTFILE " k                Tree search procedure?  ");
			if (puzzlemode == QUARTPUZ) FPRINTF(STDOUTFILE "Quartet puzzling\n");
			if (puzzlemode == USERTREE) FPRINTF(STDOUTFILE "User defined trees\n");
			if (puzzlemode == PAIRDIST) FPRINTF(STDOUTFILE "Pairwise distances only (no tree)\n");
			if (puzzlemode == QUARTPUZ) {
				FPRINTF(STDOUTFILE " v       Approximate quartet likelihood?  %s\n",
					(approxqp ? "Yes" : "No"));
				FPRINTF(STDOUTFILE " u             List unresolved quartets?  %s\n",
					(show_optn ? "Yes" : "No"));
				FPRINTF(STDOUTFILE " n             Number of puzzling steps?  %lu\n",
						Numtrial);
				FPRINTF(STDOUTFILE " j             List puzzling step trees?  ");
				switch (listqptrees) {
					case PSTOUT_NONE:      FPRINTF(STDOUTFILE "No\n"); break;
					case PSTOUT_ORDER:     FPRINTF(STDOUTFILE "Unique topologies\n"); break;
					case PSTOUT_LISTORDER: FPRINTF(STDOUTFILE "Unique topologies & Chronological list\n"); break;
					case PSTOUT_LIST:      FPRINTF(STDOUTFILE "Chronological list only\n"); break;
				}

				FPRINTF(STDOUTFILE " o                  Display as outgroup?  ");
				fputid(STDOUT, outgroup);
				FPRINTF(STDOUTFILE "\n");
			}
			if (puzzlemode == QUARTPUZ || puzzlemode == USERTREE) {
				FPRINTF(STDOUTFILE " z     Compute clocklike branch lengths?  ");
				if (compclock) FPRINTF(STDOUTFILE "Yes\n");
				else FPRINTF(STDOUTFILE "No\n");
			}
			if (compclock)
				if (puzzlemode == QUARTPUZ || puzzlemode == USERTREE) {
					FPRINTF(STDOUTFILE " l                     Location of root?  ");
					if (locroot < 0) FPRINTF(STDOUTFILE "Best place (automatic search)\n");
					else if (locroot < Maxspc) {
						FPRINTF(STDOUTFILE "Branch %d (", locroot + 1);
						fputid(STDOUT, locroot);
						FPRINTF(STDOUTFILE ")\n");
					} else FPRINTF(STDOUTFILE "Branch %d (internal branch)\n", locroot + 1);
				}
		}
		if (typ_optn == LIKMAPING_OPTN) {
			  FPRINTF(STDOUTFILE " g          Group sequences in clusters?  ");
			  if (numclust == 1) FPRINTF(STDOUTFILE "No\n");
			  else FPRINTF(STDOUTFILE "Yes (%d clusters as specified)\n", numclust);
			  FPRINTF(STDOUTFILE " n                   Number of quartets?  ");
			  if (lmqts == 0) FPRINTF(STDOUTFILE "%lu (all possible)\n", Numquartets);
			  else FPRINTF(STDOUTFILE "%lu (random choice)\n", lmqts);
		}
		FPRINTF(STDOUTFILE " e                  Parameter estimates?  ");
		if (approxp_optn) FPRINTF(STDOUTFILE "Approximate (faster)\n");
		else FPRINTF(STDOUTFILE "Exact (slow)\n");
		if (!(puzzlemode == USERTREE && typ_optn == TREERECON_OPTN)) {	
			FPRINTF(STDOUTFILE " x            Parameter estimation uses?  ");
			if (qcalg_optn) FPRINTF(STDOUTFILE "Quartet sampling + NJ tree\n");
			else FPRINTF(STDOUTFILE "Neighbor-joining tree\n");
			
		} else {
			FPRINTF(STDOUTFILE " x            Parameter estimation uses?  ");
			if (utree_optn)
				FPRINTF(STDOUTFILE "1st input tree\n");
			else if (qcalg_optn) FPRINTF(STDOUTFILE "Quartet sampling + NJ tree\n");
			else FPRINTF(STDOUTFILE "Neighbor-joining tree\n");
		}
		FPRINTF(STDOUTFILE "SUBSTITUTION PROCESS\n");
		FPRINTF(STDOUTFILE " d          Type of sequence input data?  ");
		if (auto_datatype == AUTO_GUESS) FPRINTF(STDOUTFILE "Auto: ");
		if (data_optn == NUCLEOTIDE) FPRINTF(STDOUTFILE "Nucleotides\n");
		if (data_optn == AMINOACID) FPRINTF(STDOUTFILE "Amino acids\n");
		if (data_optn == BINARY) FPRINTF(STDOUTFILE "Binary states\n");
		if (data_optn == NUCLEOTIDE && (Maxseqc % 3) == 0  && !SH_optn) {
			FPRINTF(STDOUTFILE " h             Codon positions selected?  ");
			if (codon_optn == 0) FPRINTF(STDOUTFILE "Use all positions\n");
			if (codon_optn == 1) FPRINTF(STDOUTFILE "Use only 1st positions\n");
			if (codon_optn == 2) FPRINTF(STDOUTFILE "Use only 2nd positions\n");
			if (codon_optn == 3) FPRINTF(STDOUTFILE "Use only 3rd positions\n");
			if (codon_optn == 4) FPRINTF(STDOUTFILE "Use 1st and 2nd positions\n");
		}
		FPRINTF(STDOUTFILE " m                Model of substitution?  ");
		if (data_optn == NUCLEOTIDE) { /* nucleotides */
			if (nuc_optn) {
				if(HKY_optn)
					FPRINTF(STDOUTFILE "HKY (Hasegawa et al. 1985)\n");	
				else {
					FPRINTF(STDOUTFILE "TN (Tamura-Nei 1993)\n");
					FPRINTF(STDOUTFILE " p      Constrain TN model to F84 model?  ");
					if (tstvf84 == 0.0)
						FPRINTF(STDOUTFILE "No\n");
					else FPRINTF(STDOUTFILE "Yes (Ts/Tv ratio = %.2f)\n", tstvf84);
				}
				FPRINTF(STDOUTFILE " t    Transition/transversion parameter?  ");
				if (optim_optn)
					FPRINTF(STDOUTFILE "Estimate from data set\n");
				else
					FPRINTF(STDOUTFILE "%.2f\n", TSparam);	
				if (TN_optn) {
					FPRINTF(STDOUTFILE " r             Y/R transition parameter?  ");
					if (optim_optn)
						FPRINTF(STDOUTFILE "Estimate from data set\n");
					else
						FPRINTF(STDOUTFILE "%.2f\n", YRparam);		
 				}			
			}
			if (SH_optn) {
					FPRINTF(STDOUTFILE "SH (Schoeniger-von Haeseler 1994)\n");
					FPRINTF(STDOUTFILE " t    Transition/transversion parameter?  ");
					if (optim_optn)
						FPRINTF(STDOUTFILE "Estimate from data set\n");
					else
						FPRINTF(STDOUTFILE "%.2f\n", TSparam);
			}	
		}
		if (data_optn == NUCLEOTIDE && SH_optn) {
			FPRINTF(STDOUTFILE " h                  Doublets defined by?  ");
			if (SHcodon)
				FPRINTF(STDOUTFILE "1st and 2nd codon positions\n");
			else
				FPRINTF(STDOUTFILE "1st+2nd, 3rd+4th, etc. site\n");
		}
		if (data_optn == AMINOACID) { /* amino acids */
			switch (auto_aamodel) {
				case AUTO_GUESS:
					FPRINTF(STDOUTFILE "Auto: ");
					break;
				case AUTO_DEFAULT:
					FPRINTF(STDOUTFILE "Def.: ");
					break;
			}
			if (Dayhf_optn) FPRINTF(STDOUTFILE "Dayhoff (Dayhoff et al. 1978)\n");	
			if (Jtt_optn) FPRINTF(STDOUTFILE "JTT (Jones et al. 1992)\n");
			if (mtrev_optn) FPRINTF(STDOUTFILE "mtREV24 (Adachi-Hasegawa 1996)\n");
			if (cprev_optn) FPRINTF(STDOUTFILE "cpREV45 (Adachi et al. 2000)\n");
			if (blosum62_optn) FPRINTF(STDOUTFILE "BLOSUM62 (Henikoff-Henikoff 92)\n");
			if (vtmv_optn) FPRINTF(STDOUTFILE "VT (Mueller-Vingron 2000)\n");
			if (wag_optn) FPRINTF(STDOUTFILE "WAG (Whelan-Goldman 2000)\n");
		}
		if (data_optn == BINARY) { /* binary states */
			FPRINTF(STDOUTFILE "Two-state model (Felsenstein 1981)\n");
		}
		if (data_optn == AMINOACID)
			FPRINTF(STDOUTFILE " f               Amino acid frequencies?  ");
		else if (data_optn == NUCLEOTIDE && SH_optn)
			FPRINTF(STDOUTFILE " f                  Doublet frequencies?  ");
		else if (data_optn == NUCLEOTIDE && nuc_optn)
			FPRINTF(STDOUTFILE " f               Nucleotide frequencies?  ");
		else if (data_optn == BINARY)
			FPRINTF(STDOUTFILE " f             Binary state frequencies?  ");
		FPRINTF(STDOUTFILE "%s\n", (Frequ_optn ? "Estimate from data set" :
			"Use specified values"));
		if (data_optn == NUCLEOTIDE && SH_optn)
			FPRINTF(STDOUTFILE " s       Symmetrize doublet frequencies?  %s\n",
				(sym_optn ? "Yes" : "No"));
				
		FPRINTF(STDOUTFILE "RATE HETEROGENEITY\n");
		FPRINTF(STDOUTFILE " w          Model of rate heterogeneity?  ");
		if (rhetmode == UNIFORMRATE) FPRINTF(STDOUTFILE "Uniform rate\n");
		if (rhetmode == GAMMARATE  ) FPRINTF(STDOUTFILE "Gamma distributed rates\n");
		if (rhetmode == TWORATE    ) FPRINTF(STDOUTFILE "Two rates (1 invariable + 1 variable)\n");
		if (rhetmode == MIXEDRATE  ) FPRINTF(STDOUTFILE "Mixed (1 invariable + %d Gamma rates)\n", numcats);
		
		if (rhetmode == TWORATE || rhetmode == MIXEDRATE) {
			FPRINTF(STDOUTFILE " i         Fraction of invariable sites?  ");
			if (fracinv_optim) FPRINTF(STDOUTFILE "Estimate from data set");
			else FPRINTF(STDOUTFILE "%.2f", fracinv);
			if (fracinv == 0.0 && !fracinv_optim) FPRINTF(STDOUTFILE " (all sites variable)");
			FPRINTF(STDOUTFILE "\n");
		}
		if (rhetmode == GAMMARATE || rhetmode == MIXEDRATE) {
			FPRINTF(STDOUTFILE " a   Gamma distribution parameter alpha?  ");
			if (grate_optim)
				FPRINTF(STDOUTFILE "Estimate from data set\n");
			else if (Geta > 0.5)
				FPRINTF(STDOUTFILE "%.2f (strong rate heterogeneity)\n", (1.0-Geta)/Geta);
			else FPRINTF(STDOUTFILE "%.2f (weak rate heterogeneity)\n", (1.0-Geta)/Geta);
			FPRINTF(STDOUTFILE " c      Number of Gamma rate categories?  %d\n", numcats);
		}
	
		FPRINTF(STDOUTFILE "\nQuit [q], confirm [y], or change [menu] settings: ");
		
		/* read one char */
		ch = getchar();
		if (ch != '\n') {
			do ;
			while (getchar() != '\n');
		}
		ch = (char) tolower((int) ch);
		
		/* letters in use: d m */
		/* letters not in use:  */

		switch (ch) {

			case '\n':	break;
			
			

			case 'd':	if (auto_datatype == AUTO_GUESS) {
						auto_datatype = AUTO_OFF;
						guessdata_optn = data_optn;
						data_optn = 0;
					} else {
						data_optn = data_optn + 1;
						if (data_optn == 3) {
							auto_datatype = AUTO_GUESS;
							data_optn = guessdata_optn;
						}
					}
					/* translate characters into format used by ML engine */
	                    		translatedataset();
	                    		estimatebasefreqs();
					break;
						
			
	
			case 'm':	if (data_optn == NUCLEOTIDE) { /* nucleotide data */
							if(HKY_optn && nuc_optn) {
								/* HKY -> TN */
								tstvf84       = 0.0;
								TSparam       = 2.0;
								YRparam       = 0.9;
								HKY_optn      = FALSE;
								TN_optn       = TRUE;
								optim_optn    = TRUE;
								nuc_optn      = TRUE;
								SH_optn       = FALSE;
								break;
							}
							if(TN_optn && nuc_optn) {
								if (Maxseqc % 2 == 0 || Maxseqc % 3 == 0) {
									/* number of chars needs to be a multiple 2 or 3 */
									/* TN -> SH */		
									if (Maxseqc % 2 != 0 && Maxseqc % 3 == 0)
										SHcodon = TRUE;
									else
										SHcodon = FALSE;								
									tstvf84       = 0.0;
									TSparam       = 2.0;
									YRparam       = 1.0;
									HKY_optn      = TRUE;
									TN_optn       = FALSE;
									optim_optn    = TRUE;
									nuc_optn      = FALSE;
									SH_optn       = TRUE;
									/* translate characters into format */
									/* used by ML engine */
									translatedataset();
									estimatebasefreqs();
								} else {
									FPRINTF(STDOUTFILE "\n\n\nSH model not ");
									FPRINTF(STDOUTFILE "available for the data set!\n");
									/* TN -> HKY */
									tstvf84       = 0.0;
									TSparam       = 2.0;
									YRparam       = 1.0;
									HKY_optn      = TRUE;
									TN_optn       = FALSE;
									optim_optn    = TRUE;
									nuc_optn      = TRUE;
									SH_optn       = FALSE;
								}
								break;
							}
							if(SH_optn) {
								/* SH -> HKY */
								tstvf84       = 0.0;
								TSparam       = 2.0;
								YRparam       = 1.0;
								HKY_optn      = TRUE;
								TN_optn       = FALSE;
								optim_optn    = TRUE;
								nuc_optn      = TRUE;
								SH_optn       = FALSE;
								/* translate characters into format */
								/* used by ML engine */
								translatedataset();
								estimatebasefreqs();
								break;
							}
							break;
						}
						if (data_optn == AMINOACID) { /* amino acid data */
							if (auto_aamodel) {
								/* AUTO -> Dayhoff */
								Dayhf_optn    = TRUE;
								Jtt_optn      = FALSE;
								mtrev_optn    = FALSE;
								cprev_optn    = FALSE;
								blosum62_optn = FALSE;
								vtmv_optn     = FALSE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
							if (Dayhf_optn) {
								/* Dayhoff -> JTT */
								Dayhf_optn    = FALSE;
								Jtt_optn      = TRUE;
								mtrev_optn    = FALSE;
								cprev_optn    = FALSE;
								blosum62_optn = FALSE;
								vtmv_optn     = FALSE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
							if (Jtt_optn) {
								/* JTT -> mtREV */
								Dayhf_optn    = FALSE;
								Jtt_optn      = FALSE;
								mtrev_optn    = TRUE;
								cprev_optn    = FALSE;
								blosum62_optn = FALSE;
								vtmv_optn     = FALSE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
#ifdef CPREV
							if (mtrev_optn) {
								/* mtREV -> cpREV */
								Dayhf_optn    = FALSE;
								Jtt_optn      = FALSE;
								mtrev_optn    = FALSE;
								cprev_optn    = TRUE;
								blosum62_optn = FALSE;
								vtmv_optn     = FALSE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
#else /* ! CPREV */
							if (mtrev_optn) {
								/* mtREV -> BLOSUM 62 */
								Dayhf_optn    = FALSE;
								Jtt_optn      = FALSE;
								mtrev_optn    = FALSE;
								cprev_optn    = FALSE;
								blosum62_optn = TRUE;
								vtmv_optn     = FALSE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
#endif /* ! CPREV */

#ifdef CPREV
							if (cprev_optn) {
								/* cpREV -> BLOSUM 62 */
								Dayhf_optn    = FALSE;
								Jtt_optn      = FALSE;
								mtrev_optn    = FALSE;
								cprev_optn    = FALSE;
								blosum62_optn = TRUE;
								vtmv_optn     = FALSE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
#endif
							if (blosum62_optn) {
								/* BLOSUM 62 -> VT model */
								Dayhf_optn    = FALSE;
								Jtt_optn      = FALSE;
								mtrev_optn    = FALSE;
								cprev_optn    = FALSE;
								blosum62_optn = FALSE;
								vtmv_optn     = TRUE;
								wag_optn      = FALSE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
							if (vtmv_optn) {
								/* VT model -> WAG model */
								Dayhf_optn    = FALSE;
								Jtt_optn      = FALSE;
								mtrev_optn    = FALSE;
								cprev_optn    = FALSE;
								blosum62_optn = FALSE;
								vtmv_optn     = FALSE;
								wag_optn      = TRUE;
								auto_aamodel  = AUTO_OFF;
								break;
							}
							if (wag_optn) {
								/* WAG model -> AUTO */
								Dayhf_optn    = guessDayhf_optn;
								Jtt_optn      = guessJtt_optn;
								mtrev_optn    = guessmtrev_optn;
								cprev_optn    = guesscprev_optn;
								blosum62_optn = guessblosum62_optn;
								vtmv_optn     = guessvtmv_optn;
								wag_optn      = guesswag_optn;
								auto_aamodel  = guessauto_aamodel;
								break;
							}
							break;
						}
						if (data_optn == BINARY) {
							FPRINTF(STDOUTFILE "\n\n\nNo other model available!\n");
						}
						break;
						
			

			case 'y':	break;

			default:	FPRINTF(STDOUTFILE "\n\n\nThis is not a possible option!\n");
						break;
		}
	} while (ch != 'y');

	FPRINTF(STDOUTFILE "\n\n\n");
}

/* open file for reading */
void openfiletoread(FILE **fp, char name[], char descr[])
{
	int count = 0;
	cvector str;

	if ((*fp = fopen(name, "r")) == NULL) {
		FPRINTF(STDOUTFILE "\n\n\nPlease enter a file name for the %s: ", descr);
		str = mygets();
		while ((*fp = fopen(str, "r")) == NULL)
		{
			count++;
			if (count > 10)
			{
				FPRINTF(STDOUTFILE "\n\n\nToo many trials - quitting ...\n");
				exit(1);
			}
   			FPRINTF(STDOUTFILE "File '%s' not found, ", str);
			FPRINTF(STDOUTFILE "please enter alternative name: ");
			free_cvector(str);
			str = mygets();
		}
		free_cvector(str);
		FPRINTF(STDOUTFILE "\n");
	}
} /* openfiletoread */


/* open file for writing */
void openfiletowrite(FILE **fp, char name[], char descr[])
{
	int count = 0;
	cvector str;

	if ((*fp = fopen(name, "w")) == NULL) {
   		FPRINTF(STDOUTFILE "\n\n\nPlease enter a file name for the %s: ", descr);
		str = mygets();
		while ((*fp = fopen(str, "w")) == NULL)
		{
			count++;
			if (count > 10)
			{
				FPRINTF(STDOUTFILE "\n\n\nToo many trials - quitting ...\n");
				exit(1);
			}
   			FPRINTF(STDOUTFILE "File '%s' not created, ", str);
			FPRINTF(STDOUTFILE "please enter other name: ");
			free_cvector(str);
			str = mygets();
		}
		free_cvector(str);
		FPRINTF(STDOUTFILE "\n");
	}
} /* openfiletowrite */


/* open file for appending */
void openfiletoappend(FILE **fp, char name[], char descr[])
{
	int count = 0;
	cvector str;

	if ((*fp = fopen(name, "a")) == NULL) {
   		FPRINTF(STDOUTFILE "\n\n\nPlease enter a file name for the %s: ", descr);
		str = mygets();
		while ((*fp = fopen(str, "a")) == NULL)
		{
			count++;
			if (count > 10)
			{
				FPRINTF(STDOUTFILE "\n\n\nToo many trials - quitting ...\n");
				exit(1);
			}
   			FPRINTF(STDOUTFILE "File '%s' not created, ", str);
			FPRINTF(STDOUTFILE "please enter other name: ");
			free_cvector(str);
			str = mygets();
		}
		free_cvector(str);
		FPRINTF(STDOUTFILE "\n");
	}
} /* openfiletowrite */


/* close file */
void closefile(FILE *fp)
{	
	fclose(fp);
} /* closefile */

/* symmetrize doublet frequencies */
void symdoublets()
{
	int i, imean;
	double mean;
	
	if (data_optn == NUCLEOTIDE && SH_optn && sym_optn) {
		/* ML frequencies */
		mean = (Freqtpm[1] + Freqtpm[4])/2.0; /* AC CA */
		Freqtpm[1] = mean;
		Freqtpm[4] = mean;
		mean = (Freqtpm[2] + Freqtpm[8])/2.0; /* AG GA */
		Freqtpm[2] = mean;
		Freqtpm[8] = mean;
		mean = (Freqtpm[3] + Freqtpm[12])/2.0; /* AT TA */
		Freqtpm[3] = mean;
		Freqtpm[12] = mean;
		mean = (Freqtpm[6] + Freqtpm[9])/2.0; /* CG GC */
		Freqtpm[6] = mean;
		Freqtpm[9] = mean;
		mean = (Freqtpm[7] + Freqtpm[13])/2.0; /* CT TC */
		Freqtpm[7] = mean;
		Freqtpm[13] = mean;
		mean = (Freqtpm[11] + Freqtpm[14])/2.0; /* GT TG */
		Freqtpm[11] = mean;
		Freqtpm[14] = mean;
		
		/* base composition of each taxon */
		for (i = 0; i < Maxspc; i++) {
			imean = (Basecomp[i][1] + Basecomp[i][4])/2; /* AC CA */
			Basecomp[i][1] = imean;
			Basecomp[i][4] = imean;
			imean = (Basecomp[i][2] + Basecomp[i][8])/2; /* AG GA */
			Basecomp[i][2] = imean;
			Basecomp[i][8] = imean;
			imean = (Basecomp[i][3] + Basecomp[i][12])/2; /* AT TA */
			Basecomp[i][3] = imean;
			Basecomp[i][12] = imean;
			imean = (Basecomp[i][6] + Basecomp[i][9])/2; /* CG GC */
			Basecomp[i][6] = imean;
			Basecomp[i][9] = imean;
			imean = (Basecomp[i][7] + Basecomp[i][13])/2; /* CT TC */
			Basecomp[i][7] = imean;
			Basecomp[i][13] = imean;
			imean = (Basecomp[i][11] + Basecomp[i][14])/2; /* GT TG */
			Basecomp[i][11] = imean;
			Basecomp[i][14] = imean;
		}
	}
}

/* show Ts/Tv ratio and Ts Y/R ratio */
void computeexpectations()
{
	/* CZ */
}

/* write ML distance matrix to file. Modified CZ 05/29/01 */
void putdistance(FILE *fp)
{
	/*int i;*/
    int i, j;
	
    for (i = 0; i < Maxspc - 1; i++) {
		/*fprintf(fp, "%.5f  ", Distanmat[i]/100.0);*/
        for ( j = 0; j < 26; j++ ) {
            fputc( Identif[i][j], fp ); /*CZ*/
        }
        fprintf(fp, "%.5f\n", Distanmat[i]/100.0);
	}
	fprintf(fp, "\n");
	
}




/* first lines of EPSF likelihood mapping file */
void initps(FILE *ofp)
{
	/* CZ */
}

/* plot one point of likelihood mapping analysis */
void plotlmpoint(FILE *ofp, double w1, double w2)
{
	/* CZ */
}

/* last lines of EPSF likelihood mapping file */
void finishps(FILE *ofp)
{
	/* CZ */
}

/* computes LM point from the three log-likelihood values,
   plots the point, and does some statistics */
void makelmpoint(FILE *fp, double b1, double b2, double b3)
{
	double w1, w2, w3, temp;
	unsigned char qpbranching;
	double temp1, temp2, temp3, onethird;
	unsigned char discreteweight[3], treebits[3];
	
	onethird = 1.0/3.0;
	treebits[0] = (unsigned char) 1;
	treebits[1] = (unsigned char) 2;
	treebits[2] = (unsigned char) 4;

	/* sort in descending order */
	qweight[0] = b1;
	qweight[1] = b2;
	qweight[2] = b3;
	sort3doubles(qweight, qworder);

	/* compute Bayesian weights */
	qweight[qworder[1]] = exp(qweight[qworder[1]]-qweight[qworder[0]]);
	qweight[qworder[2]] = exp(qweight[qworder[2]]-qweight[qworder[0]]);
	qweight[qworder[0]] = 1.0;
	temp = qweight[0] + qweight[1] + qweight[2];
	qweight[0] = qweight[0]/temp;
	qweight[1] = qweight[1]/temp;
	qweight[2] = qweight[2]/temp;

	/* plot one point in likelihood mapping triangle */
	w1 = qweight[0];
	w2 = qweight[1];
	w3 = qweight[2];
	plotlmpoint(fp, w1, w2);
	
	/* check areas 1,2,3 */	
	if (treebits[qworder[0]] == 1) ar1++;
	else if (treebits[qworder[0]] == 2) ar2++;
	else ar3++;				

	/* check out regions 1,2,3,4,5,6,7 */

	/* 100 distribution */
	temp1 = 1.0 - qweight[qworder[0]];
	sqdiff[0] = temp1*temp1 +
		qweight[qworder[1]]*qweight[qworder[1]] +
		qweight[qworder[2]]*qweight[qworder[2]];
	discreteweight[0] = treebits[qworder[0]];

	/* 110 distribution */
	temp1 = 0.5 - qweight[qworder[0]];
	temp2 = 0.5 - qweight[qworder[1]];
	sqdiff[1] = temp1*temp1 + temp2*temp2 +
		qweight[qworder[2]]*qweight[qworder[2]];
	discreteweight[1] = treebits[qworder[0]] + treebits[qworder[1]];

	/* 111 distribution */
	temp1 = onethird - qweight[qworder[0]];
	temp2 = onethird - qweight[qworder[1]];
	temp3 = onethird - qweight[qworder[2]];
	sqdiff[2] = temp1 * temp1 + temp2 * temp2 + temp3 * temp3;
	discreteweight[2] = (unsigned char) 7;

	/* sort in descending order */
	sort3doubles(sqdiff, sqorder);
			
	qpbranching = (unsigned char) discreteweight[sqorder[2]];
							
	if (qpbranching == 1) {
		reg1++;
		if (w2 < w3) reg1l++;
		else reg1r++;
	}
	if (qpbranching == 2) {
		reg2++;
		if (w1 < w3) reg2d++;
		else reg2u++;
	}
	if (qpbranching == 4) {
		reg3++;
		if (w1 < w2) reg3d++;
		else reg3u++;
	}
	if (qpbranching == 3) {
		reg4++;
		if (w1 < w2) reg4d++;
		else reg4u++;
	}
	if (qpbranching == 6) {
		reg5++;
		if (w2 < w3) reg5l++;
		else reg5r++;
	}
	if (qpbranching == 5) {
		reg6++;
		if (w1 < w3) reg6d++;
		else reg6u++;
	}
	if (qpbranching == 7) reg7++;
}

/* print tree statistics */
void printtreestats(FILE *ofp)
{
	int i, j, besttree;
	double bestlkl, difflkl, difflklps, temp, sum;
	
	/* find best tree */
	besttree = 0;
	bestlkl = ulkl[0];
	for (i = 1; i < numutrees; i++)
		if (ulkl[i] > bestlkl) {
			besttree = i;
			bestlkl = ulkl[i];
		}
	
	fprintf(ofp, "\n\nCOMPARISON OF USER TREES (NO CLOCK)\n\n");
	fprintf(ofp, "Tree   log L   difference     S.E.   Significantly worse\n");
	fprintf(ofp, "--------------------------------------------------------\n");
	for (i = 0; i < numutrees; i++) {
		difflkl = ulkl[besttree]-ulkl[i];
		fprintf(ofp, "%2d %10.2f %8.2f ", i+1, ulkl[i], difflkl);
		if (i == besttree) {
			fprintf(ofp, " <----------------- best tree");
		} else {
			/* compute variance of Log L differences over sites */
			difflklps = difflkl/(double)Maxsite;
			sum = 0.0;
			for (j = 0; j < Numptrn; j++) {
				temp = allsites[besttree][j] - allsites[i][j] - difflklps;
				sum += temp*temp*Weight[j];
			}
			sum = sqrt(fabs(sum/(Maxsite-1.0)*Maxsite));
			fprintf(ofp, "%11.2f         ", sum);
			if (difflkl > 1.96*sum)
				fprintf(ofp, "yes");
			else
				fprintf(ofp, "no");
		}
		fprintf(ofp, "\n");
	}
	fprintf(ofp, "\nThis test (5%% significance) follows Kishino and Hasegawa (1989).\n");
	
	if (compclock) {
	
		/* find best tree */
		besttree = 0;
		bestlkl = ulklc[0];
		for (i = 1; i < numutrees; i++)
			if (ulklc[i] > bestlkl) {
				besttree = i;
				bestlkl = ulklc[i];
			}
	
		fprintf(ofp, "\n\nCOMPARISON OF USER TREES (WITH CLOCK)\n\n");
		fprintf(ofp, "Tree   log L   difference     S.E.   Significantly worse\n");
		fprintf(ofp, "--------------------------------------------------------\n");
		for (i = 0; i < numutrees; i++) {
			difflkl = ulklc[besttree]-ulklc[i];
			fprintf(ofp, "%2d %10.2f %8.2f ", i+1, ulklc[i], difflkl);
			if (i == besttree) {
				fprintf(ofp, " <----------------- best tree");
			} else {
				/* compute variance of Log L differences over sites */
				difflklps = difflkl/(double)Maxsite;
				sum = 0.0;
				for (j = 0; j < Numptrn; j++) {
					temp = allsitesc[besttree][j] - allsitesc[i][j] - difflklps;
					sum += temp*temp*Weight[j];
				}
				sum = sqrt(fabs(sum/(Maxsite-1.0)*Maxsite));
				fprintf(ofp, "%11.2f         ", sum);
				if (difflkl > 1.96*sum)
					fprintf(ofp, "yes");
				else
					fprintf(ofp, "no");
			}
			fprintf(ofp, "\n");
		}
		fprintf(ofp, "\nThis test (5%% significance) follows Kishino and Hasegawa (1989).\n");
	}
}

/* time stamp */
void timestamp(FILE* ofp)
{
	double timespan;
	double cpuspan;
	timespan = difftime(Stoptime, Starttime);
	cpuspan  = ((double) (Stopcpu - Startcpu) / CLOCKS_PER_SEC);
	fprintf(ofp, "\n\nTIME STAMP\n\n");
	fprintf(ofp, "Date and time: %s", asctime(localtime(&Starttime)) );
	fprintf(ofp, "Runtime (excl. input) : %.0f seconds (= %.1f minutes = %.1f hours)\n",
		timespan, timespan/60., timespan/3600.);
	fprintf(ofp, "Runtime (incl. input) : %.0f seconds (= %.1f minutes = %.1f hours)\n",
		fulltime, fulltime/60., fulltime/3600.);
#ifdef TIMEDEBUG
	fprintf(ofp, "CPU time (incl. input): %.0f seconds (= %.1f minutes = %.1f hours)\n\n",
		fullcpu, fullcpu/60., fullcpu/3600.);
#endif /* TIMEDEBUG */

}

/* extern int bestrfound; */

/* write output file */
void writeoutputfile(FILE *ofp, int part)
{
	/* CZ */
}


/******************************************************************************/
/* timer routines                                                             */
/******************************************************************************/

/* start timer */
void starttimer()
{
	time(&time0);
	time1 = time0;
}

/* check remaining time and print message if necessary */
void checktimer(uli numqts)
{
	double tc2, mintogo, minutes, hours;

	time(&time2);
	if ( (time2 - time1) > 900) { /* generate message every 15 minutes */
		/* every 900 seconds */
		/* percentage of completed quartets */
		if (mflag == 0) {
			mflag = 1;
			FPRINTF(STDOUTFILE "\n");
		}
		tc2 = 100.*numqts/Numquartets;
		mintogo = (100.0-tc2) *
			(double) (time2-time0)/60.0/tc2;
		hours = floor(mintogo/60.0);
		minutes = mintogo - 60.0*hours;
		FPRINTF(STDOUTFILE "%.2f%%", tc2);
		FPRINTF(STDOUTFILE " completed (remaining");
		FPRINTF(STDOUTFILE " time: %.0f", hours);
		FPRINTF(STDOUTFILE " hours %.0f", minutes);
		FPRINTF(STDOUTFILE " minutes)\n");
		fflush(STDOUT);
		time1 = time2;
	}

}

/* check remaining time and print message if necessary */
void checktimer2(uli numqts, uli all, int flag)
{
	double tc2, mintogo, minutes, hours;

	static time_t tt1;
	static time_t tt2;

	if (flag == 1) {
		time(&tt1);
		time(&tt2);
	} else {
		time(&tt2);
		if ( (tt2 - tt1) > 900) { /* generate message every 15 minutes */
			/* every 900 seconds */
			/* percentage of completed quartets */
			if (mflag == 0) {
				mflag = 1;
				FPRINTF(STDOUTFILE "\n");
			}
			tc2 = 100.*numqts/Numquartets;
			mintogo = (100.0-tc2) *
				(double) (tt2-time0)/60.0/tc2;
			hours = floor(mintogo/60.0);
			minutes = mintogo - 60.0*hours;
			FPRINTF(STDOUTFILE "%.2f%%", tc2);
			FPRINTF(STDOUTFILE " completed (remaining");
			FPRINTF(STDOUTFILE " time: %.0f", hours);
			FPRINTF(STDOUTFILE " hours %.0f", minutes);
			FPRINTF(STDOUTFILE " minutes)\n");
			fflush(STDOUT);
			tt1 = tt2;
		}
	}
}

void resetqblocktime(timearray_t *ta)
{
	ta->quartcpu += ta->quartblockcpu;
	ta->quartblockcpu = 0.0;
	ta->quarttime += ta->quartblocktime;
	ta->quartblocktime = 0.0;
} /* resetqblocktime */


void resetpblocktime(timearray_t *ta)
{
	ta->puzzcpu += ta->puzzblockcpu;
	ta->puzzblockcpu = 0.0;
	ta->puzztime += ta->puzzblocktime;
	ta->puzzblocktime = 0.0;
} /* resetpblocktime */


#ifdef TIMEDEBUG
void printtimearr(timearray_t *ta)
{
#	if ! PARALLEL
	  int PP_Myid;
	  PP_Myid = -1;
#	endif
	printf("(%2d) MMCPU: %11ld  /  %11ld \n", PP_Myid, ta->maxcpu, ta->mincpu);
	printf("(%2d) CTick: %11.6f [tks]  /  %11.6f [s] \n", PP_Myid, ta->mincputick, ta->mincputicktime);

	printf("(%2d) MMTIM: %11ld  /  %11ld \n", PP_Myid, ta->maxtime, ta->mintime);

	printf("(%2d) Mxblk: %11.6e  /  %11.6e \n", PP_Myid, ta->maxcpublock, ta->maxtimeblock);
	printf("(%2d) Mnblk: %11.6e  /  %11.6e \n", PP_Myid, ta->mincpublock, ta->mintimeblock);

	printf("(%2d) Gnrl: %11.6e  /  %11.6e \n", PP_Myid, ta->generalcpu, ta->generaltime);
	printf("(%2d) Optn: %11.6e  /  %11.6e \n", PP_Myid, ta->optionscpu, ta->optionstime);
	printf("(%2d) Estm: %11.6e  /  %11.6e \n", PP_Myid, ta->paramestcpu, ta->paramesttime);
	printf("(%2d) Qurt: %11.6e  /  %11.6e \n", PP_Myid, ta->quartcpu, ta->quarttime);
	printf("(%2d) QBlk: %11.6e  /  %11.6e \n", PP_Myid, ta->quartblockcpu, ta->quartblocktime);
	printf("(%2d) QMax: %11.6e  /  %11.6e \n", PP_Myid, ta->quartmaxcpu, ta->quartmaxtime);
	printf("(%2d) QMin: %11.6e  /  %11.6e \n", PP_Myid, ta->quartmincpu, ta->quartmintime);

	printf("(%2d) Puzz: %11.6e  /  %11.6e \n", PP_Myid, ta->puzzcpu, ta->puzztime);
	printf("(%2d) PBlk: %11.6e  /  %11.6e \n", PP_Myid, ta->puzzblockcpu, ta->puzzblocktime);
	printf("(%2d) PMax: %11.6e  /  %11.6e \n", PP_Myid, ta->puzzmaxcpu, ta->puzzmaxtime);
	printf("(%2d) PMin: %11.6e  /  %11.6e \n", PP_Myid, ta->puzzmincpu, ta->puzzmintime);

	printf("(%2d) Tree: %11.6e  /  %11.6e \n", PP_Myid, ta->treecpu, ta->treetime);
	printf("(%2d) TBlk: %11.6e  /  %11.6e \n", PP_Myid, ta->treeblockcpu, ta->treeblocktime);
	printf("(%2d) TMax: %11.6e  /  %11.6e \n", PP_Myid, ta->treemaxcpu, ta->treemaxtime);
	printf("(%2d) TMin: %11.6e  /  %11.6e \n", PP_Myid, ta->treemincpu, ta->treemintime);

	printf("(%2d) C/T : %11.6e / %11.6e \n", PP_Myid, 
		(ta->generalcpu + ta->optionscpu + ta->paramestcpu + ta->quartblockcpu + ta->puzzblockcpu + ta->treeblockcpu),
		(ta->generaltime + ta->optionstime + ta->paramesttime + ta->quartblocktime + ta->puzzblocktime + ta->treeblocktime));
	printf("(%2d)  CPU: %11.6e /  Time: %11.6e \n", PP_Myid, ta->cpu, ta->time);
	printf("(%2d) aCPU: %11.6e / aTime: %11.6e \n", PP_Myid, ta->fullcpu, ta->fulltime);

} /* printtimearr */
#endif /* TIMEDEBUG */

char *jtype [7];

void inittimearr(timearray_t *ta)
{
	clock_t c0, c1, c2; 

	jtype[OVERALL]  = "OVERALL";
	jtype[GENERAL]  = "GENERAL";
	jtype[OPTIONS]  = "OPTIONS";
	jtype[PARAMEST] = "PARAMeter ESTimation";
	jtype[QUARTETS] = "QUARTETS";
	jtype[PUZZLING] = "PUZZLING steps";
	jtype[TREEEVAL] = "TREE EVALuation";
	ta->currentjob = GENERAL;

	c1 = clock();
	c2 = clock();
	while (c1 == c2)
	   c2 = clock(); 
	ta->mincputick = (double)(c2 - c1);
	ta->mincputicktime = ((double)(c2 - c1))/CLOCKS_PER_SEC;

	ta->tempcpu = clock();
	ta->tempcpustart = ta->tempcpu;
	ta->tempfullcpu  = ta->tempcpu;
	time(&(ta->temptime));
	ta->temptimestart = ta->temptime;
        ta->tempfulltime  = ta->temptime;

	c0=0; c1=0; c2=(clock_t)((2 * c1) + 1);;
	while (c1 < c2) {
	   c0 = c1;
	   c1 = c2;
	   c2 = (clock_t)((2 * c1) + 1);
	}
	if (c1 == c2) ta->maxcpu=c0;
	if (c1  > c2) ta->maxcpu=c1;

	c0=0; c1=0; c2=(clock_t)((2 * c1) - 1);
	while (c1 > c2) {
	   c0 = c1;
	   c1 = c2;
	   c2 = (clock_t)((2 * c1) - 1);
	}
	if (c1 == c2) ta->mincpu=c0;
	if (c1  < c2) ta->mincpu=c1;



	ta->maxtime        = 0;
	ta->mintime        = 0;

	ta->maxcpublock    = 0;
	ta->mincpublock    = DBL_MAX;
	ta->maxtimeblock   = 0;
	ta->mintimeblock   = DBL_MAX;

	ta->cpu            = 0.0;
	ta->time           = 0.0;

	ta->fullcpu        = 0.0;
	ta->fulltime       = 0.0;

	ta->generalcpu     = 0.0;
	ta->optionscpu     = 0.0;
	ta->paramestcpu    = 0.0;
	ta->quartcpu       = 0.0;
	ta->quartblockcpu  = 0.0;
	ta->quartmaxcpu    = 0.0;
	ta->quartmincpu    = ((double) ta->maxcpu)/CLOCKS_PER_SEC;
	ta->puzzcpu        = 0.0;
	ta->puzzblockcpu   = 0.0;
	ta->puzzmaxcpu     = 0.0;
	ta->puzzmincpu     = ((double) ta->maxcpu)/CLOCKS_PER_SEC;
	ta->treecpu        = 0.0;
	ta->treeblockcpu   = 0.0;
	ta->treemaxcpu     = 0.0;
	ta->treemincpu     = ((double) ta->maxcpu)/CLOCKS_PER_SEC;

	ta->generaltime    = 0.0;
	ta->optionstime    = 0.0;
	ta->paramesttime   = 0.0;
	ta->quarttime      = 0.0;
	ta->quartblocktime = 0.0;
	ta->quartmaxtime   = 0.0;
	ta->quartmintime   = DBL_MAX;
	ta->puzztime       = 0.0;
	ta->puzzblocktime  = 0.0;
	ta->puzzmaxtime    = 0.0;
	ta->puzzmintime    = DBL_MAX;
	ta->treetime       = 0.0;
	ta->treeblocktime  = 0.0;
	ta->treemaxtime    = 0.0;
	ta->treemintime    = DBL_MAX;
} /* inittimearr */


/***************/

void addup(int jobtype, clock_t c1, clock_t c2, time_t t1, time_t t2, timearray_t *ta)
{
   double  c,
           t;

   if (t2 != t1) t = difftime(t2, t1); 
   else          t = 0.0;

   if (c2 < c1)
      c = ((double)(c2 - ta->mincpu))/CLOCKS_PER_SEC +
          ((double)(ta->maxcpu - c1))/CLOCKS_PER_SEC;
   else
      c = ((double)(c2 - c1))/CLOCKS_PER_SEC;

   if (jobtype != OVERALL) {

      if (ta->mincpublock  > c) ta->mincpublock = c;
      if (ta->maxcpublock  < c) ta->maxcpublock = c;
      if (ta->mintimeblock > t) ta->mintimeblock = t;
      if (ta->maxtimeblock < t) ta->maxtimeblock = t;

      switch (jobtype) {
         case GENERAL: ta->generalcpu  += c;
                       ta->generaltime += t;
                       break;
         case OPTIONS: ta->optionscpu  += c;
                       ta->optionstime += t;
                       break;
         case PARAMEST: ta->paramestcpu  += c;
                        ta->paramesttime += t;
                        break;
         case QUARTETS: ta->quartblockcpu  += c;
                        ta->quartblocktime += t;
                        if (ta->quartmincpu  > c) ta->quartmincpu = c;
                        if (ta->quartmaxcpu  < c) ta->quartmaxcpu = c;
                        if (ta->quartmintime > t) ta->quartmintime = t;
                        if (ta->quartmaxtime < t) ta->quartmaxtime = t;
                        break;
         case PUZZLING: ta->puzzblockcpu  += c;
                        ta->puzzblocktime += t;
                        if (ta->puzzmincpu  > c) ta->puzzmincpu = c;
                        if (ta->puzzmaxcpu  < c) ta->puzzmaxcpu = c;
                        if (ta->puzzmintime > t) ta->puzzmintime = t;
                        if (ta->puzzmaxtime < t) ta->puzzmaxtime = t;
                        break;
         case TREEEVAL: ta->treeblockcpu  += c;
                        ta->treeblocktime += t;
                        if (ta->treemincpu  > c) ta->treemincpu = c;
                        if (ta->treemaxcpu  < c) ta->treemaxcpu = c;
                        if (ta->treemintime > t) ta->treemintime = t;
                        if (ta->treemaxtime < t) ta->treemaxtime = t;
                        break;
      }
      ta->cpu  += c;
      ta->time += t;
   
   } else {
         ta->fullcpu  += c;
         ta->fulltime += t;
   }

#     ifdef TIMEDEBUG
         {
#        if ! PARALLEL
            int PP_Myid =  -1;
#        endif /* !PARALLEL */
         printf("(%2d) CPU: +%10.6f / Time: +%10.6f (%s)\n", PP_Myid, c, t, jtype[jobtype]);
         printf("(%2d) CPU: %11.6f / Time: %11.6f (%s)\n", PP_Myid, ta->cpu, ta->time, jtype[jobtype]);
         printf("(%2d) CPU: %11.6f / Time: %11.6f (%s)\n", PP_Myid, ta->fullcpu, ta->fulltime, jtype[jobtype]);
         }
#     endif /* TIMEDEBUG */
} /* addup */


/***************/


void addtimes(int jobtype, timearray_t *ta)
{
   clock_t tempc;
   time_t  tempt;

   time(&tempt);
   tempc = clock();

   if ((tempc < ta->tempfullcpu) || (jobtype == OVERALL)) {   /* CPU counter overflow for overall time */
   	addup(OVERALL, ta->tempfullcpu, tempc, ta->tempfulltime, tempt, ta);
	ta->tempfullcpu  = tempc;
	ta->tempfulltime = tempt;
        if (jobtype == OVERALL) {
   	   addup(ta->currentjob, ta->tempcpustart, tempc, ta->temptimestart, tempt, ta);
	   ta->tempcpustart = ta->tempcpu;
	   ta->tempcpu      = tempc;
	   ta->temptimestart = ta->temptime;
	   ta->temptime      = tempt;
        }
   }

   if((jobtype != ta->currentjob) && (jobtype != OVERALL)) {   /* change of job type */
   	addup(ta->currentjob, ta->tempcpustart, ta->tempcpu, ta->temptimestart, ta->temptime, ta);
	ta->tempcpustart  = ta->tempcpu;
	ta->tempcpu       = tempc;
	ta->temptimestart = ta->temptime;
	ta->temptime      = tempt;
	ta->currentjob    = jobtype;
   }

   if (tempc < ta->tempcpustart) {   /* CPU counter overflow */
   	addup(jobtype, ta->tempcpustart, tempc, ta->temptimestart, tempt, ta);
	ta->tempcpustart = ta->tempcpu;
	ta->tempcpu      = tempc;
	ta->temptimestart = ta->temptime;
	ta->temptime      = tempt;
   }

} /* addtimes */



/******************************************************************************/

/* estimate parameters of substitution process and rate heterogeneity - no tree
   n-taxon tree is not needed because of quartet method or NJ tree topology */
void estimateparametersnotree()
{
	int it, nump, change;
	double TSold, YRold, FIold, GEold;

	it = 0;
	nump = 0;

	/* count number of parameters */
	if (data_optn == NUCLEOTIDE && optim_optn) nump++;
	if (fracinv_optim || grate_optim) nump++;

	do { /* repeat until nothing changes any more */
		it++;
		change = FALSE;

		/* optimize substitution parameters */
		if (data_optn == NUCLEOTIDE && optim_optn) {
		
			TSold = TSparam;
			YRold = YRparam;
			

			/*
			 * optimize
			 */
			 
			FPRINTF(STDOUTFILE "Optimizing missing substitution process parameters\n");
			fflush(STDOUT);

			if (qcalg_optn) { /* quartet sampling */
				optimseqevolparamsq();
			} else { /* NJ tree */
				tmpfp = tmpfile();
				njtree(tmpfp);				
				rewind(tmpfp);
				readusertree(tmpfp);
				closefile(tmpfp);
				optimseqevolparamst();				
			}
			
			computedistan(); /* update ML distances */
		
			/* same tolerance as 1D minimization */
			if ((fabs(TSparam - TSold) > 3.3*PEPS1) ||
				(fabs(YRparam - YRold) > 3.3*PEPS1)
			) change = TRUE;

		}

		/* optimize rate heterogeneity variables */
		if (fracinv_optim || grate_optim) {

			FIold = fracinv;
			GEold = Geta;


			/* 
			 * optimize 
			 */

			FPRINTF(STDOUTFILE "Optimizing missing rate heterogeneity parameters\n");
			fflush(STDOUT);
			/* compute NJ tree */
			tmpfp = tmpfile();
			njtree(tmpfp);
			/* use NJ tree topology to estimate parameters */
			rewind(tmpfp);
			readusertree(tmpfp);
			closefile(tmpfp);

			optimrateparams();				
			computedistan(); /* update ML distances */


			/* same tolerance as 1D minimization */
			if ((fabs(fracinv - FIold) > 3.3*PEPS2) ||
				(fabs(Geta - GEold) > 3.3*PEPS2)
			) change = TRUE;
			
		}

		if (nump == 1) return;
		
	} while (it != MAXITS && change);

	return;
}


/* estimate parameters of substitution process and rate heterogeneity - tree
   same as above but here the n-taxon tree is already in memory */
void estimateparameterstree()
{
	int it, nump, change;
	double TSold, YRold, FIold, GEold;

	it = 0;
	nump = 0;

	/* count number of parameters */
	if (data_optn == NUCLEOTIDE && optim_optn) nump++;
	if (fracinv_optim || grate_optim) nump++;

	do { /* repeat until nothing changes any more */
		it++;
		change = FALSE;

		/* optimize substitution process parameters */
		if (data_optn == NUCLEOTIDE && optim_optn) {
		
			TSold = TSparam;
			YRold = YRparam;
			

			/*
			 * optimize
			 */
			 
			FPRINTF(STDOUTFILE "Optimizing missing substitution process parameters\n");
			fflush(STDOUT);
			optimseqevolparamst();
			computedistan(); /* update ML distances */

		
			/* same tolerance as 1D minimization */
			if ((fabs(TSparam - TSold) > 3.3*PEPS1) ||
				(fabs(YRparam - YRold) > 3.3*PEPS1)
			) change = TRUE;

		}

		/* optimize rate heterogeneity variables */
		if (fracinv_optim || grate_optim) {

			FIold = fracinv;
			GEold = Geta;


			/* 
			 * optimize 
			 */

			FPRINTF(STDOUTFILE "Optimizing missing rate heterogeneity parameters\n");
			fflush(STDOUT);
			optimrateparams();				
			computedistan(); /* update ML distances */


			/* same tolerance as 1D minimization */
			if ((fabs(fracinv - FIold) > 3.3*PEPS2) ||
				(fabs(Geta - GEold) > 3.3*PEPS2)
			) change = TRUE;
			
		}

		if (nump == 1) return;
		
	} while (it != MAXITS && change);

	return;
}


/******************************************************************************/
/* exported from main                                                         */
/******************************************************************************/

void compute_quartlklhds(int a, int b, int c, int d, double *d1, double *d2, double *d3, int approx)
{
	if (approx == APPROX) {

		*d1 = quartet_alklhd(a,b, c,d); /* (a,b)-(c,d) */
		*d2 = quartet_alklhd(a,c, b,d); /* (a,c)-(b,d) */
		*d3 = quartet_alklhd(a,d, b,c); /* (a,d)-(b,c) */

	} else /* approx == EXACT */ {

		*d1 = quartet_lklhd(a,b, c,d);  /* (a,b)-(c,d) */
		*d2 = quartet_lklhd(a,c, b,d);  /* (a,c)-(b,d) */
		*d3 = quartet_lklhd(a,d, b,c);  /* (a,d)-(b,c) */

	}
}

/***************************************************************/ 

void recon_tree()
{	
	int i;
#	if ! PARALLEL
		int a, b, c;
		uli nq;
		double tc2, mintogo, minutes, hours;
#	endif

	/* allocate memory for taxon list of bad quartets */
	badtaxon = new_ulivector(Maxspc);
	for (i = 0; i < Maxspc; i++) badtaxon[i] = 0;

	/* allocate variable used for randomizing input order */
	trueID = new_ivector(Maxspc);

	/* allocate memory for quartets */
	quartetinfo = mallocquartets(Maxspc);
	
	/* prepare for consensus tree analysis */
	initconsensus();
		
	if (!(readquart_optn) || (readquart_optn && savequart_optn)) {
	  /* compute quartets */
	  FPRINTF(STDOUTFILE "Computing quartet maximum likelihood trees\n");
	  fflush(STDOUT);
	  computeallquartets();
	}

	if (savequart_optn) 
		writeallquarts(Maxspc, ALLQUART, quartetinfo);
	if (readquart_optn) {
		int xx1, xx2, xx3, xx4, count;
		readallquarts (Maxspc, ALLQUART, quartetinfo);
		if (show_optn) { /* list all unresolved quartets */
			openfiletowrite(&unresfp, UNRESOLVED, "unresolved quartet trees");
			fprintf(unresfp, "List of all completely unresolved quartets:\n\n");
		}

		/* initialize bad quartet memory */
		for (count = 0; count < Maxspc; count++) badtaxon[count] = 0;
		badqs = 0;

		for (xx4 = 3; xx4 < Maxspc; xx4++)
		  for (xx3 = 2; xx3 < xx4; xx3++)
		    for (xx2 = 1; xx2 < xx3; xx2++)
		      for (xx1 = 0; xx1 < xx2; xx1++) {
			if (readquartet(xx1, xx2, xx3, xx4) == 7) {
				badqs++;
				badtaxon[xx1]++;
				badtaxon[xx2]++;
				badtaxon[xx3]++;
				badtaxon[xx4]++;
				if (show_optn) {
					fputid10(unresfp, xx1);
					fprintf(unresfp, "  ");
					fputid10(unresfp, xx2);
					fprintf(unresfp, "  ");
					fputid10(unresfp, xx3);
					fprintf(unresfp, "  ");
					fputid  (unresfp, xx4);
					fprintf(unresfp, "\n");
				}
			}
		} /* end for xx4; for xx3; for xx2; for xx1 */
		if (show_optn)  /* list all unresolved quartets */
			fclose(unresfp);
	} /* readquart_optn */

#	if PARALLEL
		PP_SendAllQuarts(numquarts(Maxspc), quartetinfo);
#	endif /* PARALLEL */

	FPRINTF(STDOUTFILE "Computing quartet puzzling tree\n");
	fflush(STDOUT);
	
	/* start timer - percentage of completed trees */
	time(&time0);
	time1 = time0;
	mflag = 0;
			
	/* open file for chronological list of puzzling step trees */
	if((listqptrees == PSTOUT_LIST) || (listqptrees == PSTOUT_LISTORDER))
		openfiletowrite(&qptlist, OUTPTLIST, "puzzling step trees (chonological)");
							
#	if PARALLEL
		{
		PP_SendDoPermutBlock(Numtrial);
		}
#	else
		addtimes(GENERAL, &tarr);
		for (Currtrial = 0; Currtrial < Numtrial; Currtrial++) {
			
			/* randomize input order */
			chooser(Maxspc, Maxspc, trueID);
				
			/* initialize tree */
			inittree();

			/* adding all other leafs */
			for (i = 3; i < Maxspc; i++) { 
					
				/* clear all edgeinfos */
				resetedgeinfo();
				
				/* clear counter of quartets */
				nq = 0;

				/*
				 * core of quartet puzzling algorithm
			 	 */

				for (a = 0; a < nextleaf - 2; a++)
					for (b = a + 1; b < nextleaf - 1; b++)
						for (c = b + 1; c < nextleaf; c++) {

							/*  check which two _leaves_ out of a, b, c
							    are closer related to each other than
							    to leaf i according to a least squares
							    fit of the continous Baysian weights to the
							    seven trivial "attractive regions". We assign 
							    a score of 1 to all edges between these two leaves
							    chooseA and chooseB */

							checkquartet(a, b, c, i);
							incrementedgeinfo(chooseA, chooseB);

							nq++;

							/* generate message every 15 minutes */
								
							/* check timer */
							time(&time2);
							if ( (time2 - time1) > 900) {
								/* every 900 seconds */
								/* percentage of completed trees */
								if (mflag == 0) {
									FPRINTF(STDOUTFILE "\n");
									mflag = 1;
								}
								tc2 = 100.0*Currtrial/Numtrial + 
									100.0*nq/Numquartets/Numtrial;
									mintogo = (100.0-tc2) *
									(double) (time2-time0)/60.0/tc2;
									hours = floor(mintogo/60.0);
									minutes = mintogo - 60.0*hours;
									FPRINTF(STDOUTFILE "%2.2f%%", tc2);
									FPRINTF(STDOUTFILE " completed (remaining");
									FPRINTF(STDOUTFILE " time: %.0f", hours);
									FPRINTF(STDOUTFILE " hours %.0f", minutes);
									FPRINTF(STDOUTFILE " minutes)\n");
									fflush(STDOUT);
									time1 = time2;
								}
							}

				/* find out which edge has the lowest edgeinfo */
				minimumedgeinfo();

				/* add the next leaf on minedge */
				addnextleaf(minedge);
			}
	
			/* compute bipartitions of current tree */
			computebiparts();
			makenewsplitentries();

			{
				int *ctree, startnode;
				char *trstr;
				treelistitemtype *treeitem;
				ctree = initctree();
				copytree(ctree);
				startnode = sortctree(ctree);
				trstr=sprintfctree(ctree, psteptreestrlen);


				treeitem = addtree2list(&trstr, 1, &psteptreelist, &psteptreenum, &psteptreesum);

				if((listqptrees == PSTOUT_LIST) 
				  || (listqptrees == PSTOUT_LISTORDER)) {
					/* print: order no/# topol per this id/tree id/sum of unique topologies/sum of trees so far */
					fprintf(qptlist, "%ld.\t1\t%d\t%d\t%d\t%d\n", 
						Currtrial + 1, (*treeitem).count, (*treeitem).id, psteptreenum, psteptreesum);
				}

#				ifdef VERBOSE1
					printf("%s\n", trstr);
					printfsortedpstrees(psteptreelist);
#				endif
				freectree(&ctree);
			}



			/* free tree before building the next tree */
			freetree();
			
			addtimes(PUZZLING, &tarr);
		}
#	endif /* PARALLEL */
			
			/* close file for list of puzzling step trees */
			if((listqptrees == PSTOUT_LIST) || (listqptrees == PSTOUT_LISTORDER))
				closefile(qptlist);
			
			if (mflag == 1) FPRINTF(STDOUTFILE "\n");

			/* garbage collection */
			free(splitcomp);
			free_ivector(trueID);

#                       if ! PARALLEL
				free_cmatrix(biparts);
#			endif /* PARALLEL */

			freequartets();

			/* compute majority rule consensus tree */
			makeconsensus();
			
			/* write consensus tree to tmp file */
			tmpfp = tmpfile();
			writeconsensustree(tmpfp);
} /* recon_tree */

/***************************************************************/ 

void map_lklhd()
{	
	int i, a, a1, a2, b, b1, b2, c, c1, c2, d;
	uli nq;
	double logs[3], d1, d2, d3, temp;
	ivector qts, mlorder, gettwo;
		/* reset variables */
		ar1 = ar2 = ar3 = 0;
		reg1 = reg2 = reg3 = reg4 = reg5 = reg6 = reg7 = 0;
		reg1l = reg1r = reg2u = reg2d = reg3u = reg3d = reg4u =
			reg4d = reg5l = reg5r = reg6u = reg6d = 0;
			
		/* place for random quartet */
			qts = new_ivector(4);

		/* initialize output file */
		openfiletowrite(&trifp, TRIANGLE, "Postscript output");
		initps(trifp);
		FPRINTF(STDOUTFILE "Performing likelihood mapping analysis\n");
		fflush(STDOUT);
			
		/* start timer */
		starttimer();
		nq = 0;
		mflag = 0;

		addtimes(GENERAL, &tarr);
		if (lmqts == 0) { /* all possible quartets */
			
			if (numclust == 4) { /* four-cluster analysis */

				for (a = 0; a < clustA; a++)
					for (b = 0; b < clustB; b++)
						for (c = 0; c < clustC; c++)
							for (d = 0; d < clustD; d++) {
									
								nq++;
									
								/* check timer */
								checktimer(nq);

								/* maximum likelihood values */
								/* approximate ML is sufficient */
								compute_quartlklhds(clusterA[a],clusterB[b],clusterC[c],clusterD[d],&d1,&d2,&d3, APPROX);
								
								/* draw point for LM analysis */
								makelmpoint(trifp, d1, d2, d3);
								addtimes(QUARTETS, &tarr);

							}					
			}
				
			if (numclust == 3) { /* three-cluster analysis */

				gettwo = new_ivector(2);

				for (a = 0; a < clustA; a++)
					for (b = 0; b < clustB; b++)
						for (c1 = 0; c1 < clustC-1; c1++)
							for (c2 = c1+1; c2 < clustC; c2++) {

								nq++;

								/* check timer */
								checktimer(nq);
								
								/* maximum likelihood values */
								/* approximate ML is sufficient */
								compute_quartlklhds(clusterA[a],clusterB[b],clusterC[c1],clusterC[c2],&d1,&d2,&d3, APPROX);
								
								/* randomize order of d2 and d3 */
								if (randominteger(2) == 1) {
									temp = d3;
									d3 = d2;
									d2 = temp;
								}
						
								/* draw point for LM analysis */
								makelmpoint(trifp, d1, d2, d3);
								addtimes(QUARTETS, &tarr);

				}				
				free_ivector(gettwo);
			}
				
			if (numclust == 2) { /* two-cluster analysis */
					
				gettwo = new_ivector(2);

				for (a1 = 0; a1 < clustA-1; a1++)
					for (a2 = a1+1; a2 < clustA; a2++)
						for (b1 = 0; b1 < clustB-1; b1++)
							for (b2 = b1+1; b2 < clustB; b2++) {

								nq++;

								/* check timer */
								checktimer(nq);
								
								/* maximum likelihood values */
								/* approximate ML is sufficient */
								compute_quartlklhds(clusterA[a1],clusterA[a2],clusterB[b1],clusterB[b2],&d1,&d2,&d3, APPROX);
								
								/* randomize order of d2 and d3 */
								if (randominteger(2) == 1) {
									temp = d3;
									d3 = d2;
									d2 = temp;
								}
						
								/* draw point for LM analysis */
								makelmpoint(trifp, d1, d2, d3);
								addtimes(QUARTETS, &tarr);

				}
					
				free_ivector(gettwo);
			}
			
			if (numclust == 1) { /* normal likelihood mapping (one cluster) */
			
				mlorder = new_ivector(3);

#if 0
				for (i = 3; i < Maxspc; i++) 
					for (a = 0; a < i - 2; a++) 
						for (b = a + 1; b < i - 1; b++) 
							for (c = b + 1; c < i; c++)
   for (d = 3; d < Maxspc; d++) 
      for (c = 2; c < d; c++) 
         for (b = 1; b < c; b++) 
            for (a = 0; a < b; a++) 
#endif

				for (i = 3; i < Maxspc; i++) 
					for (c = 2; c < i; c++) 
						for (b = 1; b < c; b++) 
							for (a = 0; a < b; a++) {
							
								nq++;

								/* check timer */
								checktimer(nq);

								/* maximum likelihood values */
								/* approximate ML is sufficient */
								compute_quartlklhds(a,b,c,i,&logs[0],&logs[1],&logs[2], APPROX);

								/* randomize order */
								chooser(3,3,mlorder);
								d1 = logs[mlorder[0]];
								d2 = logs[mlorder[1]];
								d3 = logs[mlorder[2]];
								
								/* draw point for LM analysis */
								makelmpoint(trifp, d1, d2, d3);
								addtimes(QUARTETS, &tarr);
				
				}
				free_ivector(mlorder);
			}
			
		} else { /* randomly selected quartets */
			
			if (numclust == 4) { /* four-cluster analysis */

				for (lmqts = 0; lmqts < Numquartets; lmqts++) {
				
					nq++;

					/* check timer */
					checktimer(nq);

					/* choose random quartet */	
					qts[0] = clusterA[ randominteger(clustA) ];
					qts[1] = clusterB[ randominteger(clustB) ];
					qts[2] = clusterC[ randominteger(clustC) ];
					qts[3] = clusterD[ randominteger(clustD) ];			

					/* maximum likelihood values */
					/* approximate ML is sufficient */
					compute_quartlklhds(qts[0],qts[1],qts[2],qts[3],&d1,&d2,&d3, APPROX);
					
					/* draw point for LM analysis */
					makelmpoint(trifp, d1, d2, d3);
					addtimes(QUARTETS, &tarr);

				}
			}
				
			if (numclust == 3) { /* three-cluster analysis */

				gettwo = new_ivector(2);

				for (lmqts = 0; lmqts < Numquartets; lmqts++) {
				
					nq++;

					/* check timer */
					checktimer(nq);

					/* choose random quartet */	
					qts[0] = clusterA[ randominteger(clustA) ];
					qts[1] = clusterB[ randominteger(clustB) ];
					chooser(clustC, 2, gettwo);
					qts[2] = clusterC[gettwo[0]];
					qts[3] = clusterC[gettwo[1]];		

					/* maximum likelihood values */
					/* approximate ML is sufficient */
					compute_quartlklhds(qts[0],qts[1],qts[2],qts[3],&d1,&d2,&d3, APPROX);
					
					/* order of d2 and d3 is already randomized! */
						
					/* draw point for LM analysis */
					makelmpoint(trifp, d1, d2, d3);
					addtimes(QUARTETS, &tarr);

				}
					
				free_ivector(gettwo);
			}
				
			if (numclust == 2) { /* two-cluster analysis */

				gettwo = new_ivector(2);

				for (lmqts = 0; lmqts < Numquartets; lmqts++) {
				
					nq++;

					/* check timer */
					checktimer(nq);

					/* choose random quartet */	
					chooser(clustA, 2, gettwo);
					qts[0] = clusterA[gettwo[0]];
					qts[1] = clusterA[gettwo[1]];		
					chooser(clustB, 2, gettwo);
					qts[2] = clusterB[gettwo[0]];
					qts[3] = clusterB[gettwo[1]];		

					/* maximum likelihood values */
					/* approximate ML is sufficient */
					compute_quartlklhds(qts[0],qts[1],qts[2],qts[3],&d1,&d2,&d3, APPROX);
					
					/* order of d2 and d3 is already randomized! */
						
					/* draw point for LM analysis */
					makelmpoint(trifp, d1, d2, d3);
					addtimes(QUARTETS, &tarr);

				}
				free_ivector(gettwo);
			}
				
			if (numclust == 1) { /* normal likelihood mapping (one cluster) */
			
				for (lmqts = 0; lmqts < Numquartets; lmqts++) {
				
					nq++;

					/* check timer */
					checktimer(nq);

					/* choose random quartet */	
					chooser(Maxspc, 4, qts);				

					/* maximum likelihood values */
					/* approximate ML is sufficient */
					compute_quartlklhds(qts[0],qts[1],qts[2],qts[3],&d1,&d2,&d3, APPROX);
					
					/* order of d1, d2, and d3 is already randomized! */
				
					/* draw point for LM analysis */
					makelmpoint(trifp, d1, d2, d3);
					addtimes(QUARTETS, &tarr);

				}
			}
		}

		finishps(trifp);
		closefile(trifp);
		free_ivector(qts);

} /* map_lklhd */

/***************************************************************/ 

void setdefaults() {

  strcpy(INFILE,     INFILEDEFAULT);
  strcpy(OUTFILE,    OUTFILEDEFAULT);
  strcpy(TREEFILE,   TREEFILEDEFAULT);
  strcpy(INTREE,     INTREEDEFAULT);
  strcpy(DISTANCES,  DISTANCESDEFAULT);
  strcpy(TRIANGLE,   TRIANGLEDEFAULT);
  strcpy(UNRESOLVED, UNRESOLVEDDEFAULT);
  strcpy(ALLQUART,   ALLQUARTDEFAULT);
  strcpy(ALLQUARTLH, ALLQUARTLHDEFAULT);
  strcpy(OUTPTLIST,  OUTPTLISTDEFAULT);
  strcpy(OUTPTORDER, OUTPTORDERDEFAULT);

  usebestq_optn    = FALSE;
  savequartlh_optn = FALSE;
  savequart_optn   = FALSE;
  readquart_optn   = FALSE;

  randseed = -1;                  /* to set random random seed */

} /* setdefaults */

/***************************************************************/ 

void printversion()
{
#  if ! PARALLEL
   fprintf(stderr, "puzzle (%s) %s\n", PACKAGE, VERSION);
#else
   fprintf(stderr, "ppuzzle (%s) %s\n", PACKAGE, VERSION);
#  endif
   exit (0);
}
/***************************************************************/ 

void printusage(char *fname)
{
   fprintf(stderr, "\n\nUsage: %s [-h] [ Infilename [ UserTreeFilename ] ]\n\n", fname);
#  if PARALLEL
	PP_SendDone();
	MPI_Finalize();
#  endif
   exit (1);
}

/***************************************************************/ 

#ifdef HHH
void printusagehhh(char *fname)
{
   fprintf(stderr, "\n\nUsage: %s [options] [ Infilename [ UserTreeFilename ] ]\n\n", fname);
   fprintf(stderr, "  -h           - print usage\n");
   fprintf(stderr, "  -wqf         - write quartet file to Infilename.allquart\n");
   fprintf(stderr, "  -rqf         - read quartet file from Infilename.allquart\n");
   fprintf(stderr, "  -wqlb        - write quart lhs to Infilename.allquartlh (binary)\n");
   fprintf(stderr, "  -wqla        - write quart lhs to Infilename.allquartlh (ASCII)\n");
   fprintf(stderr, "  -bestq       - use best quart, no basian weights\n");
   fprintf(stderr, "  -randseed<#> - use <#> as random number seed, for debug purposes only\n");
#  if PARALLEL
	PP_SendDone();
	MPI_Finalize();
#  endif
   exit (2);
}
#endif /* HHH */

/***************************************************************/ 


void scancmdline(int *argc, char **argv[]) 
{
   static short infileset     = 0;
   static short intreefileset = 0;
   short flagused;
   int n;
   int count, dummyint;

   for (n = 1; n < *argc; n++) {
#     ifdef VERBOSE1
         printf("argv[%d] = %s\n", n, (*argv)[n]);
#     endif
 
      flagused = FALSE;

#     ifdef HHH
      dummyint = 0;
      count = sscanf((*argv)[n], "-wqlb%n", &dummyint);
      if (dummyint == 5) {
        savequartlh_optn = TRUE;
        saveqlhbin_optn  = TRUE;
        flagused         = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n], "-wqla%n", &dummyint);
      if (dummyint == 5) {
        savequartlh_optn = TRUE;
        saveqlhbin_optn  = FALSE;
        flagused         = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n], "-wqf%n", &dummyint);
      if (dummyint == 4) {
        savequart_optn = TRUE;
        flagused = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n],"-rqf%n", &dummyint);
      if (dummyint == 4) {
        readquart_optn = TRUE; 
        flagused = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n],"-bestq%n", &dummyint);
      if (dummyint == 6) {
        usebestq_optn = TRUE; 
        flagused = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n],"-hhh%n", &dummyint);
      if (dummyint==4) {
        printusagehhh((*argv)[0]); 
        flagused = TRUE;
      }
#     endif /* HHH */

      dummyint = 0;
      count = sscanf((*argv)[n],"-V%n", &dummyint);
      if (dummyint==2) {
        printversion((*argv)[0]); 
        flagused = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n],"-version%n", &dummyint);
      if (dummyint==8) {
        printversion((*argv)[0]); 
        flagused = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n],"--version%n", &dummyint);
      if (dummyint>=4) {
        printversion((*argv)[0]); 
        flagused = TRUE;
      }

      dummyint = 0;
      count = sscanf((*argv)[n],"-h%n", &dummyint);
      if (dummyint==2) {
        printusage((*argv)[0]); 
        flagused = TRUE;
      }

      count = sscanf((*argv)[n],"-randseed%d", &dummyint);
      if (count == 1) {
        randseed = dummyint; 
        flagused = TRUE;
      }

#if 0
      count = sscanf((*argv)[n],"-h%n", &dummyint);
      if ((count == 1) && (dummyint>=2)) printusage((*argv)[0]); 

      count = sscanf((*argv)[n],"-writequarts%n", &dummyint);
      if (count == 1) writequartstofile = 1;; 

      count = sscanf((*argv)[n],"-ws%d", &dummyint);
      if (count == 1) windowsize = dummyint; 
#endif

      if ((*argv)[n][0] != '-') {
         if (infileset == 0) {
            strcpy(INFILE, (*argv)[n]);
            infileset++;
            sprintf(OUTFILE    ,"%s.%s", INFILE, OUTFILEEXT);
            sprintf(TREEFILE   ,"%s.%s", INFILE, TREEFILEEXT);
            sprintf(DISTANCES  ,"%s.%s", INFILE, DISTANCESEXT);
            sprintf(TRIANGLE   ,"%s.%s", INFILE, TRIANGLEEXT);
            sprintf(UNRESOLVED ,"%s.%s", INFILE, UNRESOLVEDEXT);
            sprintf(ALLQUART   ,"%s.%s", INFILE, ALLQUARTEXT);
            sprintf(ALLQUARTLH ,"%s.%s", INFILE, ALLQUARTLHEXT);
            sprintf(OUTPTLIST  ,"%s.%s", INFILE, OUTPTLISTEXT);
            sprintf(OUTPTORDER ,"%s.%s", INFILE, OUTPTORDEREXT);
            FPRINTF(STDOUTFILE "Input file: %s\n", INFILE);
            flagused = TRUE;
         } else {
            if (intreefileset == 0) {
               strcpy(INTREE, (*argv)[n]);
               intreefileset++;
               sprintf(OUTFILE    ,"%s.%s", INTREE, OUTFILEEXT);
               sprintf(TREEFILE   ,"%s.%s", INTREE, TREEFILEEXT);
               sprintf(DISTANCES  ,"%s.%s", INTREE, DISTANCESEXT);
               FPRINTF(STDOUTFILE "Usertree file: %s\n", INTREE);
               flagused = TRUE;
            }
         }
      }
      if (flagused == FALSE) {
         fprintf(stderr, "WARNING: commandline parameter %d not recognized (\"%s\")\n", n, (*argv)[n]);
      }
      flagused = FALSE;
   }

} /* scancmdline */


/***************************************************************/ 

void inputandinit(int *argc, char **argv[]) {

	int ci;

	/* vectors used in QP and LM analysis */
	qweight = new_dvector(3);
	sqdiff = new_dvector(3);
	qworder = new_ivector(3);
	sqorder = new_ivector(3);
	
	/* Initialization and parsing of Commandline */
	setdefaults();
	scancmdline(argc, argv);

	/* initialize random numbers generator */
	if (randseed >= 0)
	   fprintf(stderr, "WARNING: random seed set to %d for debugging!\n", randseed);
	randseed = initrandom(randseed);

	psteptreelist = NULL;
        psteptreesum = 0;
	bestratefound = 0;

#	ifndef ALPHA
		FPRINTF(STDOUTFILE "\n\n\nWELCOME TO TREE-PUZZLE %s!\n\n\n", VERSION);
#	else
		FPRINTF(STDOUTFILE "\n\n\nWELCOME TO TREE-PUZZLE %s%s!\n\n\n", VERSION, ALPHA);
#	endif


	/* get sequences */
	openfiletoread(&seqfp, INFILE, "sequence data");
	getsizesites(seqfp);
	FPRINTF(STDOUTFILE "\nInput data set contains %d sequences of length %d\n", Maxspc, Maxseqc);
	getdataset(seqfp);
	closefile(seqfp);		
	data_optn = guessdatatype();
	
	/* translate characters into format used by ML engine */
	nuc_optn = TRUE; 
	SH_optn = FALSE;
	Seqchar = NULL;
	translatedataset();
	
	/* estimate base frequencies from data set */
	Freqtpm = NULL;
	Basecomp = NULL;
	estimatebasefreqs();
	
	/* guess model of substitution */
	guessmodel();

	/* initialize guess variables */
	auto_datatype = AUTO_GUESS;
	if (data_optn == AMINOACID) auto_aamodel = AUTO_GUESS;
	else                        auto_aamodel = AUTO_DEFAULT;
	/* save guessed amino acid options */
	guessDayhf_optn    = Dayhf_optn;
	guessJtt_optn      = Jtt_optn;
	guessmtrev_optn    = mtrev_optn;
	guesscprev_optn    = cprev_optn;
	guessblosum62_optn = blosum62_optn;
	guessvtmv_optn     = vtmv_optn;
	guesswag_optn     = wag_optn;
	guessauto_aamodel  = auto_aamodel;


	/* check for user specified tree */
	if ((utfp = fopen(INTREE, "r")) != NULL) {
		fclose(utfp);
		puzzlemode = USERTREE;
	} else {
		puzzlemode = QUARTPUZ;
	}

	/* reserve memory for cluster LM analysis */
	clusterA = new_ivector(Maxspc);
	clusterB = new_ivector(Maxspc);
	clusterC = new_ivector(Maxspc);
	clusterD = new_ivector(Maxspc);

	/* set options interactively */
	setoptions();

	/* open usertree file right after start */
	if (typ_optn == TREERECON_OPTN && puzzlemode == USERTREE) {
		openfiletoread(&utfp, INTREE, "user trees");
	}

	/* start main timer */
	time(&Starttime);
	Startcpu=clock();
	addtimes(OPTIONS, &tarr);

	/* symmetrize doublet frequencies if specified */
	symdoublets();

	/* initialise ML */
	mlstart();

	/* determine how many usertrees */
	if (typ_optn == TREERECON_OPTN && puzzlemode == USERTREE) {
		numutrees = 0;		
		do {
			ci = fgetc(utfp);
			if ((char) ci == ';') numutrees++;
		} while (ci != EOF);
		rewind(utfp);
		if (numutrees < 1) {
			FPRINTF(STDOUTFILE "Unable to proceed (no tree in input tree file)\n\n\n");
			exit(1);
		}		
	}
	
	/* check fraction of invariable sites */
	if ((rhetmode == TWORATE || rhetmode == MIXEDRATE) && !fracinv_optim)
		/* fraction of invariable site was specified manually */
		if (fracinv > MAXFI)
			fracinv = MAXFI;

	addtimes(GENERAL, &tarr);
	/* estimate parameters */
	if (!(typ_optn == TREERECON_OPTN && puzzlemode == USERTREE)) {
		/* no tree present */
		estimateparametersnotree();
	} else {
		if (utree_optn) {
			/* use 1st user tree */
			readusertree(utfp);
			rewind(utfp);
			estimateparameterstree();
		} else {
			/* don't use first user tree */
			estimateparametersnotree();
		}
	}
	addtimes(PARAMEST, &tarr);

	/* compute expected Ts/Tv ratio */
	if (data_optn == NUCLEOTIDE) computeexpectations();

} /* inputandinit */



/***************************************************************/ 

void evaluatetree(FILE *intreefp, FILE *outtreefp, int pmode, int utreenum, int maxutree, int *oldlocroot)
{

	switch (pmode) {
		case QUARTPUZ: /* read QP tree */
			readusertree(intreefp);			
			FPRINTF(STDOUTFILE "Computing maximum likelihood branch lengths (without clock)\n");
			fflush(STDOUT);
			usertree_lklhd();
			findbestratecombination();
			break;
		case USERTREE: /* read user tree */
			readusertree(intreefp);
			FPRINTF(STDOUTFILE "Computing maximum likelihood branch lengths (without clock) for tree # %d\n", utreenum+1);
			fflush(STDOUT);
			usertree_lklhd();
	 		if (maxutree > 1) {
				ulkl[utreenum] = Ctree->lklhd;
				allsitelkl(Ctree->condlkl, allsites[utreenum]);
			}
			if (utreenum==0) findbestratecombination();
			break;
	}


	if (compclock) { /* clocklike branch length */
		switch (pmode) {
			case QUARTPUZ:
				FPRINTF(STDOUTFILE "Computing maximum likelihood branch lengths (with clock)\n");
				fflush(STDOUT);
				break;
			case USERTREE:
				FPRINTF(STDOUTFILE "Computing maximum likelihood branch lengths (with clock) for tree # %d\n", utreenum+1);
				fflush(STDOUT);
				break;
		}
									
		/* find best place for root */
		rootsearch = 0;

		if (utreenum==0) locroot = *oldlocroot;
		else             *oldlocroot = locroot;

		if (locroot < 0) {
			locroot = findrootedge();
			rootsearch = 1;
		}
		/* if user-specified edge for root does not exist use displayed outgroup */
		if (!checkedge(locroot)) {
			locroot = outgroup; 
			rootsearch = 2;
		}
		/* compute likelihood */
		clock_lklhd(locroot);
		if (maxutree > 1) {
			ulklc[utreenum] = Ctree->lklhdc;
			allsitelkl(Ctree->condlkl, allsitesc[utreenum]);
		}

	}

	if (clockmode == 0)
		fprintf(outtreefp, "[ lh=%.6f ]", Ctree->lklhd);
	else
		fprintf(outtreefp, "[ lh=%.6f ]", Ctree->lklhdc);

	/* write ML branch length tree to outree file */
	clockmode = 0; /* nonclocklike branch lengths */
	fputphylogeny(outtreefp);
		
	/* clocklike branch lengths */
	if (compclock) {
		clockmode = 1;
		fputrooted(outtreefp, locroot);
	}
} /* evaluatetree */

/***************************************************************/ 

void memcleanup() { 
	if (puzzlemode == QUARTPUZ && typ_optn == TREERECON_OPTN) {
		free(splitfreqs);
		free(splitpatterns);
		free(splitsizes);
		free_ivector(consconfid);
		free_ivector(conssizes);
		free_cmatrix(consbiparts);
		free_ulivector(badtaxon);
	}
	free_cmatrix(Identif);
	free_dvector(Freqtpm);
	free_imatrix(Basecomp);
	free_ivector(clusterA);
	free_ivector(clusterB);
	free_ivector(clusterC);
	free_ivector(clusterD);
	free_dvector(qweight);
	free_dvector(sqdiff);
	free_ivector(qworder);
	free_ivector(sqorder);
	freetreelist(&psteptreelist, &psteptreenum, &psteptreesum);
} /* memcleanup */

/***************************************************************/ 


/******************************************************************************/
/* main part                                                                  */
/******************************************************************************/

int main(int argc, char *argv[])
{	
	int i, oldlocroot=0;
	
	/* start main timer */
	time(&walltimestart);
	cputimestart = clock();
	inittimearr(&tarr);


	
	inputandinit(&argc, &argv);

	

	/* write distance matrix */
	FPRINTF(STDOUTFILE "Writing pairwise distances to file %s\n", DISTANCES);
	openfiletowrite(&dfp, DISTANCES, "pairwise distances");
	putdistance(dfp);
	closefile(dfp);



	free_cmatrix(Seqchar);
	free_cmatrix(seqchars);

	


	/* write CPU/Wallclock times and parallel statistics */
	time(&walltimestop);
	cputimestop = clock();
	addtimes(OVERALL, &tarr);

        fullcpu          = tarr.fullcpu;
        fulltime         = tarr.fulltime;



	/* stop timer */
	
    time(&Stoptime);
	Stopcpu=clock();
    /*
	timestamp(ofp);
	closefile(ofp);
    CZ 05/16/01*/


     /* printbestratecombination(stderr); */
	mlfinish();

	FPRINTF(STDOUTFILE "\nAll results written to disk:\n");
	/*FPRINTF(STDOUTFILE "       Puzzle report file:         %s\n", OUTFILE);*/
	FPRINTF(STDOUTFILE "       Likelihood distances:       %s\n", DISTANCES);
	
	if (typ_optn == TREERECON_OPTN && puzzlemode != PAIRDIST) 
		FPRINTF(STDOUTFILE "       Phylip tree file:           %s\n", TREEFILE);
	if (typ_optn == TREERECON_OPTN && puzzlemode == QUARTPUZ) {
		if ((listqptrees == PSTOUT_ORDER) ||(listqptrees == PSTOUT_LISTORDER)) 
			FPRINTF(STDOUTFILE "       Unique puzzling step trees: %s\n", OUTPTORDER);	
		if ((listqptrees == PSTOUT_LIST) ||(listqptrees == PSTOUT_LISTORDER)) 
			FPRINTF(STDOUTFILE "       Puzzling step tree list:    %s\n", OUTPTLIST);	
	}
	if (show_optn && typ_optn == TREERECON_OPTN && puzzlemode == QUARTPUZ) 
		FPRINTF(STDOUTFILE "       Unresolved quartets:        %s\n", UNRESOLVED);
	if (typ_optn == LIKMAPING_OPTN) 
		FPRINTF(STDOUTFILE "       Likelihood mapping diagram: %s\n", TRIANGLE);
	FPRINTF(STDOUTFILE "\n");

	/* runtime message */
	FPRINTF(STDOUTFILE 
		"The computation took %.0f seconds (= %.1f minutes = %.1f hours)\n",
			difftime(Stoptime, Starttime), difftime(Stoptime, Starttime)/60.,
			difftime(Stoptime, Starttime)/3600.);
	FPRINTF(STDOUTFILE 
		"     including input %.0f seconds (= %.1f minutes = %.1f hours)\n",
			fulltime, fulltime/60., fulltime/3600.);


	/* free memory */
	memcleanup();



	return 0;
}


/* compare function for uli - sort largest numbers first */
int ulicmp(const void *ap, const void *bp)
{
	uli a, b;
	
	a = *((uli *) ap);
	b = *((uli *) bp);

	if (a > b) return -1;
	else if (a < b) return 1;
	else return 0;
}

/* compare function for int - sort smallest numbers first */
int intcmp(const void *ap, const void *bp)
{
	int a, b;
	
	a = *((int *) ap);
	b = *((int *) bp);

	if (a < b) return -1;
	else if (a > b) return 1;
	else return 0;
}
