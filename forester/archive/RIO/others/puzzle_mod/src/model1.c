/*
 * model1.c
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


/* definitions */
#define EXTERN extern

/* prototypes */
#include <stdio.h>
#include "util.h"
#include "ml.h"

/* number of states of the selected model */ 
int gettpmradix()
{
	if (data_optn == 0) { /* nucleotides */
		if (nuc_optn) return 4;
		if (SH_optn) return 16;
	} else if (data_optn == 1) { /* amino acids */
		return 20;
	} else { /* two-state model */
		return 2;
	}
	return 1;
}

/* relative transition frequencies */
void rtfdata(dmatrix q, double *f)
{
	double alp, alpy, alpr;
	int i, j;

	if (data_optn == 0)
	{ /* nucleotides */

		if (nuc_optn)
		{ /* 4x4 nucleotides */
			alp = 2.0*TSparam;
			alpr = (alp * 2.0) / (YRparam + 1.0);
			alpy = YRparam * alpr;
			
			q[0][1] = 1; q[0][2] = alpr; q[0][3] = 1;
			q[1][2] = 1; q[1][3] = alpy;
			q[2][3] = 1;

			f[0] = 0.25; f[1] = 0.25; f[2] = 0.25; f[3] = 0.25;
		}
		
		if (SH_optn)
		{ /* 16x16 nucleotides */
		
			alp = 2.0*TSparam;
	
			q[0][1] = 1; q[0][2] = alp; q[0][3] = 1; q[0][4] = 1;
			q[0][5] = 0; q[0][6] = 0; q[0][7] = 0; q[0][8] = alp;
			q[0][9] = 0; q[0][10] = 0; q[0][11] = 0; q[0][12] = 1;
			q[0][13] = 0; q[0][14] = 0; q[0][15] = 0;
			
			q[1][2] = 1; q[1][3] = alp; q[1][4] = 0; q[1][5] = 1;
			q[1][6] = 0; q[1][7] = 0; q[1][8] = 0; q[1][9] = alp;
			q[1][10] = 0; q[1][11] = 0; q[1][12] = 0; q[1][13] = 1;
			q[1][14] = 0; q[1][15] = 0;
			
			q[2][3] = 1; q[2][4] = 0; q[2][5] = 0; q[2][6] = 1;
			q[2][7] = 0; q[2][8] = 0; q[2][9] = 0; q[2][10] = alp;
			q[2][11] = 0; q[2][12] = 0; q[2][13] = 0; q[2][14] = 1;
			q[2][15] = 0;
			
			q[3][4] = 0; q[3][5] = 0; q[3][6] = 0; q[3][7] = 1;
			q[3][8] = 0; q[3][9] = 0; q[3][10] = 0; q[3][11] = alp;
			q[3][12] = 0; q[3][13] = 0; q[3][14] = 0; q[3][15] = 1;
			
			q[4][5] = 1; q[4][6] = alp; q[4][7] = 1; q[4][8] = 1;
			q[4][9] = 0; q[4][10] = 0; q[4][11] = 0; q[4][12] = alp;
			q[4][13] = 0; q[4][14] = 0; q[4][15] = 0;
			
			q[5][6] = 1; q[5][7] = alp; q[5][8] = 0; q[5][9] = 1;
			q[5][10] = 0; q[5][11] = 0; q[5][12] = 0; q[5][13] = alp;
			q[5][14] = 0; q[5][15] = 0;
			
			q[6][7] = 1; q[6][8] = 0; q[6][9] = 0; q[6][10] = 1;
			q[6][11] = 0; q[6][12] = 0; q[6][13] = 0; q[6][14] = alp;
			q[6][15] = 0;

			q[7][8] = 0; q[7][9] = 0; q[7][10] = 0; q[7][11] = 1;
			q[7][12] = 0; q[7][13] = 0; q[7][14] = 0; q[7][15] = alp;
			
			q[8][9] = 1; q[8][10] = alp; q[8][11] = 1; q[8][12] = 1;
			q[8][13] = 0; q[8][14] = 0; q[8][15] = 0;
			
			q[9][10] = 1; q[9][11] = alp; q[9][12] = 0; q[9][13] = 1;
			q[9][14] = 0; q[9][15] = 0;
			
			q[10][11] = 1; q[10][12] = 0; q[10][13] = 0; q[10][14] = 1;
			q[10][15] = 0;
			
			q[11][12] = 0; q[11][13] = 0; q[11][14] = 0; q[11][15] = 1;
			
			q[12][13] = 1; q[12][14] = alp; q[12][15] = 1;
			
			q[13][14] = 1; q[13][15] = alp;
			
			q[14][15] = 1;

			
			for (i = 0; i < 16; i++) f[i] = 0.0625;
		}
	}
	else if (data_optn == 1)
	{ /* amino acids */
		if (Dayhf_optn) /* Dayhoff model */
		{
			dyhfdata(q, f);			
		}
		else if (Jtt_optn) /* JTT model */
		{
			jttdata(q, f);
		}
		else if (blosum62_optn) /* BLOSUM 62 model */
		{
			blosum62data(q, f);
		}
		else if (mtrev_optn) /* mtREV model */
		{
			mtrevdata(q, f);
		} 
		else if (cprev_optn) /* cpREV model */
		{
			cprev45data(q, f);
		} 
		else if (vtmv_optn) /* VT model */
		{
			vtmvdata(q, f);
		}
		else /* if (wag_optn) */ /* WAG model */
		{
			wagdata(q, f);
		} 
		
	}
	else /* two-state model */
	{
		q[0][1] = 1.0;
		
		f[0] = 0.5; f[1] = 0.5;
	}
	
	/* fill matrix from upper triangle */
	for (i = 0; i < tpmradix; i++)
	{
		q[i][i] = 0.0;
		for (j = i+1; j < tpmradix; j++)
		{
			q[j][i] = q[i][j];
		}
	}
}

/* transform letter codes to state numbers */
int code2int(cvector c)
{	if (data_optn == 0) { /* nucleotides */
		if (nuc_optn) { /* 4x4 */
			switch (c[0]) {
				case 'A': return 0;
				case 'C': return 1;
				case 'G': return 2;
				case 'T': return 3;
				case 'U': return 3;
				default : return 4;
			}
		}
		if (SH_optn) { /* 16x16 */
			if (c[0] == 'A') {
				switch (c[1]) {
					case 'A': return 0; /* AA */
					case 'C': return 1; /* AC */
					case 'G': return 2; /* AG */
					case 'T': return 3; /* AT */
					case 'U': return 3; /* AT */
					default:  return 16;
				}
			}
			if (c[0] == 'C') {
				switch (c[1]) {
					case 'A': return 4; /* CA */
					case 'C': return 5; /* CC */
					case 'G': return 6; /* CG */
					case 'T': return 7; /* CT */
					case 'U': return 7; /* CT */
					default:  return 16;
				}
			}
			if (c[0] == 'G') {
				switch (c[1]) {
					case 'A': return 8; /* GA */
					case 'C': return 9; /* GC */
					case 'G': return 10; /* GG */
					case 'T': return 11; /* GT */
					case 'U': return 11; /* GT */
					default:  return 16;
				}
			}
			if (c[0] == 'T' || c[0] == 'U') {
				switch (c[1]) {
					case 'A': return 12; /* TA */
					case 'C': return 13; /* TC */
					case 'G': return 14; /* TG */
					case 'T': return 15; /* TT */
					case 'U': return 15; /* TT */
					default:  return 16;
				}
			}
			return 16;
		}
	} else if (data_optn == 1) { /* amino acids */
		switch (c[0]) {
			case 'A': return 0;
			case 'C': return 4;
			case 'D': return 3;
			case 'E': return 6;
			case 'F': return 13;
			case 'G': return 7;
			case 'H': return 8;
			case 'I': return 9;
			case 'K': return 11;
			case 'L': return 10;
			case 'M': return 12;
			case 'N': return 2;
			case 'P': return 14;
			case 'Q': return 5;
			case 'R': return 1;
			case 'S': return 15;
			case 'T': return 16;
			case 'V': return 19;
			case 'W': return 17;
			case 'Y': return 18;
			default : return 20;
		}
	} else { /* two-state model */
		switch (c[0]) {
			case '0': return 0;
			case '1': return 1;
			default : return 2;
		}
	}
	return 0;
}

/* return letter code belonging to state number */
char *int2code(int s)
{
	if (data_optn == 0) { /* nucleotides */
		if (nuc_optn) { /* 4x4 */
			switch (s) {
				case 0: return "A";
				case 1: return "C";
				case 2: return "G";
				case 3: return "T";
				default : return "?";
			}
		}
		if (SH_optn) { /* 16x16 */
			switch (s) {
				case 0: return "AA";
				case 1: return "AC";
				case 2: return "AG";
				case 3: return "AT";
				case 4: return "CA";
				case 5: return "CC";
				case 6: return "CG";
				case 7: return "CT";
				case 8: return "GA";
				case 9: return "GC";
				case 10: return "GG";
				case 11: return "GT";
				case 12: return "TA";
				case 13: return "TC";
				case 14: return "TG";
				case 15: return "TT";
				default : return "??";
			}
		}
	} else if (data_optn == 1) { /* amino acids */
		switch (s) {
			case 0: return "A";
			case 1: return "R";
			case 2: return "N";
			case 3: return "D";
			case 4: return "C";
			case 5: return "Q";
			case 6: return "E";
			case 7: return "G";
			case 8: return "H";
			case 9: return "I";
			case 10: return "L";
			case 11: return "K";
			case 12: return "M";
			case 13: return "F";
			case 14: return "P";
			case 15: return "S";
			case 16: return "T";
			case 17: return "W";
			case 18: return "Y";
			case 19: return "V";
			default : return "?";
		}
	} else { /* two-state model */
		switch (s) {
			case 0: return "0";
			case 1: return "1";
			default : return "?";
		}
	}
	return "?";
}
