/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1999 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* revcomp.c
 * 
 * Reverse complement of a IUPAC character string
 * RCS $Id: revcomp.c,v 1.1.1.1 2005/03/22 08:34:16 cmzmasek Exp $
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


char *
revcomp(char *comp, char *seq)
{
  long  bases;
  char *bckp, *fwdp;
  int   idx;
  long  pos;
  int   c;

  if (comp == NULL) return NULL;
  if (seq == NULL)  return NULL;
  bases = strlen(seq);

  fwdp = comp;
  bckp = seq + bases -1;
  for (pos = 0; pos < bases; pos++)
    {
      c = *bckp;
      c = sre_toupper(c);
      for (idx = 0; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);
      if (idx == IUPACSYMNUM)
	{
	  Warn("Can't reverse complement an %c, pal. Using N.", c);
	  *fwdp = 'N';
	}
      else
	*fwdp = iupac[idx].symcomp;
      if (islower((int) *bckp)) *fwdp = (char) sre_tolower((int) *fwdp);
      fwdp++;
      bckp--;
    }
  *fwdp = '\0';
  return comp;
}
  
