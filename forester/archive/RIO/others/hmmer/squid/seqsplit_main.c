/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1999 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/


/* seqsplit_main.c
 * SRE, Mon Sep 25 11:43:58 2000
 * 
 * Split sequences into smaller chunks of defined size and overlap;
 * output a FASTA file.
 *
 * Limitations:
 *   still working in 32 bits -- no sequence can be more than 2 GB
 *   in size.
 * CVS $Id: seqsplit_main.c,v 1.1.1.1 2005/03/22 08:34:26 cmzmasek Exp $
 */

#include <stdio.h>
#include <string.h>
#include "squid.h"
#include "msa.h"

static char banner[] = "seqsplit - split seqs into chunks of defined size and overlap";

static char usage[]  = "\
Usage: seqsplit [-options] <seqfile>\n\
  Available options:\n\
  -h        : help; display usage and version\n\
  -o <file> : output the new FASTA file to <file>\n\
";  

static char experts[] = "\
  --informat <s> : specify sequence file format <s>\n\
  --length <n>   : set max length of each unique seq frag to <n>\n\
  --overlap <n>  : set overlap length to <n> (total frag size = length+overlap)\n\
";

struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE },    
  { "-o", TRUE, sqdARG_STRING },    
  { "--informat", FALSE, sqdARG_STRING },
  { "--length",   FALSE, sqdARG_INT },
  { "--overlap",  FALSE, sqdARG_INT },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


int
main(int argc, char **argv)
{
  char     *seqfile;            /* name of sequence file     */
  char     *outfile;		/* name of output file       */
  SQFILE   *dbfp;		/* open sequence file        */
  FILE     *ofp;		/* open output file          */
  int       fmt;		/* format of seqfile         */
  char     *seq;		/* sequence                  */
  SQINFO    sqinfo;             /* extra info about sequence */
  char     *seqfrag;		/* space for a seq fragment  */
  int       fraglength;		/* length of unique seq per frag     */
  int       overlap;            /* length of overlap. frags are fraglength+overlap*/
  char      seqname[256];	/* renamed fragment, w/ coord info */
  int       num;		/* number of this fragment   */
  int       pos;		/* position in a sequence */
  int       len;		/* length of a fragment   */
  char     *desc;	
  
  int       nseqs;		/* total number of sequences */
  int       nsplit;		/* number of seqs that get split */
  int       nnewfrags;		/* total number of new fragments */

  char  *optname;
  char  *optarg;
  int    optind;

  /***********************************************
   * Parse command line
   ***********************************************/

  fmt       = SQFILE_UNKNOWN;	/* default: autodetect      */
  fraglength = 100000;
  overlap   = 1000;
  outfile   = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, 
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "-o")         == 0)  outfile    = optarg;
      else if (strcmp(optname, "--length")   == 0)  fraglength = atoi(optarg);
      else if (strcmp(optname, "--overlap")  == 0)  overlap    = atoi(optarg);
      else if (strcmp(optname, "--informat") == 0) {
	fmt = String2SeqfileFormat(optarg);
	if (fmt == SQFILE_UNKNOWN) 
	  Die("unrecognized sequence file format \"%s\"", optarg);
      }
      else if (strcmp(optname, "-h") == 0) {
	Banner(stdout, banner);
	puts(usage);
	puts(experts);
        exit(EXIT_SUCCESS);
      }
    }

  if (argc - optind != 1) Die("%s\n", usage);
  seqfile = argv[argc-1];

  seqfrag = MallocOrDie(sizeof(char) * (fraglength+overlap));
  seqfrag[fraglength+overlap] = '\0';

  /***********************************************
   * Read the file.
   ***********************************************/

  if (outfile == NULL)  ofp = stdout; 
  else {
    if ((ofp = fopen(outfile, "w")) == NULL)
      Die("Failed to open output sequence file %s for writing", outfile);
  }

  if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
    Die("Failed to open sequence file %s for reading", seqfile);
  
  nseqs = nsplit = nnewfrags = 0;
  while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo))
    {
      nseqs++;
      if (sqinfo.flags & SQINFO_DESC) desc = sqinfo.desc;
      else desc = NULL;

      if (sqinfo.len <= fraglength+overlap) {
	WriteSimpleFASTA(ofp, seq, sqinfo.name, desc);
	continue;
      }
      
      num = 1;
      nsplit++;
      for (pos = 0; pos < sqinfo.len; pos += fraglength)
	{
	  if (sqinfo.len - pos <= overlap) continue;
	  strncpy(seqfrag, seq+pos, fraglength+overlap);
	  len = strlen(seqfrag);
	  sprintf(seqname, "%s/frag%d/%d-%d", 
		  sqinfo.name, num, pos+1, pos+len);
	  WriteSimpleFASTA(ofp, seqfrag, seqname, desc);
	  nnewfrags++;
	  num ++;
	}
      FreeSequence(seq, &sqinfo);
    }
  SeqfileClose(dbfp);
  if (outfile != NULL) fclose(ofp);

  printf("Total # of seqs:         %d\n", nseqs);
  printf("Affected by splitting:   %d\n", nsplit);
  printf("New # of seqs:           %d\n", nseqs-nsplit + nnewfrags);

  return 0;
}
