/*
 * sched.h
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


#ifndef SCHED_H
#define SCHED_H
#ifndef SCHEDTEST
#   include "util.h" 
#else
    typedef unsigned long int uli; 
#endif


typedef struct sched_t{
   uli    truetasks;
   uli    alltasks;
   uli    numtasks;
   uli    minchunk;
   int    numprocs;
   int    delta;
   double ddelta;
   int    overhead;
   int    rest;
   int    nconst;
   double fconst;
   double lconst;
   double kconst;
   int    inited;
} schedtype;

void num2quart(uli qnum, int *a, int *b, int *c, int *d);
uli numquarts(int maxspc);
uli quart2num (int a, int b, int c, int d);

void printsched(schedtype sch);
void initsched(schedtype *sch, uli tasks, int procs, uli minchunk);
uli sc(schedtype *sch);
uli gss(schedtype *sch);
uli sgss(schedtype *sch);
uli tss(schedtype *sch);

#endif /* SCHED_H */
