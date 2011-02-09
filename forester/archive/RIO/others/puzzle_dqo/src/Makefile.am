bin_PROGRAMS = puzzle

puzzle_SOURCES = gamma.c ml1.c ml2.c ml3.c model1.c model2.c puzzle1.c puzzle2.c util.c ml.h util.h puzzle.h gamma.h 
puzzle_LDADD = sgamma.o sml1.o sml2.o sml3.o smodel1.o smodel2.o spuzzle1.o spuzzle2.o sutil.o

SDEFS =
SCFLAGS = 
SLDFLAGS = @LIBS@

SCOMPILE = $(CC) $(SDEFS) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(SCFLAGS) $(CFLAGS) 
SCCLD = $(CC)
SLINK = $(SCCLD) $(AM_CFLAGS) $(CFLAGS) $(SLDFLAGS) $(LDFLAGS)



PCC = @MPICC@
PDEFS = -DPARALLEL
PCFLAGS = 
PLDFLAGS = @LIBS@ @MPILIBS@



 
puzzle: $(puzzle_LDADD) $(puzzle_SOURCES)
	$(SLINK) $(puzzle_LDADD) -o $@

sml1.o: ml1.c ml.h util.h
	$(SCOMPILE) -c ml1.c && mv ml1.o $@
sml2.o: ml2.c ml.h util.h
	$(SCOMPILE) -c ml2.c && mv ml2.o $@
sml3.o: ml3.c ml.h util.h gamma.h
	$(SCOMPILE) -c ml3.c && mv ml3.o $@
smodel1.o: model1.c ml.h util.h
	$(SCOMPILE) -c model1.c && mv model1.o $@
smodel2.o: model2.c ml.h util.h
	$(SCOMPILE) -c model2.c && mv model2.o $@
spuzzle1.o: puzzle1.c ml.h util.h puzzle.h gamma.h ppuzzle.h
	$(SCOMPILE) -c puzzle1.c && mv puzzle1.o $@
spuzzle2.o: puzzle2.c ml.h util.h puzzle.h ppuzzle.h
	$(SCOMPILE) -c puzzle2.c && mv puzzle2.o $@
sutil.o: util.c util.h
	$(SCOMPILE) -c util.c && mv util.o $@
sgamma.o: gamma.c gamma.h util.h
	$(SCOMPILE) -c gamma.c && mv gamma.o $@





