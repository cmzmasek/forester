# Modified by Christian Zmasek. Use at your own risk.
#
# Generic Linux/Unix Makefile for PHYLIP 3.6.
#
# You should not need to change anything, though if you want you could
# change the first (noncomment) statement to some directory location
# that might be more useful.
# The compressed tar archive phylip.tar.Z when uncompressed and extracted
# puts the source code into a directory ./src, and also makes two other
# directories ./exe and ./doc for the final executables and the documentation
# files.  Only change the EXEDIR settings if you want something different
# from that structure.
# If it causes trouble in compiling, the CFLAGS statement below may also need
# to be changed.
#
#  To use the PHYLIP v3.6 Makefile, type
#      make install       to compile the whole package and install
#                           the executables in $(EXEDIR), and then
#                           remove the object files to save space
#      make all           to compile the whole package but not install it
#                           or remove the object files
#      make put           to move the executables into $(EXEDIR)
#      make clean         to remove all object files and executables from the
#                           current directory
#      make dnaml         to compile and link one program, (in this example,
#                           DnaML) and leave the executable and object files
#                           in the current directory (where the source code is).
#                           You will have to move the executable in to the
#                           executables directory (e.g. "mv dnaml ../exe")
#                           Note that the program name should be lower case.
# 
# ----------------------------------------------------------------------------
#  (Starting here is the section where you may want to change things)
# ----------------------------------------------------------------------------
#
# these are the statements we have been talking about:
# one of the reasons for changing them would be to put the executables
# on a different file system.
# The default configuration is to have within the overall PHYLIP
# directory three subdirectories:  "src" for source code, "exe" for the
# executables, and "doc" for the documentation files.
#
# the following specifies the directory where the executables will be placed
EXEDIR  = ../exe
#
# ----------------------------------------------------------------------------
#
# In the following statements (the ones that set CFLAGS, DFLAGS, LIBS
# and DLIBS, CC and DC) you should make sure each is set properly.
# Usually this will simply involve making sure that the proper statement
# has no "#" as its first character and that all other possibilities
# have "#" for their first character.
#
# ----------------------------------------------------------------------------
#
# This is the CFLAGS statement:
#
# if these statements say "-g"  and that causes trouble, remove the "-g"
# if the "-lX11" does not work (X compiling not working), you may want to
# remove that switch or fix it.
#
# Here are some possible CFLAGS statements:
#
#
#A minimal one
#CFLAGS =
#
# A basic one for debugging
#CFLAGS  = -g 
#
# An optimized one for gcc
CFLAGS  = -O3 -fomit-frame-pointer
#
# For some serious debugging using Gnu gcc
#
#CFLAGS=-g -Wall -Wmain -Wmissing-prototypes -Wreturn-type -Wstrict-prototypes  -Wunused -Werror -Wredundant-decls -Waggregate-return -Wcast-align -Wcomment
#
# For Digital Alpha systems with Compaq Tru64 Unix
# (however, be aware that this may cause floating-point problems in programs
#  like Dnaml owing to not using IEEE floating point standards).
#CFLAGS = -fast
#
# ----------------------------------------------------------------------------
#
# and here are some possible DFLAGS statements:
#
# A minimal one
#DFLAGS  = -DX $(CFLAGS)
#
# A basic one for debugging
#DFLAGS  = -g -DX
#
#
# For Gnu C++ for runs
#DFLAGS  = -I/usr/X11R6/include -O3 -DX -fomit-frame-pointer
#
# For Digital Alpha systems with Compaq Tru64 Unix
#DFLAGS = -DX -fast 
#
# for Linux with X Windows development packages installed
# or for MacOS X with X Windows installed
DFLAGS = $(CFLAGS) -DX -I/usr/X11R6/include
#
# ----------------------------------------------------------------------------
#
# These are the libraries for the CC and DC compiles, respectively
#
LIBS    = -lm
#
# if the Xlib library for the X windowing system is somewhere
# unexpected, you may have to change the path /usr/X11R6/lib in this one
#
# For gcc for Linux with X windows development packages installed
# or for MacOS X with X windows installed
DLIBS= -L/usr/X11R6/lib/ -lX11 -lXaw -lXt
#
# ----------------------------------------------------------------------------
#
# The next two assignments are the invocations of the compiler for the
# ordinary compiles and the tree-drawing programs, CC and DC
#
# This one specifies the "cc" C compiler
CC        = cc $(CFLAGS)
#
#  To use GCC instead, if it is not the compiler that "cc" invokes
#CC        = gcc $(CFLAGS)
#
# This one specifies the "cc" C compiler for the Draw programs
DC        = cc $(DFLAGS)
#
# To use GCC instead, if it is not the compiler that "cc" invokes
#DC        = gcc $(DFLAGS)
#
# ----------------------------------------------------------------------------
#  (After this point there should not be any reason to change anything)
# ----------------------------------------------------------------------------
#
#
# the list of programs
#
#PROGS   =          clique consense contml contrast dnacomp dnadist \
#                   dnainvar dnaml dnamlk dnamove dnapars dnapenny \
#                   dolmove dollop dolpenny factor fitch gendist kitsch \
#                   mix move neighbor pars penny proml promlk protdist \
#                   protpars restdist restml retree seqboot treedist \
#                   drawgram drawtree
                   
PROGS   =  consense fitch neighbor proml promlk protdist protpars seqboot                   
                   
#
# general commands
#

#
#  The first uses a symbol you are unlikely to type.  It is the one that
#  is executed if you just type "make".  It tells you how to use the
#  Makefile.
#
a1b2c3d4:
	@echo ""
	@echo " To use the PHYLIP v3.6 Makefile, type"
	@echo "     make all           to compile the whole package but not install it"
	@echo "                          or remove the object files"
	@echo " "

introduce:
	@echo "Building PHYLIP based on version 3.6 - Modified by Christian Zmasek. Use at your own risk."

all:        introduce $(PROGS)
	@echo "Finished compiling."
	@echo ""

#install:        all put clean
#	@echo "Done."
#	@echo ""

#put:
#	@echo "Installing PHYLIP v3.6 binaries in $(EXEDIR)"
#	@mkdir -p $(EXEDIR)
#	@cp $(PROGS) $(EXEDIR)
#	@echo "Installing font files in $(EXEDIR)"
#	@cp font* $(EXEDIR)
#	@echo "Finished installation."
#	@echo ""

#clean:
#	@echo "Removing object files to save space"
#	@rm -f *.o
#	@echo "Finished removing object files.  Now will remove"
#	@echo "executable files from the current directory, but not from the"
#	@echo "executables directory.  (If some are not here, the makefile"
#	@echo "will terminate with an error message but this is not a problem)"
#	@echo ""
#	@echo "Removing executables from this directory"
#	@rm -f $(PROGS)
#	@echo "Finished cleanup."
#	@echo ""

#
# compile the shared stuff
#

phylip.o:     phylip.c phylip.h
	$(CC) -c phylip.c

seq.o:        seq.c phylip.h seq.h
	$(CC) -c seq.c

#disc.o:   disc.c phylip.h disc.h
#	$(CC) -c disc.c

#discrete.o:        discrete.c discrete.h phylip.h
#	$(CC) -c discrete.c

#dollo.o:        dollo.c phylip.h dollo.h
#	$(CC) -c dollo.c

#wagner.o:        wagner.c phylip.h wagner.h
#	$(CC) -c wagner.c

dist.o:   dist.c phylip.h dist.h
	$(CC) -c dist.c

#cont.o:       cont.c cont.h phylip.h
#	$(CC) -c cont.c

#moves.o:   moves.c phylip.h moves.h
#	$(CC) -c moves.c

#
# compile the individual programs
#

#clique.o:       clique.c disc.h phylip.h
#	$(CC) -c clique.c

#clique:       clique.o clique.c disc.o disc.c phylip.o phylip.c disc.h phylip.h
#	$(CC) clique.o disc.o phylip.o $(LIBS) -o clique

cons.o:          cons.c cons.h phylip.h 
	$(CC) -c cons.c

consense.o:     consense.c cons.h phylip.h
	$(CC) -c consense.c

consense:     consense.o consense.c phylip.o phylip.c cons.o cons.c cons.h phylip.h
	$(CC) consense.o cons.o phylip.o $(LIBS) -o consense

#contml.o:       contml.c cont.h phylip.h
#	$(CC) -c contml.c

#contml:       contml.o contml.c cont.o cont.c phylip.o phylip.c cont.h phylip.h
#	$(CC) contml.o cont.o phylip.o $(LIBS) -o contml

#contrast.o:       contrast.c cont.h phylip.h
#	$(CC) -c contrast.c

#contrast:       contrast.o contrast.c cont.o cont.c phylip.o phylip.c cont.h phylip.h
#	$(CC) contrast.o cont.o phylip.o $(LIBS) -o contrast

#dnacomp.o:      dnacomp.c seq.h phylip.h
#	$(CC) -c dnacomp.c

#dnacomp:      dnacomp.o seq.o phylip.o dnacomp.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnacomp.o seq.o phylip.o $(LIBS) -o dnacomp

#dnadist.o:      dnadist.c seq.h phylip.h
#	$(CC) -c dnadist.c

#dnadist:      dnadist.o seq.o phylip.o dnadist.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnadist.o seq.o phylip.o $(LIBS) -o dnadist

#dnainvar.o:      dnainvar.c seq.h phylip.h
#	$(CC) -c dnainvar.c

#dnainvar:      dnainvar.o seq.o phylip.o dnainvar.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnainvar.o seq.o phylip.o $(LIBS) -o dnainvar

#dnaml.o:      dnaml.c seq.h phylip.h
#	$(CC) -c dnaml.c

#dnaml:      dnaml.o seq.o phylip.o dnaml.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnaml.o seq.o phylip.o $(LIBS) -o dnaml

#dnamlk.o:      dnamlk.c seq.h phylip.h
#	$(CC) -c dnamlk.c

#dnamlk:      dnamlk.o seq.o phylip.o dnamlk.c seq.c phylip.c
#	$(CC) dnamlk.o seq.o phylip.o $(LIBS) -o dnamlk

#dnamove.o:      dnamove.c seq.h moves.h phylip.h
#	$(CC) -c dnamove.c

#dnamove:      dnamove.o seq.o moves.o phylip.o dnamove.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnamove.o seq.o moves.o phylip.o $(LIBS) -o dnamove

#dnapenny.o:      dnapenny.c seq.h phylip.h
#	$(CC) -c dnapenny.c

#dnapenny:      dnapenny.o seq.o phylip.o dnapenny.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnapenny.o seq.o phylip.o $(LIBS) -o dnapenny

#dnapars.o:      dnapars.c seq.h phylip.h
#	$(CC) -c dnapars.c

#dnapars:      dnapars.o seq.o phylip.o dnapars.c seq.c phylip.c seq.h phylip.h
#	$(CC) dnapars.o seq.o phylip.o $(LIBS) -o dnapars

#dolmove.o:       dolmove.c disc.h moves.h dollo.h phylip.h
#	$(CC) -c dolmove.c

#dolmove:       dolmove.o disc.o moves.o dollo.o phylip.o dolmove.c disc.c moves.c dollo.c phylip.c disc.h moves.h dollo.h phylip.h
#	$(CC) dolmove.o disc.o moves.o dollo.o phylip.o $(LIBS) -o dolmove

#dollop.o:       dollop.c disc.h dollo.h phylip.h
#	$(CC) -c dollop.c

#dollop:       dollop.o disc.o dollo.o phylip.o dollop.c disc.c dollo.c phylip.c disc.h dollo.h phylip.h
#	$(CC) dollop.o disc.o dollo.o phylip.o $(LIBS) -o dollop

#dolpenny.o:       dolpenny.c disc.h dollo.h phylip.h
#	$(CC) -c dolpenny.c

#dolpenny:       dolpenny.o disc.o dollo.o phylip.o dolpenny.c disc.c dollo.c phylip.c disc.h dollo.h phylip.h
#	$(CC) dolpenny.o disc.o dollo.o phylip.o $(LIBS) -o dolpenny

#draw.o:   draw.c draw.h phylip.h 
#	$(DC) -c draw.c

#draw2.o:   draw2.c draw.h phylip.h 
#	$(DC) -c draw2.c

#drawgram.o:     drawgram.c draw.h phylip.h
#	$(DC) -c drawgram.c

#drawgram:     drawgram.o draw.o draw2.o phylip.o drawgram.c draw.c draw2.c draw.h phylip.h
#	$(DC) $(DLIBS) draw.o draw2.o drawgram.o phylip.o $(LIBS) -o drawgram

#drawtree.o:     drawtree.c draw.h phylip.h
#	$(DC) -c drawtree.c

#drawtree:     drawtree.o draw.o draw2.o phylip.o drawtree.c draw.c draw2.c draw.h phylip.h
#	$(DC) $(DLIBS) draw.o draw2.o drawtree.o phylip.o $(LIBS) -o drawtree

#factor.o:       factor.c phylip.h
#	$(CC) -c factor.c

#factor:       factor.o phylip.o factor.c phylip.c phylip.h
#	$(CC) factor.o phylip.o $(LIBS) -o factor

fitch.o:        fitch.c dist.h phylip.h
	$(CC) -c fitch.c

fitch:        fitch.o dist.o phylip.o fitch.c dist.c phylip.c dist.h phylip.h
	$(CC) fitch.o dist.o phylip.o $(LIBS) -o fitch

#gendist.o:      gendist.c phylip.h
#	$(CC) -c gendist.c

#gendist:      gendist.o phylip.o gendist.c phylip.c phylip.h
#	$(CC) gendist.o phylip.o $(LIBS) -o gendist

#kitsch.o:        kitsch.c dist.h phylip.h
#	$(CC) -c kitsch.c

#kitsch:        kitsch.o dist.o phylip.o kitsch.c dist.c phylip.c dist.h phylip.h
#	$(CC) kitsch.o dist.o phylip.o $(LIBS) -o kitsch

#mix.o:        mix.c disc.h wagner.h phylip.h
#	$(CC) -c mix.c

#mix:        mix.o disc.o wagner.o phylip.o mix.c disc.c wagner.c phylip.c disc.h wagner.h phylip.h
#	$(CC) mix.o disc.o wagner.o phylip.o $(LIBS) -o mix

#move.o:        move.c disc.h moves.h wagner.h phylip.h
#	$(CC) -c move.c

#move:        move.o disc.o moves.o wagner.o phylip.o move.c disc.c moves.c wagner.c phylip.c disc.h moves.h wagner.h phylip.h
#	$(CC) move.o disc.o moves.o wagner.o phylip.o $(LIBS) -o move

neighbor.o:        neighbor.c dist.h phylip.h
	$(CC) -c neighbor.c

neighbor:        neighbor.o dist.o phylip.o neighbor.c dist.c phylip.c dist.h phylip.h
	$(CC) neighbor.o dist.o phylip.o $(LIBS) -o neighbor

#pars.o:   pars.c discrete.h phylip.h
#	$(CC) -c pars.c

#pars: pars.o pars.c discrete.o discrete.c phylip.o phylip.c discrete.h phylip.h
#	$(CC) pars.o discrete.o phylip.o $(LIBS) -o pars

#penny.o:  penny.c disc.h wagner.h phylip.h
#	$(CC) -c penny.c

#penny:  penny.o disc.o wagner.o phylip.o penny.c disc.c wagner.c disc.h wagner.h phylip.h
#	$(CC) penny.o disc.o wagner.o  phylip.o $(LIBS) -o penny

proml.o:      proml.c seq.h phylip.h
	$(CC) -c proml.c

proml:      proml.o seq.o phylip.o proml.c seq.c phylip.c seq.h phylip.h
	$(CC) proml.o seq.o phylip.o $(LIBS) -o proml

promlk.o:      promlk.c seq.h phylip.h
	$(CC) -c promlk.c

promlk:      promlk.o seq.o phylip.o promlk.c seq.c phylip.c
	$(CC) promlk.o seq.o phylip.o $(LIBS) -o promlk

protdist.o:      protdist.c seq.h phylip.h
	$(CC) -c protdist.c

protdist:      protdist.o seq.o phylip.o protdist.c seq.c phylip.c seq.h phylip.h
	$(CC) protdist.o seq.o phylip.o $(LIBS) -o protdist

protpars.o: protpars.c seq.h phylip.h
	$(CC) -c protpars.c

protpars: protpars.o seq.o phylip.o protpars.c seq.c phylip.c seq.h phylip.h
	$(CC) protpars.o seq.o phylip.o $(LIBS) -o protpars

#restdist.o: restdist.c seq.h phylip.h
#	 $(CC) -c restdist.c

#restdist: restdist.o seq.o phylip.o restdist.c seq.c phylip.c seq.h phylip.h
#	$(CC) restdist.o seq.o phylip.o $(LIBS) -o restdist

#restml.o: restml.c seq.h phylip.h
#	 $(CC) -c restml.c

#restml: restml.o seq.o phylip.o restml.c seq.c phylip.c seq.h phylip.h
#	$(CC) restml.o seq.o phylip.o $(LIBS) -o restml

#retree.o:       retree.c moves.h phylip.h
#	$(CC) -c retree.c

#retree:       retree.o moves.o phylip.o retree.c moves.c phylip.c moves.h phylip.h
#	$(CC) retree.o moves.o phylip.o $(LIBS) -o retree

seqboot.o:      seqboot.c phylip.h
	$(CC) -c seqboot.c

seqboot:      seqboot.o seq.o phylip.o seqboot.c seq.c phylip.c seq.h phylip.h
	$(CC) seqboot.o seq.o phylip.o $(LIBS) -o seqboot

#treedist.o:     treedist.c cons.h phylip.h
#	$(CC) -c treedist.c

#treedist:     treedist.o phylip.o cons.o
#	$(CC) treedist.o cons.o phylip.o $(LIBS) -o treedist
