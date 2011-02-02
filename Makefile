#
# Makefile for pregraphs
#

SHELL = /bin/sh

# Compiling executing program with DWORDSIZE=32 is slightly faster, 
# but limits the order of the graphs to 32.
CC32 = gcc -DWORDSIZE=32 -DMAXN=WORDSIZE -DSNARKHUNTERMAIN=sh_nomain
CC64 = gcc -DWORDSIZE=64 -DMAXN=WORDSIZE -DSNARKHUNTERMAIN=sh_nomain
CFLAGS = -O4 -Wall

all :
	${CC32} $(CFLAGS) pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs

64bit :
	${CC64} $(CFLAGS) pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs-64

profile :
	${CC32} -Wall -pg -g pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs-profile

debug :
	${CC32} -Wall -rdynamic -g pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs-debug

translator: 3regpregraphtranslator.c 3regpregraphtranslator.h
	$(CC) $(CFLAGS) -o translator 3regpregraphtranslator.c

multi2simple: multi2simple.c multi2simple.h
	$(CC) $(CFLAGS) -o multi2simple multi2simple.c

printmulticode: printmulticode.c printmulticode.h
	$(CC) $(CFLAGS) -o printmulticode printmulticode.c

multicode2dreadnaut: multicode2dreadnaut.c multicode2dreadnaut.h
	$(CC) $(CFLAGS) -o multicode2dreadnaut multicode2dreadnaut.c

admissable_c: admissable_c.c
	$(CC) $(CFLAGS) -o admissable_c admissable_c.c

c4cover: c4cover.c
	$(CC) $(CFLAGS) -o c4cover c4cover.c

bipartite: bipartite.c
	$(CC) $(CFLAGS) -o bipartite bipartite.c

pgfilter: pgfilter.c
	$(CC) $(CFLAGS) -o pgfilter pgfilter.c

