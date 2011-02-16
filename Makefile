#
# Makefile for pregraphs
#

SHELL = /bin/sh

# Compiling executing program with DWORDSIZE=32 is slightly faster, 
# but limits the order of the graphs to 32.
CC32 = gcc -DWORDSIZE=32 -DMAXN=WORDSIZE -DSNARKHUNTERMAIN=sh_nomain
CC64 = gcc -DWORDSIZE=64 -DMAXN=WORDSIZE -DSNARKHUNTERMAIN=sh_nomain
CFLAGS = -O4 -Wall
COMPLETE = bipartite has3edgecolouring pgfilter

all : 32bit admissable_c c4cover

complete: all $(COMPLETE)

32bit: pregraphs

64bit : pregraphs-64

profile : pregraphs-profile

debug : pregraphs-debug

pregraphs: pregraphs.c snarkhunter.c  nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c
	${CC32} $(CFLAGS) pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs

pregraphs-64: pregraphs.c snarkhunter.c  nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c
	${CC64} $(CFLAGS) pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs-64

pregraphs-profile: pregraphs.c snarkhunter.c  nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c
	${CC32} -Wall -pg -g pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs-profile

pregraphs-debug: pregraphs.c snarkhunter.c  nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c
	${CC32} -Wall -rdynamic -g pregraphs.c snarkhunter.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c -o pregraphs-debug

admissable_c: admissable_c.c
	$(CC) $(CFLAGS) -o admissable_c admissable_c.c

c4cover: c4cover.c
	$(CC) $(CFLAGS) -o c4cover c4cover.c

bipartite: bipartite.c
	$(CC) $(CFLAGS) -o bipartite bipartite.c

has3edgecolouring: has3edgecolouring.c
	$(CC) $(CFLAGS) -o has3edgecolouring has3edgecolouring.c

pgfilter: pgfilter.c
	$(CC) $(CFLAGS) -o pgfilter pgfilter.c

clean:
	rm -f pregraphs pregraphs-* admissable_c c4cover $(COMPLETE)
