#
# Makefile for pregraphs
#

SHELL = /bin/sh

# Compiling executing program with DWORDSIZE=32 is slightly faster, 
# but limits the order of the graphs to 32.
CC32 = gcc -DWORDSIZE=32 -DMAXN=WORDSIZE -DSNARKHUNTERMAIN=sh_nomain -DMINIBAUM_NO_MAIN
CC64 = gcc -DWORDSIZE=64 -DMAXN=WORDSIZE -DSNARKHUNTERMAIN=sh_nomain -DMINIBAUM_NO_MAIN
CFLAGS = -O4 -Wall
COMPLETE = pregraphs pregraphs-64 pregraphs-profile pregraphs-debug admissable_c c4cover bipartite has3edgecolouring pgfilter
SOURCES = pregraphs.c pregraphs.h util.h snarkhunter.c snarkhunter.h admissable_c.c\
          admissable_c.h c4cover.c c4cover.h bipartite.c bipartite.h has3edgecolouring.c\
          has3edgecolouring.h pgfilter.c pgfilter.h Makefile COPYRIGHT.txt LICENSE.txt
PREGRAPHS_SOURCES = pregraphs.c snarkhunter.c minibaum5.c nauty/nautil.c nauty/nausparse.c nauty/naugraph.c nauty/nauty.c

all : 32bit

complete: $(COMPLETE)

32bit: pregraphs

64bit : pregraphs-64

profile : pregraphs-profile

debug : pregraphs-debug

pregraphs: $(PREGRAPHS_SOURCES)
	${CC32} $(CFLAGS) $(PREGRAPHS_SOURCES) -o pregraphs

pregraphs-64: $(PREGRAPHS_SOURCES)
	${CC64} $(CFLAGS) $(PREGRAPHS_SOURCES) -o pregraphs-64

pregraphs-profile: $(PREGRAPHS_SOURCES)
	${CC32} -Wall -pg -g $(PREGRAPHS_SOURCES) -o pregraphs-profile

pregraphs-debug: $(PREGRAPHS_SOURCES)
	${CC32} -Wall -rdynamic -g $(PREGRAPHS_SOURCES) -o pregraphs-debug

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

sources: pregraphs-sources.zip pregraphs-sources.tar.gz

pregraphs-sources.zip: $(SOURCES)
	zip pregraphs-sources $(SOURCES)

pregraphs-sources.tar.gz: $(SOURCES)
	tar czf pregraphs-sources.tar.gz $(SOURCES)

clean:
	rm -f $(COMPLETE)
