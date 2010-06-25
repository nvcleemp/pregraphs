CFLAGS = -Wall -g -DMAXN=50

TARG = pregraphs$(SUFFIX)
OBJS = pregraphs.o snarkhunter.o nauty/nauty.o nauty/nautil.o nauty/naugraph.o nauty/nausparse.o nauty/nautinv.o

all: $(TARG) translator multi2simple printmulticode multicode2dreadnaut admissable_c

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

$(TARG): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARG) $(OBJS)

clean:
	rm -f *.o $(TARG) translator multi2simple printmulticode multicode2dreadnaut
