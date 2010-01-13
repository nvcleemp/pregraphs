CFLAGS = -Wall -g -DMAXN=50

TARG = pregraphs$(SUFFIX)
OBJS = pregraphs.o nauty/nauty.o nauty/nautil.o nauty/naugraph.o

all: $(TARG)

$(TARG): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARG) $(OBJS)

clean:
	rm -f *.o $(TARG)
