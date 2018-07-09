CFLAGS   = -O -g -Wall -fPIC -ansi -lm
ROOT = $(shell $(ROOTSYS)/bin/root-config --glibs) $(shell $(ROOTSYS)/bin/root-config --cflags)

all: a2a4

a2a4: a2a4.c a2a4.h
	g++ a2a4.c $(CFLAGS) $(ROOT) -o a2a4
clean:
	rm -rf *~ *.o a2a4 *tmpdatafile*
