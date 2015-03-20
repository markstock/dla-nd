# Makefile for dla-nd, an arbitrary-dimensional particle system
#
# Eventually support building win, mac, or Linux executables
# cc -O2 -o stickkit stickkit.c -lm
# /usr/local/bin/i386-mingw32-gcc -O2 -o stickkit.exe stickkit.c -lm


CC=gcc

# the neccessary c files
#CFILES = main.c ndtree.c density.c setup.c writeout.c readin.c
#CFILES = main.c ndtree.c setup.c writeout.c readin.c
CFILES = main.c ndtree.c setup.c inout.c
HFILES = structs.h

# To make a normal version 
#CFLAGS = -lm -lpng -O1 -funroll-loops
#CFLAGS = -lm -lpng -pg
#CFLAGS = -lm -lpng -pg -ggdb -Wall
CFLAGS=-O3 -funroll-loops -ffast-math -Wall -std=c99
LIBS=-lm -lpng
#LIBS=-I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include -L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib -lm -lpng


all: dla-2d dla-3d dla-4d
all: dla-3d

dla-1d: $(CFILES) $(HFILES) Makefile
	$(CC) -o dla-1d -DDIM=1 -DNCHILD=2 $(CFLAGS) $(CFILES) $(LIBS)
	@echo "dla-1d made"

dla-2d: $(CFILES) $(HFILES) Makefile
	$(CC) -o dla-2d -DDIM=2 -DNCHILD=4 $(CFLAGS) $(CFILES) $(LIBS)
	@echo "dla-2d made"

dla-3d: $(CFILES) $(HFILES) Makefile
	$(CC) -o dla-3d -DDIM=3 -DNCHILD=8 $(CFLAGS) $(CFILES) $(LIBS)
	@echo "dla-3d made"

dla-4d: $(CFILES) $(HFILES) Makefile
	$(CC) -o dla-4d -DDIM=4 -DNCHILD=16 $(CFLAGS) $(CFILES) $(LIBS)
	@echo "dla-4d made"

dla-5d: $(CFILES) $(HFILES) Makefile
	$(CC) -o dla-5d -DDIM=5 -DNCHILD=32 $(CFLAGS) $(CFILES) $(LIBS)
	@echo "dla-5d made"

lint:
	lint -abchp $(CFILES)

clean:
	rm -f *.o *.a gmon* a.out dla-?d
