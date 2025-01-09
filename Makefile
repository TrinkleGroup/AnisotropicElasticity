# This is the Makefile for the anisotropic codes

CCPP = g++
CC = gcc

INCLUDE = .
LIBS = .

CFLAGS = -O5
CPPFLAGS = -O5

# bcc map removed, as well as nnpair.H drawfig.H
TARGET = make-slab anisotropic anisotropic-xyz anisotropic-xyz-ref \
         anisotropic-xyz-ref-outputstrain anisotropic-xyz-strain
INCLUDES = cell.H dcomp.H drawfig.H elastic.H integrate.H io.H matrix.H nnpair.H slab.H

all: ${TARGET}

make-slab: make-slab.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

anisotropic: anisotropic.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

anisotropic-xyz: anisotropic-xyz.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm
	
anisotropic-xyz-ref: anisotropic-xyz-ref.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

anisotropic-xyz-ref-outputstrain: anisotropic-xyz-ref-outputstrain.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

anisotropic-xyz-strain: anisotropic-xyz-strain.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

map: map.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

map-edge: map-edge.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

bcc: bcc.C ${INCLUDES}
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) $< -o $@ -lm

.SUFFIXES: .C .c .o

.c.o:
	$(CC) $(CFLAGS) -I$(INCLUDE) -c $< -o $@ -lm

.C.o:
	$(CCPP) $(CPPFLAGS) -I$(INCLUDE) -c $< -o $@ -lm

clean:
	rm -f *.o *.a
