INCLUDES = -I.
AR = ar rv
RANLIB = ranlib
# add path to c compiler 
CC = icc 
# add c compiler optimization flags
OPTFLAGS = -O2
# Include custom compiler and other configuration in local.mk, which is
# ignored by git
ifneq ("$(wildcard local.mk))","")
        include local.mk
endif

CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
LD = $(CC) -L. 


OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o memory.o \
       pmetis.o pqueue.o refine.o util.o timing.o debug.o \
       bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
       kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
       mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
       myqsort.o compress.o parmetis.o estmem.o \
       mpmetis.o mcoarsen.o mmatch.o minitpart.o mbalance.o \
       mrefine.o mutil.o mfm.o mkmetis.o mkwayrefine.o mkwayfmh.o \
       mrefine2.o minitpart2.o mbalance2.o mfm2.o \
       kvmetis.o kwayvolrefine.o kwayvolfm.o subdomains.o partition.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

libmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f libmetis.a
