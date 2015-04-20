BIN    = .
CC     = gcc
CFLAGS = -O3 -Wall -DTHREADS -fstrict-aliasing -ffast-math -msse2
LIBS   = -lm -lpthread
DEPS   = defs.h
PROGS  = $(BIN)/HighFCM
OBJS   = mem.o ac.o common.o context.o
all:
	$(MAKE) progs
progs: $(PROGS)
$(BIN)/HighFCM: HighFCM.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/HighFCM HighFCM.c $(OBJS) $(LIBS)
ac.o: ac.c ac.h $(DEPS)
	$(CC) -c $(CFLAGS) ac.c
mem.o: mem.c mem.h $(DEPS)
	$(CC) -c $(CFLAGS) mem.c
common.o: common.c common.h $(DEPS)
	$(CC) -c $(CFLAGS) common.c
context.o: context.c context.h $(DEPS)
	$(CC) -c $(CFLAGS) context.c
clean:
	/bin/rm -f *.o
cleanall:
	/bin/rm -f *.o $(PROGS)

