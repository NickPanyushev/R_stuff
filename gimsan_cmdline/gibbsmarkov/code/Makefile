#CC = g++ -O3 -pg -g #gprof
#CC = g++ -O0 #valgrind
CC = g++ -O3
DEBUG = 0
VERBOSE = 0

INCDIRS = -I. 
LIBS = -lm

ifeq (${DEBUG}, 1)
	GDB = -ggdb
else
	GDB =
endif

#-O3		For optimization
#-Wall		To turn on "all warnings"
#-m32		The 32-bit environment sets int, long and pointer to 32 bits. 
#
CFLAGS = -Wall -m32 ${GDB} -D DEBUG=${DEBUG} -D VERBOSE=${VERBOSE} ${INCDIRS}

OBJS = main_entsamp_ilr.o entsamp_ilr.o gibbs_util.o dataset.o profile.o transprob.o\
print_results.o em_alg.o entsamp.o markov.o markov_score.o mt19937ar.o symbols.o\
zoops.o psprior.o

OBJS_RAND = main_rand.o mt19937ar.o

all: gibbsmarkov rand

.c.o .cpp.o: 
	${CC} ${CFLAGS} -c $<

gibbsmarkov: ${OBJS}
	${CC} ${CFLAGS} -o gibbsmarkov.out ${OBJS} ${LIBS}
	mv gibbsmarkov.out ../

rand: ${OBJS_RAND}
	${CC} ${CFLAGS} -o rand.out ${OBJS_RAND} ${LIBS}
	mv rand.out ../

clean: 
	@ \rm -f *.o depend

depend: 
	${CC} ${CFLAGS} -MM *.c *.cpp > depend

include depend
