#CC = g++ -O3 -pg -g #gprof
#CC = g++ -O0 -ggdb #valgrind
CC = g++ -O3
DEBUG = 0
VERBOSE = 0

#INCDIRS = -I. -I/home/ppn3/boost_1_35_0 
#INCDIRS = -I. -I/Users/png/boost_1_35_0 
LIBS = -lm

ifeq (${DEBUG}, 1)
	GDB = -ggdb
else
	GDB =
endif

#-O3		For optimization
#-Wall		To turn on "all warnings"
#-m32		The 32-bit environment sets int, long and pointer to 32 bits. 
CFLAGS = -Wall -m32 ${GDB} -D DEBUG=${DEBUG} -D VERBOSE=${VERBOSE} ${INCDIRS}

OBJS = main_col_depend.o binomial_distribution.o col_depend.o matrix.o random.o mt19937ar.o dataset.o params.o results.o symbols.o siteloc.o dcdflib.o
TEST_OBJS =  testmain.o binomial_distribution.o random.o mt19937ar.o dcdflib.o

all: column_dependency testmain

.c.o: 
	${CC} ${CFLAGS} -c $*.c

column_dependency: ${OBJS}
	${CC} ${CFLAGS} -o column_dependency.out ${OBJS} ${LIBS}
	mv column_dependency.out ../

testmain: ${TEST_OBJS}
	${CC} ${CFLAGS} -o testmain.out ${TEST_OBJS} ${LIBS}
	mv testmain.out ../

clean: 
	@ \rm *.o 

random.o: random.h
col_depend.o: col_depend.h random.h binomial_distribution.h
main_col_depend.o: col_depend.h
results.o: results.h
params.o: params.h
symbols.o: symbols.h
siteloc.o: siteloc.h
matrix.o: matrix.h
testmain.o: binomial_distribution.h random.h 
binomial_distribution.o: binomial_distribution.h dcdflib.h

