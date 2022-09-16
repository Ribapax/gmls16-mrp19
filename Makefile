compiler=gcc -g -O3 -mavx -march=native -ftree-vectorize -fopt-info-vec-optimized -fopenmp -pthread -Ofast -ffast-math -fassociative-math
rm=rm -rf
FLAGS=-lm -L${LIKWID_LIB} -llikwid

all: main.c common.o io.o linear_system.o lu_factorization.o Makefile
	$(compiler) -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -o invmat main.c common.o io.o linear_system.o lu_factorization.o $(FLAGS)

linear_system.o: linear_system/linear_system.c linear_system/linear_system.h Makefile 
	$(compiler) -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c linear_system/linear_system.c linear_system/linear_system.h

lu_factorization.o: lu_factorization/lu_factorization.c lu_factorization/lu_factorization.h
	${compiler} -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c lu_factorization/lu_factorization.c lu_factorization/lu_factorization.h

common.o: common/common.c common/common.h Makefile
	$(compiler) -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c common/common.c common/common.h

io.o: io/io.c io/io.h Makefile
	$(compiler) -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c io/io.c io/io.h

clean:
	$(rm) invmat *.o linear_system/*.o linear_system/*.h.gch io/*.o io/*.h.gch common/*.o common/*.h.gch lu_factorization/*.o lu_factorization/*.h.gch