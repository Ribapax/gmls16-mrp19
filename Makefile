compiler=gcc
rm=rm -rf
FLAGS=-lm

all: main.c common.o io.o linear_system.o lu_factorization.o Makefile
	$(compiler) -o invmat main.c common.o io.o linear_system.o lu_factorization.o $(FLAGS)

linear_system.o: linear_system/linear_system.c linear_system/linear_system.h Makefile 
	$(compiler) -c linear_system/linear_system.c linear_system/linear_system.h 

lu_factorization.o: linear_system/lu_factorization.c linear_system/lu_factorization.h
	${compiler} -c linear_system/lu_factorization.c linear_system/lu_factorization.h 

common.o: common/common.c common/common.h Makefile
	$(compiler) -c common/common.c common/common.h 

io.o: io/io.c io/io.h Makefile
	$(compiler) -c io/io.c io/io.h 

clean:
	$(rm) invmat *.o linear_system/*.o linear_system/*.h.gch io/*.o io/*.h.gch common/*.o common/*.h.gch