HH.o:	HeavyHittersBenchmark.cpp EfficientCount.h SimpleCount.h
	mpicxx -O3 -std=c++14 -c -o HH.o HeavyHittersBenchmark.cpp 

hash.o:
	mpicxx -O3 -std=c++14 -c -o hash.o hash.cpp

HHbench: HH.o hash.o 
	mpicxx -O3 -std=c++14 -o HHbench HH.o hash.o

clean:
	rm *.o
	rm -rf HHbench

